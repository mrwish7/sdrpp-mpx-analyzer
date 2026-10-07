#pragma once
#include <dsp/demod/quadrature.h>
#include <dsp/multirate/rational_resampler.h>
#include <dsp/sink/handler_sink.h>
#include <dsp/window/blackman_harris.h>
#include <utils/wav.h>
#include <utils/flog.h>
#include <fftw3.h>
#include <atomic>
#include <functional>
#include <mutex>
#include <string>
#include <vector>
#include <string.h>
#include <math.h>

// IQ rate delivered by the VFO. 384k decimates by exactly 2 to the 192k MPX rate
#define MPX_IF_SAMPLERATE   384000.0
#define MPX_SAMPLERATE      192000.0

// Deviation that maps to +/-1.0 at the demodulator output (broadcast FM full deviation)
#define MPX_DEVIATION       75000.0

#define MPX_FFT_SIZE        4096
#define MPX_FFT_BINS        ((MPX_FFT_SIZE / 2) + 1)
#define MPX_FFT_HOP         (MPX_FFT_SIZE / 2)

// Split recordings before the 4GiB limit of the 32bit WAV size fields
#define MPX_MAX_WAV_BYTES   4000000000ULL

// FM demodulator -> 192kHz MPX -> (WAV recorder, spectrum analyzer)
class MPXChain {
public:
    MPXChain() {
        // Window, normalized so that a full scale (+/-1.0 = 75kHz deviation) sine reads 0dB
        window.resize(MPX_FFT_SIZE);
        double sum = 0.0;
        for (int i = 0; i < MPX_FFT_SIZE; i++) {
            window[i] = dsp::window::blackmanHarris(i, MPX_FFT_SIZE);
            sum += window[i];
        }
        float norm = 2.0f / sum;
        for (auto& w : window) { w *= norm; }

        frame.resize(MPX_FFT_SIZE, 0.0f);
        avgSpec.resize(MPX_FFT_BINS, 0.0f);
        peakSpec.resize(MPX_FFT_BINS, 0.0f);
        outAvg.resize(MPX_FFT_BINS, 0.0f);
        outPeak.resize(MPX_FFT_BINS, 0.0f);

        fftIn = (float*)fftwf_malloc(MPX_FFT_SIZE * sizeof(float));
        fftOut = (fftwf_complex*)fftwf_malloc(MPX_FFT_BINS * sizeof(fftwf_complex));
        plan = fftwf_plan_dft_r2c_1d(MPX_FFT_SIZE, fftIn, fftOut, FFTW_ESTIMATE);
    }

    ~MPXChain() {
        stop();
        stopRecording();
        fftwf_destroy_plan(plan);
        fftwf_free(fftIn);
        fftwf_free(fftOut);
    }

    void init(dsp::stream<dsp::complex_t>* in) {
        demod.init(in, MPX_DEVIATION, MPX_IF_SAMPLERATE);
        resamp.init(&demod.out, MPX_IF_SAMPLERATE, MPX_SAMPLERATE);
        sink.init(&resamp.out, handler, this);
    }

    void setInput(dsp::stream<dsp::complex_t>* in) {
        demod.setInput(in);
    }

    void start() {
        if (running) { return; }
        demod.start();
        resamp.start();
        sink.start();
        running = true;
    }

    void stop() {
        if (!running) { return; }
        demod.stop();
        resamp.stop();
        sink.stop();
        running = false;
    }

    // ---- Spectrum ----

    // Exponential averaging over roughly `frames` FFT frames (1 = no averaging)
    void setAveraging(int frames) {
        avgAlpha = 1.0f / (float)std::max<int>(frames, 1);
    }

    void resetPeak() {
        peakReset = true;
    }

    // Copies the latest linear power spectra (MPX_FFT_BINS each). Returns false if nothing computed yet
    bool getSpectrum(std::vector<float>& avg, std::vector<float>& peak) {
        std::lock_guard<std::mutex> lck(specMtx);
        if (!specValid) { return false; }
        avg = outAvg;
        peak = outPeak;
        return true;
    }

    // ---- Recording ----

    // `nextPath` is called to name the first file and each split file when the 4GiB limit is reached
    bool startRecording(std::function<std::string()> nextPath, wav::SampleType type) {
        std::lock_guard<std::mutex> lck(recMtx);
        if (recording) { return true; }
        pathGen = nextPath;
        writer.setSampleType(type);
        bytesPerSample = (type == wav::SAMP_TYPE_INT16) ? 2 : ((type == wav::SAMP_TYPE_UINT8) ? 1 : 4);
        currentPath = pathGen();
        if (!writer.open(currentPath)) {
            flog::error("MPX Analyzer: Failed to open file for recording: {0}", currentPath);
            return false;
        }
        flog::info("MPX Analyzer: Recording to {0}", currentPath);
        totalSamples = 0;
        recording = true;
        return true;
    }

    void stopRecording() {
        std::lock_guard<std::mutex> lck(recMtx);
        if (!recording) { return; }
        writer.close();
        recording = false;
    }

    bool isRecording() { return recording; }

    uint64_t getRecordedSamples() { return totalSamples; }

    std::string getRecordingPath() {
        std::lock_guard<std::mutex> lck(recMtx);
        return currentPath;
    }

private:
    static void handler(float* data, int count, void* ctx) {
        MPXChain* _this = (MPXChain*)ctx;
        _this->record(data, count);
        _this->analyze(data, count);
    }

    void record(float* data, int count) {
        std::lock_guard<std::mutex> lck(recMtx);
        if (!recording) { return; }

        // Roll over to a new file before exceeding the WAV size limit
        if ((writer.getSamplesWritten() + count) * bytesPerSample > MPX_MAX_WAV_BYTES) {
            writer.close();
            currentPath = pathGen();
            if (!writer.open(currentPath)) {
                flog::error("MPX Analyzer: Failed to open file for recording: {0}", currentPath);
                recording = false;
                return;
            }
            flog::info("MPX Analyzer: Continuing recording in {0}", currentPath);
        }

        writer.write(data, count);
        totalSamples += count;
    }

    void analyze(float* data, int count) {
        while (count > 0) {
            int n = std::min<int>(count, MPX_FFT_SIZE - framePos);
            memcpy(&frame[framePos], data, n * sizeof(float));
            framePos += n;
            data += n;
            count -= n;

            if (framePos < MPX_FFT_SIZE) { break; }
            computeFrame();

            // 50% overlap
            memmove(&frame[0], &frame[MPX_FFT_HOP], (MPX_FFT_SIZE - MPX_FFT_HOP) * sizeof(float));
            framePos = MPX_FFT_SIZE - MPX_FFT_HOP;
        }
    }

    void computeFrame() {
        for (int i = 0; i < MPX_FFT_SIZE; i++) { fftIn[i] = frame[i] * window[i]; }
        fftwf_execute(plan);

        float alpha = avgAlpha;
        bool resetPk = peakReset.exchange(false);
        for (int i = 0; i < MPX_FFT_BINS; i++) {
            float p = (fftOut[i][0] * fftOut[i][0]) + (fftOut[i][1] * fftOut[i][1]);
            avgSpec[i] = firstFrame ? p : (avgSpec[i] + alpha * (p - avgSpec[i]));
            peakSpec[i] = (firstFrame || resetPk) ? p : std::max<float>(peakSpec[i], p);
        }
        firstFrame = false;

        std::lock_guard<std::mutex> lck(specMtx);
        memcpy(outAvg.data(), avgSpec.data(), MPX_FFT_BINS * sizeof(float));
        memcpy(outPeak.data(), peakSpec.data(), MPX_FFT_BINS * sizeof(float));
        specValid = true;
    }

    // DSP
    dsp::demod::Quadrature demod;
    dsp::multirate::RationalResampler<float> resamp;
    dsp::sink::Handler<float> sink;
    bool running = false;

    // Recording
    std::mutex recMtx;
    wav::Writer writer{ 1, (uint64_t)MPX_SAMPLERATE, wav::FORMAT_WAV, wav::SAMP_TYPE_FLOAT32 };
    std::function<std::string()> pathGen;
    std::string currentPath;
    std::atomic<bool> recording = false;
    std::atomic<uint64_t> totalSamples = 0;
    int bytesPerSample = 4;

    // Spectrum (DSP thread only)
    std::vector<float> window;
    std::vector<float> frame;
    int framePos = 0;
    float* fftIn;
    fftwf_complex* fftOut;
    fftwf_plan plan;
    std::vector<float> avgSpec;
    std::vector<float> peakSpec;
    bool firstFrame = true;
    std::atomic<float> avgAlpha = 0.1f;
    std::atomic<bool> peakReset = false;

    // Spectrum (shared with the GUI)
    std::mutex specMtx;
    std::vector<float> outAvg;
    std::vector<float> outPeak;
    bool specValid = false;
};
