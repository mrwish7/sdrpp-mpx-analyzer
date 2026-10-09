#pragma once
#include <dsp/demod/quadrature.h>
#include <dsp/filter/decimating_fir.h>
#include <dsp/taps/low_pass.h>
#include <dsp/sink/handler_sink.h>
#include <dsp/window/blackman_harris.h>
#include <dsp/window/hann.h>
#include <dsp/buffer/buffer.h>
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

// Anti-alias filter for the /2 decimation. Windowed sinc so the passband is flat for level measurements
#define MPX_DECIM_CUTOFF    91000.0
#define MPX_DECIM_TRANS     10000.0

#define MPX_FFT_SIZE        4096
#define MPX_FFT_BINS        ((MPX_FFT_SIZE / 2) + 1)
#define MPX_FFT_HOP         (MPX_FFT_SIZE / 2)

// Split recordings before the 4GiB limit of the 32bit WAV size fields
#define MPX_MAX_WAV_BYTES   4000000000ULL

// Level metering is done in 100ms blocks
#define MPX_METER_BLOCK         19200
#define MPX_POWER_LONG_BLOCKS   600     // 60s integration for MPX power (ITU-R BS.412)
#define MPX_POWER_SHORT_BLOCKS  10      // 1s short term MPX power
#define MPX_PEAK_BLOCKS         10      // Peak deviation over the last 1s
#define MPX_SETTLE_BLOCKS       5       // Blocks ignored after a reset while the DC filter settles

// 0dBr = power of a sine giving +/-19kHz deviation (ITU-R BS.412), in demodulator units (1.0 = 75kHz)
#define MPX_POWER_REF           (((19000.0 / MPX_DEVIATION) * (19000.0 / MPX_DEVIATION)) / 2.0)

// One pole DC removal (~3Hz) so a slightly off-centre carrier doesn't affect the measurements
#define MPX_DC_ALPHA            1e-4f

// Subcarrier levels. 19k and 57k both have a whole number of cycles in this many samples at 192kHz
#define MPX_TONE_TABLE          384
#define MPX_PILOT_FREQ          19000.0
#define MPX_RDS_FREQ            57000.0

// Pilot amplitude from 10ms Hann windowed correlations (~ +/-100Hz detector)
#define MPX_PILOT_BLOCK         1920
#define MPX_PILOT_MIN_KHZ       0.5f    // Below this the station is considered to have no pilot

// RDS is shifted to 0Hz, filtered to its bandwidth and decimated to 12kHz
#define MPX_RDS_DECIM           16
#define MPX_RDS_CUTOFF          3000.0  // Flat to 2.4kHz (RDS bandwidth), rejects L-R audio from 4kHz below 57kHz
#define MPX_RDS_TRANS           1200.0
#define MPX_RDS_CHUNK           4096
#define MPX_RDS_SUB_BLOCK       120     // 10ms at 12kHz
#define MPX_RDS_MIN_KHZ         0.3f
// Biphase coded RDS has almost no energy at exactly 57kHz while noise does. Noise power in the band is
// estimated from the power of 10ms means (~ +/-50Hz around 57kHz), see computeNoiseDcFraction()
// RDS is considered present when its power is at least this many times the noise power in its band
#define MPX_RDS_MIN_SNR         2.0
// Peak/RMS ratio of a standard RDS signal (EN 50067 biphase symbol shaping), to report peak deviation
#define MPX_RDS_CREST_FACTOR    1.444

struct MPXMeasurements {
    float powerDBr;         // MPX power over the last 60s (or since reset, if shorter)
    float powerShortDBr;    // MPX power over the last 1s
    float powerSeconds;     // Seconds integrated in powerDBr (up to 60)
    float peakDevKHz;       // Peak deviation over the last 1s
    float maxDevKHz;        // Peak deviation since reset
    float pilotKHz;         // 19kHz pilot level (deviation) over the last 1s
    bool pilotPresent;
    float rdsKHz;           // RDS level (peak deviation, from noise corrected RMS) over the last 1s
    bool rdsPresent;
};

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

        // Decimation filter, normalized to exactly unity gain
        decimTaps = dsp::taps::lowPass(MPX_DECIM_CUTOFF, MPX_DECIM_TRANS, MPX_IF_SAMPLERATE, true);
        double tapSum = 0.0;
        for (int i = 0; i < decimTaps.size; i++) { tapSum += decimTaps.taps[i]; }
        for (int i = 0; i < decimTaps.size; i++) { decimTaps.taps[i] /= tapSum; }

        // Subcarrier reference tables
        for (int i = 0; i < MPX_TONE_TABLE; i++) {
            double t = (double)i / MPX_SAMPLERATE;
            pilotCos[i] = cos(2.0 * FL_M_PI * MPX_PILOT_FREQ * t);
            pilotSin[i] = sin(2.0 * FL_M_PI * MPX_PILOT_FREQ * t);
            rdsCos[i] = cos(2.0 * FL_M_PI * MPX_RDS_FREQ * t);
            rdsSin[i] = sin(2.0 * FL_M_PI * MPX_RDS_FREQ * t);
        }
        pilotWin.resize(MPX_PILOT_BLOCK);
        pilotWinSum = 0.0;
        for (int i = 0; i < MPX_PILOT_BLOCK; i++) {
            pilotWin[i] = dsp::window::hann(i, MPX_PILOT_BLOCK);
            pilotWinSum += pilotWin[i];
        }

        // RDS channel filter, only used through process() (not run as a block)
        rdsTaps = dsp::taps::lowPass(MPX_RDS_CUTOFF, MPX_RDS_TRANS, MPX_SAMPLERATE);
        rdsFir.init(NULL, rdsTaps, MPX_RDS_DECIM);
        rdsFir.out.free();
        rdsMix = dsp::buffer::alloc<dsp::complex_t>(MPX_RDS_CHUNK);
        rdsBase = dsp::buffer::alloc<dsp::complex_t>(MPX_RDS_CHUNK);
        rdsNoiseDcFraction = computeDcFraction(false);
        rdsSignalDcFraction = computeDcFraction(true);
    }

    ~MPXChain() {
        stop();
        stopRecording();
        fftwf_destroy_plan(plan);
        fftwf_free(fftIn);
        fftwf_free(fftOut);
        dsp::taps::free(decimTaps);
        dsp::taps::free(rdsTaps);
        dsp::buffer::free(rdsMix);
        dsp::buffer::free(rdsBase);
    }

    void init(dsp::stream<dsp::complex_t>* in) {
        demod.init(in, MPX_DEVIATION, MPX_IF_SAMPLERATE);
        decim.init(&demod.out, decimTaps, (int)(MPX_IF_SAMPLERATE / MPX_SAMPLERATE));
        sink.init(&decim.out, handler, this);
    }

    void setInput(dsp::stream<dsp::complex_t>* in) {
        demod.setInput(in);
    }

    void start() {
        if (running) { return; }
        demod.start();
        decim.start();
        sink.start();
        running = true;
    }

    void stop() {
        if (!running) { return; }
        demod.stop();
        decim.stop();
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

    // ---- Measurements ----

    // Restart the MPX power integration and peak deviation hold (e.g. after retuning)
    void resetMeasurements() {
        measReset = true;
    }

    // Returns false while there isn't any measurement yet (just started or reset)
    bool getMeasurements(MPXMeasurements& meas) {
        std::lock_guard<std::mutex> lck(measMtx);
        if (!measValid) { return false; }
        meas = measOut;
        return true;
    }

    // ---- Audio output ----

    // Sends the 192kHz MPX (copied to both channels) to `out`, or stops when NULL
    void setAudioOutput(dsp::stream<dsp::stereo_t>* out) {
        // Release the DSP thread if it's waiting on the old output
        dsp::stream<dsp::stereo_t>* old = audioOut;
        if (old) { old->stopWriter(); }

        std::lock_guard<std::mutex> lck(audioMtx);
        if (old) { old->clearWriteStop(); }
        audioOut = out;
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
        _this->output(data, count);
        _this->analyze(data, count);
        _this->measure(data, count);
    }

    void output(float* data, int count) {
        std::lock_guard<std::mutex> lck(audioMtx);
        if (!audioOut) { return; }
        for (int i = 0; i < count; i++) { audioOut->writeBuf[i] = { data[i], data[i] }; }
        audioOut->swap(count);
    }

    void measure(float* data, int count) {
        if (measReset.exchange(false)) {
            powerFill = 0;
            powerPos = 0;
            peakFill = 0;
            peakPos = 0;
            maxPeak = 0.0f;
            blkSumSq = 0.0;
            blkPeak = 0.0f;
            blkCount = 0;
            settleBlocks = MPX_SETTLE_BLOCKS;
            pilotI = pilotQ = 0.0;
            pilotCount = 0;
            pilotAmpSum = 0.0;
            pilotAmpCount = 0;
            rdsSubSum = { 0.0f, 0.0f };
            rdsSubPow = 0.0;
            rdsSubCount = 0;
            rdsDcPow = rdsTotPow = 0.0;
            rdsSubBlocks = 0;
            std::lock_guard<std::mutex> lck(measMtx);
            measValid = false;
        }

        measureRDS(data, count);

        for (int i = 0; i < count; i++) {
            dcEst += MPX_DC_ALPHA * (data[i] - dcEst);
            float x = data[i] - dcEst;
            blkSumSq += x * x;
            blkPeak = std::max<float>(blkPeak, fabsf(x));

            // Pilot: correlate with 19kHz over a windowed block
            float wx = pilotWin[pilotCount] * x;
            pilotI += wx * pilotCos[pilotPos];
            pilotQ -= wx * pilotSin[pilotPos];
            pilotPos = (pilotPos + 1) % MPX_TONE_TABLE;
            if (++pilotCount == MPX_PILOT_BLOCK) {
                pilotAmpSum += 2.0 * sqrt((pilotI * pilotI) + (pilotQ * pilotQ)) / pilotWinSum;
                pilotAmpCount++;
                pilotI = pilotQ = 0.0;
                pilotCount = 0;
            }

            if (++blkCount == MPX_METER_BLOCK) { finishBlock(); }
        }
    }

    // Fraction of the power (after the RDS filter) that ends up in the 10ms means, for white noise (rds = false)
    // or for a standard RDS signal (rds = true)
    double computeDcFraction(bool rds) {
        const double outRate = MPX_SAMPLERATE / MPX_RDS_DECIM;
        const double td = 1.0 / 1187.5;
        const int points = 4096;
        double total = 0.0, dc = 0.0;
        for (int p = 0; p < points; p++) {
            // Frequencies across the decimated band (-outRate/2 .. outRate/2)
            double f = ((p + 0.5) / points - 0.5) * outRate;

            // Filter response
            double re = 0.0, im = 0.0;
            for (int i = 0; i < rdsTaps.size; i++) {
                double ph = -2.0 * FL_M_PI * f * i / MPX_SAMPLERATE;
                re += rdsTaps.taps[i] * cos(ph);
                im += rdsTaps.taps[i] * sin(ph);
            }
            double h2 = (re * re) + (im * im);

            // EN 50067 RDS spectrum: biphase symbol pair times the cos shaping, zero above 2/td
            if (rds) {
                double af = fabs(f);
                double shaping = (af <= 2.0 / td) ? cos(FL_M_PI * af * td / 4.0) : 0.0;
                double biphase = sin(FL_M_PI * af * td / 2.0);
                h2 *= (shaping * shaping) * (biphase * biphase);
            }

            // Response of the mean over MPX_RDS_SUB_BLOCK samples
            double x = FL_M_PI * f / outRate;
            double d = (fabs(sin(x)) < 1e-12) ? 1.0 : sin(MPX_RDS_SUB_BLOCK * x) / (MPX_RDS_SUB_BLOCK * sin(x));

            total += h2;
            dc += h2 * d * d;
        }
        return dc / total;
    }

    void measureRDS(float* data, int count) {
        for (int off = 0; off < count; off += MPX_RDS_CHUNK) {
            int n = std::min<int>(MPX_RDS_CHUNK, count - off);

            // Shift 57kHz to 0Hz, then filter and decimate
            for (int i = 0; i < n; i++) {
                float x = data[off + i];
                rdsMix[i] = { x * rdsCos[rdsPos], -x * rdsSin[rdsPos] };
                rdsPos = (rdsPos + 1) % MPX_TONE_TABLE;
            }
            int outCount = rdsFir.process(n, rdsMix, rdsBase);

            for (int i = 0; i < outCount; i++) {
                dsp::complex_t z = rdsBase[i];
                float pow = (z.re * z.re) + (z.im * z.im);

                // Power close to 57kHz (10ms means) vs total power, to tell RDS from noise
                rdsSubSum.re += z.re;
                rdsSubSum.im += z.im;
                rdsSubPow += pow;
                if (++rdsSubCount == MPX_RDS_SUB_BLOCK) {
                    float mre = rdsSubSum.re / MPX_RDS_SUB_BLOCK;
                    float mim = rdsSubSum.im / MPX_RDS_SUB_BLOCK;
                    rdsDcPow += (mre * mre) + (mim * mim);
                    rdsTotPow += rdsSubPow / MPX_RDS_SUB_BLOCK;
                    rdsSubBlocks++;
                    rdsSubSum = { 0.0f, 0.0f };
                    rdsSubPow = 0.0;
                    rdsSubCount = 0;
                }
            }
        }
    }

    void finishBlock() {
        double meanSq = blkSumSq / (double)MPX_METER_BLOCK;
        float peak = blkPeak;
        blkSumSq = 0.0;
        blkPeak = 0.0f;
        blkCount = 0;

        float pilotAmp = pilotAmpCount ? (pilotAmpSum / pilotAmpCount) : 0.0f;
        pilotAmpSum = 0.0;
        pilotAmpCount = 0;
        double rdsDc = rdsSubBlocks ? (rdsDcPow / rdsSubBlocks) : 0.0;
        double rdsTot = rdsSubBlocks ? (rdsTotPow / rdsSubBlocks) : 0.0;
        rdsDcPow = rdsTotPow = 0.0;
        rdsSubBlocks = 0;

        // Let the DC filter settle after a reset (e.g. a retune) before measuring
        if (settleBlocks > 0) {
            settleBlocks--;
            return;
        }

        powerRing[powerPos] = meanSq;
        powerPos = (powerPos + 1) % MPX_POWER_LONG_BLOCKS;
        powerFill = std::min<int>(powerFill + 1, MPX_POWER_LONG_BLOCKS);

        peakRing[peakPos] = peak;
        pilotRing[peakPos] = pilotAmp;
        rdsDcRing[peakPos] = rdsDc;
        rdsTotRing[peakPos] = rdsTot;
        peakPos = (peakPos + 1) % MPX_PEAK_BLOCKS;
        peakFill = std::min<int>(peakFill + 1, MPX_PEAK_BLOCKS);
        maxPeak = std::max<float>(maxPeak, peak);

        // Average over the most recent blocks
        double longSum = 0.0, shortSum = 0.0;
        int shortCount = std::min<int>(powerFill, MPX_POWER_SHORT_BLOCKS);
        for (int k = 1; k <= powerFill; k++) {
            double p = powerRing[(powerPos - k + MPX_POWER_LONG_BLOCKS) % MPX_POWER_LONG_BLOCKS];
            longSum += p;
            if (k <= shortCount) { shortSum += p; }
        }
        float recentPeak = 0.0f;
        double pilotSum = 0.0, rdsDcSum = 0.0, rdsTotSum = 0.0;
        for (int k = 0; k < peakFill; k++) {
            recentPeak = std::max<float>(recentPeak, peakRing[k]);
            pilotSum += pilotRing[k];
            rdsDcSum += rdsDcRing[k];
            rdsTotSum += rdsTotRing[k];
        }

        MPXMeasurements meas;
        meas.powerDBr = 10.0 * log10(std::max<double>(longSum / powerFill, 1e-20) / MPX_POWER_REF);
        meas.powerShortDBr = 10.0 * log10(std::max<double>(shortSum / shortCount, 1e-20) / MPX_POWER_REF);
        meas.powerSeconds = (powerFill * MPX_METER_BLOCK) / MPX_SAMPLERATE;
        meas.peakDevKHz = recentPeak * (MPX_DEVIATION / 1000.0);
        meas.maxDevKHz = maxPeak * (MPX_DEVIATION / 1000.0);
        meas.pilotKHz = (pilotSum / peakFill) * (MPX_DEVIATION / 1000.0);
        meas.pilotPresent = (meas.pilotKHz >= MPX_PILOT_MIN_KHZ);

        // RDS: split the band power into signal and noise using how much of each falls close to 57kHz:
        //   total = signal + noise,  dc = rdsSignalDcFraction * signal + rdsNoiseDcFraction * noise
        // then RMS -> peak. A DSB signal a(t)cos(wt) becomes a(t)/2 at 0Hz
        double rdsTotal = rdsTotSum / peakFill;
        double rdsDcPower = rdsDcSum / peakFill;
        double rdsSignal = std::clamp<double>(((rdsNoiseDcFraction * rdsTotal) - rdsDcPower) / (rdsNoiseDcFraction - rdsSignalDcFraction), 0.0, rdsTotal);
        double rdsNoise = rdsTotal - rdsSignal;
        meas.rdsKHz = 2.0 * sqrt(rdsSignal) * MPX_RDS_CREST_FACTOR * (MPX_DEVIATION / 1000.0);
        meas.rdsPresent = (meas.rdsKHz >= MPX_RDS_MIN_KHZ) && (rdsSignal >= MPX_RDS_MIN_SNR * rdsNoise);

        std::lock_guard<std::mutex> lck(measMtx);
        measOut = meas;
        measValid = true;
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
    dsp::tap<float> decimTaps;
    dsp::filter::DecimatingFIR<float, float> decim;
    dsp::sink::Handler<float> sink;
    bool running = false;

    // Audio output
    std::mutex audioMtx;
    dsp::stream<dsp::stereo_t>* audioOut = NULL;

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

    // Measurements (DSP thread only)
    float dcEst = 0.0f;
    double blkSumSq = 0.0;
    float blkPeak = 0.0f;
    int blkCount = 0;
    int settleBlocks = MPX_SETTLE_BLOCKS;
    double powerRing[MPX_POWER_LONG_BLOCKS];
    int powerPos = 0;
    int powerFill = 0;
    float peakRing[MPX_PEAK_BLOCKS];
    int peakPos = 0;
    int peakFill = 0;
    float maxPeak = 0.0f;
    std::atomic<bool> measReset = false;

    // Pilot and RDS level (DSP thread only)
    float pilotCos[MPX_TONE_TABLE];
    float pilotSin[MPX_TONE_TABLE];
    float rdsCos[MPX_TONE_TABLE];
    float rdsSin[MPX_TONE_TABLE];
    std::vector<float> pilotWin;
    double pilotWinSum;
    int pilotPos = 0;
    double pilotI = 0.0, pilotQ = 0.0;
    int pilotCount = 0;
    double pilotAmpSum = 0.0;
    int pilotAmpCount = 0;
    float pilotRing[MPX_PEAK_BLOCKS];

    dsp::tap<float> rdsTaps;
    dsp::filter::DecimatingFIR<dsp::complex_t, float> rdsFir;
    dsp::complex_t* rdsMix;
    dsp::complex_t* rdsBase;
    int rdsPos = 0;
    dsp::complex_t rdsSubSum = { 0.0f, 0.0f };
    double rdsSubPow = 0.0;
    int rdsSubCount = 0;
    double rdsDcPow = 0.0, rdsTotPow = 0.0;
    int rdsSubBlocks = 0;
    double rdsNoiseDcFraction;
    double rdsSignalDcFraction;
    double rdsDcRing[MPX_PEAK_BLOCKS];
    double rdsTotRing[MPX_PEAK_BLOCKS];

    // Measurements (shared with the GUI)
    std::mutex measMtx;
    MPXMeasurements measOut;
    bool measValid = false;

    // Spectrum (shared with the GUI)
    std::mutex specMtx;
    std::vector<float> outAvg;
    std::vector<float> outPeak;
    bool specValid = false;
};
