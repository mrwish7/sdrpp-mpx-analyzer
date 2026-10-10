#pragma once
#include <dsp/demod/quadrature.h>
#include <dsp/filter/decimating_fir.h>
#include <dsp/taps/low_pass.h>
#include <dsp/sink/handler_sink.h>
#include <dsp/window/blackman_harris.h>
#include <dsp/window/hann.h>
#include <dsp/window/nuttall.h>
#include <dsp/window/blackman.h>
#include <dsp/channel/rx_vfo.h>
#include <dsp/multirate/power_decimator.h>
#include <dsp/filter/fir.h>
#include <dsp/taps/estimate_tap_count.h>
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
#include <complex>

// IQ rate delivered by the VFO. 384k decimates by exactly 2 to the 192k MPX rate
#define MPX_IF_SAMPLERATE   384000.0
#define MPX_SAMPLERATE      192000.0

// Deviation that maps to +/-1.0 at the demodulator output (broadcast FM full deviation)
#define MPX_DEVIATION       75000.0

// Anti-alias filter for the /2 decimation. Windowed sinc so the passband is flat for level measurements
#define MPX_DECIM_CUTOFF    91000.0
#define MPX_DECIM_TRANS     10000.0

// VFO equalizer: probe tones every 2kHz to +/-160kHz, corrected up to 150kHz (the FM signal's extent)
#define MPX_EQ_PROBE_STEP   2000.0
#define MPX_EQ_PROBE_TONES  80
#define MPX_EQ_PROBE_LEN    3072    // Output samples analysed (tones fall exactly on bins)
#define MPX_EQ_MAX_FREQ     150000.0
#define MPX_EQ_TAPS         401

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
// Peak/RMS ratio used to report the RDS level as peak deviation. The theoretical EN 50067 waveform gives
// 1.444, but real encoders and reference instruments read ~4.5% higher: a calibrated transmitter checked
// with MPX Tool implies 1.51-1.52, and a Pira P175-calibrated receiver (VibeSDR, 21 UK stations) 1.507
#define MPX_RDS_CREST_FACTOR    1.51

// RDS to pilot phase lock: the 57kHz carrier phase is compared to 3x the pilot phase every 10ms and the
// result averaged over this many 100ms blocks. Locked when the phase difference is stable enough
#define MPX_LOCK_BLOCKS         50      // 5s
#define MPX_LOCK_MIN_BLOCKS     10      // Don't report before 1s of data
#define MPX_LOCK_MIN_COHERENCE  0.9

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
    bool lockValid;         // Pilot and RDS both present and enough data for the lock measurement
    bool locked;            // RDS 57kHz carrier locked to the 3rd harmonic of the pilot
    float lockCoherence;    // 0..1, stability of the RDS/pilot phase difference
    float lockPhaseDeg;     // RDS carrier phase relative to the 3rd pilot harmonic, -90..90 (0 = in phase, +/-90 = quadrature)
    float lockSeconds;      // Seconds of data in the lock measurement (up to 5)
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

        decimTaps = designDecimTaps();

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

        // Pilot at 0Hz through an identical filter, so its phase lines up exactly with the RDS samples
        lockPilotFir.init(NULL, rdsTaps, MPX_RDS_DECIM);
        lockPilotFir.out.free();
        lockPilotMix = dsp::buffer::alloc<dsp::complex_t>(MPX_RDS_CHUNK);
        lockPilotBase = dsp::buffer::alloc<dsp::complex_t>(MPX_RDS_CHUNK);
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
        dsp::taps::free(eqTaps);
        dsp::taps::free(rdsTaps);
        dsp::buffer::free(rdsMix);
        dsp::buffer::free(rdsBase);
        dsp::buffer::free(lockPilotMix);
        dsp::buffer::free(lockPilotBase);
    }

    // Low pass for the /2 decimation that also corrects the quadrature demodulator's droop. The demodulator
    // measures the phase change over one sample, i.e. the average frequency over 1/384000 s, which scales
    // modulation at frequency f by sinc(f / 384kHz): -0.4% at 19kHz, -3.6% at 57kHz. The passband is
    // shaped by the inverse of that (windowed frequency sampling design), and normalized to unity at DC
    static dsp::tap<float> designDecimTaps() {
        int count = dsp::taps::estimateTapCount(MPX_DECIM_TRANS, MPX_IF_SAMPLERATE);
        if (!(count % 2)) { count++; }
        dsp::tap<float> taps = dsp::taps::alloc<float>(count);

        const int steps = 4000;
        const double df = MPX_DECIM_CUTOFF / steps;
        const double half = (count - 1) / 2.0;
        double sum = 0.0;
        for (int n = 0; n < count; n++) {
            double t = n - half;
            double acc = 0.0;
            for (int k = 0; k < steps; k++) {
                double f = (k + 0.5) * df;
                double x = FL_M_PI * f / MPX_IF_SAMPLERATE;
                acc += (x / sin(x)) * cos(2.0 * FL_M_PI * f * t / MPX_IF_SAMPLERATE);
            }
            taps.taps[n] = (2.0 * acc * df / MPX_IF_SAMPLERATE) * dsp::window::nuttall(n, count - 1);
            sum += taps.taps[n];
        }
        for (int n = 0; n < count; n++) { taps.taps[n] /= sum; }
        return taps;
    }

    // ---- VFO equalizer ----
    //
    // SDR++'s VFO resampler (and front end decimator) has up to ~0.4dB of passband ripple, depending on
    // the SDR sample rate. FM sidebands are scaled by it, so e.g. a -0.4dB dip at +/-19kHz reads a pilot
    // 4.5% low. The response is measured by passing a comb of tones through a private copy of the same
    // SDR++ blocks, and the IQ is filtered by its inverse before the demodulator.

    // Measures the gain of SDR++'s IQ path (front end decimation + VFO resampling to 384kHz) at
    // k * MPX_EQ_PROBE_STEP Hz from the channel centre, relative to the centre. Returns +/- pairs averaged
    static std::vector<double> measureVfoResponse(double vfoInRate, int frontDecim) {
        const int tones = MPX_EQ_PROBE_TONES;
        const double rawRate = vfoInRate * frontDecim;

        dsp::multirate::PowerDecimator<dsp::complex_t> front;
        if (frontDecim > 1) { front.init(NULL, frontDecim); front.out.free(); }
        dsp::channel::RxVFO vfo;
        vfo.init(NULL, vfoInRate, MPX_IF_SAMPLERATE, MPX_IF_SAMPLERATE, 0.0);
        vfo.out.free();

        // Comb of equal tones with fixed pseudo random phases
        std::vector<double> toneFreq, tonePhase;
        uint32_t lcg = 12345;
        for (int k = -tones; k <= tones; k++) {
            toneFreq.push_back(k * MPX_EQ_PROBE_STEP);
            lcg = lcg * 1664525u + 1013904223u;
            tonePhase.push_back(2.0 * FL_M_PI * (lcg / 4294967296.0));
        }

        // The comb repeats every 1/MPX_EQ_PROBE_STEP seconds: compute one period when that is a whole
        // number of samples (all usual SDR rates), otherwise rotate one phasor per tone
        double amp = 1.0 / (2 * tones + 1);
        double periodExact = rawRate / MPX_EQ_PROBE_STEP;
        bool usePeriod = fabs(periodExact - round(periodExact)) < 1e-6;
        std::vector<dsp::complex_t> period;
        std::vector<std::complex<double>> phasor, step;
        if (usePeriod) {
            period.resize((size_t)round(periodExact));
            for (size_t i = 0; i < period.size(); i++) {
                double t = i / rawRate, re = 0.0, im = 0.0;
                for (size_t k = 0; k < toneFreq.size(); k++) {
                    double ph = 2.0 * FL_M_PI * toneFreq[k] * t + tonePhase[k];
                    re += cos(ph);
                    im += sin(ph);
                }
                period[i] = { (float)(re * amp), (float)(im * amp) };
            }
        }
        else {
            for (size_t k = 0; k < toneFreq.size(); k++) {
                phasor.push_back(std::polar(1.0, tonePhase[k]));
                step.push_back(std::polar(1.0, 2.0 * FL_M_PI * toneFreq[k] / rawRate));
            }
        }

        // Run long enough for every filter to settle, keep the last MPX_EQ_PROBE_LEN output samples
        const int chunk = 65536;
        long rawCount = (long)(rawRate * 0.25);
        std::vector<dsp::complex_t> in(chunk), mid(chunk), out(chunk);
        std::vector<dsp::complex_t> tail;
        for (long pos = 0; pos < rawCount; pos += chunk) {
            int n = (int)std::min<long>(chunk, rawCount - pos);
            for (int i = 0; i < n; i++) {
                if (usePeriod) {
                    in[i] = period[(pos + i) % period.size()];
                    continue;
                }
                std::complex<double> s = 0.0;
                for (size_t k = 0; k < phasor.size(); k++) {
                    s += phasor[k];
                    phasor[k] *= step[k];
                }
                in[i] = { (float)(s.real() * amp), (float)(s.imag() * amp) };
            }
            int m = n;
            dsp::complex_t* vin = in.data();
            if (frontDecim > 1) { m = front.process(n, in.data(), mid.data()); vin = mid.data(); }
            int outCount = vfo.process(m, vin, out.data());
            tail.insert(tail.end(), out.begin(), out.begin() + outCount);
            if ((int)tail.size() > 4 * MPX_EQ_PROBE_LEN) { tail.erase(tail.begin(), tail.end() - MPX_EQ_PROBE_LEN); }
        }
        if ((int)tail.size() > MPX_EQ_PROBE_LEN) { tail.erase(tail.begin(), tail.end() - MPX_EQ_PROBE_LEN); }

        // Amplitude of every tone in the output (tones sit exactly on DFT bins of the window)
        std::vector<double> gain(toneFreq.size());
        for (size_t k = 0; k < toneFreq.size(); k++) {
            double re = 0.0, im = 0.0;
            for (size_t i = 0; i < tail.size(); i++) {
                double ph = -2.0 * FL_M_PI * toneFreq[k] * i / MPX_IF_SAMPLERATE;
                re += tail[i].re * cos(ph) - tail[i].im * sin(ph);
                im += tail[i].re * sin(ph) + tail[i].im * cos(ph);
            }
            gain[k] = sqrt(re * re + im * im) / (tail.size() * amp);
        }

        // Average +f and -f, relative to the centre
        std::vector<double> resp(tones + 1);
        for (int k = 0; k <= tones; k++) { resp[k] = 0.5 * (gain[tones + k] + gain[tones - k]) / gain[tones]; }
        return resp;
    }

    // Linear phase FIR (real taps, for the complex IQ) with the inverse of a measured response
    static dsp::tap<float> designEqualizer(const std::vector<double>& resp) {
        const int count = MPX_EQ_TAPS;
        const double half = (count - 1) / 2.0;
        auto inverseAt = [&](double f) {
            double x = std::min<double>(fabs(f), MPX_EQ_MAX_FREQ) / MPX_EQ_PROBE_STEP;
            int i = std::min<int>((int)x, (int)resp.size() - 2);
            double g = resp[i] + (resp[i + 1] - resp[i]) * (x - i);
            return 1.0 / std::max<double>(g, 0.1);
        };

        dsp::tap<float> taps = dsp::taps::alloc<float>(count);
        const int steps = 1536;                         // 125Hz grid up to fs/2
        const double df = (MPX_IF_SAMPLERATE / 2.0) / steps;
        double sum = 0.0;
        for (int n = 0; n < count; n++) {
            double t = n - half;
            double acc = inverseAt(0.0);
            for (int k = 1; k < steps; k++) {
                double f = k * df;
                acc += 2.0 * inverseAt(f) * cos(2.0 * FL_M_PI * f * t / MPX_IF_SAMPLERATE);
            }
            taps.taps[n] = (acc * df / MPX_IF_SAMPLERATE) * dsp::window::blackman(n, count - 1);
            sum += taps.taps[n];
        }
        for (int n = 0; n < count; n++) { taps.taps[n] /= sum; }
        return taps;
    }

    // Replaces the equalizer (GUI thread). The chain starts with a pass-through
    void setEqualizer(dsp::tap<float> taps) {
        dsp::tap<float> old = eqTaps;
        eqTaps = taps;
        eq.setTaps(eqTaps);
        dsp::taps::free(old);
    }

    void init(dsp::stream<dsp::complex_t>* in) {
        eqTaps = dsp::taps::alloc<float>(1);
        eqTaps.taps[0] = 1.0f;
        eq.init(in, eqTaps);
        demod.init(&eq.out, MPX_DEVIATION, MPX_IF_SAMPLERATE);
        decim.init(&demod.out, decimTaps, (int)(MPX_IF_SAMPLERATE / MPX_SAMPLERATE));
        sink.init(&decim.out, handler, this);
    }

    void setInput(dsp::stream<dsp::complex_t>* in) {
        eq.setInput(in);
    }

    void start() {
        if (running) { return; }
        eq.start();
        demod.start();
        decim.start();
        sink.start();
        running = true;
    }

    void stop() {
        if (!running) { return; }
        eq.stop();
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

    // ---- Analysis on/off ----

    // Spectrum and measurements are only needed while they're displayed. Recording and audio output
    // are not affected. Called from the GUI thread
    void setAnalysis(bool enable) {
        if (enable && !analysisOn) {
            // Start fresh rather than mixing in data from before the pause
            measReset = true;
            spectrumReset = true;
        }
        analysisOn = enable;
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
        if (_this->analysisOn) {
            _this->analyze(data, count);
            _this->measure(data, count);
        }
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
            lockPilotSum = lockRdsSqSum = 0.0;
            lockVec = 0.0;
            lockWeight = 0.0;
            lockPos = 0;
            lockFill = 0;
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

            // Shift 57kHz (and the pilot, for the lock check) to 0Hz, then filter and decimate
            for (int i = 0; i < n; i++) {
                float x = data[off + i];
                rdsMix[i] = { x * rdsCos[rdsPos], -x * rdsSin[rdsPos] };
                lockPilotMix[i] = { x * pilotCos[rdsPos], -x * pilotSin[rdsPos] };
                rdsPos = (rdsPos + 1) % MPX_TONE_TABLE;
            }
            int outCount = rdsFir.process(n, rdsMix, rdsBase);
            lockPilotFir.process(n, lockPilotMix, lockPilotBase);

            for (int i = 0; i < outCount; i++) {
                dsp::complex_t z = rdsBase[i];
                float pow = (z.re * z.re) + (z.im * z.im);

                // Squaring removes the BPSK data, leaving twice the 57kHz carrier phase
                std::complex<double> zc(z.re, z.im);
                lockRdsSqSum += zc * zc;
                lockPilotSum += std::complex<double>(lockPilotBase[i].re, lockPilotBase[i].im);

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

                    // Phase of (RDS carrier)^2 relative to (3rd pilot harmonic)^2 = 6x pilot phase,
                    // weighted by the RDS strength in this 10ms
                    double pilotMag = std::abs(lockPilotSum);
                    if (pilotMag > 0.0) {
                        std::complex<double> u = lockPilotSum / pilotMag;
                        std::complex<double> u3 = u * u * u;
                        lockVec += lockRdsSqSum * std::conj(u3 * u3);
                        lockWeight += std::abs(lockRdsSqSum);
                    }
                    lockPilotSum = lockRdsSqSum = 0.0;
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
        std::complex<double> blkLockVec = lockVec;
        double blkLockWeight = lockWeight;
        lockVec = 0.0;
        lockWeight = 0.0;
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

        lockVecRing[lockPos] = blkLockVec;
        lockWeightRing[lockPos] = blkLockWeight;
        lockPos = (lockPos + 1) % MPX_LOCK_BLOCKS;
        lockFill = std::min<int>(lockFill + 1, MPX_LOCK_BLOCKS);
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

        // Lock: a stable phase difference gives a long average vector, a drifting one averages out
        std::complex<double> lockSum = 0.0;
        double lockWeightSum = 0.0;
        for (int k = 0; k < lockFill; k++) {
            lockSum += lockVecRing[k];
            lockWeightSum += lockWeightRing[k];
        }
        meas.lockSeconds = (lockFill * MPX_METER_BLOCK) / MPX_SAMPLERATE;
        meas.lockCoherence = (lockWeightSum > 0.0) ? (std::abs(lockSum) / lockWeightSum) : 0.0f;
        meas.lockValid = meas.pilotPresent && meas.rdsPresent && (lockFill >= MPX_LOCK_MIN_BLOCKS);
        meas.locked = meas.lockValid && (meas.lockCoherence >= MPX_LOCK_MIN_COHERENCE);

        // Halve the doubled phase. BPSK leaves a 180 degree ambiguity, so fold into -90..90
        double phase = (std::arg(lockSum) / 2.0) * (180.0 / FL_M_PI);
        if (phase >= 90.0) { phase -= 180.0; }
        if (phase < -90.0) { phase += 180.0; }
        meas.lockPhaseDeg = phase;

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
        if (spectrumReset.exchange(false)) {
            framePos = 0;
            firstFrame = true;
        }

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
    dsp::tap<float> eqTaps;
    dsp::filter::FIR<dsp::complex_t, float> eq;
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
    std::atomic<bool> spectrumReset = false;
    std::atomic<bool> analysisOn = true;

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

    // RDS to pilot phase lock (DSP thread only)
    dsp::filter::DecimatingFIR<dsp::complex_t, float> lockPilotFir;
    dsp::complex_t* lockPilotMix;
    dsp::complex_t* lockPilotBase;
    std::complex<double> lockPilotSum = 0.0;
    std::complex<double> lockRdsSqSum = 0.0;
    std::complex<double> lockVec = 0.0;
    double lockWeight = 0.0;
    std::complex<double> lockVecRing[MPX_LOCK_BLOCKS];
    double lockWeightRing[MPX_LOCK_BLOCKS];
    int lockPos = 0;
    int lockFill = 0;
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
