#include <imgui.h>
#include <module.h>
#include <config.h>
#include <core.h>
#include <gui/gui.h>
#include <gui/style.h>
#include <gui/widgets/folder_select.h>
#include <signal_path/signal_path.h>
#include <utils/optionlist.h>
#include <utils/flog.h>
#include <filesystem>
#include <regex>
#include <ctime>
#include "mpx_dsp.h"

#define CONCAT(a, b) ((std::string(a) + b).c_str())

#define MPX_MIN_BANDWIDTH       50000.0
#define MPX_MAX_BANDWIDTH       350000.0
#define MPX_DEFAULT_BANDWIDTH   250000.0

#define MPX_METER_FONT_SCALE    1.6f

SDRPP_MOD_INFO{
    /* Name:            */ "mpx_analyzer",
    /* Description:     */ "FM MPX spectrum analyzer and recorder for SDR++",
    /* Author:          */ "Adam Wisher",
    /* Version:         */ 0, 1, 0,
    /* Max instances    */ -1
};

ConfigManager config;

struct MPXMarker {
    double freq;
    const char* label;
};

const MPXMarker MPX_MARKERS[] = {
    { 19000.0, "Pilot" },
    { 38000.0, "38k" },
    { 57000.0, "RDS" },
    { 67000.0, "SCA" },
    { 92000.0, "SCA" }
};

class MPXAnalyzerModule : public ModuleManager::Instance {
public:
    MPXAnalyzerModule(std::string name) : folderSelect("%ROOT%/recordings") {
        this->name = name;
        root = (std::string)core::args["root"];

        // First entry is the default
        sampleTypes.define("int16", "Int16", wav::SAMP_TYPE_INT16);
        sampleTypes.define("int32", "Int32", wav::SAMP_TYPE_INT32);
        sampleTypes.define("float32", "Float32", wav::SAMP_TYPE_FLOAT32);

        loadConfig();

        // Create the VFO and DSP chain
        vfo = createVFO();
        chain.init(vfo->output);
        chain.setAveraging(avgFrames);
        chain.start();

        // Draw the window every frame (not only when the menu is open) and block waterfall input under it
        fftRedrawHandler.ctx = this;
        fftRedrawHandler.handler = fftRedraw;
        inputHandler.ctx = this;
        inputHandler.handler = inputProcess;
        gui::waterfall.onFFTRedraw.bindHandler(&fftRedrawHandler);
        gui::waterfall.onInputProcess.bindHandler(&inputHandler);

        gui::menu.registerEntry(name, menuHandler, this, this);
    }

    ~MPXAnalyzerModule() {
        gui::menu.removeEntry(name);
        gui::waterfall.onFFTRedraw.unbindHandler(&fftRedrawHandler);
        gui::waterfall.onInputProcess.unbindHandler(&inputHandler);
        chain.stopRecording();
        chain.stop();
        if (vfo) { sigpath::vfoManager.deleteVFO(vfo); }
    }

    void postInit() {}

    void enable() {
        vfo = createVFO();
        chain.setInput(vfo->output);
        chain.resetMeasurements();
        chain.start();
        enabled = true;
    }

    void disable() {
        chain.stopRecording();
        chain.stop();
        sigpath::vfoManager.deleteVFO(vfo);
        vfo = NULL;
        enabled = false;
    }

    bool isEnabled() {
        return enabled;
    }

private:
    void loadConfig() {
        config.acquire();
        bool created = false;
        if (!config.conf.contains(name)) {
            config.conf[name] = json::object();
            created = true;
        }
        json& conf = config.conf[name];
        if (conf.contains("followVfo")) { followName = conf["followVfo"]; }
        if (conf.contains("bandwidth")) { bandwidth = conf["bandwidth"]; }
        if (conf.contains("showWindow")) { showWindow = conf["showWindow"]; }
        if (conf.contains("avgFrames")) { avgFrames = conf["avgFrames"]; }
        if (conf.contains("maxFreq")) { maxFreq = conf["maxFreq"]; }
        if (conf.contains("dbMin")) { dbMin = conf["dbMin"]; }
        if (conf.contains("dbMax")) { dbMax = conf["dbMax"]; }
        if (conf.contains("peakHold")) { peakHold = conf["peakHold"]; }
        if (conf.contains("recPath")) { folderSelect.setPath(conf["recPath"]); }
        if (conf.contains("sampleType") && sampleTypes.keyExists(conf["sampleType"])) {
            sampleTypeId = sampleTypes.keyId(conf["sampleType"]);
        }
        config.release(created);

        bandwidth = std::clamp<double>(bandwidth, MPX_MIN_BANDWIDTH, MPX_MAX_BANDWIDTH);
        maxFreq = std::clamp<float>(maxFreq, 20.0f, MPX_SAMPLERATE / 2000.0f);
    }

    template <class T>
    void saveSetting(const char* key, T value) {
        config.acquire();
        config.conf[name][key] = value;
        config.release(true);
    }

    // Returns the followed VFO, or NULL when not following or it doesn't exist
    ImGui::WaterfallVFO* getFollowedVFO() {
        if (followName.empty() || followName == name) { return NULL; }
        auto it = gui::waterfall.vfos.find(followName);
        return (it != gui::waterfall.vfos.end()) ? it->second : NULL;
    }

    VFOManager::VFO* createVFO() {
        double offset = 0.0;
        double bw = bandwidth;
        ImGui::WaterfallVFO* followed = getFollowedVFO();
        if (followed) {
            offset = followed->centerOffset;
            bw = std::clamp<double>(followed->bandwidth, MPX_MIN_BANDWIDTH, MPX_MAX_BANDWIDTH);
        }
        VFOManager::VFO* v = sigpath::vfoManager.createVFO(name, ImGui::WaterfallVFO::REF_CENTER, offset, bw, MPX_IF_SAMPLERATE, MPX_MIN_BANDWIDTH, MPX_MAX_BANDWIDTH, true);
        v->setColor(IM_COL32(255, 160, 0, 40));
        return v;
    }

    // Keep our VFO on top of the followed one, with the same filter bandwidth
    void followVFO() {
        if (!enabled || !vfo) { return; }
        ImGui::WaterfallVFO* followed = getFollowedVFO();
        if (!followed) { return; }

        double offset = followed->centerOffset;
        if (vfo->wtfVFO->centerOffset != offset) { vfo->setCenterOffset(offset); }

        double bw = std::clamp<double>(followed->bandwidth, MPX_MIN_BANDWIDTH, MPX_MAX_BANDWIDTH);
        if (vfo->getBandwidth() != bw) { vfo->setBandwidth(bw); }

        // Our VFO can't be moved while following, so hand the selection back to the followed VFO
        if (gui::waterfall.selectedVFO == name) {
            gui::waterfall.selectedVFO = followName;
            gui::waterfall.selectedVFOChanged = true;
        }
    }

    // Measurements from a different station or filter setting shouldn't be mixed
    void checkRetune() {
        if (!enabled || !vfo) { return; }
        double freq = getFrequency();
        double bw = vfo->getBandwidth();
        if (freq != measFreq || bw != measBandwidth) {
            chain.resetMeasurements();
            measFreq = freq;
            measBandwidth = bw;
        }
    }

    void updateFollowList() {
        followList.clear();
        followList.define("", "None", "");
        for (auto const& [vfoName, wtfVfo] : gui::waterfall.vfos) {
            if (vfoName == name) { continue; }
            followList.define(vfoName, vfoName, vfoName);
        }
        // Keep a configured VFO selectable even if it doesn't exist (yet)
        if (!followList.keyExists(followName)) {
            followList.define(followName, followName + " (not found)", followName);
        }
        followId = followList.keyId(followName);
    }

    double getFrequency() {
        double freq = gui::waterfall.getCenterFrequency();
        if (vfo) { freq += vfo->wtfVFO->generalOffset; }
        return freq;
    }

    std::string expandString(std::string input) {
        input = std::regex_replace(input, std::regex("%ROOT%"), root);
        return std::regex_replace(input, std::regex("//"), "/");
    }

    std::string genFileName() {
        time_t now = time(0);
        tm* ltm = localtime(&now);
        char buf[256];
        snprintf(buf, sizeof(buf), "mpx_%.0lfHz_%04d%02d%02d_%02d%02d%02d", getFrequency(),
                 ltm->tm_year + 1900, ltm->tm_mon + 1, ltm->tm_mday, ltm->tm_hour, ltm->tm_min, ltm->tm_sec);
        std::string base = expandString(folderSelect.path + "/" + buf);

        // Don't overwrite an existing file (e.g. two splits within the same second)
        std::string path = base + ".wav";
        for (int i = 1; std::filesystem::exists(path); i++) {
            path = base + "_" + std::to_string(i) + ".wav";
        }
        return path;
    }

    void startRecording() {
        if (!folderSelect.pathIsValid()) { return; }
        chain.startRecording([this]() { return genFileName(); }, sampleTypes.value(sampleTypeId));
    }

    static void menuHandler(void* ctx) {
        MPXAnalyzerModule* _this = (MPXAnalyzerModule*)ctx;
        float menuWidth = ImGui::GetContentRegionAvail().x;

        if (!_this->enabled) { style::beginDisabled(); }

        // Tuning
        _this->updateFollowList();
        ImGui::LeftLabel("Follow VFO");
        ImGui::FillWidth();
        if (ImGui::Combo(CONCAT("##_mpx_follow_", _this->name), &_this->followId, _this->followList.txt)) {
            _this->followName = _this->followList.key(_this->followId);
            _this->saveSetting("followVfo", _this->followName);

            // Back to our own bandwidth when no longer following
            if (_this->vfo && !_this->getFollowedVFO()) { _this->vfo->setBandwidth(_this->bandwidth); }
        }

        // While following, the bandwidth comes from the followed VFO and is only displayed here
        bool following = (_this->getFollowedVFO() != NULL);
        float bwKHz = ((following && _this->vfo) ? _this->vfo->getBandwidth() : _this->bandwidth) / 1000.0;
        ImGui::LeftLabel("IF Bandwidth");
        ImGui::FillWidth();
        if (following) { style::beginDisabled(); }
        if (ImGui::SliderFloat(CONCAT("##_mpx_bw_", _this->name), &bwKHz, MPX_MIN_BANDWIDTH / 1000.0, MPX_MAX_BANDWIDTH / 1000.0, following ? "%.0f kHz (followed)" : "%.0f kHz")) {
            _this->bandwidth = std::round(bwKHz) * 1000.0;
            if (_this->vfo) { _this->vfo->setBandwidth(_this->bandwidth); }
            _this->saveSetting("bandwidth", _this->bandwidth);
        }
        if (following) { style::endDisabled(); }

        if (_this->enabled && sigpath::iqFrontEnd.getSampleRate() < MPX_IF_SAMPLERATE) {
            ImGui::TextColored(ImVec4(1.0f, 1.0f, 0.0f, 1.0f), "Source sample rate is below 384 kS/s");
        }

        // Display
        if (ImGui::Checkbox(CONCAT("Show MPX window##_mpx_show_", _this->name), &_this->showWindow)) {
            _this->saveSetting("showWindow", _this->showWindow);
        }

        ImGui::LeftLabel("Averaging");
        ImGui::FillWidth();
        if (ImGui::SliderInt(CONCAT("##_mpx_avg_", _this->name), &_this->avgFrames, 1, 100)) {
            _this->chain.setAveraging(_this->avgFrames);
            _this->saveSetting("avgFrames", _this->avgFrames);
        }

        ImGui::LeftLabel("Max Freq");
        ImGui::FillWidth();
        if (ImGui::SliderFloat(CONCAT("##_mpx_maxfreq_", _this->name), &_this->maxFreq, 20.0f, MPX_SAMPLERATE / 2000.0f, "%.0f kHz")) {
            _this->saveSetting("maxFreq", _this->maxFreq);
        }

        ImGui::LeftLabel("dB Max");
        ImGui::FillWidth();
        if (ImGui::SliderFloat(CONCAT("##_mpx_dbmax_", _this->name), &_this->dbMax, -60.0f, 20.0f, "%.0f dB")) {
            _this->dbMax = std::max<float>(_this->dbMax, _this->dbMin + 10.0f);
            _this->saveSetting("dbMax", _this->dbMax);
        }

        ImGui::LeftLabel("dB Min");
        ImGui::FillWidth();
        if (ImGui::SliderFloat(CONCAT("##_mpx_dbmin_", _this->name), &_this->dbMin, -160.0f, -20.0f, "%.0f dB")) {
            _this->dbMin = std::min<float>(_this->dbMin, _this->dbMax - 10.0f);
            _this->saveSetting("dbMin", _this->dbMin);
        }

        if (ImGui::Checkbox(CONCAT("Peak hold##_mpx_peak_", _this->name), &_this->peakHold)) {
            _this->chain.resetPeak();
            _this->saveSetting("peakHold", _this->peakHold);
        }
        ImGui::SameLine();
        if (ImGui::Button(CONCAT("Reset peak##_mpx_peak_reset_", _this->name))) {
            _this->chain.resetPeak();
        }

        // Recording
        bool recording = _this->chain.isRecording();
        if (recording) { style::beginDisabled(); }
        if (_this->folderSelect.render("##_mpx_rec_fold_" + _this->name)) {
            if (_this->folderSelect.pathIsValid()) {
                _this->saveSetting("recPath", _this->folderSelect.path);
            }
        }

        ImGui::LeftLabel("Sample type");
        ImGui::FillWidth();
        if (ImGui::Combo(CONCAT("##_mpx_st_", _this->name), &_this->sampleTypeId, _this->sampleTypes.txt)) {
            _this->saveSetting("sampleType", _this->sampleTypes.key(_this->sampleTypeId));
        }
        if (recording) { style::endDisabled(); }

        if (!recording) {
            bool canRecord = _this->folderSelect.pathIsValid();
            if (!canRecord) { style::beginDisabled(); }
            if (ImGui::Button(CONCAT("Record MPX##_mpx_rec_", _this->name), ImVec2(menuWidth, 0))) {
                _this->startRecording();
            }
            if (!canRecord) { style::endDisabled(); }
            ImGui::TextUnformatted("Idle --:--:--");
        }
        else {
            if (ImGui::Button(CONCAT("Stop##_mpx_rec_", _this->name), ImVec2(menuWidth, 0))) {
                _this->chain.stopRecording();
            }
            uint64_t seconds = _this->chain.getRecordedSamples() / (uint64_t)MPX_SAMPLERATE;
            ImGui::TextColored(ImVec4(1.0f, 0.0f, 0.0f, 1.0f), "Recording %02d:%02d:%02d", (int)(seconds / 3600), (int)((seconds / 60) % 60), (int)(seconds % 60));
            std::string file = std::filesystem::path(_this->chain.getRecordingPath()).filename().string();
            ImGui::TextWrapped("%s", file.c_str());
        }

        if (!_this->enabled) { style::endDisabled(); }
    }

    static void fftRedraw(ImGui::WaterFall::FFTRedrawArgs args, void* ctx) {
        MPXAnalyzerModule* _this = (MPXAnalyzerModule*)ctx;
        _this->followVFO();
        _this->checkRetune();
        _this->drawWindow();
    }

    static void inputProcess(ImGui::WaterFall::InputHandlerArgs args, void* ctx) {
        MPXAnalyzerModule* _this = (MPXAnalyzerModule*)ctx;
        ImVec2 mouse = ImGui::GetMousePos();
        bool inside = _this->winVisible && mouse.x >= _this->winMin.x && mouse.y >= _this->winMin.y && mouse.x < _this->winMax.x && mouse.y < _this->winMax.y;

        // A click that starts on the window keeps the waterfall locked until released (e.g. while dragging the window)
        if (ImGui::IsMouseClicked(ImGuiMouseButton_Left) || ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
            _this->mouseCaptured = inside;
        }
        else if (!ImGui::IsMouseDown(ImGuiMouseButton_Left) && !ImGui::IsMouseDown(ImGuiMouseButton_Right)) {
            _this->mouseCaptured = false;
        }

        if (inside || _this->mouseCaptured) { gui::waterfall.inputHandled = true; }
    }

    void drawWindow() {
        winVisible = false;
        if (!showWindow) { return; }

        float s = style::uiScale;
        ImGui::SetNextWindowSize(ImVec2(720.0f * s, 340.0f * s), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSizeConstraints(ImVec2(440.0f * s, 280.0f * s), ImVec2(FLT_MAX, FLT_MAX));

        bool open = true;
        std::string title = "MPX Spectrum (" + name + ")###_mpx_win_" + name;
        if (ImGui::Begin(title.c_str(), &open, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse)) {
            drawToolbar();

            // Spectrum on the left, level readouts on the right
            ImVec2 avail = ImGui::GetContentRegionAvail();
            float metersWidth = getMetersWidth();
            drawPlot(ImVec2(avail.x - metersWidth - ImGui::GetStyle().ItemSpacing.x, avail.y));
            ImGui::SameLine();
            drawMeters(metersWidth);
        }
        winMin = ImGui::GetWindowPos();
        winMax = ImVec2(winMin.x + ImGui::GetWindowWidth(), winMin.y + ImGui::GetWindowHeight());
        winVisible = true;
        ImGui::End();

        if (!open) {
            showWindow = false;
            saveSetting("showWindow", showWindow);
        }
    }

    void drawToolbar() {
        if (!enabled) {
            ImGui::TextDisabled("Module disabled");
            return;
        }

        ImGui::Text("%.4f MHz", getFrequency() / 1e6);
        if (!followName.empty() && gui::waterfall.vfos.find(followName) != gui::waterfall.vfos.end()) {
            ImGui::SameLine();
            ImGui::TextDisabled("(following %s)", followName.c_str());
        }

        ImGui::SameLine();
        if (ImGui::Checkbox(CONCAT("Peak hold##_mpx_win_peak_", name), &peakHold)) {
            chain.resetPeak();
            saveSetting("peakHold", peakHold);
        }
        if (peakHold) {
            ImGui::SameLine();
            if (ImGui::SmallButton(CONCAT("Reset##_mpx_win_peak_reset_", name))) { chain.resetPeak(); }
        }

        if (chain.isRecording()) {
            uint64_t seconds = chain.getRecordedSamples() / (uint64_t)MPX_SAMPLERATE;
            ImGui::SameLine();
            ImGui::TextColored(ImVec4(1.0f, 0.2f, 0.2f, 1.0f), "REC %02d:%02d:%02d", (int)(seconds / 3600), (int)((seconds / 60) % 60), (int)(seconds % 60));
        }

        ImGui::SameLine();
        ImGui::TextDisabled("0 dB = %.0f kHz dev.", MPX_DEVIATION / 1000.0);
    }

    float getMetersWidth() {
        ImGui::SetWindowFontScale(MPX_METER_FONT_SCALE);
        float valueWidth = ImGui::CalcTextSize("-00.0 dBr").x;
        ImGui::SetWindowFontScale(1.0f);
        float labelWidth = ImGui::CalcTextSize("1 s: -00.0 dBr").x;
        return std::max<float>(valueWidth, labelWidth) + 8.0f * style::uiScale;
    }

    void drawMeters(float width) {
        MPXMeasurements meas;
        bool valid = enabled && chain.getMeasurements(meas);
        ImVec4 normalCol = ImGui::GetStyleColorVec4(ImGuiCol_Text);
        ImVec4 warnCol = ImVec4(1.0f, 0.65f, 0.0f, 1.0f);
        ImVec4 overCol = ImVec4(1.0f, 0.2f, 0.2f, 1.0f);

        ImGui::BeginGroup();
        ImGui::PushItemWidth(width);

        // MPX power (ITU-R BS.412)
        ImGui::TextDisabled("MPX power");
        ImGui::SetWindowFontScale(MPX_METER_FONT_SCALE);
        if (valid) {
            ImGui::TextColored((meas.powerDBr > 0.0f) ? warnCol : normalCol, "%+.1f dBr", meas.powerDBr);
        }
        else {
            ImGui::TextDisabled("--.- dBr");
        }
        ImGui::SetWindowFontScale(1.0f);
        if (valid) {
            if (meas.powerSeconds < 59.95f) {
                ImGui::TextDisabled("%.0f s avg", meas.powerSeconds);
            }
            else {
                ImGui::TextDisabled("60 s avg");
            }
            ImGui::TextDisabled("1 s: %+.1f dBr", meas.powerShortDBr);
        }
        else {
            ImGui::TextDisabled("60 s avg");
            ImGui::TextDisabled("1 s: --.- dBr");
        }

        ImGui::Spacing();
        ImGui::Spacing();

        // Peak deviation
        ImGui::TextDisabled("Deviation");
        ImGui::SetWindowFontScale(MPX_METER_FONT_SCALE);
        if (valid) {
            ImGui::TextColored((meas.peakDevKHz > MPX_DEVIATION / 1000.0) ? overCol : normalCol, "%.1f kHz", meas.peakDevKHz);
        }
        else {
            ImGui::TextDisabled("--.- kHz");
        }
        ImGui::SetWindowFontScale(1.0f);
        ImGui::TextDisabled("1 s peak");
        if (valid) {
            ImGui::TextColored((meas.maxDevKHz > MPX_DEVIATION / 1000.0) ? overCol : ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled), "Max: %.1f kHz", meas.maxDevKHz);
        }
        else {
            ImGui::TextDisabled("Max: --.- kHz");
        }

        ImGui::Spacing();
        if (ImGui::SmallButton(CONCAT("Reset##_mpx_meas_reset_", name))) { chain.resetMeasurements(); }

        ImGui::PopItemWidth();
        ImGui::EndGroup();
    }

    void drawPlot(ImVec2 size) {
        float s = style::uiScale;
        ImVec2 pos = ImGui::GetCursorScreenPos();
        if (size.x < 50.0f * s || size.y < 50.0f * s) { return; }

        ImGui::InvisibleButton(CONCAT("##_mpx_plot_", name), size);
        bool hovered = ImGui::IsItemHovered();
        ImDrawList* dl = ImGui::GetWindowDrawList();

        // Plot area with room for the axis labels
        float textH = ImGui::GetTextLineHeight();
        ImVec2 pMin(pos.x + ImGui::CalcTextSize("-100").x + 6.0f * s, pos.y + textH * 0.5f);
        ImVec2 pMax(pos.x + size.x - 10.0f * s, pos.y + size.y - textH - 4.0f * s);
        float pw = pMax.x - pMin.x;
        float ph = pMax.y - pMin.y;
        if (pw < 10.0f || ph < 10.0f) { return; }

        ImU32 textCol = ImGui::GetColorU32(ImGuiCol_Text);
        ImU32 gridCol = ImGui::GetColorU32(ImGuiCol_Text, 0.15f);
        ImU32 markerCol = IM_COL32(255, 200, 0, 150);
        ImU32 avgCol = IM_COL32(0, 200, 255, 255);
        ImU32 avgFillCol = IM_COL32(0, 200, 255, 51);
        ImU32 peakCol = IM_COL32(255, 120, 60, 200);

        dl->AddRectFilled(pMin, pMax, ImGui::GetColorU32(ImGuiCol_FrameBg));

        double span = maxFreq * 1000.0;
        float dbRange = std::max<float>(dbMax - dbMin, 1.0f);
        auto fToX = [&](double f) { return pMin.x + (float)(f / span) * pw; };
        auto dbToY = [&](float db) { return std::clamp<float>(pMax.y - ((db - dbMin) / dbRange) * ph, pMin.y, pMax.y); };

        // Level grid
        char buf[64];
        float dbStep = (dbRange > 100.0f) ? 20.0f : 10.0f;
        for (float db = ceilf(dbMin / dbStep) * dbStep; db <= dbMax + 0.01f; db += dbStep) {
            float y = dbToY(db);
            dl->AddLine(ImVec2(pMin.x, y), ImVec2(pMax.x, y), gridCol);
            snprintf(buf, sizeof(buf), "%.0f", db);
            ImVec2 ts = ImGui::CalcTextSize(buf);
            dl->AddText(ImVec2(pMin.x - ts.x - 4.0f * s, y - ts.y * 0.5f), textCol, buf);
        }

        // Frequency grid
        double fStep = (span > 50000.0) ? 10000.0 : 5000.0;
        for (double f = 0.0; f <= span + 1.0; f += fStep) {
            float x = fToX(f);
            dl->AddLine(ImVec2(x, pMin.y), ImVec2(x, pMax.y), gridCol);
            snprintf(buf, sizeof(buf), "%.0fk", f / 1000.0);
            ImVec2 ts = ImGui::CalcTextSize(buf);
            dl->AddText(ImVec2(std::min<float>(x - ts.x * 0.5f, pMax.x - ts.x), pMax.y + 2.0f * s), textCol, buf);
        }

        // Markers for the well known MPX components
        dl->PushClipRect(pMin, pMax, true);
        for (auto const& m : MPX_MARKERS) {
            if (m.freq > span) { continue; }
            float x = fToX(m.freq);
            for (float y = pMin.y; y < pMax.y; y += 6.0f * s) {
                dl->AddLine(ImVec2(x, y), ImVec2(x, std::min<float>(y + 3.0f * s, pMax.y)), markerCol);
            }
            dl->AddText(ImVec2(x + 3.0f * s, pMin.y + 2.0f * s), markerCol, m.label);
        }

        // Traces
        bool haveSpectrum = enabled && chain.getSpectrum(avgSpec, peakSpec);
        if (haveSpectrum) {
            // Shade below the average trace, like the main SDR++ spectrum, then draw the lines on top
            computeTrace(avgSpec, pMin, pw, span, dbToY);
            for (auto const& p : tracePoints) {
                dl->AddLine(ImVec2(p.x, p.y), ImVec2(p.x, pMax.y), avgFillCol);
            }
            if (peakHold) {
                computeTrace(peakSpec, pMin, pw, span, dbToY);
                dl->AddPolyline(tracePoints.data(), (int)tracePoints.size(), peakCol, 0, s);
                computeTrace(avgSpec, pMin, pw, span, dbToY);
            }
            dl->AddPolyline(tracePoints.data(), (int)tracePoints.size(), avgCol, 0, s);
        }
        else {
            const char* msg = enabled ? "No signal" : "Module disabled";
            ImVec2 ts = ImGui::CalcTextSize(msg);
            dl->AddText(ImVec2(pMin.x + (pw - ts.x) * 0.5f, pMin.y + (ph - ts.y) * 0.5f), textCol, msg);
        }

        // Cursor readout
        ImVec2 mouse = ImGui::GetMousePos();
        if (hovered && haveSpectrum && mouse.x >= pMin.x && mouse.x <= pMax.x && mouse.y >= pMin.y && mouse.y <= pMax.y) {
            double binHz = MPX_SAMPLERATE / (double)MPX_FFT_SIZE;
            double f = ((mouse.x - pMin.x) / pw) * span;
            int bin = std::clamp<int>((int)std::round(f / binHz), 0, MPX_FFT_BINS - 1);
            dl->AddLine(ImVec2(mouse.x, pMin.y), ImVec2(mouse.x, pMax.y), ImGui::GetColorU32(ImGuiCol_Text, 0.5f));

            ImGui::BeginTooltip();
            ImGui::Text("%.2f kHz", (bin * binHz) / 1000.0);
            ImGui::Text("Avg:  %.1f dB", toDB(avgSpec[bin]));
            if (peakHold) { ImGui::Text("Peak: %.1f dB", toDB(peakSpec[bin])); }
            ImGui::EndTooltip();
        }
        dl->PopClipRect();

        dl->AddRect(pMin, pMax, gridCol);
    }

    // Fills tracePoints with one point per pixel column
    template <class F>
    void computeTrace(const std::vector<float>& spec, ImVec2 pMin, float pw, double span, F dbToY) {
        double binHz = MPX_SAMPLERATE / (double)MPX_FFT_SIZE;
        int width = (int)pw;
        tracePoints.resize(width + 1);

        // Show the highest bin within each pixel column so narrow carriers (pilot) never disappear
        for (int px = 0; px <= width; px++) {
            int b0 = std::clamp<int>((int)((span * px / width) / binHz), 0, MPX_FFT_BINS - 1);
            int b1 = std::clamp<int>((int)((span * (px + 1) / width) / binHz), b0, MPX_FFT_BINS - 1);
            float p = spec[b0];
            for (int b = b0 + 1; b < b1; b++) { p = std::max<float>(p, spec[b]); }
            tracePoints[px] = ImVec2(pMin.x + px, dbToY(toDB(p)));
        }
    }

    static float toDB(float p) {
        return 10.0f * log10f(std::max<float>(p, 1e-20f));
    }

    std::string name;
    std::string root;
    bool enabled = true;

    // DSP
    VFOManager::VFO* vfo = NULL;
    MPXChain chain;

    // Settings
    std::string followName = "Radio";
    double bandwidth = MPX_DEFAULT_BANDWIDTH;
    bool showWindow = true;
    int avgFrames = 10;
    float maxFreq = 96.0f;
    float dbMin = -100.0f;
    float dbMax = 0.0f;
    bool peakHold = false;
    FolderSelect folderSelect;
    OptionList<std::string, wav::SampleType> sampleTypes;
    int sampleTypeId = 0;

    // GUI state
    OptionList<std::string, std::string> followList;
    int followId = 0;
    std::vector<float> avgSpec;
    std::vector<float> peakSpec;
    std::vector<ImVec2> tracePoints;
    bool winVisible = false;
    ImVec2 winMin;
    ImVec2 winMax;
    bool mouseCaptured = false;
    double measFreq = 0.0;
    double measBandwidth = 0.0;

    EventHandler<ImGui::WaterFall::FFTRedrawArgs> fftRedrawHandler;
    EventHandler<ImGui::WaterFall::InputHandlerArgs> inputHandler;
};

MOD_EXPORT void _INIT_() {
    // Create default recording directory
    std::string root = (std::string)core::args["root"];
    if (!std::filesystem::exists(root + "/recordings")) {
        flog::warn("Recordings directory does not exist, creating it");
        if (!std::filesystem::create_directory(root + "/recordings")) {
            flog::error("Could not create recordings directory");
        }
    }
    json def = json({});
    config.setPath(root + "/mpx_analyzer_config.json");
    config.load(def);
    config.enableAutoSave();
}

MOD_EXPORT ModuleManager::Instance* _CREATE_INSTANCE_(std::string name) {
    return new MPXAnalyzerModule(name);
}

MOD_EXPORT void _DELETE_INSTANCE_(ModuleManager::Instance* instance) {
    delete (MPXAnalyzerModule*)instance;
}

MOD_EXPORT void _END_() {
    config.disableAutoSave();
    config.save();
}
