# MPX Analyzer for SDR++

An [SDR++](https://github.com/AlexandreRouma/SDRPlusPlus) plugin for broadcast FM that:

- **Shows the MPX (composite) spectrum** from 0 to 96 kHz in its own window, so you can see the 19 kHz stereo pilot, the 38 kHz L−R subcarrier, 57 kHz RDS and any SCA subcarriers.
- **Records the raw MPX signal** as a mono, 192 kHz WAV file. This is the demodulated composite signal, not de-multiplexed audio.

## How it works

The plugin creates its own VFO (384 kS/s IQ, 250 kHz IF bandwidth by default) and FM-demodulates it with SDR++'s quadrature demodulator. It then decimates the result by 2 to 192 kHz. That 192 kHz signal feeds two places:

- the spectrum display (4096-point FFT, Blackman-Harris window, 50% overlap)
- the WAV recorder

By default the plugin VFO **follows the `Radio` VFO**, copying both its frequency and its filter bandwidth. You tune and set the bandwidth in the main window as normal, and the MPX view tracks the station you're listening to. Pick a different VFO, or `None` to tune independently and set the IF bandwidth in the plugin menu, with *Follow VFO*.

The Radio's WFM default bandwidth of 150 kHz is narrow for MPX work. It cuts into the FM sidebands of the higher subcarriers, which raises noise and distortion around 38–60 kHz. For a cleaner MPX view, widen the Radio bandwidth to about 200–250 kHz on a station without strong neighbours.

**Levels:** the demodulator output is scaled so that ±1.0 = ±75 kHz deviation. On the spectrum, **0 dB is a sine at 75 kHz deviation**. A 9% (6.75 kHz) pilot therefore reads about −21 dB, and RDS at 2–4 kHz injection reads roughly −30 to −40 dB, spread over its ~5 kHz bandwidth.

**Measurements** (shown to the right of the spectrum):
- **MPX power (dBr):** the power of the complete MPX signal (audio, pilot, RDS and everything else), as defined in ITU-R BS.412.
  - **Reference:** 0 dBr is the power of a sine tone giving ±19 kHz deviation.
  - **Main reading:** the RMS over a sliding 60-second window, as BS.412 specifies. For the first minute after tuning it shows the average so far, with the time covered (e.g. "12 s avg").
  - **Short-term reading:** a 1-second value is shown underneath for spotting changes quickly.
  - **Highlight:** the reading turns amber above 0 dBr, the most common limit. Some regulators allow +3 or +6 dBr.
- **Deviation (kHz):** the peak frequency deviation over the last second, plus a max hold. Both turn red above ±75 kHz.
  - **What "peak" means here:** this is the true sample peak of the demodulated signal, so noise on weak signals pushes it up. Professional modulation meters may apply extra filtering or a peak-counting rule, so expect small differences from them.
- **Reset:** the measurements restart automatically when the frequency or bandwidth changes. *Reset* restarts them by hand.
- **Accuracy:** a slow DC filter removes any carrier offset first, so slight mistuning doesn't affect the readings. The IF filter does affect them: a narrow bandwidth (e.g. the Radio's 150 kHz WFM default) distorts the MPX and changes both readings. For measurements, use a bandwidth of about 200–250 kHz on a clean signal.

**WAV files:**
- **Sample format:** Int16 (default) is enough for most uses and keeps files small. Int32 gives more dynamic range. Both clip at ±75 kHz deviation; choose Float32 if you need to capture over-deviating stations without clipping.
- **Bandwidth:** 192 kHz sampling gives a 96 kHz Nyquist limit. Content up to about 90 kHz is clean.
- **File size:** recordings are split into a new file before reaching the 4 GiB WAV limit, which is about 3 hours at Int16 (about 93 minutes at Int32/Float32). Int16 uses 384 kB/s, about 1.4 GB per hour.
- **File names:** `mpx_<frequency>Hz_<YYYYMMDD>_<HHMMSS>.wav`, in the folder you choose (default `<SDR++ root>/recordings`).

## Requirements

- An SDR source running at **≥ 384 kS/s**. The menu shows a warning otherwise.
- The **SDR++ source tree for the exact version you run**. SDR++ has no stable plugin ABI, so a plugin built against one version may crash another. Distribution `sdrpp` packages are known to be incompatible with out-of-tree modules; build SDR++ from source or use an official nightly matching your source checkout.

## Building

### Linux (standalone)

```sh
# Dependencies of SDR++'s core (Debian/Ubuntu names)
sudo apt install build-essential cmake pkg-config libfftw3-dev libvolk-dev libglfw3-dev libzstd-dev libgl-dev

git clone https://github.com/AlexandreRouma/SDRPlusPlus.git "$HOME/SDRPlusPlus"   # same commit as the SDR++ you run

cd sdrpp-mpx-analyzer
cmake -B build -DSDRPP_SOURCE_DIR="$HOME/SDRPlusPlus"
cmake --build build -j
sudo cmake --install build      # installs mpx_analyzer.so to /usr/lib/sdrpp/plugins
```

This also compiles `sdrpp_core`, but only to link against; it is not installed. If SDR++ was installed to another prefix, pass the same `-DCMAKE_INSTALL_PREFIX` (default `/usr`). You can also copy `build/mpx_analyzer.so` by hand into the `modulesDirectory` set in SDR++'s `config.json`.

### As part of an SDR++ build (any platform)

Add one line near the end of SDR++'s top-level `CMakeLists.txt`, after the other modules:

```cmake
add_subdirectory("/path/to/sdrpp-mpx-analyzer" "${CMAKE_BINARY_DIR}/mpx_analyzer")
```

Then build SDR++ as usual. The plugin is built and installed with the other modules.

### Windows

**Pre-built (GitHub Actions):** every push builds `mpx_analyzer.dll` on GitHub's Windows runners, using the same toolchain and dependencies as the official SDR++ Windows builds. Download it from the run's *Artifacts* (`mpx_analyzer_windows_x64`). `SDRPP_VERSION.txt` in the artifact shows which SDR++ commit it was built against. That must match the SDR++ you run, so use an SDR++ nightly from around the same date, or run the workflow manually (*Actions → Build → Run workflow*) with `sdrpp_ref` set to your SDR++ build's commit.

**Building locally:** Set up the build environment described in SDR++'s readme (Visual Studio, vcpkg, PothosSDR). Then use either method above, adding `-DCMAKE_TOOLCHAIN_FILE=<vcpkg>/scripts/buildsystems/vcpkg.cmake` for the standalone build. Copy `mpx_analyzer.dll` into SDR++'s `modules` folder.

### macOS

Use the standalone build and copy `mpx_analyzer.dylib` into the app bundle's `Plugins` folder, or into the configured `modulesDirectory`.

## Usage

1. Start SDR++. Open **Module Manager**, choose `mpx_analyzer` in the type list, give it a name (e.g. `MPX`) and press **+**.
2. Tune the Radio VFO to an FM broadcast station. The plugin's orange VFO follows it, and the **MPX Spectrum** window appears. It stays open when the menu is collapsed; reopen it with *Show MPX window* after closing it.
3. Hover over the spectrum to read frequency and level. Enable **Peak hold** to catch transients.
4. Choose a folder and sample type, then press **Record MPX**.

### Menu options

| Option | Description |
|---|---|
| Follow VFO | VFO whose frequency and bandwidth the plugin tracks (`None` = set them yourself) |
| IF Bandwidth | Channel filter before demodulation (50–350 kHz). Copied from the followed VFO and read-only while following; adjustable when Follow VFO is `None`. |
| Averaging | Number of FFT frames averaged (about 94 frames/s) |
| Max Freq | Upper edge of the display (20–96 kHz) |
| dB Max / dB Min | Vertical scale |
| Peak hold | Overlay the maximum level seen since the last reset |
| Sample type | Int16 (default) / Int32 / Float32 WAV samples |

## License

GPL-3.0, the same as SDR++, since the plugin links against SDR++'s core.
