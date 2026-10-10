# MPX Analyzer for SDR++

An [SDR++](https://github.com/AlexandreRouma/SDRPlusPlus) plugin for broadcast FM that:

- **Shows the MPX (composite) spectrum** from 0 to 96 kHz in its own window, so you can see the 19 kHz stereo pilot, the 38 kHz L−R subcarrier, 57 kHz RDS and any SCA subcarriers.
- **Records the raw MPX signal** as a mono, 192 kHz WAV file. This is the demodulated composite signal, not de-multiplexed audio.

## How it works

The plugin creates its own VFO (384 kS/s IQ, 250 kHz IF bandwidth by default) and FM-demodulates it with SDR++'s quadrature demodulator. It then decimates the result by 2 to 192 kHz. The decimation filter also corrects the demodulator's high-frequency roll-off (it measures the phase change over one 384 kHz sample, which would otherwise read 0.4% low at 19 kHz and 3.6% low at 57 kHz), so the MPX is flat within 0.01% up to about 76 kHz. That 192 kHz signal feeds two places:

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
- **Pilot (kHz):** the level of the 19 kHz stereo pilot as deviation, averaged over 1 second. A typical 9% pilot reads 6.75 kHz; hover over it to see it as a percentage.
  - **How:** the MPX is correlated with a 19 kHz reference over 10 ms windows, so stereo audio and noise hardly affect it.
  - **No pilot:** shows "none" below 0.5 kHz, e.g. on mono stations.
- **RDS (kHz):** the RDS injection level as peak deviation, averaged over 1 second. Typical values are 2–4 kHz; hover over it to see it as a percentage.
  - **How:** the 57 kHz band is shifted to 0 Hz and filtered to the RDS bandwidth. RDS has almost no energy exactly at 57 kHz, while noise does, so the noise in the band is estimated and removed. The remaining RMS level is converted to peak deviation with a crest factor of 1.51. That's the value real encoders show against reference instruments: a calibrated transmitter measured with MPX Tool, and a Pira P175-calibrated receiver across 21 stations, agree within 0.5%. The theoretical EN 50067 waveform would give 1.444, about 4.5% lower.
  - **Result:** readings stay accurate on noisy signals, matching a peak reading on a clean one.
  - **No RDS:** shows "none" when there's no RDS, i.e. when the band holds less than twice as much RDS power as noise, or the level is below 0.3 kHz.
- **RDS lock:** whether the RDS 57 kHz subcarrier is phase-locked to the third harmonic of the 19 kHz pilot, as EN 50067 / IEC 62106 require for stereo stations. The standard allows in phase (0°) or quadrature (90°), within ±10°.
  - **Readings:** *Yes* with the phase offset (e.g. "Yes +2°"), green when within tolerance and amber when not. *No* (amber) means the RDS encoder isn't locked to the pilot. *n/a* means there is no pilot or no RDS.
  - **How:** squaring the RDS signal removes its data and leaves the carrier phase. This is compared with 3× the pilot phase every 10 ms and averaged over 5 seconds.
  - **Lock decision:** a stable difference means locked; a drifting one averages out. Even a 0.2 Hz frequency difference is detected as "No".
  - **Hover:** the tooltip shows the exact phase and the stability figure (locked above 0.90).
  - **180° ambiguity:** RDS's data modulation hides 180° phase steps, so the offset is shown in the range −90° to +90°.
- **CPU use:** the spectrum and all measurements are only computed while the MPX window is open. Recording and audio output carry on when it's closed. Reopening the window restarts the measurements, so the 60 s MPX power and 5 s RDS lock windows fill up again. The analysis is light anyway, about 0.3% of one desktop CPU core. Most of the plugin's CPU use is the channel filter that extracts its 384 kHz channel from the SDR's bandwidth, which recording and audio output also need.
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

**Pre-built (GitHub Actions):** every push builds `mpx_analyzer.dylib` for Intel (`mpx_analyzer_macos_intel`) and Apple Silicon (`mpx_analyzer_macos_arm`) Macs. The build uses the same runners, dependencies and macOS 10.15 deployment target as the official SDR++ macOS builds.
- **Prepared for the app:** the dylib's library references point into `SDR++.app/Contents/Frameworks`, like the app's own plugins, so it doesn't need Homebrew on the user's Mac.
- **Checked against the nightly:** the build checks that every library the plugin needs is in the current SDR++ macOS nightly, and shows a warning in the build log if not.
- **Version match:** as on Windows, `SDRPP_VERSION.txt` shows which SDR++ commit it was built against, which must match the SDR++ you run.

To install a downloaded build:
```sh
# Downloaded files are quarantined by macOS and would be refused when SDR++ loads them
xattr -d com.apple.quarantine mpx_analyzer.dylib
cp mpx_analyzer.dylib /Applications/SDR++.app/Contents/Plugins/
```
Then start SDR++ and add an `mpx_analyzer` instance in Module Manager. If macOS complains that SDR++ is damaged after the plugin is added, re-sign the app locally: `codesign --force --deep -s - /Applications/SDR++.app`.

**Building locally:**
```sh
brew install cmake pkg-config fftw volk glfw zstd
cmake -B build -DSDRPP_SOURCE_DIR="$HOME/SDRPlusPlus" -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
cp build/mpx_analyzer.dylib /Applications/SDR++.app/Contents/Plugins/
```
A local build links to your Homebrew libraries, so it works on your own Mac but isn't suitable for giving to others. Use the GitHub Actions build for that.

## Usage

1. Start SDR++. Open **Module Manager**, choose `mpx_analyzer` in the type list, give it a name (e.g. `MPX`) and press **+**.
2. Tune the Radio VFO to an FM broadcast station. The plugin's orange VFO follows it, and the **MPX Spectrum** window appears. It stays open when the menu is collapsed; reopen it with *Show MPX window* after closing it.
3. Hover over the spectrum to read frequency and level. Enable **Peak hold** to catch transients.
4. Choose a folder and sample type in the plugin's menu section, then press **Record MPX**. You can also use the **red record button** at the top right of the MPX window; it turns into a square stop button, with the elapsed time beside it, while recording.

### Sending the MPX to a sound device (Windows, Linux, macOS)

Tick **Audio output** in the plugin menu; it's off by default. The plugin then adds an audio stream named after the plugin instance (e.g. `MPX`) to SDR++'s **Sinks** menu, separate from the Radio's decoded audio.
1. **Pick the output:** in the Sinks menu, choose the `MPX` stream's sink (*Audio*) and its output device.
2. **Pick the rate:** choose **192000 Hz**, so the MPX passes through unchanged.
   - **Lower rates:** if the device only offers lower rates, the MPX is resampled to the chosen rate and the plugin menu shows a warning, because content above half that rate is lost. At 48 kHz, for example, nothing above 24 kHz survives, so the RDS is lost and so are most of the stereo subcarriers.
3. **Level:** full scale (±1.0) is 75 kHz deviation, the same as the WAV recordings. The stream's volume slider in the Sinks menu can reduce it if over-deviating stations clip.
4. **Other sinks:** the same stream can also go to the *Network* sink, sending MPX at 192 kHz over TCP or UDP.

SDR++ remembers the device and rate selection for the stream. The option isn't available on Android.

### Menu options

| Option | Description |
|---|---|
| Follow VFO | VFO whose frequency and bandwidth the plugin tracks (`None` = set them yourself) |
| IF Bandwidth | Channel filter before demodulation (50–350 kHz). Copied from the followed VFO and read-only while following; adjustable when Follow VFO is `None`. |
| Averaging | Number of FFT frames averaged (about 94 frames/s) |
| Max Freq | Upper edge of the display (20–96 kHz) |
| dB Max / dB Min | Vertical scale |
| Peak hold | Overlay the maximum level seen since the last reset |
| Level units | Show deviation, pilot and RDS levels in kHz (default) or percent, where 100% = 75 kHz. Hover over Pilot/RDS to see the other unit |
| Sample type | Int16 (default) / Int32 / Float32 WAV samples |
| Audio output | Send the MPX to an SDR++ sink (sound device / network), configured in the Sinks menu. Off by default; not on Android |

## Calibration test bench (HackRF)

`tools/mpx_testgen.py` generates FM broadcast test signals with exactly known levels, for checking the plugin end to end with a HackRF looped back into a receiver (e.g. an RTL-SDR). It writes a seamlessly looping 8-bit IQ file and can play it with `hackrf_transfer`. It needs Python 3, numpy and `hackrf-tools`.

**Safety first:**
- **Coax only:** connect the HackRF to the receiver by coax through **30–40 dB of fixed attenuation**. An RTL-SDR can be damaged above about +10 dBm.
- **Low output:** keep the HackRF's TX amplifier off (the script uses `-a 0`) and its TX gain low (`--txvga 0..10`).
- **Don't radiate:** nothing should be transmitted over the air. Follow your local regulations.

The deviation of the generated signal is exact by construction. The only error is the HackRF's clock, about 0.002%. The FM is generated with continuous phase, like a real transmitter, so a loopback test is realistic.

```sh
# Pilot only at 10%, on 100.0 MHz (receiver/plugin tuned to 100.3 MHz)
python3 tools/mpx_testgen.py --pilot 10 --out pilot10.cs8 --tx 100.0e6
# Tones: 1 kHz at 50% plus 57 kHz at 5%
python3 tools/mpx_testgen.py --tone 1000:50 --tone 57000:5 --out tones.cs8 --tx 100.0e6
# Stereo station: L 400 Hz / R 1.5 kHz at 40%, pilot 9%, RDS 4% at +30 deg to the pilot
python3 tools/mpx_testgen.py --left 400:40 --right 1500:40 --pilot 9 --rds 4 --rds-phase 30 --out station.cs8 --tx 100.0e6
# Unlocked RDS (0.5 Hz off), and a Bessel null check (carrier vanishes at 31.185 kHz, 100%)
python3 tools/mpx_testgen.py --pilot 9 --rds 4 --rds-offset 0.5 --out unlocked.cs8 --tx 100.0e6
python3 tools/mpx_testgen.py --tone 31185:100 --out bessel.cs8 --tx 100.0e6
```

**What to expect:**
- **Pilot, tones and MPX power:** within about 0.5%. The generator prints the expected peak deviation and MPX power for each signal.
- **RDS:** about 4.6% above the `--rds` setting. The generator uses the theoretical EN 50067 waveform, while the plugin's RDS reading follows real encoders and reference instruments (see Measurements).
- **Deviation:** reads a little high, because the receiver's noise and the 8-bit IQ add to the true sample peak.

**SDR++ VFO ripple correction:**
- **The problem:** this bench found that SDR++'s own channel resampler has up to about 0.45 dB of passband ripple, depending on the SDR sample rate. On a low-deviation signal that alone read a 10% pilot as 9.55%.
- **The fix:** the plugin measures this response by passing a comb of tones through a private copy of SDR++'s resampler, and cancels it with an equalising filter before the demodulator. This happens whenever the sample rate or decimation changes.
- **In the menu:** the correction applied is shown as "VFO ripple corrected: x.xx dB".

## License

GPL-3.0, the same as SDR++, since the plugin links against SDR++'s core.
