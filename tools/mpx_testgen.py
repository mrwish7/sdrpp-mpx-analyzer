#!/usr/bin/env python3
"""
MPX test signal generator for calibrating the MPX Analyzer with a HackRF looped back into a receiver.

Builds an FM broadcast test signal (tones, pilot, stereo, RDS, noise) as a seamlessly looping 8-bit IQ
file and optionally plays it with hackrf_transfer. The deviation is exact by construction: levels are
in percent of 75 kHz, and the HackRF's clock error (~20 ppm) is the only scale error.

The FM phase is computed by integrating the MPX in the frequency domain, so the IQ behaves like a real
continuous-phase FM signal. (A sample-by-sample modulator would boost high modulation frequencies by
exactly the amount a phase-difference demodulator loses, hiding that error in a loopback test.)

The signal is placed --offset Hz from the HackRF centre frequency, away from the DC / LO spikes of
both radios. Tune the receiver to centre + offset.

SAFETY: connect the HackRF to the receiver with coax and at least 30-40 dB of fixed attenuation,
keep the HackRF's TX amp OFF (-a 0) and TX VGA low (-x 0..10). An RTL-SDR can be damaged above
about +10 dBm. Nothing is meant to be radiated: use a dummy-load style setup and observe your local
regulations on transmitting.

Examples:
  # Pilot only at 10% (7.5 kHz)
  mpx_testgen.py --pilot 10 --out pilot10.cs8
  # 1 kHz tone at 50%, 57 kHz tone at 5%
  mpx_testgen.py --tone 1000:50 --tone 57000:5 --out tones.cs8
  # Stereo station: L 400 Hz / R 1.5 kHz at 40%, pilot 9%, RDS 4% at +30 deg to the pilot
  mpx_testgen.py --left 400:40 --right 1500:40 --pilot 9 --rds 4 --rds-phase 30 --out station.cs8
  # Bessel null check: 31.18 kHz tone at 100% makes the carrier vanish (modulation index 2.405)
  mpx_testgen.py --tone 31185:100 --out bessel.cs8
  # Generate and transmit on 100.0 MHz (receiver tuned to 100.3 MHz)
  mpx_testgen.py --pilot 10 --out pilot10.cs8 --tx 100.0e6
"""
import argparse
import subprocess
import sys

import numpy as np

DEVIATION = 75000.0     # Hz for 100%
LOOP_SECONDS = 2.0      # All frequencies must be multiples of 1/LOOP_SECONDS for a seamless loop
PILOT_HZ = 19000.0
RDS_BITRATE = 1187.5


def parse_tone(text):
    try:
        freq, level = text.split(":")
        return float(freq), float(level)
    except ValueError:
        raise argparse.ArgumentTypeError("expected FREQ:PERCENT, e.g. 1000:50")


def check_periodic(name, freq):
    cycles = freq * LOOP_SECONDS
    if abs(cycles - round(cycles)) > 1e-9:
        sys.exit(f"{name} {freq} Hz doesn't loop seamlessly: use a multiple of {1 / LOOP_SECONDS} Hz")


def rds_baseband(n, fs, rng):
    """Standard EN 50067 RDS data signal: differentially coded biphase symbols with cos shaping.
    Built where a bit is a whole number of samples, then resampled exactly (circularly) to fs.
    Normalised to a peak of 1."""
    bits = int(round(RDS_BITRATE * LOOP_SECONDS))      # 2375 bits in 2 s
    spb = 320                                          # samples per bit at the build rate
    fs_build = RDS_BITRATE * spb                       # 380 kHz
    n_build = bits * spb
    data = rng.integers(0, 2, bits)
    sym = 2 * (np.cumsum(data) % 2) - 1                # differential coding
    x = np.zeros(n_build)
    x[np.arange(bits) * spb] += sym                    # biphase: +impulse then -impulse half a bit later
    x[np.arange(bits) * spb + spb // 2] -= sym
    td = 1.0 / RDS_BITRATE
    spec = np.fft.rfft(x)
    f = np.fft.rfftfreq(n_build, 1.0 / fs_build)
    spec *= np.where(f <= 2.0 / td, np.cos(np.pi * f * td / 4.0), 0.0)
    # Exact periodic resampling to fs by zero padding the spectrum
    out_spec = np.zeros(n // 2 + 1, dtype=complex)
    out_spec[:len(spec)] = spec
    a = np.fft.irfft(out_spec, n) * (n / n_build)
    return a / np.max(np.abs(a))


def build_mpx(args, fs, n):
    t = np.arange(n) / fs
    mpx = np.zeros(n)

    for freq, level in args.tone:
        check_periodic("Tone", freq)
        mpx += (level / 100.0) * np.cos(2 * np.pi * freq * t)

    pilot_phase = 2 * np.pi * PILOT_HZ * t
    if args.left or args.right:
        left = np.zeros(n)
        right = np.zeros(n)
        for freq, level in args.left:
            check_periodic("Left tone", freq)
            left += (level / 100.0) * np.cos(2 * np.pi * freq * t)
        for freq, level in args.right:
            check_periodic("Right tone", freq)
            right += (level / 100.0) * np.cos(2 * np.pi * freq * t)
        # Levels are the deviation each channel alone would cause (L+R and L-R share it)
        mpx += (left + right) / 2 + (left - right) / 2 * np.cos(2 * pilot_phase)

    if args.pilot:
        mpx += (args.pilot / 100.0) * np.cos(pilot_phase)

    if args.rds:
        check_periodic("RDS carrier", 3 * PILOT_HZ + args.rds_offset)
        a = rds_baseband(n, fs, np.random.default_rng(args.seed))
        carrier = 3 * pilot_phase + 2 * np.pi * args.rds_offset * t + np.deg2rad(args.rds_phase)
        mpx += (args.rds / 100.0) * a * np.cos(carrier)

    return mpx


def fm_modulate(mpx, fs, offset):
    """Continuous-phase FM: integrate 2*pi*75kHz*mpx in the frequency domain (exact for a periodic signal)."""
    n = len(mpx)
    spec = np.fft.rfft(mpx)
    f = np.fft.rfftfreq(n, 1.0 / fs)
    spec[0] = 0.0                                      # DC would be a frequency offset, not modulation
    phase_spec = np.zeros_like(spec)
    phase_spec[1:] = spec[1:] * (2 * np.pi * DEVIATION) / (1j * 2 * np.pi * f[1:])
    phase = np.fft.irfft(phase_spec, n)
    t = np.arange(n) / fs
    return np.exp(1j * (phase + 2 * np.pi * offset * t))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tone", type=parse_tone, action="append", default=[], metavar="FREQ:PERCENT",
                   help="MPX tone, e.g. 1000:50 (repeatable)")
    p.add_argument("--left", type=parse_tone, action="append", default=[], metavar="FREQ:PERCENT",
                   help="Left channel tone (stereo, use with --pilot)")
    p.add_argument("--right", type=parse_tone, action="append", default=[], metavar="FREQ:PERCENT",
                   help="Right channel tone (stereo, use with --pilot)")
    p.add_argument("--pilot", type=float, default=0.0, metavar="PERCENT", help="19 kHz pilot level")
    p.add_argument("--rds", type=float, default=0.0, metavar="PERCENT", help="RDS level (peak of the EN 50067 waveform)")
    p.add_argument("--rds-phase", type=float, default=0.0, metavar="DEG", help="RDS carrier phase vs the pilot's 3rd harmonic")
    p.add_argument("--rds-offset", type=float, default=0.0, metavar="HZ", help="RDS carrier frequency offset (unlocked RDS)")
    p.add_argument("--snr", type=float, default=None, metavar="DB", help="Add noise at this carrier-to-noise ratio (over the sample rate)")
    p.add_argument("--rate", type=float, default=2.4e6, help="IQ sample rate (default 2.4 MS/s)")
    p.add_argument("--offset", type=float, default=300e3, help="Signal offset from the HackRF centre (default 300 kHz)")
    p.add_argument("--seed", type=int, default=1, help="Random seed for RDS data and noise")
    p.add_argument("--out", required=True, help="Output file (8-bit signed IQ, hackrf_transfer format)")
    p.add_argument("--tx", type=float, default=None, metavar="HZ", help="Transmit with hackrf_transfer at this centre frequency")
    p.add_argument("--txvga", type=int, default=0, help="HackRF TX VGA gain in dB (default 0, keep it low)")
    args = p.parse_args()

    fs = args.rate
    n = int(round(fs * LOOP_SECONDS))
    check_periodic("Offset", args.offset)

    mpx = build_mpx(args, fs, n)
    peak = np.max(np.abs(mpx)) * 100.0
    rms_power = np.mean(mpx ** 2)
    ref = (19000.0 / DEVIATION) ** 2 / 2
    print(f"MPX peak deviation {peak:.2f}% ({peak * DEVIATION / 100e3:.2f} kHz), "
          f"MPX power {10 * np.log10(max(rms_power, 1e-20) / ref):+.2f} dBr")

    iq = fm_modulate(mpx, fs, args.offset)
    if args.snr is not None:
        rng = np.random.default_rng(args.seed + 1)
        sigma = np.sqrt(10 ** (-args.snr / 10) / 2)
        iq = iq + sigma * (rng.standard_normal(n) + 1j * rng.standard_normal(n))

    # 8-bit interleaved I/Q with headroom
    scale = 100.0 / max(1.0, np.max(np.abs(np.concatenate([iq.real, iq.imag]))))
    out = np.empty(2 * n, dtype=np.int8)
    out[0::2] = np.clip(np.round(iq.real * scale), -127, 127)
    out[1::2] = np.clip(np.round(iq.imag * scale), -127, 127)
    out.tofile(args.out)
    print(f"Wrote {args.out}: {LOOP_SECONDS:.0f} s loop at {fs / 1e6:g} MS/s, signal at centre {args.offset / 1e3:+g} kHz")

    if args.tx is not None:
        cmd = ["hackrf_transfer", "-t", args.out, "-f", str(int(args.tx)), "-s", str(int(fs)),
               "-x", str(args.txvga), "-a", "0", "-R"]
        print(f"Transmitting: {' '.join(cmd)}   (tune the receiver to {(args.tx + args.offset) / 1e6:.4f} MHz, Ctrl+C to stop)")
        subprocess.run(cmd)


if __name__ == "__main__":
    main()
