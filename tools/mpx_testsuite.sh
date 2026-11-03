#!/bin/bash
# Runs the MPX Analyzer calibration suite with a HackRF looped back (through attenuators!) into a receiver.
# Each test transmits until you press Enter, so you can read the plugin's values.
#
# Usage: tools/mpx_testsuite.sh [centre frequency Hz] [TX VGA dB]
#   The test signal is 300 kHz above the centre: with the default 100.0 MHz, tune the plugin to 100.3 MHz.
#
# SAFETY: coax only, 30-40 dB of fixed attenuation between HackRF and receiver, HackRF amp off (-a 0).

FREQ=${1:-100000000}
TXVGA=${2:-0}
DIR=$(cd "$(dirname "$0")" && pwd)
WORK=$(mktemp -d)
trap 'kill $TXPID 2>/dev/null; rm -rf "$WORK"' EXIT

# A typo like 40 instead of 4 shouldn't reach a receiver connected without an attenuator
if [ "$TXVGA" -gt 20 ] 2> /dev/null && [ "$FORCE" != "1" ]; then
    echo "TX VGA $TXVGA dB is high for a direct loopback (0-10 is normally plenty). Run with FORCE=1 to allow it."
    exit 1
fi

command -v hackrf_transfer > /dev/null || { echo "hackrf_transfer not found (install hackrf / hackrf-tools)"; exit 1; }
python3 -c "import numpy" 2> /dev/null || { echo "numpy not found (install python3-numpy)"; exit 1; }

# name | expected readings | generator arguments
TESTS=(
"Pilot 10%|Pilot 10.00%, MPX power -8.07 dBr, Deviation ~10%|--pilot 10"
"Tone 1 kHz 10%|MPX power (1 s) -8.07 dBr, Deviation ~10%|--tone 1000:10"
"Tone 5 kHz 10%|MPX power (1 s) -8.07 dBr, Deviation ~10%|--tone 5000:10"
"Tone 10 kHz 10%|MPX power (1 s) -8.07 dBr, Deviation ~10%|--tone 10000:10"
"Tone 15 kHz 10%|MPX power (1 s) -8.07 dBr, Deviation ~10%|--tone 15000:10"
"Tone 19 kHz 10%|MPX power (1 s) -8.07 dBr, Pilot 10.00%, Deviation ~10%|--tone 19000:10"
"Tone 30 kHz 10%|MPX power (1 s) -8.07 dBr, Deviation ~10%|--tone 30000:10"
"Tone 45 kHz 10%|MPX power (1 s) -8.07 dBr, Deviation ~10%|--tone 45000:10"
"Tone 57 kHz 10%|MPX power (1 s) -8.07 dBr, Deviation ~10%|--tone 57000:10"
"Tone 76 kHz 10%|MPX power (1 s) -8.07 dBr, Deviation ~10%|--tone 76000:10"
"Tone 19 kHz 50%|MPX power (1 s) +5.91 dBr, Pilot 50.0%, Deviation 50%|--tone 19000:50"
"Tone 1 kHz 100% (75 kHz)|MPX power (1 s) +11.93 dBr, Deviation 100% = 75 kHz (note how much it over-reads)|--tone 1000:100"
"Bessel null, 31.185 kHz at 100%|In SDR++'s main spectrum, zoomed in on the carrier: the carrier line should vanish|--tone 31185:100"
"RDS 5% only|RDS ~5.23% (theoretical waveform, see README)|--rds 5"
"Stereo station + RDS at +30 deg|Pilot 9.00%, RDS ~4.18%, RDS lock Yes +30 deg|--left 400:40 --right 1500:40 --pilot 9 --rds 4 --rds-phase 30"
"Unlocked RDS (0.5 Hz off)|Pilot 9.00%, RDS ~4.18%, RDS lock No|--pilot 9 --rds 4 --rds-offset 0.5"
"Pilot 10% with noise (30 dB)|Pilot ~10.0%, Deviation reads high from the noise|--pilot 10 --snr 30"
)

echo "Generating ${#TESTS[@]} test signals..."
i=0
for T in "${TESTS[@]}"; do
    IFS='|' read -r NAME EXPECT ARGS <<< "$T"
    python3 "$DIR/mpx_testgen.py" $ARGS --out "$WORK/test$i.cs8" > "$WORK/test$i.txt" || exit 1
    i=$((i + 1))
done

echo
echo "HackRF centre $(awk "BEGIN { printf \"%.4f\", $FREQ / 1e6 }") MHz, TX VGA $TXVGA dB, amp off."
echo "Tune the plugin (or the Radio VFO it follows) to $(awk "BEGIN { printf \"%.4f\", ($FREQ + 300000) / 1e6 }") MHz."
echo "Allow ~5 s per test for the readings to settle (RDS lock takes 5 s; MPX power uses the 1 s value)."

i=0
for T in "${TESTS[@]}"; do
    IFS='|' read -r NAME EXPECT ARGS <<< "$T"
    echo
    echo "=== Test $((i + 1))/${#TESTS[@]}: $NAME"
    echo "    Generator: $(head -n 1 "$WORK/test$i.txt")"
    echo "    Expect:    $EXPECT"
    hackrf_transfer -t "$WORK/test$i.cs8" -f "$FREQ" -s 2400000 -x "$TXVGA" -a 0 -R > /dev/null 2>&1 &
    TXPID=$!
    read -r -p "    Transmitting. Note the readings, then press Enter for the next test (Ctrl+C to quit)... "
    kill $TXPID 2> /dev/null
    wait $TXPID 2> /dev/null
    sleep 1    # let the HackRF release before the next test
    i=$((i + 1))
done
echo
echo "Suite complete."
