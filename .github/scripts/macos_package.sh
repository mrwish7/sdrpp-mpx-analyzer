#!/bin/bash
# Prepares mpx_analyzer.dylib for use inside the official SDR++.app bundle, the same way SDR++'s
# macos/bundle_utils.sh treats its own plugins: libraries are found through @rpath in the bundle's
# Contents/Frameworks folder instead of at the paths they had on the build machine.
#
# Usage: macos_package.sh <built dylib> <output dir> [SDR++ nightly zip to check against]
set -e

SRC=$1
OUT=$2
NIGHTLY_ZIP=$3

mkdir -p "$OUT"
DYLIB="$OUT/$(basename "$SRC")"
cp "$SRC" "$DYLIB"

echo "Dependencies as built:"
otool -L "$DYLIB"

# Point every non-system library at @rpath/<name>
otool -L "$DYLIB" | tail -n +2 | awk '{print $1}' | while read -r DEP; do
    case "$DEP" in
        /usr/lib/*|/System/*) continue ;;
    esac
    NAME=$(basename "$DEP")
    if [ "$NAME" = "$(basename "$DYLIB")" ]; then continue; fi
    install_name_tool -change "$DEP" "@rpath/$NAME" "$DYLIB"
done

# Replace the build machine's rpaths with the bundle's Frameworks folder (plugins live in Contents/Plugins)
otool -l "$DYLIB" | grep -A2 LC_RPATH | grep " path " | awk '{print $2}' | while read -r RPATH; do
    install_name_tool -delete_rpath "$RPATH" "$DYLIB"
done
install_name_tool -add_rpath @loader_path/../Frameworks "$DYLIB"

# install_name_tool invalidates the signature; Apple Silicon requires at least an ad-hoc one
codesign --force -s - "$DYLIB"

echo "Dependencies after packaging:"
otool -L "$DYLIB"

# Optionally check every @rpath library exists in the official SDR++ bundle
if [ -n "$NIGHTLY_ZIP" ] && [ -f "$NIGHTLY_ZIP" ]; then
    CHECK_DIR=$(mktemp -d)
    unzip -q "$NIGHTLY_ZIP" -d "$CHECK_DIR"
    FRAMEWORKS=$(find "$CHECK_DIR" -type d -path "*SDR++.app/Contents/Frameworks" | head -n 1)
    MISSING=0
    for NAME in $(otool -L "$DYLIB" | tail -n +2 | awk '{print $1}' | grep "^@rpath/" | sed 's#@rpath/##'); do
        if [ -f "$FRAMEWORKS/$NAME" ]; then
            echo "OK: $NAME is in the SDR++ nightly bundle"
        else
            echo "::warning::$NAME is not in the SDR++ nightly bundle ($(ls "$FRAMEWORKS" | grep -i "${NAME%%.*}" | tr '\n' ' ')). The plugin may fail to load."
            MISSING=1
        fi
    done
    rm -rf "$CHECK_DIR"
    if [ $MISSING -eq 0 ]; then echo "All plugin libraries are present in the SDR++ nightly bundle"; fi
fi
