#!/usr/bin/env bash
# Stage the libevogen artifact content (patched libraries + headers) for
# tools/build_wrapper.jl and tools/install_artifact.jl.
#
# libevogen.so and libplinkio.so are built without a SONAME and without an
# rpath, so out of the box they only find their sibling libraries (the
# Intel oneAPI runtimes, each other) via LD_LIBRARY_PATH or oneAPI env
# vars, and a full-path dlopen does not register their names for NEEDED
# lookup, so anything that links them (like the libcxxwrap wrapper) would
# silently load a second copy of the library. This script copies the
# pristine binaries into artifact_staging/ and patches them:
#   * RUNPATH = $ORIGIN, so they resolve their siblings wherever they are
#     installed (no environment setup needed), and
#   * SONAME = the file name, so that loading them by full path (as the
#     JLL does) registers them for NEEDED lookup by name.
#
# Usage:
#   tools/rebuild_libs.sh [source_root]
#
#   source_root  directory holding the pristine binaries and headers
#                (default: the package root itself; it must contain
#                lib/ and include/)
#
# Requirements: bash, plus either patchelf or python3 (the patch is done
# with patchelf if available, otherwise with the stdlib-only python3 ELF
# patcher embedded below).
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
pkgroot="$(cd "$here/.." && pwd)"
default_srcroot="$pkgroot"
srcroot="${1:-$default_srcroot}"
srclib="$srcroot/lib"
srcinc="$srcroot/include"
stage="$pkgroot/artifact_staging"

# Libraries that must be present, in the order used by the JLL wrapper
# (dependencies first).
libs=(libintlc.so.5 libimf.so libsvml.so libirng.so libiomp5.so libplinkio.so libevogen.so)

if [[ -z "$srcroot" ]]; then
    echo "error: no source root given" >&2
    echo "usage: tools/rebuild_libs.sh [source_root]" >&2
    exit 1
fi
for lib in "${libs[@]}"; do
    [[ -f "$srclib/$lib" ]] || { echo "error: $srclib/$lib not found" >&2; exit 1; }
done
[[ -f "$srcinc/Population.hpp" ]] || { echo "error: $srcinc does not look like the evogen include dir" >&2; exit 1; }

echo "Staging artifact content in $stage"
rm -rf "$stage"
mkdir -p "$stage/lib" "$stage/include"
cp "${libs[@]/#/"$srclib"/}" "$stage/lib/"
cp "$srcinc"/*.hpp "$stage/include/"

# ---------------------------------------------------------------------------
# ELF patch: RUNPATH ($ORIGIN) + SONAME (= file name) on the two libraries
# built by the evogen project. The Intel runtime libraries are left as-is.
# ---------------------------------------------------------------------------
patch_lib() {
    local file="$1"
    local soname
    soname="$(basename "$file")"
    if command -v patchelf >/dev/null 2>&1; then
        patchelf --set-rpath '$ORIGIN' "$file"
        patchelf --set-soname "$soname" "$file"
    elif command -v python3 >/dev/null 2>&1; then
        python3 - "$file" "$soname" <<'PYEOF'
# Minimal ELF64 RUNPATH/SONAME patcher (no third-party tools).
#
# Adds DT_RUNPATH = "$ORIGIN" and DT_SONAME = <file name>, using only the
# unused slots that follow the DT_NULL terminator inside .dynamic:
#
#   [ ... entries ... ][NULL -> RUNPATH][SONAME][NULL]["$ORIGIN\0"][soname\0"]
#
# The string values are offsets relative to DT_STRTAB; since .dynamic lies
# at a higher address than .dynstr in these files, strings stored inside
# .dynamic are reachable with plain positive offsets. No other section,
# string, or symbol is touched.
import struct
import sys

DT_NULL, DT_STRTAB, DT_RUNPATH, DT_RPATH, DT_SONAME = 0, 5, 29, 15, 14

def die(msg):
    sys.exit(f"elf patch: {msg}")

def patch(path, soname):
    with open(path, "r+b") as f:
        data = f.read()
        if data[:4] != b"\x7fELF" or data[4] != 2 or data[5] != 1:
            die(f"{path}: not a 64-bit little-endian ELF")
        e_phoff, = struct.unpack_from("<Q", data, 0x20)
        e_phentsize, e_phnum = struct.unpack_from("<HH", data, 0x36)

        dynamic = None  # (file offset, vaddr, filesz)
        for i in range(e_phnum):
            off = e_phoff + i * e_phentsize
            p_type, = struct.unpack_from("<I", data, off)
            if p_type == 2:  # PT_DYNAMIC
                p_offset, p_vaddr = struct.unpack_from("<QQ", data, off + 8)
                p_filesz, = struct.unpack_from("<Q", data, off + 32)
                dynamic = (p_offset, p_vaddr, p_filesz)
                break
        if dynamic is None:
            die(f"{path}: no PT_DYNAMIC segment")

        p_offset, p_vaddr, p_filesz = dynamic
        nslots = p_filesz // 16

        null_idx = None
        strtab_vaddr = None
        for k in range(nslots):
            tag, val = struct.unpack_from("<QQ", data, p_offset + 16 * k)
            if tag == DT_NULL:
                null_idx = k
                break
            if tag == DT_STRTAB:
                strtab_vaddr = val
            if tag in (DT_RUNPATH, DT_RPATH, DT_SONAME):
                die(f"{path}: already has a RUNPATH/RPATH/SONAME; refusing to patch")
        if null_idx is None or strtab_vaddr is None:
            die(f"{path}: DT_NULL or DT_STRTAB not found")

        # slots: null_idx = RUNPATH, +1 = SONAME, +2 = new NULL,
        #        +3 = "$ORIGIN\0", +4 = soname + "\0" (one 16-byte slot each)
        if null_idx + 5 > nslots:
            die(f"{path}: no slack after DT_NULL in .dynamic "
                f"(need 4 free slots, have {nslots - null_idx - 1}); "
                "install patchelf and re-run")
        for s in (b"$ORIGIN\x00", soname.encode() + b"\x00"):
            if len(s) > 16:
                die(f"{path}: string too long for one slot: {s!r}")

        p_dyn_vaddr = p_vaddr
        def stroff(slot):
            return (p_dyn_vaddr + 16 * slot - strtab_vaddr) & 0xFFFFFFFFFFFFFFFF

        f.seek(p_offset + 16 * null_idx)
        f.write(struct.pack("<QQ", DT_RUNPATH, stroff(null_idx + 3)))
        f.seek(p_offset + 16 * (null_idx + 1))
        f.write(struct.pack("<QQ", DT_SONAME, stroff(null_idx + 4)))
        f.seek(p_offset + 16 * (null_idx + 2))
        f.write(struct.pack("<QQ", DT_NULL, 0))
        f.seek(p_offset + 16 * (null_idx + 3))
        f.write(b"$ORIGIN\x00".ljust(16, b"\x00"))
        f.seek(p_offset + 16 * (null_idx + 4))
        f.write((soname.encode() + b"\x00").ljust(16, b"\x00"))
    print(f"elf patch: {path}: RUNPATH=$ORIGIN SONAME={soname} (python fallback)")

patch(sys.argv[1], sys.argv[2])
PYEOF
    else
        echo "error: need patchelf or python3 to patch the libraries" >&2
        exit 1
    fi
}

patch_lib "$stage/lib/libevogen.so"
patch_lib "$stage/lib/libplinkio.so"

if command -v readelf >/dev/null 2>&1; then
    echo "Patched ELF metadata:"
    for lib in libevogen.so libplinkio.so; do
        readelf -d "$stage/lib/$lib" | grep -E 'SONAME|RPATH|RUNPATH' | sed "s|^|  $lib |"
    done
fi
echo "Done. Next step: build the wrapper (tools/README) and run tools/install_artifact.jl"
