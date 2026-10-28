#!/bin/bash
# Quick check that a fresh build works:
#   regression/check_install.sh [BUILD_DIR] [DUMP_FILE]
#
# 1. main_cluster --check validates every regression/bench config.
# 2. Test_exactness evaluates energies and gradients of the same configs.
#
# The configs name the checkout of whoever wrote them (init_file), so the
# script runs on copies that point at this checkout. BUILD_DIR defaults to
# ../build. With DUMP_FILE the Test_exactness dump is kept there; dumps of two
# builds on the same machine can be compared with cmp.
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../../.." && pwd)
build=${1:-$here/../build}
dump=${2:-}

for exe in main_cluster Test_exactness; do
    if [ ! -x "$build/bin/$exe" ]; then
        echo "Missing $build/bin/$exe: build it first (see the README)" >&2
        exit 1
    fi
done

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

configs=()
for f in "$here"/configs/*.json "$here"/bench/*.json; do
    out="$tmp/$(basename "$f")"
    sed "s#/home/mrojasve/Documents/DDG/MiraculousMembranes#$root#g" "$f" > "$out"
    configs+=("$out")
done

for c in "${configs[@]}"; do
    printf '%-40s ' "$(basename "$c")"
    "$build/bin/main_cluster" --check "$c"
done

"$build/bin/Test_exactness" "$tmp/dump.txt" "${configs[@]}" > "$tmp/exactness.log" 2>&1
echo "Test_exactness: $(wc -c < "$tmp/dump.txt") bytes written"
if [ -n "$dump" ]; then
    # the dump names each config by its path: keep only the file name
    sed -E 's#^(config ).*/#\1#' "$tmp/dump.txt" > "$dump"
    echo "Dump kept in $dump"
fi
echo "Install check passed."
