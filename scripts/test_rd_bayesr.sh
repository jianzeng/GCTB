#!/usr/bin/env bash
# Run after building scr/gctb; uses the same compiler and include paths as Makefile.
set -euo pipefail
cd "$(dirname "$0")/.."
compiler=${CXX:-c++}
flags=(-std=c++11 -O1 -I "${EIGEN3_INCLUDE_DIR:?}" -I "${BOOST_LIB:?}" -DUNIX -Dfopen64=fopen)
if [[ $(uname -s) == Darwin ]]; then
    compiler=${CXX:-/opt/homebrew/opt/llvm/bin/clang++}
    flags+=(-stdlib=libc++ -fopenmp)
else
    flags+=(-fopenmp)
fi
objects=()
for object in scr/*.o; do
    [[ $object == scr/main.o ]] || objects+=("$object")
done
test_dir=$(mktemp -d "${TMPDIR:-/tmp}/gctb-rd-bayesr.XXXXXX")
for test in test_rd_bayesr test_annotation_integration; do
    "$compiler" "${flags[@]}" "scripts/$test.cpp" "${objects[@]}" -lz -o "$test_dir/$test"
    OMP_NUM_THREADS=2 "$test_dir/$test"
done
