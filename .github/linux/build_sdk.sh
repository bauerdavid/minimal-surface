#!/bin/bash
# Builds ITK and SimpleITK as static libraries inside the manylinux_2_28 container, and leaves
# the parts setup.py needs in $WORK/sdk (its SuperBuild layout: include/, lib/, ITK-build/lib/),
# together with the license files of everything that was built.
#
#   docker run --rm -v "$PWD:/work" -e SITK_TAG=v2.5.6 quay.io/pypa/manylinux_2_28_x86_64 /work/.github/linux/build_sdk.sh
set -euxo pipefail

WORK=${WORK:-/work}
SITK_TAG=${SITK_TAG:?set SITK_TAG to the SimpleITK release to build}
BUILD=$WORK/sdk-build
OUT=$WORK/sdk

export PATH=/opt/python/cp312-cp312/bin:$PATH
python -m pip install cmake ninja

# ccache, if it is available, makes a rerun after a failure or a cancelled job cheap
if yum install -y epel-release ccache >/dev/null 2>&1 || command -v ccache; then
    export CCACHE_DIR=$WORK/ccache CCACHE_MAXSIZE=4G
    export CMAKE_C_COMPILER_LAUNCHER=ccache CMAKE_CXX_COMPILER_LAUNCHER=ccache
    ccache --zero-stats
fi

rm -rf "$BUILD" "$OUT"
git clone --depth 1 --branch "$SITK_TAG" https://github.com/SimpleITK/SimpleITK "$BUILD/src"

# Position independent code, because these static libraries are linked into a shared library.
cmake -S "$BUILD/src/SuperBuild" -B "$BUILD/bin" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_FLAGS=-fPIC -DCMAKE_CXX_FLAGS=-fPIC \
    -DBUILD_SHARED_LIBS=OFF \
    -DWRAP_DEFAULT=OFF \
    -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF
cmake --build "$BUILD/bin" -j"$(nproc)"

command -v ccache >/dev/null && ccache --show-stats || true

# Where the SuperBuild put things differs between SimpleITK versions, so look
ls "$BUILD/bin"
sitk_lib=$(dirname "$(find "$BUILD/bin" -maxdepth 3 -name 'libSimpleITKCommon-*.a' | head -1)")
sitk_include=$(find "$BUILD/bin" -maxdepth 3 -type d -name 'SimpleITK-*' -path '*include*' | head -1)
itk_lib=$(dirname "$(find "$BUILD/bin/ITK-build" -maxdepth 2 -name 'libITKCommon-*.a' | head -1)")
echo "SimpleITK libraries: $sitk_lib; headers: $sitk_include; ITK libraries: $itk_lib"
test -n "$sitk_lib" && test -n "$sitk_include" && test -n "$itk_lib"

mkdir -p "$OUT/lib" "$OUT/include" "$OUT/ITK-build/lib"
cp -L "$sitk_lib"/*.a "$OUT/lib/"
cp -R "$sitk_include" "$OUT/include/"
cp -L "$itk_lib"/*.a "$OUT/ITK-build/lib/"

# The licenses of SimpleITK, ITK and the third-party libraries ITK builds into itself
mkdir -p "$OUT/licenses/SimpleITK" "$OUT/licenses/ITK"
cp "$BUILD"/src/LICENSE* "$BUILD"/src/NOTICE* "$OUT/licenses/SimpleITK/" 2>/dev/null || true
cp "$BUILD"/bin/ITK/LICENSE* "$BUILD"/bin/ITK/NOTICE* "$OUT/licenses/ITK/" 2>/dev/null || true
find "$BUILD/bin/ITK/Modules/ThirdParty" -maxdepth 6 -type f \
    \( -iname 'LICENSE*' -o -iname 'LICENCE*' -o -iname 'COPYING*' -o -iname 'COPYRIGHT*' -o -iname 'NOTICE*' \) \
    -not -path '*/test/*' -not -path '*/doc/*' -not -path '*/examples/*' | sort | while read -r f; do
    rel=${f#"$BUILD"/bin/ITK/Modules/ThirdParty/}
    mkdir -p "$OUT/licenses/ITK-third-party/$(dirname "$rel")"
    cp "$f" "$OUT/licenses/ITK-third-party/$rel"
done
ls -R "$OUT/licenses" | head -100
test -f "$OUT/licenses/SimpleITK/LICENSE" && test -f "$OUT/licenses/ITK/LICENSE"

rm -rf "$BUILD"
du -sh "$OUT"
chown -R "${HOST_UID:-0}:${HOST_GID:-0}" "$OUT" "${CCACHE_DIR:-/nonexistent}" 2>/dev/null || true
