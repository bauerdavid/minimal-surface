#!/bin/bash
# Builds the wheel for one Python version inside the manylinux_2_28 container, against the SDK
# that build_sdk.sh made, and repairs it (auditwheel) and adds the third-party licenses.
#
#   docker run --rm -v "$PWD:/work" -e PYTAG=cp313-cp313 quay.io/pypa/manylinux_2_28_x86_64 /work/.github/linux/build_wheel.sh
set -euxo pipefail

WORK=${WORK:-/work}
PYTAG=${PYTAG:?set PYTAG, for example cp313-cp313}
PY=/opt/python/$PYTAG/bin/python

cd "$WORK"
export SIMPLEITK_PATH=$WORK/sdk
rm -rf /tmp/dist "$WORK/wheelhouse"
"$PY" -m pip wheel . --no-deps --wheel-dir /tmp/dist -v
auditwheel show /tmp/dist/*.whl
auditwheel repair /tmp/dist/*.whl --plat manylinux_2_28_x86_64 --wheel-dir "$WORK/wheelhouse"

"$PY" .github/scripts/licenses.py linux --wheel "$WORK"/wheelhouse/*.whl \
    --sdk-licenses "$WORK/sdk/licenses" --eigen "$WORK/src/minimal-surface/code/eigen-3.4.0"

ls -l "$WORK/wheelhouse"
chown -R "${HOST_UID:-0}:${HOST_GID:-0}" "$WORK/wheelhouse" 2>/dev/null || true
