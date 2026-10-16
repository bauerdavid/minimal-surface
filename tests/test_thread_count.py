"""The extension has to work, and give the same result, whatever number of OpenMP threads there is.

The evolution used to run the two sides in parallel, each with an inner team of half the
available threads. With a single thread that half was zero, and the process crashed (SIGFPE or a
segmentation fault, depending on the platform). The two sides now share one team, but a small
number of threads is still the case most likely to break, and the result must not depend on the
number of threads. OMP_NUM_THREADS is read when the OpenMP runtime is loaded, so each case runs
in a fresh interpreter.
"""
import os
import subprocess
import sys

import pytest

pytest.importorskip("minimal_surface")

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))

SCRIPT = """
import sys
sys.path.insert(0, {tests_dir!r})
import test_calculate_inputs as t
image, phi = t.make_volumes()
for meeting in (False, True):
    out = t.calculate(phi, image, t.ALONG_X, meeting)
    print("RESULT", meeting, int((out >= 0).sum()), repr(float(out.sum())))
""".format(tests_dir=TESTS_DIR)


def run_with_threads(threads):
    """Run the calculation with OMP_NUM_THREADS=threads, or with it unset when threads is None."""
    env = dict(os.environ)
    env.pop("OMP_NUM_THREADS", None)
    if threads is not None:
        env["OMP_NUM_THREADS"] = str(threads)
    done = subprocess.run([sys.executable, "-c", SCRIPT], capture_output=True, text=True, env=env, timeout=300)
    assert done.returncode == 0, "exit status %s with OMP_NUM_THREADS=%s\n%s" % (
        done.returncode, threads, done.stderr[-2000:])
    results = {}
    for line in done.stdout.splitlines():
        if line.startswith("RESULT"):
            _, meeting, inside, total = line.split()
            results[meeting] = (int(inside), float(total))
    assert set(results) == {"False", "True"}, done.stdout[-2000:]
    return results


@pytest.mark.parametrize("threads", [1, 2, 3])
def test_a_small_number_of_threads_does_not_crash(threads):
    results = run_with_threads(threads)

    assert all(inside > 0 for inside, _ in results.values())


def test_the_result_does_not_depend_on_the_number_of_threads():
    single = run_with_threads(1)
    # None is the default: no OMP_NUM_THREADS, so the extension picks its own limit. 16 is more
    # threads than most test machines have cores.
    for threads in (2, 4, 16, None):
        other = run_with_threads(threads)
        for meeting in single:
            assert other[meeting][0] == single[meeting][0], (threads, meeting)
            assert other[meeting][1] == pytest.approx(single[meeting][1], rel=1e-9), (threads, meeting)
