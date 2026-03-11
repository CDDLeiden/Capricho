#!/usr/bin/env python3
"""Run a command and report wall time and peak memory usage.

Outputs a machine-parseable line:
    BENCHMARK_RESULT: <wall_seconds> <peak_rss_mb>

Uses resource.getrusage(RUSAGE_CHILDREN) for peak RSS measurement.
On macOS, ru_maxrss is reported in bytes; on Linux, in kilobytes.
"""

import resource
import subprocess
import sys
import time

cmd = sys.argv[1:]
if not cmd:
    print("Usage: measure_command.py <command> [args...]", file=sys.stderr)
    sys.exit(1)

start = time.time()
proc = subprocess.run(cmd)
elapsed = time.time() - start

usage = resource.getrusage(resource.RUSAGE_CHILDREN)
peak_rss_raw = usage.ru_maxrss
if sys.platform == "darwin":
    peak_rss_mb = peak_rss_raw / (1024 * 1024)
else:
    peak_rss_mb = peak_rss_raw / 1024

print(f"BENCHMARK_RESULT: {elapsed:.2f} {peak_rss_mb:.1f}")
sys.exit(proc.returncode)
