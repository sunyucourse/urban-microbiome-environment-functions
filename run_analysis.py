"""Reproduce the included post-profiling statistics without altering reference outputs."""
from pathlib import Path
import argparse
import os
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results",
                        help="Fresh output directory, relative to current working directory")
    parser.add_argument("--pilot", action="store_true",
                        help="Run only factor/category screening and factor two-stage tests")
    args = parser.parse_args()
    out = Path(args.output).resolve()
    if out.exists():
        raise SystemExit("Output already exists; use a fresh --output directory to preserve prior runs.")
    out.mkdir(parents=True)
    for p in (ROOT/"data/stats_seed").glob("*.csv"):
        shutil.copyfile(p, out/p.name)
    env = dict(os.environ)
    env["URBAN_RESULTS_DIR"] = str(out)
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["OMP_NUM_THREADS"] = "1"
    scripts = sorted((ROOT/"scripts").glob("*.py"))
    if args.pilot:
        scripts = scripts[:1]
    for script in scripts:
        start = time.perf_counter()
        print(f"\n=== {script.name} ===", flush=True)
        subprocess.run([sys.executable, str(script)], env=env, check=True)
        print(f"Elapsed: {time.perf_counter()-start:.2f} seconds", flush=True)
    if not args.pilot:
        subprocess.run([sys.executable, str(ROOT/"validate_results.py"),
                        "--output", str(out)], check=True)

if __name__ == "__main__":
    main()
