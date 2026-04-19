"""
Integration tests for `wgbstools segment`.

Runs the full Python orchestrator + compiled `segmentor` C++ binary on a
small region using real beta files, and golden-matches the block output.

These tests exist specifically to lock byte-identity across the segmentor
performance refactor (ring-buffer DP, pool hoisting, -O2). A regression
here means the refactor changed numerics — not a routine test failure.

Skipped when the beta fixtures are absent from the cluster.

Invocation:
    pytest -m integration                        # run all integration
    pytest tests/integration/test_segment*       # run this file only

Updating golden files:
    pytest -m integration --update-golden
"""
import subprocess
import sys
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Paths / fixtures
# ---------------------------------------------------------------------------

REPO        = Path(__file__).parents[2]
WGBSTOOLS   = REPO / "src" / "python" / "wgbs_tools.py"
SEGMENTOR   = REPO / "src" / "segment_betas" / "segmentor"
DATA_DIR    = Path(__file__).parents[1] / "data" / "segment"

BETA_DIR    = Path("/sci/labs/tommy/nloyfer/Brain/other_attempts/nanopore/from_Amir")
BETAS       = sorted(BETA_DIR.glob("*.beta"))[:10]

# Small region: deterministic, fast, exercises the DP + stitching at a
# chunk boundary, and includes the 294042/294175 adjacent-CpG pair that
# was the canary for the float/double regression in the ring buffer.
REGION      = "chr19:100000-500000"
REGION_TAG  = "chr19_100k_500k"

# Larger region: exercises multi-chunk stitching (chunk_size=60000 default
# puts this at ~7 chunks → 3 merge rounds).
LARGE_REGION = "chr19:100000-2000000"
LARGE_TAG    = "chr19_100k_2M"


# ---------------------------------------------------------------------------
# Decorators
# ---------------------------------------------------------------------------

integration = pytest.mark.integration

_skip_no_data = pytest.mark.skipif(
    len(BETAS) < 10 or not SEGMENTOR.exists(),
    reason=(
        f"Missing prerequisites: {len(BETAS)} beta files found (need 10) "
        f"or segmentor binary not built at {SEGMENTOR}"
    ),
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run_segment(tmp_path: Path, region: str, threads: int = 4) -> Path:
    out = tmp_path / "blocks.bed"
    cmd = [
        sys.executable, str(WGBSTOOLS), "segment",
        "--betas", *[str(b) for b in BETAS],
        "-r", region,
        "-o", str(out),
        "-@", str(threads),
    ]
    subprocess.run(cmd, check=True, capture_output=True)
    return out


def golden_match(actual: Path, golden: Path, update: bool) -> None:
    actual_bytes = actual.read_bytes()
    if update:
        golden.parent.mkdir(parents=True, exist_ok=True)
        golden.write_bytes(actual_bytes)
        pytest.skip(f"Golden updated: {golden}. Re-run without --update-golden.")
    assert golden.exists(), (
        f"Golden file missing: {golden}\n"
        f"Run with --update-golden to create it."
    )
    expected = golden.read_bytes()
    if actual_bytes != expected:
        # Produce a short diff-style message for debugging.
        a_lines = actual_bytes.decode().splitlines()
        e_lines = expected.decode().splitlines()
        msg = [f"segment output differs from golden {golden.name}:"]
        msg.append(f"  actual={len(a_lines)} lines, golden={len(e_lines)} lines")
        for i, (a, e) in enumerate(zip(a_lines, e_lines), 1):
            if a != e:
                msg.append(f"  line {i}: expected {e!r}")
                msg.append(f"           got      {a!r}")
                if len(msg) > 10:
                    msg.append("  ...")
                    break
        pytest.fail("\n".join(msg))


# ===========================================================================
# Small-region golden match — byte-identity at the DP/stitching boundary
# ===========================================================================

@integration
@_skip_no_data
class TestSegmentSmallRegion:
    def test_exit_code_and_non_empty(self, tmp_path):
        out = run_segment(tmp_path, REGION)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_valid_bed_format(self, tmp_path):
        out = run_segment(tmp_path, REGION)
        for line in out.read_text().splitlines():
            fields = line.split("\t")
            assert len(fields) == 5, f"expected 5 columns, got {fields!r}"
            chrom, start, end, start_cpg, end_cpg = fields
            assert chrom == "chr19"
            assert int(start) < int(end)
            assert int(start_cpg) < int(end_cpg)

    def test_matches_golden(self, tmp_path, request):
        out = run_segment(tmp_path, REGION)
        golden_match(
            out,
            DATA_DIR / f"{REGION_TAG}.blocks.bed",
            request.config.getoption("--update-golden"),
        )


# ===========================================================================
# Multi-chunk stitching golden — locks cross-chunk merge numerics
# ===========================================================================

@integration
@_skip_no_data
@pytest.mark.heavy
class TestSegmentMultiChunk:
    def test_matches_golden(self, tmp_path, request):
        out = run_segment(tmp_path, LARGE_REGION, threads=8)
        golden_match(
            out,
            DATA_DIR / f"{LARGE_TAG}.blocks.bed",
            request.config.getoption("--update-golden"),
        )
