import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def invoke(*args):
    return subprocess.run(
        [sys.executable, "-m", "mar_certified", *args], capture_output=True, text=True
    )


def test_query_cli_writes_certificate(tmp_path):
    output = tmp_path / "query.json"
    result = invoke(
        "query",
        str(ROOT / "examples/switching.json"),
        "--alpha",
        "1|2",
        "--tol",
        "1e-12",
        "--output",
        str(output),
    )
    assert result.returncode == 0, result.stderr
    data = json.loads(output.read_text())
    assert data["tolerance_met"]
    assert data["unit"] == "bit"


def test_budget_limited_cli_writes_interval_and_exits_nonzero(tmp_path):
    output = tmp_path / "query.json"
    result = invoke(
        "query",
        str(ROOT / "examples/switching.json"),
        "--alpha",
        "1|2",
        "--method",
        "endpoint",
        "--max-cells",
        "0",
        "--output",
        str(output),
    )
    assert result.returncode != 0
    data = json.loads(output.read_text())
    assert not data["tolerance_met"]
    assert data["interval"] is not None


def test_cli_rejects_comparable_antichain():
    result = invoke("query", str(ROOT / "examples/switching.json"), "--alpha", "1|1,2")
    assert result.returncode != 0
    assert "Traceback" not in result.stderr
