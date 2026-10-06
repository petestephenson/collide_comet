"""Execute the notebooks end to end (on synthetic demo data).

Slow, so only run on request:  COLLIDE_COMET_TEST_NOTEBOOKS=1 pytest tests/test_notebooks.py
"""

import os
import shutil
from pathlib import Path

import pytest

NOTEBOOK_DIR = Path(__file__).resolve().parents[1] / "notebooks"
NOTEBOOKS = sorted(p.name for p in NOTEBOOK_DIR.glob("*.ipynb"))

pytestmark = pytest.mark.skipif(
    os.environ.get("COLLIDE_COMET_TEST_NOTEBOOKS") != "1",
    reason="set COLLIDE_COMET_TEST_NOTEBOOKS=1 to execute the notebooks",
)


@pytest.mark.parametrize("name", NOTEBOOKS)
def test_notebook_runs(name, tmp_path, monkeypatch):
    nbformat = pytest.importorskip("nbformat")
    nbclient = pytest.importorskip("nbclient")
    work = tmp_path / "notebooks"
    work.mkdir()
    shutil.copy(NOTEBOOK_DIR / name, work / name)
    shutil.copy(NOTEBOOK_DIR / "run_config.py", work / "run_config.py")
    monkeypatch.setenv("COLLIDE_COMET_FIGURES", str(tmp_path / "Figures"))
    monkeypatch.setenv("MPLBACKEND", "Agg")
    nb = nbformat.read(work / name, as_version=4)
    nbclient.NotebookClient(nb, timeout=1200, kernel_name="python3",
                            resources={"metadata": {"path": str(work)}}).execute()
    assert list((tmp_path / "Figures").rglob("*.png")), "notebook saved no figures"
