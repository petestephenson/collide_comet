import matplotlib

matplotlib.use("Agg")

import pytest

import collide_comet as cc


@pytest.fixture(scope="session")
def cases(tmp_path_factory):
    """A tiny synthetic five-case set (12^3 cells)."""
    folder = tmp_path_factory.mktemp("demo_runs")
    pe, sw, names = cc.synthetic.write_demo_cases(folder, n_cells=12)
    return cc.CaseSet(folder, pe, sw, names, max_loaded=2)


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    import matplotlib.pyplot as plt

    plt.close("all")
