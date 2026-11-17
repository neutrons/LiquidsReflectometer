"""
NR_Reduction reduces one run without reduce() having been called first.

reduce_from_template builds an NR_Reduction and calls _reduce_single_run directly. What reduce() sets up for the
figures (whether to draw them, where they are kept) must therefore exist from construction.
"""

import pytest

from lr_reduction.nr_reduction_calc import NR_Reduction
from lr_reduction.nr_reduction_config import NRReductionConfig


def reduction(plot):
    config = NRReductionConfig()
    config.RBnum = [1]
    config.DBname = ["direct_beam.dat"]
    config.RB_Ymin = [140]
    config.RB_Ymax = [160]
    config.plotON = plot
    return NR_Reduction(config)


@pytest.mark.parametrize("plot", [True, False])
def test_a_new_reduction_knows_whether_to_draw_figures(plot):
    assert reduction(plot).create_figures is plot


def test_a_new_reduction_has_nowhere_else_to_keep_its_figures_than_its_own_list():
    first, second = reduction(True), reduction(True)
    assert first.fig_store == []
    assert first.fig_store is not second.fig_store
