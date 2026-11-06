"""Unit tests for the SlitScanCorrection Mantid algorithm."""

import os

import numpy as np
import pytest
from mantid.api import IEventWorkspace
from mantid.kernel import FloatTimeSeriesProperty, Int32TimeSeriesProperty
from mantid.simpleapi import ConvertUnits, CreateSampleWorkspace, DeleteWorkspace, LoadEventNexus, mtd

from lr_reduction.exceptions import LrValidationError
from lr_reduction.operations.algorithms.cd_attenuation_correction import CdAttenuationCorrection
from lr_reduction.operations.algorithms.slit_scan_correction import (
    FACTOR_LOG,
    SCALE_MULTIPLIER_LOG,
    SlitScanCorrection,
)
from lr_reduction.properties.cd_attenuation import MENU_LOG

# Event weights are stored as float32.
RTOL = 1e-6


@pytest.fixture
def make_direct_beam():
    """Factory for a small four-pixel event workspace in wavelength, with `n` recorded as the
    ScaleMultiplier series unless it is None; deletes each workspace on teardown."""
    names = []

    def make(*n, menu=None):
        name = mtd.unique_hidden_name()
        CreateSampleWorkspace(
            WorkspaceType="Event",
            Function="Flat background",
            NumBanks=1,
            BankPixelWidth=2,
            XMin=8000.0,
            XMax=56000.0,
            BinWidth=48000.0,
            OutputWorkspace=name,
        )
        ConvertUnits(InputWorkspace=name, Target="Wavelength", OutputWorkspace=name)
        run = mtd[name].getRun()
        if n:
            series = FloatTimeSeriesProperty(SCALE_MULTIPLIER_LOG)
            for second, value in enumerate(n):
                series.addValue(f"2020-01-01T00:00:0{second}", value)
            run.addProperty(SCALE_MULTIPLIER_LOG, series, True)
        if menu is not None:
            series = Int32TimeSeriesProperty(MENU_LOG)
            series.addValue("2020-01-01T00:00:00", menu)
            run.addProperty(MENU_LOG, series, True)
        names.append(name)
        return name

    yield make
    for name in names + [n + "_out" for n in names] + [n + "_mid" for n in names]:
        if mtd.doesExist(name):
            DeleteWorkspace(name)


def run_algorithm(algorithm_class, input_workspace, output_workspace=None):
    algorithm = algorithm_class()
    algorithm.PyInit()
    algorithm.setProperty("InputWorkspace", input_workspace)
    algorithm.setProperty("OutputWorkspace", output_workspace or f"{input_workspace}_out")
    algorithm.PyExec()
    return algorithm.getProperty("OutputWorkspace").value


def weights(ws):
    return np.concatenate([ws.getSpectrum(i).getWeights() for i in range(ws.getNumberHistograms())])


def errors(ws):
    return np.concatenate([ws.getSpectrum(i).getWeightErrors() for i in range(ws.getNumberHistograms())])


def test_algorithm_metadata():
    algorithm = SlitScanCorrection()

    assert algorithm.name() == "SlitScanCorrection"
    assert algorithm.category() == "Reflectometry\\SNS"
    assert algorithm.version() == 1


def test_weights_and_errors_are_multiplied_by_n_squared(make_direct_beam):
    direct_beam = make_direct_beam(7.0, 7.0)
    before = weights(mtd[direct_beam])

    corrected = run_algorithm(SlitScanCorrection, direct_beam)

    assert isinstance(corrected, IEventWorkspace)
    assert corrected.getNumberEvents() == mtd[direct_beam].getNumberEvents()
    np.testing.assert_allclose(weights(corrected), 49 * before, rtol=RTOL)
    np.testing.assert_allclose(errors(corrected), 49 * errors(mtd[direct_beam]), rtol=RTOL)
    assert corrected.getRun().getProperty(FACTOR_LOG).value == 49
    # The input is left as it was.
    np.testing.assert_array_equal(weights(mtd[direct_beam]), before)
    assert mtd[direct_beam].getSpectrum(0).getEventType().name == "TOF"


@pytest.mark.parametrize("n", [(1.0,), ()], ids=["n_is_one", "no_log"])
def test_unit_factor_leaves_events_unweighted(make_direct_beam, n):
    direct_beam = make_direct_beam(*n)

    corrected = run_algorithm(SlitScanCorrection, direct_beam)

    assert corrected.name() != direct_beam
    assert corrected.getSpectrum(0).getEventType().name == "TOF"
    assert corrected.getNumberEvents() == mtd[direct_beam].getNumberEvents()
    np.testing.assert_array_equal(weights(corrected), weights(mtd[direct_beam]))
    assert corrected.getRun().getProperty(FACTOR_LOG).value == 1
    assert not mtd[direct_beam].getRun().hasProperty(FACTOR_LOG)


@pytest.mark.parametrize(
    "n, message",
    [((7.0, 5.0), "varies"), ((0.0,), "positive"), ((-2.0,), "positive"), ((float("nan"),), "NaN")],
    ids=["changing", "zero", "negative", "nan"],
)
def test_invalid_log_is_rejected(make_direct_beam, n, message):
    with pytest.raises(LrValidationError, match=message):
        run_algorithm(SlitScanCorrection, make_direct_beam(*n))


def test_order_with_cd_correction_does_not_matter(make_direct_beam):
    direct_beam = make_direct_beam(3.0, menu=5)

    cd_first = run_algorithm(
        SlitScanCorrection, run_algorithm(CdAttenuationCorrection, direct_beam, f"{direct_beam}_mid").name()
    )
    slit_first = run_algorithm(
        CdAttenuationCorrection, run_algorithm(SlitScanCorrection, direct_beam, f"{direct_beam}_mid").name()
    )

    np.testing.assert_allclose(weights(cd_first), weights(slit_first), rtol=RTOL)


@pytest.mark.datarepo
def test_measured_slit_scan_run(nexus_dir):
    """REF_L_220063 was measured with N = 10."""
    name = "slit_scan_input_220063"
    LoadEventNexus(Filename=os.path.join(nexus_dir, "REF_L_220063.nxs.h5"), OutputWorkspace=name)
    try:
        corrected = run_algorithm(SlitScanCorrection, name)
        assert corrected.getNumberEvents() == mtd[name].getNumberEvents()
        np.testing.assert_allclose(weights(corrected).sum(), 100 * weights(mtd[name]).sum(), rtol=1e-5)
    finally:
        for ws in (name, f"{name}_out"):
            if mtd.doesExist(ws):
                DeleteWorkspace(ws)
