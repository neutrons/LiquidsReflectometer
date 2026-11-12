import pytest
from mantid.simpleapi import AddSampleLog, CreateSampleWorkspace, DeleteWorkspace, mtd

from lr_reduction.exceptions import ConfigError, SampleLogsError
from lr_reduction.models.config import RunParameters
from lr_reduction.models.run_data import RunData
from lr_reduction.operations.direct_beam import cd_order, direct_beam_roi

MENU_LOG = "BL4B:Actuator:Menu"

# Attenuator menu codes (LSB-first stage flags) and the Cd thickness they put in the beam, in cm.
NO_CD = (0, 0.0)  # no stage in the beam
THIN_CD = (1, 57.5e-4)  # stage 1
THICK_CD = (5, (57.5 + 126.5 + 123.0) * 1e-4)  # stages 1 and 3


@pytest.fixture
def make_run():
    """Factory for direct-beam runs carrying an attenuator menu log, deleted afterwards."""
    names = []

    def _make(run_number, menu_code=None):
        name = mtd.unique_hidden_name()
        CreateSampleWorkspace(WorkspaceType="Event", NumBanks=1, BankPixelWidth=1, NumEvents=1, OutputWorkspace=name)
        if menu_code is not None:
            AddSampleLog(Workspace=name, LogName=MENU_LOG, LogText=str(menu_code), LogType="Number")
        names.append(name)
        return RunData(workspace=name, run_numbers=(run_number,))

    yield _make
    for name in names:
        DeleteWorkspace(name)


# --- cd_order ----------------------------------------------------------------------------


def test_cd_order_sorts_by_ascending_thickness(make_run):
    thick = make_run(1, THICK_CD[0])
    none = make_run(2, NO_CD[0])
    thin = make_run(3, THIN_CD[0])

    ordered = cd_order([thick, none, thin], flip_atten=False)

    assert [run for run, _ in ordered] == [none, thin, thick]
    assert [thickness for _, thickness in ordered] == pytest.approx([NO_CD[1], THIN_CD[1], THICK_CD[1]])


def test_cd_order_keeps_input_order_for_equal_thickness(make_run):
    first = make_run(1, THIN_CD[0])
    second = make_run(2, THIN_CD[0])
    assert [run for run, _ in cd_order([first, second], flip_atten=False)] == [first, second]
    assert [run for run, _ in cd_order([second, first], flip_atten=False)] == [second, first]


def test_cd_order_applies_flip_atten(make_run):
    # Inverted polarity: menu code 15 (every flag set) means no stage is in the beam.
    all_flags = make_run(1, 15)
    no_flags = make_run(2, 0)
    ordered = cd_order([no_flags, all_flags], flip_atten=True)
    assert [run for run, _ in ordered] == [all_flags, no_flags]
    assert ordered[0][1] == 0.0


def test_cd_order_without_attenuator_log_raises(make_run):
    with pytest.raises(SampleLogsError):
        cd_order([make_run(1)], flip_atten=False)


# --- direct_beam_roi ---------------------------------------------------------------------


def test_direct_beam_roi_pads_peak_and_uses_low_res():
    params = RunParameters(peak={"min": 140, "max": 150}, low_res={"min": 90, "max": 160}, peak_fit={"pad": 3})
    assert direct_beam_roi(params) == ((90, 160), (137, 153))


def test_direct_beam_roi_default_pad_without_peak_fit():
    params = RunParameters(peak={"min": 140, "max": 150}, low_res={"min": 90, "max": 160})
    assert direct_beam_roi(params) == ((90, 160), (139, 151))


@pytest.mark.parametrize(
    ("params", "match"),
    [
        (RunParameters(low_res={"min": 90, "max": 160}), "no `peak` range"),
        (RunParameters(peak={"min": 140, "max": 150}), "no `low_res` range"),
        (
            RunParameters(peak={"min": 1, "max": 10}, low_res={"min": 90, "max": 160}, peak_fit={"pad": 2}),
            "start the box at pixel -1",
        ),
    ],
    ids=["no peak", "no low_res", "pad below pixel 0"],
)
def test_direct_beam_roi_rejects_incomplete_or_out_of_range_box(params, match):
    with pytest.raises(ConfigError, match=match):
        direct_beam_roi(params)
