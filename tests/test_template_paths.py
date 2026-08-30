"""Template path resolution — the scaling-factor file must not depend on cwd.

These guard a defect that the campaign gate structurally cannot see: the
sanctioned command is `cd tests/ && pytest`, and the fixture's
`<scaling_factor_file>data/sf_197912_Si_auto.cfg</scaling_factor_file>` resolves
against that cwd. Any other invocation — `pixi run pytest` from the repo root,
an IDE, `pytest tests/test_reduction.py -k …` — resolved it to nothing and then
failed four frames later with `TypeError: cannot unpack non-iterable
EventWorkspace object`. Both tests below chdir away from `tests/` deliberately,
so they keep failing if the anchoring is reverted even though CI never runs
from the repo root.
"""

import os

import pytest

from lr_reduction import template


def test_scaling_factor_missing_file_raises(tmp_path):
    """A missing scaling-factor file must fail where it is detected.

    The single call site unpacks a 4-tuple, so the old `return workspace`
    branch could never succeed — it turned a clear file-not-found into a
    TypeError in the caller. The workspace argument is deliberately None: the
    file check happens before the workspace is touched, so this needs no data.
    """
    missing = str(tmp_path / "nonexistent_sf.cfg")
    with pytest.raises(FileNotFoundError) as excinfo:
        template.scaling_factor(missing, None)
    assert missing in str(excinfo.value)


def test_template_scaling_factor_path_is_cwd_independent(monkeypatch, tmp_path, template_dir):
    """Reading a template from an unrelated cwd must still locate its
    scaling-factor file."""
    template_path = os.path.join(template_dir, "template.xml")
    monkeypatch.chdir(tmp_path)

    template_data = template.read_template(template_path, 7)

    assert template_data.scaling_factor_file, "template declares no scaling factor file"
    assert os.path.isfile(template_data.scaling_factor_file), (
        f"scaling factor file not resolved from cwd={os.getcwd()}: "
        f"{template_data.scaling_factor_file}"
    )
