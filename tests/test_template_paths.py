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
import re
from pathlib import Path

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


def _declared_scaling_factor_file(template_path):
    """The raw <scaling_factor_file> string, straight out of the XML."""
    match = re.search(r"<scaling_factor_file>([^<]*)</scaling_factor_file>", Path(template_path).read_text())
    assert match, f"no scaling_factor_file element in {template_path}"
    return match.group(1)


def test_template_scaling_factor_path_is_cwd_independent(monkeypatch, tmp_path, template_dir):
    """Reading a template from an unrelated cwd must resolve to the intended file.

    Asserts identity, not mere existence: `isfile` alone would pass on any
    same-named file that happened to be reachable, which is the failure mode
    anchoring is supposed to remove rather than introduce. The precondition
    assertion keeps the test honest if someone "fixes" the fixture by
    absolutising it — then this test is no longer exercising anything and says
    so.
    """
    template_path = os.path.join(template_dir, "template.xml")
    declared = _declared_scaling_factor_file(template_path)
    assert declared, "template declares no scaling factor file"
    assert not os.path.isabs(declared), (
        f"precondition lost: {template_path} declares an absolute scaling factor path "
        f"({declared}), so this test no longer exercises cwd-independence"
    )

    monkeypatch.chdir(tmp_path)
    template_data = template.read_template(template_path, 7)

    expected = os.path.join(template_dir, "sf_197912_Si_auto.cfg")
    assert os.path.isfile(template_data.scaling_factor_file), (
        f"not resolved from cwd={os.getcwd()}: {template_data.scaling_factor_file}"
    )
    assert os.path.samefile(template_data.scaling_factor_file, expected), (
        f"resolved to the wrong file: {template_data.scaling_factor_file} is not {expected}"
    )


def test_as_given_path_wins_over_template_anchor(monkeypatch, tmp_path):
    """Precedence: a path that resolves against the cwd must win.

    The facility and the `cd tests/` gate both depend on as-given resolution;
    anchoring is a fallback. Decoys sit at *both* locations so that deleting
    the as-given candidate — or reordering the two — turns this red rather than
    silently switching which file a reduction reads.
    """
    template_home = tmp_path / "templates"
    template_home.mkdir()
    workdir = tmp_path / "work"
    workdir.mkdir()

    (workdir / "sf.cfg").write_text("cwd-relative decoy\n")
    (template_home / "sf.cfg").write_text("template-dir decoy\n")
    template_path = template_home / "t.xml"
    template_path.write_text(
        '<Reduction><DataSeries><RefLData>'
        '<scaling_factor_file>sf.cfg</scaling_factor_file>'
        '</RefLData></DataSeries></Reduction>'
    )

    monkeypatch.chdir(workdir)
    resolved = template._resolve_scaling_factor_file("sf.cfg", str(template_path))

    assert os.path.samefile(resolved, workdir / "sf.cfg"), (
        f"as-given (cwd) candidate must win; got {resolved}"
    )


def test_unresolvable_path_is_returned_unchanged(tmp_path):
    """Nothing resolves: hand back what the template declared, so the error the
    caller raises names the author's path rather than an invented one."""
    template_path = tmp_path / "t.xml"
    template_path.write_text("<Reduction/>")
    assert template._resolve_scaling_factor_file("nowhere/sf.cfg", str(template_path)) == "nowhere/sf.cfg"
