"""Template path resolution — the scaling-factor file must not depend on cwd.

A template's `<scaling_factor_file>` used to be resolved against whatever
directory the process started in. The sanctioned gate command is
`cd tests/ && pytest`, so it resolved there and nowhere else: any other
invocation failed four frames later with `TypeError: cannot unpack
non-iterable EventWorkspace object`. Resolution is now anchored to the
template's own directory, and the fixtures declare their scaling-factor file as
a sibling — which means the gate exercises the anchoring too, so deleting it
reds the gate rather than only a repo-root run.

The four tests below still chdir away deliberately: the anchoring is what makes
the location irrelevant, and a guard that only runs where the answer is already
right proves nothing.
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


def _declared_scaling_factor_files(template_path):
    """Every raw <scaling_factor_file> string in the XML.

    All of them, not the first: the template carries one per sequence and the
    test below reads sequence 7, so checking only sequence 1's would assert
    about a value the test never exercises.
    """
    declared = re.findall(r"<scaling_factor_file>([^<]*)</scaling_factor_file>", Path(template_path).read_text())
    assert declared, f"no scaling_factor_file element in {template_path}"
    return declared


def test_template_scaling_factor_path_is_cwd_independent(monkeypatch, template_dir, tmp_path):
    """Reading a template from an unrelated cwd must resolve to the intended file.

    Asserts identity, not mere existence: `isfile` alone would pass on any
    same-named file that happened to be reachable, which is the failure mode
    anchoring is supposed to remove rather than introduce. The precondition
    assertion keeps the test honest if someone "fixes" the fixture by
    absolutising it — then this test is no longer exercising anything and says
    so.
    """
    template_path = os.path.join(template_dir, "template.xml")
    declared = _declared_scaling_factor_files(template_path)
    absolute = [d for d in declared if os.path.isabs(d)]
    assert not absolute, (
        f"precondition lost: {template_path} declares absolute scaling factor path(s) "
        f"({absolute}), so this test no longer exercises cwd-independence"
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
    # Through read_template, so what is pinned is the seam callers actually use.
    resolved = template.read_template(str(template_path), 1).scaling_factor_file

    assert os.path.samefile(resolved, workdir / "sf.cfg"), (
        f"as-given (cwd) candidate must win; got {resolved}"
    )


def test_unresolvable_path_is_returned_unchanged(tmp_path):
    """Nothing resolves: hand back what the template declared, so the error the
    caller raises names the author's path rather than an invented one."""
    template_path = tmp_path / "t.xml"
    template_path.write_text(
        "<Reduction><DataSeries><RefLData>"
        "<scaling_factor_file>nowhere/sf.cfg</scaling_factor_file>"
        "</RefLData></DataSeries></Reduction>"
    )
    assert template.read_template(str(template_path), 1).scaling_factor_file == "nowhere/sf.cfg"
