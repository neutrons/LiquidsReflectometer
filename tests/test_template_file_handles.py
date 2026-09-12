"""`read_template` must not leak the file handle it opens.

Why this needs care: in CPython the local `fd` refcount drops at return, so a
naive "is it closed afterwards?" assertion passes even against the buggy code.
The leak only shows when something keeps the frame alive — which is exactly what
a failing test does, via the traceback holding the frame. That is why the buggy
version produced dozens of unclosed-file ResourceWarnings in a suite with
failures and only a handful in a passing one.

So the test retains the frame deliberately, reproducing the condition rather
than the average case.
"""

import builtins
import gc
import os
import sys

import pytest

from lr_reduction import template


@pytest.fixture
def captured_open(monkeypatch):
    """Record every file object `open` hands out during the test."""
    handles = []
    real_open = builtins.open

    def _recording_open(*args, **kwargs):
        handle = real_open(*args, **kwargs)
        handles.append(handle)
        return handle

    monkeypatch.setattr(builtins, "open", _recording_open)
    return handles


def test_read_template_closes_its_file_even_when_the_frame_is_retained(captured_open, template_dir):
    """The frame is kept alive via an exception traceback, which is what defeats
    refcount cleanup and made the leak visible in the first place."""
    template_path = os.path.join(template_dir, "template.xml")

    try:
        template.read_template(template_path, 7)
        raise RuntimeError("retain the frame")
    except RuntimeError:
        # sys.exc_info()[2] holds frames from this call chain; keeping it bound
        # is what a failing test's traceback does implicitly.
        traceback_ref = sys.exc_info()[2]
        assert traceback_ref is not None
        gc.collect()

        opened = [h for h in captured_open if getattr(h, "name", None) == template_path]
        assert opened, "read_template did not open the template through builtins.open"
        unclosed = [h for h in opened if not h.closed]
        assert not unclosed, (
            f"{len(unclosed)} handle(s) on {template_path} still open while the frame is retained; "
            f"read_template must use `with open(...)`"
        )


def test_forked_read_template_closes_its_file_too(captured_open, template_dir):
    """The fork copy carries the same code and must not drift from it."""
    from lr_reduction import new_reduction_from_template

    template_path = os.path.join(template_dir, "template.xml")

    try:
        new_reduction_from_template.read_template(template_path, 7)
        raise RuntimeError("retain the frame")
    except RuntimeError:
        traceback_ref = sys.exc_info()[2]
        assert traceback_ref is not None
        gc.collect()

        opened = [h for h in captured_open if getattr(h, "name", None) == template_path]
        assert opened, "the forked read_template did not open the template through builtins.open"
        assert not [h for h in opened if not h.closed], (
            f"the fork copy at new_reduction_from_template leaks a handle on {template_path}"
        )
