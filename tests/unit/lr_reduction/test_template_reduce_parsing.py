"""Tests for run number parsing logic used in the Template Reduce tab."""

import pytest

from apps.template_reduce import _parse_run_numbers


class TestParseRunNumbers:
    """Test _parse_run_numbers utility function."""

    def test_single_run(self):
        assert _parse_run_numbers("211029") == [211029]

    def test_csv_runs(self):
        assert _parse_run_numbers("211029,211030,211031") == [211029, 211030, 211031]

    def test_range_runs(self):
        assert _parse_run_numbers("211029-211031") == [211029, 211030, 211031]

    def test_mixed_csv_and_range(self):
        result = _parse_run_numbers("211029,211033-211035")
        assert result == [211029, 211033, 211034, 211035]

    def test_whitespace_is_stripped(self):
        assert _parse_run_numbers(" 211029 , 211030 ") == [211029, 211030]

    def test_empty_string_raises(self):
        with pytest.raises(ValueError):
            _parse_run_numbers("")

    def test_invalid_input_raises(self):
        with pytest.raises(ValueError):
            _parse_run_numbers("abc")

    def test_reversed_range_raises(self):
        with pytest.raises(ValueError):
            _parse_run_numbers("211031-211029")
