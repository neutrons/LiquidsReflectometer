"""Unit tests for the StitchingType enum"""
import pytest

from lr_reduction.stitching import StitchingType


class TestStitchingType:
    """Test the StitchingType enum functionality"""

    def test_stitching_type_values(self):
        """Test that StitchingType enum has correct values"""
        assert StitchingType.NONE.value == "None"
        assert StitchingType.AUTOMATIC_AVERAGE.value == "AutomaticAverage"
        assert StitchingType.ABSOLUTE_NORMALIZATION.value == "AbsoluteNormalization"

    def test_from_value_none(self):
        """Test from_value with None input"""
        result = StitchingType.from_value(None)
        assert result == StitchingType.NONE

    def test_from_value_none_string(self):
        """Test from_value with 'none' string"""
        result = StitchingType.from_value("none")
        assert result == StitchingType.NONE

    def test_from_value_automatic_average(self):
        """Test from_value with 'automaticaverage' string"""
        result = StitchingType.from_value("automaticaverage")
        assert result == StitchingType.AUTOMATIC_AVERAGE

    def test_from_value_absolute_normalization(self):
        """Test from_value with 'automaticaverage' string"""
        result = StitchingType.from_value("AbsoluteNormalization")
        assert result == StitchingType.ABSOLUTE_NORMALIZATION

    def test_from_value_case_insensitive(self):
        """Test that from_value is case insensitive"""
        assert StitchingType.from_value("NONE") == StitchingType.NONE
        assert StitchingType.from_value("None") == StitchingType.NONE
        assert StitchingType.from_value("AutomaticAverage") == StitchingType.AUTOMATIC_AVERAGE
        assert StitchingType.from_value("AUTOMATICAVERAGE") == StitchingType.AUTOMATIC_AVERAGE
        assert StitchingType.from_value("absolutenormalization") == StitchingType.AUTOMATIC_AVERAGE
        assert StitchingType.from_value("ABSOLUTENORMALIZATION") == StitchingType.ABSOLUTE_NORMALIZATION

    def test_from_value_invalid(self):
        """Test that from_value raises ValueError for invalid input"""
        with pytest.raises(ValueError, match="Invalid StitchingType value"):
            StitchingType.from_value("invalid")

        with pytest.raises(ValueError, match="Invalid StitchingType value"):
            StitchingType.from_value("manual")


if __name__ == "__main__":
    pytest.main([__file__])
