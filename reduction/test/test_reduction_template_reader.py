# third-party imports
import pytest

# lr_reduction imports
from lr_reduction.reduction_template_reader import ReductionParameters


class TestReductionParameters:

    def test_two_backgrounds(self):
        r"""verify the xml dump writes what we want"""
        redparms = ReductionParameters()
        redparms.two_backgrounds = True
        assert "<two_backgrounds_flag>True</two_backgrounds_flag>" in redparms.to_xml()

    def test_emission_delay(self):
        r"""verify the xml dump writes the emission delay option"""
        redparms = ReductionParameters()

        # Default should be True
        assert redparms.use_emission_time == True
        assert "<use_emission_time>True</use_emission_time>" in redparms.to_xml()

        redparms.use_emission_time = False
        assert "<use_emission_time>False</use_emission_time>" in redparms.to_xml()

    def test_from_dict(self):
        r"""verify method from_dict raises when passed some nonsense"""

        redparms = ReductionParameters()
        # valid data dictionary
        redparms.from_dict(dict(two_backgrounds=True))
        assert redparms.two_backgrounds
        # invalid data dictionary and not permissible
        with pytest.raises(ValueError) as excinfo:
            redparms.from_dict(dict(nonsense=True), permissible=False)
        assert "data_dir contains invalid entries" == str(excinfo.value)
        # invalid data dictionary and permissible
        redparms.from_dict(dict(nonsense=True))


if __name__ == "__main__":
    pytest.main(__file__)
