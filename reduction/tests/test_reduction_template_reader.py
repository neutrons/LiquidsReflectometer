import pytest
from lr_reduction import template
from lr_reduction.instrument_settings import InstrumentSettings
from lr_reduction.reduction_template_reader import ReductionParameters


class TestReductionParameters:
    def test_from_dict(self):
        r"""Verify method from_dict raises when passed some nonsense"""

        redparms = ReductionParameters()
        # valid data dictionary
        redparms.from_dict(dict(two_backgrounds=True))
        assert redparms.two_backgrounds
        # invalid data dictionary and not permissible
        with pytest.raises(ValueError, match="data_dir contains invalid entries") as excinfo:
            redparms.from_dict(dict(nonsense=True), permissible=False)
        assert "data_dir contains invalid entries" == str(excinfo.value)
        # invalid data dictionary and permissible
        redparms.from_dict(dict(nonsense=True))

    def test_two_backgrounds(self):
        r"""Verify the xml dump writes what we want"""
        redparms = ReductionParameters()
        redparms.two_backgrounds = True
        assert "<two_backgrounds>True</two_backgrounds>" in redparms.to_xml()

    def test_emission_delay(self):
        r"""Verify the xml dump writes the emission delay option"""
        redparms = ReductionParameters()

        # Default should be True
        assert redparms.use_emission_time is True
        assert "<use_emission_time>True</use_emission_time>" in redparms.to_xml()

        redparms.use_emission_time = False
        assert "<use_emission_time>False</use_emission_time>" in redparms.to_xml()

    def test_instrument_settings_default(self):
        r"""Verify the xml dump writes the instrument settings"""
        redparms = ReductionParameters()
        for key, val in InstrumentSettings().__dict__.items():
            assert f"<{key}>{val}</{key}>" in redparms.to_xml()

    def test_instrument_settings_custom(self, template_dir):
        """Verify loading custom instrument settings from template.xml"""
        template_path = f"{template_dir}/template_with_instrument_settings.xml"
        sequence_number = 188298
        template_data = template.read_template(template_path, sequence_number)
        assert template_data.apply_instrument_settings is True
        assert template_data.source_detector_distance == 2.0
        assert template_data.sample_detector_distance == 3.0
        assert template_data.num_x_pixels == 4
        assert template_data.num_y_pixels == 5
        assert template_data.pixel_width == 6.0
        assert template_data.xi_reference == 7.0
        assert template_data.s1_sample_distance == 8.0


if __name__ == "__main__":
    pytest.main(__file__)
