"""
NOTE:
EDITED TO REMOVE MANTID PART THAT WAS THROWING AN ERROR!!
Need to decide whether to build something different for the new workflow or just use the existing templates.
Gravity Direction needs reconnecting if stick with this version...
Stripped back RefRed template reader - has some parts removed that aren't used in new workflow (might need to put them back)
and some new parts added.
Parts that are commented out are currently not linked up but could be useful.
"""

"""
LIST OF PARAMS TO ADD:
- Norm flag
- DB file
- q mthod
- Autoscale
- useCalcTheta
- Qline_threshold
- Scale Factor
- 
"""

import time
import xml.dom.minidom
from typing import Optional

from lr_reduction import __version__ as VERSION
#from lr_reduction.gravity_correction import GravityDirection
from instrument_settings import InstrumentSettings

# Get the mantid version being used, if available
#try:
#    import mantid

#    MANTID_VERSION = mantid.__version__
#except:  # noqa: E722
#    MANTID_VERSION = "None"


class ReductionParameters:
    """
    Class that hold the parameters for the reduction of a single data set.
    """

    def __init__(self):
        # Signal selection
        self.data_peak_range = [140, 150]
        self.subtract_background = True
        #self.two_backgrounds: bool = False
        self.background_roi = [137, 153, 0, 0]
        self.tof_range = [9600.0, 21600.0]
        #self.select_tof_range = True

        self.data_x_range = [115, 210]

        # Data files
        self.data_files = [0]
        #self.norm_file = 0

        # Q range
        self.q_min = 0.001
        self.q_step = 0.001
        self.const_q = False

        # Scattering angle
        self.angle_offset = 0.0
        #self.angle_offset_error = 0.0

        # Dead time correction
        self.dead_time: bool = False
        #self.paralyzable: bool = True
        self.dead_time_value = 4.2
        self.dead_time_tof_step = 100
        #self.use_dead_time_threshold = False
        #self.dead_time_threshold: Optional[Float] = 1.5

        # TODO: Compare the defaults to defaults in nr_config.
        # New parts added from new reduction config scheme
        self.norm_scale = False # Ensure name is different from prior normalization flag on DB
        self.DB_file = None
        self.q_method = None
        self.autoscale = True   # Might need to check for new flags from recent change and make sure this doesn't conflict.
        self.use_calc_theta = True
        self.qline_threshold = 0.66
        self.scale_factor = 1.0

        # Instrument geometry parameters
        instrument_settings = InstrumentSettings()
        for key, value in instrument_settings.__dict__.items():
            setattr(self, key, value)

        # Calculate emission time delay instead of using an effective distance for all wavelengths
        #self.use_emission_time: bool = True

        # Gravity correction
        #self.gravity_direction = None

    def from_dict(self, data_dict, permissible=True):
        """
        Update object's attributes with a dictionary with entries of the type  attribute_name: attribute_value.

        Parameters
        ----------
        permissible : bool
            allow keys in data_dict that are not attribute names of ReductionParameters instances. Reading from
            `data_dict` will result in this instance having new attributes not defined in `__init__()`

        Raises
        ------
        ValueError
            when `permissible=False` and one entry (or more) of the dictionary is not an attribute of this object
        """

        # check all keys are data_dict are attributes of object `self`
        attribute_names = list(vars(self))
        if permissible is False and all(key in attribute_names for key in data_dict) is False:
            raise ValueError("data_dir contains invalid entries")
        # update attribute values
        for k, v in data_dict.items():
            setattr(self, k, v)

    def to_xml(self):
        """
        Create XML from the current data.
        """
        _xml = "<RefLData>\n"
        _xml += "<peak_selection_type>narrow</peak_selection_type>\n"   # Not sure what this means.
        _xml += "<from_peak_pixels>%s</from_peak_pixels>\n" % str(self.data_peak_range[0])
        _xml += "<to_peak_pixels>%s</to_peak_pixels>\n" % str(self.data_peak_range[1])
        _xml += "<peak_discrete_selection>N/A</peak_discrete_selection>\n"  # Not sure what this means.
        _xml += "<background_flag>%s</background_flag>\n" % str(self.subtract_background)
        #_xml += "<two_backgrounds>%s</two_backgrounds>\n" % str(self.two_backgrounds)
        _xml += "<back_roi1_from>%s</back_roi1_from>\n" % str(self.background_roi[0])
        _xml += "<back_roi1_to>%s</back_roi1_to>\n" % str(self.background_roi[1])
        _xml += "<back_roi2_from>%s</back_roi2_from>\n" % str(self.background_roi[2])
        _xml += "<back_roi2_to>%s</back_roi2_to>\n" % str(self.background_roi[3])
        #_xml += "<tof_range_flag>%s</tof_range_flag>\n" % str(self.select_tof_range)
        _xml += "<from_tof_range>%s</from_tof_range>\n" % str(self.tof_range[0])
        _xml += "<to_tof_range>%s</to_tof_range>\n" % str(self.tof_range[1])
        _xml += "<data_sets>%s</data_sets>\n" % ",".join([str(i) for i in self.data_files])
        _xml += "<x_min_pixel>%s</x_min_pixel>\n" % str(self.data_x_range[0])
        _xml += "<x_max_pixel>%s</x_max_pixel>\n" % str(self.data_x_range[1])
        #_xml += "<x_range_flag>%s</x_range_flag>\n" % str(self.data_x_range_flag)

        #_xml += "<norm_dataset>%s</norm_dataset>\n" % str(self.norm_file)

        _xml += "<q_min>%s</q_min>\n" % str(self.q_min)
        _xml += "<q_step>%s</q_step>\n" % str(self.q_step)
        _xml += "<const_q>%s</const_q>\n" % str(self.const_q)

        # Angle offset
        _xml += "<angle_offset>%s</angle_offset>\n" % str(self.angle_offset)
        #_xml += "<angle_offset_error>%s</angle_offset_error>\n" % str(self.angle_offset_error)

        # Dead time correction
        _xml += "<dead_time_correction>%s</dead_time_correction>\n" % str(self.dead_time)
        #_xml += "<dead_time_paralyzable>%s</dead_time_paralyzable>\n" % str(self.paralyzable)
        _xml += "<dead_time_value>%s</dead_time_value>\n" % str(self.dead_time_value)
        _xml += "<dead_time_tof_step>%s</dead_time_tof_step>\n" % str(self.dead_time_tof_step)
        #_xml += "<use_dead_time_threshold>%s</use_dead_time_threshold>\n" % str(self.use_dead_time_threshold)
        #_xml += "<dead_time_threshold>%s</dead_time_threshold>\n" % str(self.dead_time_threshold)

        # Instrument settings
        _xml += "<apply_instrument_settings>%s</apply_instrument_settings>\n" % str(self.apply_instrument_settings)
        _xml += "<source_detector_distance>%s</source_detector_distance>\n" % str(self.source_detector_distance)
        _xml += "<sample_detector_distance>%s</sample_detector_distance>\n" % str(self.sample_detector_distance)
        _xml += "<num_x_pixels>%s</num_x_pixels>\n" % str(self.num_x_pixels)
        _xml += "<num_y_pixels>%s</num_y_pixels>\n" % str(self.num_y_pixels)
        _xml += "<pixel_width>%s</pixel_width>\n" % str(self.pixel_width)
        _xml += "<xi_reference>%s</xi_reference>\n" % str(self.xi_reference)
        _xml += "<s1_sample_distance>%s</s1_sample_distance>\n" % str(self.s1_sample_distance)
        _xml += "<wavelength_resolution_function>%s</wavelength_resolution_function>\n" % str(
            self.wavelength_resolution_function
        )

        # Gravity correction
        #if self.gravity_direction is not None:
        #    _xml += "<gravity_direction>%s</gravity_direction>\n" % str(self.gravity_direction)  # -1, 0, 1

        # Emission time correction
        #_xml += "<use_emission_time>%s</use_emission_time>\n" % str(self.use_emission_time)

        # New parts added from new reduction config scheme
        _xml += "<norm_scale>%s</norm_scale>\n" % str(self.norm_scale)
        _xml += "<DB_file>%s</DB_file>\n" % str(self.DB_file)
        _xml += "<q_method>%s</q_method>\n" % str(self.q_method)
        _xml += "<autoscale>%s</autoscale>\n" % str(self.autoscale)
        _xml += "<use_calc_theta>%s</use_calc_theta>\n" % str(self.use_calc_theta)
        _xml += "<qline_threshold>%s</qline_threshold>\n" % str(self.qline_threshold)
        _xml += "<scale_factor>%s</scale_factor>\n" % str(self.scale_factor)

        _xml += "</RefLData>\n"

        return _xml

    def from_xml_element(self, instrument_dom):
        """
        Read in data from XML

        Parameters
        ----------
        instrument_dom : xml.dom.Document
        """
        # Peak from/to pixels
        self.data_peak_range = [
            getIntElement(instrument_dom, "from_peak_pixels"),
            getIntElement(instrument_dom, "to_peak_pixels"),
        ]

        # low resolution range
        self.data_x_range = [getIntElement(instrument_dom, "x_min_pixel"), getIntElement(instrument_dom, "x_max_pixel")]

        # background flag
        self.subtract_background = getBoolElement(instrument_dom, "background_flag", default=self.subtract_background)

        # use two backgrounds flag
        #self.two_backgrounds = getBoolElement(instrument_dom, "two_backgrounds", default=self.two_backgrounds)

        # background from/to pixels
        self.background_roi = [
            getIntElement(instrument_dom, "back_roi1_from"),
            getIntElement(instrument_dom, "back_roi1_to"),
            getIntElement(instrument_dom, "back_roi2_from"),
            getIntElement(instrument_dom, "back_roi2_to"),
        ]

        # TOF range
        #self.select_tof_range = getBoolElement(instrument_dom, "tof_range_flag", default=self.select_tof_range)
        self.tof_range = [
            getFloatElement(instrument_dom, "from_tof_range"),
            getFloatElement(instrument_dom, "to_tof_range"),
        ]

        self.data_files = getIntList(instrument_dom, "data_sets")

        #self.norm_file = getIntElement(instrument_dom, "norm_dataset")

        # Q cut
        self.q_min = getFloatElement(instrument_dom, "q_min", default=self.q_min)
        self.q_step = getFloatElement(instrument_dom, "q_step", default=self.q_step)
        self.const_q = getBoolElement(instrument_dom, "const_q", default=False)

        # Angle offset
        self.angle_offset = getFloatElement(instrument_dom, "angle_offset", default=self.angle_offset)
        #self.angle_offset_error = getFloatElement(instrument_dom, "angle_offset_error", default=self.angle_offset_error)

        # Dead time correction
        self.dead_time = getBoolElement(instrument_dom, "dead_time_correction", default=self.dead_time)
        #self.paralyzable = getBoolElement(instrument_dom, "dead_time_paralyzable", default=self.paralyzable)
        self.dead_time_value = getFloatElement(instrument_dom, "dead_time_value", default=self.dead_time_value)
        self.dead_time_tof_step = getFloatElement(instrument_dom, "dead_time_tof_step", default=self.dead_time_tof_step)
        #self.use_dead_time_threshold = getBoolElement(
        #    instrument_dom, "use_dead_time_threshold", default=self.use_dead_time_threshold
        #)
        #self.dead_time_threshold = getFloatElement(
        #    instrument_dom, "dead_time_threshold", default=self.dead_time_threshold
        #)

        # New parts added from new reduction config scheme
        self.norm_scale = getBoolElement(instrument_dom, "norm_scale", default=self.norm_scale)
        self.DB_file = getBoolElement(instrument_dom, "DB_file", default=self.DB_file)
        self.q_method = getBoolElement(instrument_dom, "q_method", default=self.q_method)
        self.autoscale = getBoolElement(instrument_dom, "autoscale", default=self.autoscale)
        self.use_calc_theta = getBoolElement(instrument_dom, "use_calc_theta", default=self.use_calc_theta)
        self.qline_threshold = getBoolElement(instrument_dom, "qline_threshold", default=self.qline_threshold)
        self.scale_factor = getBoolElement(instrument_dom, "scale_factor", default=self.scale_factor)


        # Instrument settings
        self.apply_instrument_settings = getBoolElement(instrument_dom, "apply_instrument_settings", default=False)
        self.source_detector_distance = getFloatElement(
            instrument_dom, "source_detector_distance", default=self.source_detector_distance
        )
        self.sample_detector_distance = getFloatElement(
            instrument_dom, "sample_detector_distance", default=self.sample_detector_distance
        )
        self.num_x_pixels = getIntElement(instrument_dom, "num_x_pixels", default=self.num_x_pixels)
        self.num_y_pixels = getIntElement(instrument_dom, "num_y_pixels", default=self.num_y_pixels)
        self.pixel_width = getFloatElement(instrument_dom, "pixel_width", default=self.pixel_width)
        self.xi_reference = getFloatElement(instrument_dom, "xi_reference", default=self.xi_reference)
        self.s1_sample_distance = getFloatElement(instrument_dom, "s1_sample_distance", default=self.s1_sample_distance)
        self.wavelength_resolution_function = getStringElement(
            instrument_dom, "wavelength_resolution_function", default=self.wavelength_resolution_function
        )

        #self.gravity_direction = GravityDirection.from_value(
        #    getIntElement(instrument_dom, "gravity_direction", default=self.gravity_direction)
        #)

        # Emission time
        # Defaults to True, but will be skipped if the necessary meta data is not found
        #self.use_emission_time = getBoolElement(instrument_dom, "use_emission_time", default=True)


#############################################
### Utility functions to read XML content ###
#############################################


def getText(nodelist):
    """Utility method to extract text out of an XML node"""
    rc = ""
    for node in nodelist:
        if node.nodeType == node.TEXT_NODE:
            rc = rc + node.data
    return rc


def getContent(dom, tag):
    """Returns the content of a tag within a dom object"""
    element_list = dom.getElementsByTagName(tag)
    return getText(element_list[0].childNodes) if len(element_list) > 0 else None


def getIntElement(dom, tag, default=None):
    """Parse an integer element from the dom object"""
    value = getContent(dom, tag)
    return int(value) if value is not None else default


def getIntList(dom, tag, default=[]):
    """Parse a list of integers from the dom object"""
    value = getContent(dom, tag)
    if value is not None and len(value.strip()) > 0:
        return list(map(int, value.split(",")))
    else:
        return default


def getFloatElement(dom, tag, default=None):
    """
    Parse a float element from the DOM object.

    Parameters
    ----------
    dom : xml.dom.minidom.Document
        The DOM object to parse the element from.
    tag : str
        The name of the tag to search for in the DOM.
    default : float, optional
        The default value to return if the tag is not found or its content is None.

    Returns
    -------
    float
        The float value of the tag's content if it exists, otherwise the default value.
    """
    value = getContent(dom, tag)
    return float(value) if value is not None else default


def getFloatList(dom, tag, default=[]):
    """Parse a list of floats from the dom object"""
    value = getContent(dom, tag)
    if value is not None and len(value.strip()) > 0:
        return list(map(float, value.split(",")))
    else:
        return default


def getStringElement(dom, tag, default=""):
    """Parse a string element from the dom object"""
    value = getContent(dom, tag)
    return value if value is not None else default


def getStringList(dom, tag, _default=[]):
    """Parse a list of strings from the dom object"""
    elem_list = []
    element_list = dom.getElementsByTagName(tag)
    if len(element_list) > 0:
        for l_ in element_list:
            elem_list.append(getText(l_.childNodes).strip())
    return elem_list


def getBoolElement(dom, tag, true_tag="true", default=False):
    """Parse a boolean element from the dom object"""
    value = getContent(dom, tag)
    return value.lower() == true_tag.lower() if value is not None else default


###############################################
### Functions to read/write a template file ###
###############################################


def to_xml(data_sets):
    """
    Create XML from the current data.

    Parameters
    ----------
    data_sets : list
        List of ReductionParameters instances

    Returns
    -------
    str
        XML string
    """
    _xml = "<Reduction>\n"
    _xml += "    <instrument_name>REFL</instrument_name>\n"
    _xml += "    <timestamp>%s</timestamp>\n" % time.ctime()
    _xml += "    <version>%s</version>\n" % VERSION
    #_xml += "    <mantid_version>%s</mantid_version>\n" % MANTID_VERSION
    _xml += "    <generator>lr_reduction-%s</generator>\n" % VERSION
    _xml += "<DataSeries>\n"
    for item in data_sets:
        _xml += item.to_xml()
    _xml += "</DataSeries>\n"
    _xml += "</Reduction>\n"
    return _xml


def from_xml(xml_str):
    """
    Read in data from XML string

    Parameters
    ----------
    xml_str : str
        String representation of a list of ReductionParameters instances

    Returns
    -------
    list
        List of ReductionParameters instances
    """
    data_sets = []
    dom = xml.dom.minidom.parseString(xml_str)

    element_list = dom.getElementsByTagName("Data")
    if len(element_list) == 0:
        element_list = dom.getElementsByTagName("RefLData")

    if len(element_list) > 0:
        for item in element_list:
            if item is not None:
                data_set = ReductionParameters()
                data_set.from_xml_element(item)
                data_sets.append(data_set)

    if len(data_sets) == 0:
        data_sets = [ReductionParameters()]

    return data_sets
