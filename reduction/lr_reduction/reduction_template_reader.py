"""
    RefRed template reader.
    Adapted from Mantid code.
"""
import time
import xml.dom.minidom

from . import __version__ as VERSION

# Get the mantid version being used, if available
try:
    import mantid
    MANTID_VERSION = mantid.__version__
except:
    MANTID_VERSION = "None"


class ReductionParameters(object):

    def __init__(self):
        # Signal selection
        self.data_peak_range = [140, 150]
        self.subtract_background = True
        self.two_backgrounds: bool = False
        self.background_roi = [137, 153, 0, 0]
        self.tof_range = [9600., 21600.]
        self.select_tof_range = True

        self.data_x_range_flag = True
        self.data_x_range = [115,210]

        # Normalization
        self.apply_normalization = True
        self.norm_peak_range = [140, 150]
        self.subtract_norm_background = True
        self.norm_background_roi = [137, 153]
        self.norm_x_range_flag = True
        self.norm_x_range = [115,210]

        # Data files
        self.data_files = [0]
        self.norm_file = 0

        # Clean up options: cut first and last points as needed
        self.pre_cut = 1
        self.post_cut = 1

        # Q range
        self.q_min = 0.001
        self.q_step = 0.001
        self.auto_q_binning = False

        # Scattering angle
        self.tthd_value = 0
        self.ths_value = 0
        self.angle_offset = 0.0
        self.angle_offset_error = 0.0

        # Scaling factor file
        self.scaling_factor_file = ''
        self.scaling_factor_flag = True
        self.slits_width_flag = True

        # Incident medium list and selected value
        self.incident_medium_list = ['air']
        self.incident_medium_index_selected = 0

        # Dead time correction
        self.dead_time:bool = False
        self.paralyzable:bool = True
        self.dead_time_value = 4.2
        self.dead_time_tof_step = 100

        # Calculate emission time delay instead of using an effective distance for all wavelengths
        self.use_emission_time:bool = True

    def from_dict(self, data_dict, permissible=True):
        r"""
        Update object's attributes with a dictionary with entries of the type  attribute_name: attribute_value.

        Parameters
        ----------
        permissible: bool
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
        _xml  = "<RefLData>\n"
        _xml += "<peak_selection_type>narrow</peak_selection_type>\n"
        _xml += "<from_peak_pixels>%s</from_peak_pixels>\n" % str(self.data_peak_range[0])
        _xml += "<to_peak_pixels>%s</to_peak_pixels>\n" % str(self.data_peak_range[1])
        _xml += "<peak_discrete_selection>N/A</peak_discrete_selection>\n"
        _xml += "<background_flag>%s</background_flag>\n" % str(self.subtract_background)
        _xml += "<two_backgrounds>%s</two_backgrounds>\n" % str(self.two_backgrounds)
        _xml += "<back_roi1_from>%s</back_roi1_from>\n" % str(self.background_roi[0])
        _xml += "<back_roi1_to>%s</back_roi1_to>\n" % str(self.background_roi[1])
        _xml += "<back_roi2_from>%s</back_roi2_from>\n" % str(self.background_roi[2])
        _xml += "<back_roi2_to>%s</back_roi2_to>\n" % str(self.background_roi[3])
        _xml += "<tof_range_flag>%s</tof_range_flag>\n" % str(self.select_tof_range)
        _xml += "<from_tof_range>%s</from_tof_range>\n" % str(self.tof_range[0])
        _xml += "<to_tof_range>%s</to_tof_range>\n" % str(self.tof_range[1])
        _xml += "<data_sets>%s</data_sets>\n" % ','.join([str(i) for i in self.data_files])
        _xml += "<x_min_pixel>%s</x_min_pixel>\n" % str(self.data_x_range[0])
        _xml += "<x_max_pixel>%s</x_max_pixel>\n" % str(self.data_x_range[1])
        _xml += "<x_range_flag>%s</x_range_flag>\n" % str(self.data_x_range_flag)

        _xml += "<tthd_value>%s</tthd_value>\n" % str(self.tthd_value)
        _xml += "<ths_value>%s</ths_value>\n" % str(self.ths_value)

        _xml += "<norm_flag>%s</norm_flag>\n" % str(self.apply_normalization)
        _xml += "<norm_x_range_flag>%s</norm_x_range_flag>\n" % str(self.norm_x_range_flag)
        _xml += "<norm_x_max>%s</norm_x_max>\n" % str(self.norm_x_range[1])
        _xml += "<norm_x_min>%s</norm_x_min>\n" % str(self.norm_x_range[0])

        _xml += "<norm_from_peak_pixels>%s</norm_from_peak_pixels>\n" % str(self.norm_peak_range[0])
        _xml += "<norm_to_peak_pixels>%s</norm_to_peak_pixels>\n" % str(self.norm_peak_range[1])
        _xml += "<norm_background_flag>%s</norm_background_flag>\n" % str(self.subtract_norm_background)
        _xml += "<norm_from_back_pixels>%s</norm_from_back_pixels>\n" % str(self.norm_background_roi[0])
        _xml += "<norm_to_back_pixels>%s</norm_to_back_pixels>\n" % str(self.norm_background_roi[1])
        _xml += "<norm_dataset>%s</norm_dataset>\n" % str(self.norm_file)

        # Q cut
        _xml += "<pre_cut>%s</pre_cut>\n" % str(self.pre_cut)
        _xml += "<post_cut>%s</post_cut>\n" % str(self.post_cut)
        _xml += "<q_min>%s</q_min>\n" % str(self.q_min)
        _xml += "<q_step>%s</q_step>\n" % str(self.q_step)
        _xml += "<auto_q_binning>%s</auto_q_binning>\n" % str(self.auto_q_binning)

        # Angle offset
        _xml += "<angle_offset>%s</angle_offset>\n" % str(self.angle_offset)
        _xml += "<angle_offset_error>%s</angle_offset_error>\n" % str(self.angle_offset_error)

        # scaling factor file name
        _xml += "<scaling_factor_flag>%s</scaling_factor_flag>\n" % str(self.scaling_factor_flag)
        _xml += "<scaling_factor_file>%s</scaling_factor_file>\n" % str(self.scaling_factor_file)
        _xml += "<slits_width_flag>%s</slits_width_flag>\n" % str(self.slits_width_flag)

        # Incident medium
        _xml += "<incident_medium_list>%s</incident_medium_list>\n" % str(self.incident_medium_list[0])
        _xml += "<incident_medium_index_selected>%s</incident_medium_index_selected>\n" % str(self.incident_medium_index_selected)

        # Dead time correction
        _xml += "<dead_time_correction>%s</dead_time_correction>\n" % str(self.dead_time)
        _xml += "<dead_time_paralyzable>%s</dead_time_paralyzable>\n" % str(self.paralyzable)
        _xml += "<dead_time_value>%s</dead_time_value>\n" % str(self.dead_time_value)
        _xml += "<dead_time_tof_step>%s</dead_time_tof_step>\n" % str(self.dead_time_tof_step)

        # Emission time correction
        _xml += "<use_emission_time>%s</use_emission_time>\n" % str(self.use_emission_time)
        _xml += "</RefLData>\n"

        return _xml

    def from_xml_element(self, instrument_dom):
        """
            Read in data from XML
            @param xml_str: text to read the data from
        """
        #Peak from/to pixels
        self.data_peak_range = [getIntElement(instrument_dom, "from_peak_pixels"),
                               getIntElement(instrument_dom, "to_peak_pixels")]

        #data metadata
        _tthd_value = getStringElement(instrument_dom, "tthd_value")
        if _tthd_value == '':
            _tthd_value = 'N/A'
        self.tthd_value = _tthd_value

        _ths_value = getStringElement(instrument_dom, "ths_value")
        if _ths_value == '':
            _ths_value = 'N/A'
        self.ths_value = _ths_value

        #low resolution range
        self.data_x_range_flag = getBoolElement(instrument_dom, "x_range_flag",
                                               default=self.data_x_range_flag)

        self.data_x_range = [getIntElement(instrument_dom, "x_min_pixel"),
                             getIntElement(instrument_dom, "x_max_pixel")]

        self.norm_x_range_flag = getBoolElement(instrument_dom, "norm_x_range_flag",
                                                default=self.norm_x_range_flag)

        self.norm_x_range = [getIntElement(instrument_dom, "norm_x_min"),
                             getIntElement(instrument_dom, "norm_x_max")]

        # background flag
        self.subtract_background = getBoolElement(instrument_dom, "background_flag",
                                                 default=self.subtract_background)

        # use two backgrounds flag
        self.two_backgrounds = getBoolElement(instrument_dom, "two_backgrounds",
                                              default=self.two_backgrounds)

        # background from/to pixels
        self.background_roi = [getIntElement(instrument_dom, "back_roi1_from"),
                               getIntElement(instrument_dom, "back_roi1_to"),
                               getIntElement(instrument_dom, "back_roi2_from"),
                               getIntElement(instrument_dom, "back_roi2_to")]

        # TOF range
        self.select_tof_range = getBoolElement(instrument_dom, "tof_range_flag",
                                               default=self.select_tof_range)
        self.tof_range = [getFloatElement(instrument_dom, "from_tof_range"),
                          getFloatElement(instrument_dom, "to_tof_range")]

        self.data_files = getIntList(instrument_dom, "data_sets")

        #with or without norm
        self.apply_normalization = getBoolElement(instrument_dom, "norm_flag",
                                                  default=self.apply_normalization)

        #Peak from/to pixels
        self.norm_peak_range = [getIntElement(instrument_dom, "norm_from_peak_pixels"),
                                getIntElement(instrument_dom, "norm_to_peak_pixels")]

        # Background subtraction option
        self.subtract_norm_background = getBoolElement(instrument_dom, "norm_background_flag",
                                                       default=self.subtract_norm_background)

        self.norm_background_roi = [getIntElement(instrument_dom, "norm_from_back_pixels"),
                                    getIntElement(instrument_dom, "norm_to_back_pixels")]

        self.norm_file = getIntElement(instrument_dom, "norm_dataset")

        # Q cut
        self.pre_cut = getIntElement(instrument_dom, "pre_cut", default=self.pre_cut)
        self.post_cut = getIntElement(instrument_dom, "post_cut", default=self.post_cut)
        self.q_min = getFloatElement(instrument_dom, "q_min", default=self.q_min)
        self.q_step = getFloatElement(instrument_dom, "q_step", default=self.q_step)
        self.auto_q_binning = getBoolElement(instrument_dom, "auto_q_binning", default=False)

        # Angle offset
        self.angle_offset = getFloatElement(instrument_dom, "angle_offset", default=self.angle_offset)
        self.angle_offset_error = getFloatElement(instrument_dom, "angle_offset_error",
                                                  default=self.angle_offset_error)

        # Scaling factor file and options
        self.scaling_factor_file = getStringElement(instrument_dom, "scaling_factor_file")
        self.slits_width_flag = getBoolElement(instrument_dom, "slits_width_flag")
        self.scaling_factor_flag = getBoolElement(instrument_dom, "scaling_factor_flag")

        # Incident medium selected
        if getStringList(instrument_dom, "incident_medium_list") != []:
            self.incident_medium_list = getStringList(instrument_dom, "incident_medium_list")
            self.incident_medium_index_selected = getIntElement(instrument_dom, "incident_medium_index_selected")
        else:
            self.incident_medium_list = ['H2O']
            self.incident_medium_index_selected = 0

        # Dead time correction
        self.dead_time = getBoolElement(instrument_dom, "dead_time_correction",
                                        default=self.dead_time)
        self.paralyzable = getBoolElement(instrument_dom, "dead_time_paralyzable",
                                          default=self.paralyzable)
        self.dead_time_value = getFloatElement(instrument_dom, "dead_time_value",
                                               default=self.dead_time_value)
        self.dead_time_tof_step = getFloatElement(instrument_dom, "dead_time_tof_step",
                                                  default=self.dead_time_tof_step)

        # Emission time
        # Defaults to True, but will be skipped if the necessary meta data is not found
        self.use_emission_time = getBoolElement(instrument_dom, "use_emission_time",
                                                default=True)


###### Utility functions to read XML content ########################
def getText(nodelist):
    """
        Utility method to extract text out of an XML node
    """
    rc = ""
    for node in nodelist:
        if node.nodeType == node.TEXT_NODE:
            rc = rc + node.data
    return rc

def getContent(dom, tag):
    element_list = dom.getElementsByTagName(tag)
    return getText(element_list[0].childNodes) if len(element_list) > 0 else None

def getIntElement(dom, tag, default=None):
    value = getContent(dom, tag)
    return int(value) if value is not None else default

def getIntList(dom, tag, default=[]):
    value = getContent(dom, tag)
    if value is not None and len(value.strip()) > 0:
        return list(map(int, value.split(',')))
    else:
        return default

def getFloatElement(dom, tag, default=None):
    value = getContent(dom, tag)
    return float(value) if value is not None else default

def getFloatList(dom, tag, default=[]):
    value = getContent(dom, tag)
    if value is not None and len(value.strip()) > 0:
        return list(map(float, value.split(',')))
    else:
        return default

def getStringElement(dom, tag, default=''):
    value = getContent(dom, tag)
    return value if value is not None else default

def getStringList(dom, tag, _default=[]):
    elem_list = []
    element_list = dom.getElementsByTagName(tag)
    if len(element_list) > 0:
        for l in element_list:
            elem_list.append(getText(l.childNodes).strip())
    return elem_list

def getBoolElement(dom, tag, true_tag='true', default=False):
    value = getContent(dom, tag)
    return value.lower() == true_tag.lower() if value is not None else default


###### Functions to read/write a template file ######################
def to_xml(data_sets):
    """
        Create XML from the current data.
    """
    _xml = "<Reduction>\n"
    _xml += "    <instrument_name>REFL</instrument_name>\n"
    _xml += "    <timestamp>%s</timestamp>\n" % time.ctime()
    _xml += "    <version>%s</version>\n" % VERSION
    _xml += "    <mantid_version>%s</mantid_version>\n" % MANTID_VERSION
    _xml += "    <generator>lr_reduction-%s</generator>\n" % VERSION
    _xml  += "<DataSeries>\n"
    for item in data_sets:
        _xml += item.to_xml()
    _xml += "</DataSeries>\n"
    _xml += "</Reduction>\n"
    return _xml

def from_xml(xml_str):
    """
        Read in data from XML string
    """
    data_sets = []
    dom = xml.dom.minidom.parseString(xml_str)

    element_list = dom.getElementsByTagName("Data")
    if len(element_list)==0:
        element_list = dom.getElementsByTagName("RefLData")

    if len(element_list)>0:
        for item in element_list:
            if item is not None:
                data_set = ReductionParameters()
                data_set.from_xml_element(item)
                data_sets.append(data_set)

    if len(data_sets) == 0:
        data_sets = [ReductionParameters()]

    return data_sets
