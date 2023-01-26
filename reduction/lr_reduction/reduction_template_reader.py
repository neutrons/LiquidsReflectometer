"""
    RefRed template reader.
    Adapted from Mantid code.
"""
import xml.dom.minidom


class ReductionParameters(object):

    def __init__(self):
        # Signal selection
        self.data_peak_range = [140, 150]
        self.subtract_background = True
        self.background_roi = [137, 153,100, 200]
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
        self.scaling_factor_file_flag = True
        self.slits_width_flag = True

        # Incident medium list and selected value
        self.incident_medium_list = ['air']
        self.incident_medium_index_selected = 0

        # 4th column of ASCII file (q resolution)
        self.fourth_column_dq0 = 0.000
        self.fourth_column_dq_over_q = 0.027

        # How to treat overlap values
        self.overlap_lowest_error = True
        self.overlap_mean_value = False

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
        _xml += "<q_min>%s</q_min>\n" % str(self.q_min)
        _xml += "<q_step>%s</q_step>\n" % str(self.q_step)
        _xml += "<auto_q_binning>%s</auto_q_binning>\n" % str(self.auto_q_binning)
        _xml += "<overlap_lowest_error>%s</overlap_lowest_error>\n" % str(self.overlap_lowest_error)
        _xml += "<overlap_mean_value>%s</overlap_mean_value>\n" % str(self.overlap_mean_value)

        # Angle offset
        _xml += "<angle_offset>%s</angle_offset>\n" % str(self.angle_offset)
        _xml += "<angle_offset_error>%s</angle_offset_error>\n" % str(self.angle_offset_error)

        # scaling factor file name
        _xml += "<scaling_factor_flag>%s</scaling_factor_flag>\n" % str(self.scaling_factor_file_flag)
        _xml += "<scaling_factor_file>%s</scaling_factor_file>\n" % str(self.scaling_factor_file)
        _xml += "<slits_width_flag>%s</slits_width_flag>\n" % str(self.slits_width_flag)

        # Incident medium
        _xml += "<incident_medium_list>%s</incident_medium_list>\n" % str(self.incident_medium_list[0])
        _xml += "<incident_medium_index_selected>%s</incident_medium_index_selected>\n" % str(self.incident_medium_index_selected)

        # Fourth column (q resolution)
        _xml += "<fourth_column_flag>True</fourth_column_flag>\n"
        _xml += "<fourth_column_dq0>%s</fourth_column_dq0>\n" % str(self.fourth_column_dq0)
        _xml += "<fourth_column_dq_over_q>%s</fourth_column_dq_over_q>\n" % str(self.fourth_column_dq_over_q)

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

        #background flag
        self.subtract_background = getBoolElement(instrument_dom, "background_flag",
                                                 default=self.subtract_background)

        #background from/to pixels
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
        self.q_min = getFloatElement(instrument_dom, "q_min", default=self.q_min)
        self.q_step = getFloatElement(instrument_dom, "q_step", default=self.q_step)
        self.auto_q_binning = getBoolElement(instrument_dom, "auto_q_binning", default=False)

        # overlap_lowest_error
        self.overlap_lowest_error = getBoolElement(instrument_dom, "overlap_lowest_error", default=True)
        self.overlap_mean_value = getBoolElement(instrument_dom, "overlap_mean_value", default=False)

        # Angle offset
        self.angle_offset = getFloatElement(instrument_dom, "angle_offset", default=self.angle_offset)
        self.angle_offset_error = getFloatElement(instrument_dom, "angle_offset_error",
                                                  default=self.angle_offset_error)

        # Scaling factor file and options
        self.scaling_factor_file = getStringElement(instrument_dom, "scaling_factor_file")
        self.slits_width_flag = getBoolElement(instrument_dom, "slits_width_flag")
        self.scaling_factor_file_flag = getBoolElement(instrument_dom, "scaling_factor_flag")

        # Incident medium selected
        if getStringList(instrument_dom, "incident_medium_list") != []:
            self.incident_medium_list = getStringList(instrument_dom, "incident_medium_list")
            self.incident_medium_index_selected = getIntElement(instrument_dom, "incident_medium_index_selected")
        else:
            self.incident_medium_list = ['H2O']
            self.incident_medium_index_selected = 0

        # Fourth column (q resolution)
        self.fourth_column_dq0 = getFloatElement(instrument_dom, "fourth_column_dq0")
        self.fourth_column_dq_over_q = getFloatElement(instrument_dom, "fourth_column_dq_over_q")


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
    _xml  = "<DataSeries>\n"
    for item in data_sets:
        _xml += item.to_xml()
    _xml += "</DataSeries>\n"
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

    if len(data_sets)==0:
        data_sets = [ReductionParameters()]

    return data_sets