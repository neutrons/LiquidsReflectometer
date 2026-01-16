.. _workflow:

Specular reflectivity reduction workflow
========================================

The specular reflectivity data reduction is build around the ``event_reduction.EventReflectivity`` class, which
performs the reduction. A number of useful modules are available to handle parts of the workflow around the actual reduction.

Data sets
---------

Specular reflectivity measurements at BL4B are done by combining several runs, taken at different
scattering angle and wavelength band. To allow for the automation of the reduction process, several
meta data entries are stored in the data files. To be able to know which data files belong together
in a single reflectivity measurement, two important log entries are used:

- ``sequence_id``: The sequence ID identifies a unique reflectivity curve. All data runs with a matching
  sequence_id are put together to create a single reflectivity curve.
- ``sequence_number``: The sequence number identifies the location of a given run in the list of runs that define a full sequence. All sequences start at 1. For instance, a sequence number of 3 means that
  this run is the third of the complete set. This becomes important for storing reduction parameters.

Reduction parameters and templates
----------------------------------

The reduction parameters are managed using the ``reduction_template_reader.ReductionParameters`` class. This class allows users to define and store the parameters required for the reduction process. By using this class, you can easily save, load, and modify the parameters, ensuring consistency and reproducibility in your data reduction workflow.

Compatibility with RefRed
.........................

[Refred](https://github.com/neutrons/RefRed) is the user interface that helps users define reduction
parameters by selecting the data to process, peak and background regions, etc. A complete reflectivity
curve is generally comprised of multiple runs, and RefRed allows one to save a so-called template file
that contains all the information needs to reduce each run in the set. The reduction backend (this package) has utilities to read and write such templates, which are stored in XML format. A template
consist of an ordered list of ``ReductionParameters`` objects, which corresponding to a specific ``sequence_number``.

To read a templates and obtains a list of ``ReductionParameters`` objects:

.. code-block:: python

    from lr_reduction import reduction_template_reader

    with open(template_file, "r") as fd:
        xml_str = fd.read()
        data_sets = reduction_template_reader.from_xml(xml_str)

From a list of ``ReductionParameters`` objects:

.. code-block:: python

    xml_str = reduction_template_reader.to_xml(data_sets)
    with open(os.path.join(output_dir, "template.xml"), "w") as fd:
        fd.write(xml_str)

Reduction workflow
------------------

The main reduction workflow, which will extract specular reflectivity from a data file given a
reduction template, is found in the `workflow` module. This workflow will is the one performed
by the automated reduction system at BL4B:

- It will extract the correct reduction parameters from the provided template
- Perform the reduction and compute the reflectivity curve for that data
- Combine the reflectivity curve segment with other runs belonging to the same set
- Write out the complete reflectivity curve in an output file
- Write out a copy of the template by replacing the run numbers in the template by those that were used

Once you have a template, you can simply do:

.. code-block:: python

    from lr_reduction import workflow
    from mantid.simpleapi import LoadEventNexus

    # Load the data from disk
    ws = LoadEventNexus(Filename='/SNS/REF_L/IPTS-XXXX/nexus/REFL_YYYY.h5')

    # The template file you want to use
    template_file = '/SNS/REF_L/IPTS-XXXX/autoreduce/template.xml'

    # The folder where you want your output
    output_dir = '/tmp'

    workflow.reduce(ws, template_file, output_dir)

This will produce output files in the specified output directory.
