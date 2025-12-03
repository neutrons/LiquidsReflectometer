# standard imports
import os

# third-party imports
import numpy as np

# lr_reduction imports
from lr_reduction import time_resolved
from lr_reduction.utils import amend_config


def test_reduce_workflow(nexus_dir):
    """
        Test the time-resolved reduction that uses a measured reference.
        It is generally used at 30 Hz but it also works at 60 Hz.
    """
    template_path = 'data/template.xml'
    output_dir = 'data/'
    reduced_path = 'data/reference_rq_avg_overlap.txt'
    ref_data = np.loadtxt(reduced_path).T
    with amend_config(data_dir=nexus_dir):
        reduced = time_resolved.reduce_30Hz_slices(198413, 198413, ref_data_60Hz=reduced_path,
                                                   template_30Hz=template_path,
                                                   time_interval=300, output_dir=output_dir,
                                                   scan_index=5, create_plot=False)

    q_long = len(ref_data[0])
    q_short = len(reduced[0][0])
    n_match = 0
    n_pts = 0
    for i in range(q_long):
        if ref_data[0][i] > 0.03 and ref_data[0][i] < 0.047:
            n_pts += 1
            for k in range(q_short):
                if np.fabs(reduced[0][0][k] - ref_data[0][i]) < 0.0001:
                    assert(np.fabs(reduced[0][1][k] - ref_data[1][i]) < 1e-10)
                    n_match += 1
    assert(n_pts == n_match)

    # Plot data
    time_resolved.plot_slices(reduced, 'Test', 300,
                              os.path.join(output_dir, 'reduced.png'), show=False)


def test_reduce_template_workflow(nexus_dir):
    """
        Test the time-resolved reduction that uses a template.
    """
    template_path = 'data/template.xml'
    output_dir = 'data/'
    reduced_path = 'data/reference_rq_avg_overlap.txt'
    ref_data = np.loadtxt(reduced_path).T
    with amend_config(data_dir=nexus_dir):
        reduced = time_resolved.reduce_slices(198413,
                                              template_file=template_path,
                                              time_interval=300,
                                              output_dir=output_dir,
                                              scan_index=5,
                                              create_plot=False)[0]  # there's only one element in the returned list
    q_values = reduced[0]
    intensities = reduced[1]

    n_match = 0
    n_pts = 0
    for i, q_ref in enumerate(ref_data[0]):
        if 0.03 < q_ref < 0.045:
            n_pts += 1
            for k, q_red in enumerate(q_values):
                if np.fabs(q_ref - q_red) < 0.0001:
                    relative_difference = np.fabs(intensities[k] - ref_data[1][i]) / ref_data[1][i]
                    assert(relative_difference < 0.07)
                    n_match += 1
    assert(n_pts == n_match)

    # check data can be plotted to a PNG image, but don't show it on the screen
    time_resolved.plot_slices([reduced], 'Test', 300,
                              os.path.join(output_dir, 'reduced_template.png'), show=False)
