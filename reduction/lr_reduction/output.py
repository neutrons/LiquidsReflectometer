"""
    Write R(q) output
"""
import json
import numpy as np


class RunCollection():
    """
        A collection of runs to assemble into a single R(Q)
    """
    def __init__(self, average_overlap=False):
        self.collection = []
        self.average_overlap = average_overlap
        self.qz_all = []
        self.refl_all = []
        self.d_refl_all = []
        self.d_qz_all = []

    def add(self, q, r, dr, meta_data, dq=None):
        if dq is None:
            resolution = meta_data['dq_over_q']
            dq = resolution * q
        self.collection.append(dict(q=q, r=r, dr=dr, dq=dq, info=meta_data))

    def merge(self):
        """
            Merge the collection of runs
        """
        qz_all = []
        refl_all = []
        d_refl_all = []
        d_qz_all = []

        for item in self.collection:
            
            idx = np.fabs(item['r']) > 0
            qz_mid = item['q'][idx]
            refl = item['r'][idx]
            d_refl = item['dr'][idx]
            d_qz = item['dq'][idx]

            for i in range(len(qz_mid)-1, -1, -1):
                qz_all.append(qz_mid[i])
                refl_all.append(refl[i])
                d_refl_all.append(d_refl[i])
                d_qz_all.append(d_qz[i])

        qz_all = np.asarray(qz_all)
        refl_all = np.asarray(refl_all)
        d_refl_all = np.asarray(d_refl_all)
        d_qz_all = np.asarray(d_qz_all)
        idx = np.argsort(qz_all)

        self.qz_all = np.take_along_axis(qz_all, idx, axis=None)
        self.refl_all = np.take_along_axis(refl_all, idx, axis=None)
        self.d_refl_all = np.take_along_axis(d_refl_all, idx, axis=None)
        self.d_qz_all = np.take_along_axis(d_qz_all, idx, axis=None)

    def save_ascii(self, file_path, meta_as_json=False):
        """
            Save R(Q) in ascii format
        """
        self.merge()

        with open(file_path, 'w') as fd:
            # Write meta data
            initial_entry_written = False
            for item in self.collection:
                _meta = item['info']
                if not initial_entry_written:
                    fd.write("# Experiment %s Run %s\n" % (_meta['experiment'], _meta['run_number']))
                    fd.write("# Run title: %s\n" % _meta['run_title'])
                    fd.write("# Run start time: %s\n" % _meta['start_time'])
                    fd.write("# Reduction time: %s\n" % _meta['time'])
                    if meta_as_json:
                        for k in _meta.keys():
                            if type(_meta[k]) == np.int32:
                                print(k)
                                
                        fd.write("# Meta:%s\n" % json.dumps(_meta))
                    fd.write("# DataRun   NormRun   TwoTheta(deg)  LambdaMin(A)   ")
                    fd.write("LambdaMax(A) Qmin(1/A)    Qmax(1/A)    SF_A         SF_B\n")
                    fd.write("")
                if 'sf' in _meta:
                    a = _meta['sf']['a']
                    b = _meta['sf']['b']
                else:
                    a = 1
                    b = 0
                value_list = (_meta['run_number'], _meta['norm_run'], _meta['theta']*2.0,
                              _meta['wl_min'], _meta['wl_max'], _meta['q_min'], _meta['q_max'],
                              a, b)
                fd.write("# %-9s %-9s %-14.6g %-14.6g %-12.6g %-12.6s %-12.6s %-12.6s %-12.6s\n" % value_list)

                initial_entry_written = True

            # Write R(q)
            fd.write('# dQ0[1/Angstrom] = %g\n' % _meta['dq0'])
            fd.write('# dQ/Q = %g\n' % _meta['dq_over_q'])
            fd.write('# %-21s %-21s %-21s %-21s\n' % ('Q [1/Angstrom]', 'R', 'dR', 'dQ [FWHM]'))
            for i in range(len(self.qz_all)):
                fd.write('%20.16f  %20.16f  %20.16f  %20.16f\n' % (self.qz_all[i], self.refl_all[i], self.d_refl_all[i], self.d_qz_all[i]))

    def add_from_file(self, file_path):
        _q, _r, _dr, _dq, _meta = read_file(file_path)
        self.add(_q, _r, _dr, _meta, dq=_dq)


def read_file(file_path):
    """
        Read a data file and extract meta data
    """
    _meta = dict()
    with open(file_path, 'r') as fd:
        for l in fd.readlines():
            if l.startswith("# Meta:"):
                _meta = json.loads(l[len("# Meta:"):-1])
    _q, _r, _dr, _dq = np.loadtxt(file_path).T
    return _q, _r, _dr, _dq, _meta
