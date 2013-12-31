"""
IO operations, these are mostly specific to in-house formats used by Hugh Robinson's lab at Cambridge.
"""

import scipy.io
import numpy as np

def data_from_sweep(path):
    """Extract the trace from a MATLAB sweep file
    
    :param path: the full path to the sweep file

    :return: raw trace

    
    """
    mat = scipy.io.loadmat(path)
    recording = mat['sweep_struct']['data']
    recording = np.array(recording)
    recording = recording[0][0][0]
    
    return recording
