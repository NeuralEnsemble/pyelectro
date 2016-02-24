"""
IO operations, these are mostly specific to in-house formats used by Hugh Robinson's lab at Cambridge.
"""

import scipy.io
import numpy as np

from pyelectro.analysis import print_comment_v

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

def load_csv_data(file_path,plot=False):
    """Extracts time and voltage data from a csv file
    
    Data must be in a csv and in two columns, first time and second 
    voltage. Units should be SI (Volts and Seconds).

    :param file_path: full file path to file e.g /home/mike/test.csv
        
    :return: two lists - time and voltage

    """
    import csv

    csv_file=file(file_path,'r')
    csv_reader=csv.reader(csv_file)

    v=[]
    t=[]

    i=0
    for row in csv_reader:

        try:

            t_value=float(row[0])*1000 #convert to ms
            v_value=float(row[1])*1000 #convert to mV

            t.append(t_value)
            v.append(v_value)

        except:
            print_comment_v('row ',i,' invalid')

        i+=1

    if plot:
        from matplotlib import pyplot        
        pyplot.plot(t,v)
        pyplot.title('Raw data')
        pyplot.xlabel('Time (ms)')
        pyplot.ylabel('Voltage (mV)')
        pyplot.show()

    return np.array(t), np.array(v)

def summary(data, label=""):
    if len(label)> 0: label = " (%s)"%label
    info = "Data%s: "%label
    if len(data)>4:
        info+="%i points: [%f, %f, ..., %f]"%(len(data), data[0], data[1], data[-1])
    elif len(data)>0:
        info+="%i points: %s"%(len(data), data)
    else:
        info+="Empty data set"
        
    return info