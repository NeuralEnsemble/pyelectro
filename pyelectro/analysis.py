# -*- coding: utf-8 -*-
"""
Module for mathematical analysis of voltage traces from electrophysiology.

AUTHOR: Mike Vella vellamike@gmail.com

"""

import scipy.stats
import numpy as np
import math
import logging
import sys
from scipy import interpolate
import operator

import pprint
    
pp = pprint.PrettyPrinter(indent=4)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def print_comment_v(text, warning=False):
    print_comment(text, True, warning)
    
    
def print_comment(text, print_it=False, warning=False):
    
    prefix = "pyelectro >>> "
    if warning:
        prefix += "WARNING "
    if not isinstance(text, str): text = text.decode('ascii')
    if print_it:
        
        print("%s%s"%(prefix, text.replace("\n", "\n"+prefix)))
        

def voltage_plot(t,v,title=None):
    """
    Plot electrophysiology recording.
    """

    from matplotlib import pyplot as plt
    plt.xlabel('Time (ms)')
    plt.ylabel('Voltage (mV)')
    plt.title(title)
    plt.grid()
    plt.plot(t,v)
    plt.show()


def smooth(x,window_len=11,window='hanning'):
    """Smooth the data using a window with requested size.

    This function is useful for smoothing out experimental data.
    This method utilises the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    :param x: the input signal 
    :param window_len: the dimension of the smoothing window; should be an odd integer
    :param window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman', flat window will produce a moving average smoothing.

    :return: smoothed signal

    example:

    .. code-block:: python

       t=linspace(-2,2,0.1)
       x=sin(t)+randn(len(t))*0.1
       y=smooth(x)

    .. seealso::

       numpy.hanning
       numpy.hamming
       numpy.bartlett
       numpy.blackman
       numpy.convolve
       scipy.signal.lfilter
    """

    if x.ndim != 1:
        raise(ValueError, "smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise(ValueError, "Input vector needs to be bigger than window size.")

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise(ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')

    edge=int(window_len/2)
    return y[edge:-edge]

def linear_fit(t, y):
    """ Fits data to a line

    :param t: time vector
    :param y: variable which varies with time (such as voltage)
    :returns: Gradient M for a formula of the type y=C+M*x
    """

    vals=np.array(y)
    m,C = np.polyfit(t, vals, 1)
    return m


def three_spike_adaptation(t,y):
    """ Linear fit of amplitude vs time of first three AP spikes

    Initial action potential amplitudes may very substaintially in amplitude
    and then settle down.

    :param t: time vector (AP times)
    :param y: corresponding AP amplitude
    :returns: Gradient M for a formula of the type y=C+M*x for first three action potentials
    """

    t = np.array(t)
    y = np.array(y)

    t = t[0:3]
    y = y[0:3]

    m = linear_fit(t,y)

    return m


def exp_fit(t, y):
    """
    Fits data to an exponential.

        Returns K for a formula of the type y=A*exp(K*x)

        :param t: time vector
        :param y: variable which varies with time (such as voltage)

    """

    vals = np.array(y)
    C = np.min(vals)
    vals = vals-C+1e-9 #make sure the data is all positive
    vals = np.log(vals)
    K, A_log = np.polyfit(t, vals, 1)

    return K


def window_peak_detector(v,delta = 0.01):
    """
    Detects peak by comparing mean of either side of
    peak and deciding whether it exceeds some threshold.

    :return: Boolean, True if a peak is detected in that window
    """

    if len(v) % 2 == 0:
        raise Exception("Window length must be odd")

    middle_index = len(v) // 2
    middle_value = v[middle_index]

    left_mean = np.mean(v[0:middle_index])
    right_mean = np.mean(v[middle_index+1:])

    left_elevation = middle_value - left_mean
    right_elevation = middle_value - right_mean

    left_exceeds_threhold = left_elevation > delta
    right_exceeds_threshold = right_elevation > delta

    return left_exceeds_threhold and right_exceeds_threshold

def centered_slice(v, index, length=5):
    """
    Retruns slice of given length centred on index.
    """

    if length % 2 == 0:
        raise Exception("Window length must be odd")

    if len(v) < index + length // 2:
        raise Exception("Index too close to edge or window too big")

    start_index = index - length // 2
    slice = v[start_index:start_index + length]

    return slice


def max_min_simple(a,
                   times,
                   delta = 0,
                   peak_threshold = 0.0,
                   verbose = False):
    
    print_comment("Calculating max_min_simple of a: (%s,...,%s)#%i, t: (%s,...,%s)#%i; thresh %s, delta %s"%(a[0],a[-1],len(a),times[0],times[-1],len(times), peak_threshold, delta), verbose)
    
    maxima_locations = []
    maxima_number = 0
    maxima_times = []
    maxima_values = []

    minima_locations = []
    minima_number = 0
    minima_times = []
    minima_values = []
    
    spiking  = False
    has_spiked  = False
    
    last_max_loc = -1
    last_max_t = -1
    last_max_v = -1*sys.float_info.max
    
    last_min_loc = -1
    last_min_t = -1
    last_min_v = sys.float_info.max
    
    for i in range(len(a)):
        t = times[i]
        v = a[i]
        
        if not spiking and v>=peak_threshold:
            print_comment('Spike of %s at %s'%(v,t),verbose)
            spiking = True
            has_spiked = True
            if last_min_loc >0:
                minima_locations.append(last_min_loc)
                minima_times.append(last_min_t)
                minima_values.append(last_min_v)
                minima_number+=1

            last_min_loc = -1
            last_min_t = -1
            last_min_v = sys.float_info.max
                
        elif spiking and v<peak_threshold:
            
            spiking = False
            if last_max_loc >0:
                maxima_locations.append(last_max_loc)
                maxima_times.append(last_max_t)
                maxima_values.append(last_max_v)
                maxima_number+=1
            
            last_max_loc = -1
            last_max_t = -1
            last_max_v = -1*sys.float_info.max
            
        if spiking:
            
            if v >= last_max_v:
                last_max_loc = i
                last_max_t = t
                last_max_v = v
                
        elif has_spiked:
        
            if v <= last_min_v:
                last_min_loc = i
                last_min_t = t
                last_min_v = v
      
    #need to construct the dictionary here:
    turning_points = {'maxima_locations':maxima_locations,
                      'minima_locations':minima_locations,
                      'maxima_number':maxima_number,
                      'minima_number':minima_number,
                      'maxima_times':maxima_times,
                      'minima_times':minima_times, 
                      'maxima_values':maxima_values,
                      'minima_values':minima_values}
                      

    return turning_points

def max_min(a, 
            t,
            delta=0,
            peak_threshold=0.0,
            verbose = False):
    """
    Find the maxima and minima of a voltage trace.

    :note This method does not appear to be very robust when comparing to experimental data

    :param a: time-dependent variable (usually voltage)
    :param t: time-vector
    :param delta: the value by which a peak or trough has to exceed its
        neighbours to be considered outside of the noise
    :param peak_threshold: peaks below this value are discarded

    :return: turning_points, dictionary containing number of max, min and 
        their locations

    .. note::

       minimum value between two peaks is in some ways a better way
       of obtaining a minimum since it guarantees an answer, this may be
       something which should be implemented.

    """
    if peak_threshold==None: 
        import sys
        peak_threshold = -1*sys.float_info.max
    print_comment("Calculating max_min of a: (%s,...,%s)#%i, t: (%s,...,%s)#%i; thresh %s, delta %s"%(a[0],a[-1],len(a),t[0],t[-1],len(t), peak_threshold, delta),verbose)
    gradients = np.diff(a)

    maxima_info = []
    minima_info = []

    count = 0

    for i in gradients[:-1]:
        count+=1

        if i > 0 and gradients[count] < 0 and i != gradients[count]:
            #found a maximum
            maximum_value=a[count]
            maximum_location=count
            maximum_time=t[count]
            preceding_point_value=a[maximum_location-1]
            succeeding_point_value=a[maximum_location+1]

            #filter:
            maximum_valid=False #logically consistent but not very pythonic..
            if ((maximum_value-preceding_point_value)>delta)*((maximum_value-succeeding_point_value)>delta):
                maximum_valid=True
            if maximum_value<peak_threshold:
                maximum_valid=False
            if maximum_valid:
                maxima_info.append((maximum_value,maximum_location,maximum_time))

    maxima_num=len(maxima_info)

    if maxima_num>0:
        minima_num=maxima_num-1
    else:
        minima_num=0



    values_getter=operator.itemgetter(0)
    location_getter=operator.itemgetter(1)
    time_getter=operator.itemgetter(2)

    maxima_locations=list(map(location_getter,maxima_info))
    maxima_times=list(map(time_getter,maxima_info))
    maxima_values=list(map(values_getter,maxima_info))

    for i in range(maxima_num-1):
        maximum_0_location=maxima_locations[i]
        maximum_1_location=maxima_locations[i+1]

        interspike_slice=a[maximum_0_location:maximum_1_location]
        minimum_value=min(interspike_slice)
        minimum_location=list(interspike_slice).index(minimum_value)+maximum_0_location
        minimum_time=t[minimum_location]

        minima_info.append((minimum_value,minimum_location,minimum_time))

    minima_locations=list(map(location_getter,minima_info))
    minima_times=list(map(time_getter,minima_info))
    minima_values=list(map(values_getter,minima_info))

    #need to construct the dictionary here:
    turning_points = {'maxima_locations':maxima_locations,'minima_locations':minima_locations,'maxima_number':maxima_num,'minima_number':minima_num,'maxima_times':maxima_times,'minima_times':minima_times, 'maxima_values':maxima_values,'minima_values':minima_values}

    return turning_points

'''  PG removing this...
def max_min2(v,t,delta=0.1,peak_threshold=0.0,window_length=11):
    """
    Uses the max_min function but then does a second pass with
    window peak detector to discard peaks.

    This is being prepared as an enhancement to the old
    peak detector.
    """

    max_min_dict = max_min(v,t,delta=0.0,peak_threshold=peak_threshold)

    maxima_locations = max_min_dict['maxima_locations']

    peak_mask = []

    for location in maxima_locations:
        slice = centered_slice(v,location,window_length)
        peak_flag = window_peak_detector(slice, delta=delta)
        peak_mask.append(peak_flag)

    #this anonymous function strips a list of all corresponding
    #non-zero elements in the mask:
    print("peak_mask: "+peak_mask)

    mask_filter = lambda l, mask : list(itertools.compress(l,mask))

    max_min_dict.pop('maxima_number',None)
    max_min_dict.pop('minima_number',None)    

    dict_keys = max_min_dict.keys()

    for key in dict_keys:
        max_min_dict[key] = mask_filter(max_min_dict[key],peak_mask)

    max_min_dict['maxima_number'] = len(max_min_dict['maxima_locations'])
    max_min_dict['minima_number'] = max_min_dict['maxima_number'] - 1

    return max_min_dict'''

def spike_frequencies(t):
    """
    Calculate frequencies associated with interspike times

    :param t: a list of spike times in ms

    :return: list of frequencies in Hz associated with interspike times and
        times associated with the frequency (time of first spike in pair)

    """
    spike_times=np.array(t)
    interspike_times=np.diff(spike_times)
    interspike_frequencies=1000/interspike_times

    return [t[:-1],interspike_frequencies]


def max_min_interspike_time(t):
    """
    Calculate the maximum & minimum interspike interval from the list of maxima times
    
    :param t: a list of spike times in ms
    
    :return: (max, min) interspike time
    """
    
    spike_times=np.array(t)
    interspike_times=np.diff(spike_times)
    return max(interspike_times), min(interspike_times)
    


def mean_spike_frequency(t):
    """
    Find the average frequency of spikes

    :param t: a list of spike times in ms

    :return: mean spike frequency in Hz, calculated from mean interspike time

    """
    interspike_times=np.diff(t)
    mean_interspike_time=np.mean(interspike_times)
    mean_frequency=1000.0/(mean_interspike_time) #factor of 1000 to give frequency in Hz

    if (math.isnan(mean_frequency)):
        mean_frequency=0
    return mean_frequency


def y_from_x(y,x,y_to_find):
    """
    Returns list of x values corresponding to a y after a doing a 
    univariate spline interpolation

    :param x: x-axis numerical data
    :param y: corresponding y-axis numerical data
    :param y_to_find: x value for desired y-value,
        interpolated from nearest two measured x/y value pairs

    :return: interpolated y value

    """

    #TODO:should have the ability to return indices, this should be a flag

    yreduced = np.array(y) - y_to_find
    freduced = interpolate.UnivariateSpline(x, yreduced, s=None)

    return freduced.roots()


def single_spike_width(y,t,baseline):
    """ Find the width of a spike at a fixed height

    calculates the width of the spike at height baseline. If the spike shape
    does not intersect the height at both sides of the peak the method
    will return value 0. If the peak is below the baseline 0 will also 
    be returned.

    The input must be a single spike or nonsense may be returned.
    Multiple-spike data can be handled by the interspike_widths method.

    :param y: voltage trace (array) corresponding to the spike
    :param t: time value array corresponding to y
    :param baseline: the height (voltage) where the width is to be measured.        

    :return: width of spike at height defined by baseline

    """

    logger.debug('Baseline: %f' %baseline)

    try:
        y = np.array(y)
        t = np.array(t)

        value = np.max(y)
        location = np.argmax(y)

        logger.debug('Max voltage: %f' %value)
        logger.debug('Index of max: %f' %location)

        #moving left:
        while value > baseline:
            location -= 1
            value = y[location]
            undershoot_value = y[location + 1]
            overshoot_time = t[location]
            undershoot_time = t[location + 1]
            interpolated_left_time = np.interp(baseline, [value, undershoot_value], [overshoot_time, undershoot_time])

            if location < 0:
                raise MathsError('Baseline does not intersect spike')

        #now go right
        value = np.max(y)
        location = np.argmax(y)

        while value > baseline :
            location += 1
            value = y[location]
            undershoot_value = y[location - 1]
            overshoot_time = t[location]
            undershoot_time = t[location - 1]
            interpolated_right_time = np.interp(baseline, [value, undershoot_value], [overshoot_time, undershoot_time])

            if location > len(y) - 1:
                raise MathsError('Baseline does not intersect spike')

        width = interpolated_right_time - interpolated_left_time

    except:
        logger.warning('Single spike width algorithm failure - setting to 0')
        width = 0.0

    return width


def spike_widths(y,t,max_min_dictionary,baseline=0,delta=0):
    """
    Find the widths of each spike at a fixed height in a train of spikes.

    Returns the width of the spike of each spike in a spike train at height 
    baseline. If the spike shapes do not intersect the height at both sides
    of the peak the method will return value 0 for that spike.
    If the peak is below the baseline 0 will also be returned for that spike.

    :param y: voltage trace (array) corresponding to the spike train
    :param t: time value array corresponding to y
    :param max_min_dictionary: precalculated max_min_dictionary
    :param baseline: the height (voltage) where the width is to be measured.

    :return: width of spike at height defined by baseline

    """

    max_num=max_min_dictionary['maxima_number']
    maxima_times=max_min_dictionary['maxima_times']
    minima_locations=max_min_dictionary['minima_locations']

    spike_widths=[]
    for i in range(max_num):
        #need to splice down the y:
        if i==0:
            left_min_location=0
            right_min_location=minima_locations[i]+1
        elif i==max_num-1:
            left_min_location=minima_locations[i-1]
            right_min_location=len(y)
        else:
            left_min_location=minima_locations[i-1]
            right_min_location=minima_locations[i]+1

        spike_shape=y[left_min_location:right_min_location]
        spike_t=t[left_min_location:right_min_location]

        try:
            width=single_spike_width(spike_shape,spike_t,baseline)
            logger.debug('Spike width: %f' %width)
        except:
            logger.warning('Spike width set to 0, this indicates a problem')
            width=0

        spike_widths.append(width)

    maxima_times_widths=[maxima_times,spike_widths]
    return maxima_times_widths

def burst_analyser(t):
    """ Pearson's correlation coefficient applied to interspike times

    :param t: Rank-1 array containing spike times

    :return: pearson's correlation coefficient of interspike times 
    """

    x=np.arange(len(t))
    pearsonr=scipy.stats.pearsonr(x,t)[0]
    return pearsonr

def spike_covar(t):
    """ Calculates the coefficient of variation of interspike times 

    :param t: Rank-1 array containing spike times

    :return: coefficient of variation of interspike times 
    """

    interspike_times=np.diff(t)
    covar=scipy.stats.variation(interspike_times)
    return covar

def inflexion_spike_detector(v,t,threshold=0.4,indices=False,max_data_points=2000,voltage_threshold = -30):
    """
    Computes spike start and stop times based on extent of
    voltage deflection.

    This function requires some familiarity with Python to understand.

    :param indices: whether to return tuples of indices for each spike or times

    :return list of tuples with start and end indices of every AP
    """

    v = smooth(v)

    voltage_derivative = np.diff(v)

    voltage_above_threshold= np.where(v > voltage_threshold)

    voltage_derivative_above_threshold = np.where(voltage_derivative > threshold)

    voltage_derivative_above_threshold = np.intersect1d(voltage_derivative_above_threshold[0],
                                                        voltage_above_threshold[0])

    voltage_derivative_above_threshold = (np.array(voltage_derivative_above_threshold),)

    logging.debug('Indices where voltage derivative exceeds\
                  threshold: %s' %voltage_derivative_above_threshold)

    #this method actually sucks, we want the indices where a gap > 1,
    #use a reduce?
    diff_te = np.diff(voltage_derivative_above_threshold)

    initial_deflection_indices = np.where(diff_te>1.0)[1]

    ap_initiation_indices = [voltage_derivative_above_threshold[0][i+1] for i in initial_deflection_indices]

    ap_initiation_indices = np.append(voltage_derivative_above_threshold[0][0],ap_initiation_indices)

    logging.debug('Indices where initial deflection occurs: %s'\
                  %ap_initiation_indices)

    ap_initiation_times = t[ap_initiation_indices]

    logging.debug('Times where initial deflection occurs: %s'\
                  %ap_initiation_times)

    #we now have the times and indices of all the AP initiations, need
    #to find the corresponding end indices

    nearest_index = lambda value,arr : np.abs(arr - value).argmin()

    ap_indices = []
    ap_times = []

    for ap_initiation_index in ap_initiation_indices:

        ap_start_time = t[ap_initiation_index]
        ap_start_voltage = v[ap_initiation_index]

        offset = 10 # offset prevents corresponding time from being

        v_slice = v[ap_initiation_index + offset:ap_initiation_index+max_data_points]
        t_slice = t[ap_initiation_index+ offset:ap_initiation_index+max_data_points]

        corresponding_times = y_from_x(v_slice,
                                       t_slice,
                                       ap_start_voltage)

        logger.debug('Corresponding times: %s' %corresponding_times)

        try:
            ap_end_time = corresponding_times[nearest_index(ap_start_time,corresponding_times)]

        except:
            logger.critical('AP end time not found, AP start time: %f' %ap_start_time)
            ap_end_time = ap_start_time + 0.002 # TODO: this fix is nonsense

#            plt.plot(t,v)
#            plt.show()

            logger.critical('Corresponding times: %s' %corresponding_times)
            logger.critical('AP start time: %f' %ap_start_time)

#            voltage_plot(t,v,title='Error during spike detection')
#            voltage_plot(t[:-1],np.diff(v)*10)

        ap_end_index = nearest_index(ap_end_time,t)

        ap_times.append((ap_start_time,ap_end_time))
        ap_indices.append((ap_initiation_index,ap_end_index))

        logger.debug('Action potential start and end time: %f %f' %(ap_start_time,ap_end_time))

    if indices:
        return_value = ap_indices
    else:
        return_value = ap_times

    return return_value


def ap_integrals(v,t):
    """
    TODO:explain this fn
    """

    logger.info('Estimating AP indices')
    ap_indices = inflexion_spike_detector(v,t,indices=True)
    logger.info('AP indices found')

    integrals = []

    for ap_index_tuple in ap_indices:
        ap = v[ap_index_tuple[0]:ap_index_tuple[1]]
        ap_zeroed = ap - ap.min()

        #assume constant timestep:
        dt = t[1] - t[0]
        #estimate integral using trapezoidal rule
        integral = np.trapz(ap_zeroed,dx=dt)
        integrals.append(integral)
        logger.debug('AP integral calculated: %f' %integral)

    return np.array(integrals)

def broadening_index(v,t):
    """
    TODO:explain this fn
    TODO:add logging to this module
    """

    logger.info('Estimating integral values of spike train')

    integrals = ap_integrals(v,t)
    logger.info('AP integrals calcuated')
    integral_0 = integrals[0]
    mean_remaining_integrals = np.mean(integrals[1:])
    bi = integral_0/mean_remaining_integrals

    logger.debug('Broadening index: %f' %bi)

    return bi


def elburg_bursting(spike_times):
    """ bursting measure B as described by Elburg & Ooyen 2004

    :param spike_times: sequence of spike times

    :return: bursting measure B as described by Elburg & Ooyen 2004
    """

    interspikes_1=np.diff(spike_times)

    num_interspikes=len(spike_times)-1

    interspikes_2=[]
    for i in range(num_interspikes-1):
        interspike=interspikes_1[i]+interspikes_1[i+1]
        interspikes_2.append(interspike)

    mean_interspike=np.mean(interspikes_1)    

    var_i_1=np.var(interspikes_1)
    var_i_2=np.var(interspikes_2)

    B=(2*var_i_1-var_i_2)/(2*mean_interspike**2)

    return B

 

def load_csv_data(file_path, delimiter=',',plot=False):
    """Extracts time and voltage data from a csv file

    Data must be in a csv and in two columns, first time and second 
    voltage. Units should be SI (Volts and Seconds).

    :param file_path: full file path to file e.g /home/mike/test.csv

    :return: two lists - time and voltage

    """
    import csv

    csv_file= open(file_path, 'r')
    csv_reader=csv.reader(csv_file, delimiter=delimiter)

    v=[]
    t=[]

    i=0
    warnings_left = 5
    for row in csv_reader:

        try:

            t_value=float(row[0])*1000 #convert to ms
            v_value=float(row[1])*1000 #convert to mV

            t.append(t_value)
            v.append(v_value)

        except:
            if warnings_left >0:
                print_comment_v('Row %i invalid in %s: %s, delimiter = [%s]'%(i, file_path, row, delimiter))
                warnings_left-=1
            elif warnings_left == 0:
                print_comment_v('Supressing further warnings about %s'%(file_path))
                warnings_left-=1

        i+=1

    if plot:
        from matplotlib import pyplot        
        pyplot.plot(t,v)
        pyplot.title('Raw data')
        pyplot.xlabel('Time (ms)')
        pyplot.ylabel('Voltage (mV)')
        pyplot.show()

    return t,v


def phase_plane(t,y,plot=False): #plot should be here really
    """
    Return a tuple with two vectors corresponding to the phase plane of
    the tracetarget
    """
    dv=np.diff(y)
    dt=np.diff(t)
    dy_dt=dv/dt

    y=list(y)
    y=y[:-1]

    if plot:
        from matplotlib import pyplot
        pyplot.title('Phase Plot')
        pyplot.ylabel('dV/dt')
        pyplot.xlabel('Voltage (mV)')
        pyplot.plot(y,dy_dt)
        pyplot.show()

    return [y,dy_dt]

def filter(t,v): #still experimental

    import scipy

    fft=scipy.fft(v) # (G) and (H)  
    bp=fft[:]  
    for i in range(len(bp)): # (H-red)  
     if i>=500:bp[i]=0  
    ibp=scipy.ifft(bp) # (I), (J), (K) and (L) 

    return ibp

def pptd(t,y,bins=10,xyrange=None,dvdt_threshold=None,plot=False):
    """
    Returns a 2D map of x vs y data and the xedges and yedges. 
    in the form of a vector (H,xedges,yedges) Useful for the 
    PPTD method described by Van Geit 2007.
    """

    phase_space=phase_plane(t,y)

    #filter the phase space data
    phase_dvdt_new=[]
    phase_v_new=[]
    if dvdt_threshold!=None:
        i=0
        for dvdt in phase_space[1]:
           if dvdt>dvdt_threshold:
               phase_dvdt_new.append(phase_space[1][i])
               phase_v_new.append(phase_space[0][i])
           i+=1
        phase_space[1]=phase_dvdt_new
        phase_space[0]=phase_v_new

    if xyrange!=None:
        density_map=np.histogram2d(phase_space[1], phase_space[0], bins=bins, 
                                normed=False, weights=None)
    elif xyrange==None:
        density_map=np.histogram2d(phase_space[1], phase_space[0], bins=bins, range=xyrange, 
                                normed=False, weights=None)

    #Reverse the density map (probably not necessary as
    #it's being done because imshow has a funny origin):
    density=density_map[0][::-1]
    xedges=density_map[1]
    yedges=density_map[2]

    if plot:
        from matplotlib import pyplot
        extent = [yedges[0], yedges[-1],xedges[0], xedges[-1]]
        imgplot=pyplot.imshow(density, extent=extent)
        imgplot.set_interpolation('nearest') #makes image pixilated
        pyplot.title('Phase Plane Trajectory Density')
        pyplot.ylabel('dV/dt')
        pyplot.xlabel('Voltage (mV)')
        pyplot.colorbar()
        pyplot.show()

    return [density,xedges,yedges]

def spike_broadening(spike_width_list):
    """
    Returns the value of the width of the first AP over
    the mean value of the following APs.
    """

    first_spike = spike_width_list[0]

    if first_spike < 1e-6:
        logger.warning('First spike width <1e-6s, this indicates a problem')

    mean_following_spikes = np.mean(spike_width_list[1:])
    broadening = first_spike/mean_following_spikes


    logger.debug('Spike widths: %s' %spike_width_list)
    logger.debug('First spike: %f, Mean of following spikes: %f' %(first_spike,mean_following_spikes))
    logger.debug('Spike broadening estimate: %f' %broadening)

    return broadening

def pptd_error(t_model,v_model,t_target,v_target,dvdt_threshold=None):
    """
    Returns error function value from comparison of two phase
    pptd maps as described by Van Geit 2007.
    """

    pptd_data=pptd(t_target,v_target,dvdt_threshold=dvdt_threshold)
    target_density_map=pptd_data[0]

    xedges=pptd_data[1]
    xmin=xedges[0]
    xmax=xedges[-1]
    yedges=pptd_data[1]
    ymin=yedges[0]
    ymax=yedges[-1]
    xyrng=[[xmin, xmax], [ymin, ymax]]

    model_density_map=pptd(t_model,v_model,xyrange=xyrng,
                           dvdt_threshold=dvdt_threshold)[0]

    #calculate number of data points for the model and target:
    N_target=sum(sum(target_density_map))
    N_model=sum(sum(model_density_map))

    #normalise each map:
    normalised_target_density_map=target_density_map/float(N_target)
    normalised_model_density_map=model_density_map/float(N_model)

    #calculate the differences and calculate the mod
    difference_matrix=normalised_target_density_map-normalised_model_density_map
    difference_matrix=abs(difference_matrix)

    #root each value:
    root_matrix=difference_matrix**0.5

    #sum each element:
    summed_matrix=sum(sum(root_matrix))

    #calculate the error:
    error=summed_matrix**2

    print_comment_v('pptd error:'+ error)

    return error

def minima_phases(max_min_dictionary):
    """
    Find the phases of minima.

    Minima are found by finding the minimum value between sets of two peaks.
    The phase of the minimum relative to the two peaks is then returned.
    i.e the fraction of time elapsed between the two peaks when the minimum
    occurs is returned.

    It is very important to make sure the correct delta is specified for
    peak discrimination, otherwise unexpected results may be returned.

    :param max_min_dictionary: max_min_dictionary

    :return: phase of minimum relative to peaks.

    """

    minima_num=max_min_dictionary['minima_number']
    maxima_num=max_min_dictionary['maxima_number']
    maxima_times=max_min_dictionary['maxima_times']
    minima_times=max_min_dictionary['minima_times']

    minima_phases=[]

    for i in range(min(minima_num,maxima_num-1)):
        maximum_0_t=maxima_times[i]
        maximum_1_t=maxima_times[i+1]
        minimum_time=minima_times[i]
        phase=(minimum_time-maximum_0_t)/(maximum_1_t-maximum_0_t)
        minima_phases.append(phase)

    phase_list=[minima_times,minima_phases]

    return phase_list


class TraceAnalysis(object):
    """
    Base class for analysis of electrophysiology data

    Constructor for TraceAnalysis base class takes the following arguments:

    :param v: time-dependent variable (usually voltage)
    :param t: time-array (1-to-1 correspondence with v_array)
    :param start_analysis: time in v,t where analysis is to start
    :param end_analysis: time in v,t where analysis is to end
    """

    def __init__(self,v,t,start_analysis=0,end_analysis=None):

        self.v = np.array(v)
        self.t = np.array(t)

        if end_analysis is None:
            end_analysis = t[-1]

        start_index=self.__nearest_index(self.t,start_analysis)
        end_index=self.__nearest_index(self.t,end_analysis)

        if end_analysis!=None or start_analysis!=0: 
            self.v=v[start_index:end_index]
            self.t=t[start_index:end_index]

    def __nearest_index(self,
            array,
            target_value):

        """Finds index of first nearest value to target_value in array"""
        nparray=np.array(array)
        differences=np.abs(nparray-target_value)
        min_difference=differences.min()
        index=np.nonzero(differences==min_difference)[0][0]
        return index

    def plot_trace(self,
           save_fig=False,
           trace_name='voltage_trace.png',
           show_plot=True):
        """
        Plot the trace and save it if requested by user.
        """

        if save_fig or show_plot:
            import matplotlib.pyplot as plt

            plt.plot(self.t,self.v)
            plt.xlabel('Time (ms)')
            plt.ylabel('Votage(mV)')

            if save_fig:
                plt.savefig(trace_name)

            if show_plot:
                plt.show()


class IClampAnalysis(TraceAnalysis):
    """Analysis class for data from whole cell current injection experiments

    This is designed to work with simulations of spiking cells or
    current clamp experimental data.

    A lot of the logic here is hardcoded to work well with Cortical Layer II/III
    Pyramidal cells in Rats.

    :param v: time-dependent variable (usually voltage)
    :param t: time-vector
    :param analysis_var: dictionary containing parameters to be used
        in analysis such as delta for peak detection
    :param start_analysis: time t where analysis is to start
    :param end_analysis: time in t where analysis is to end

    """

    def __init__(self,
                 v,
                 t,
                 analysis_var,
                 start_analysis=0,
                 end_analysis=None,
                 target_data_path=None,
                 smooth_data=False,
                 show_smoothed_data=False,
                 smoothing_window_len=11,
                 max_min_method=max_min,
                 verbose=False):

        #call the parent constructor to prepare the v,t vectors:
        super(IClampAnalysis,self).__init__(v,t,start_analysis,end_analysis)
        
        self.verbose = verbose

        if smooth_data == True:
            self.v = smooth(self.v,window_len=smoothing_window_len)

        if show_smoothed_data == True:
            from matplotlib import pyplot as plt
            plt.plot(self.t,self.v)
            plt.show()

        self.delta = analysis_var['peak_delta']
        self.baseline = analysis_var['baseline']
        self.dvdt_threshold = analysis_var['dvdt_threshold']

        self.target_data_path=target_data_path

        if "peak_threshold" in analysis_var.keys():
            peak_threshold = analysis_var["peak_threshold"]
        else:
            peak_threshold = None

        
        self.max_min_dictionary = max_min_method(self.v,
                                          self.t,
                                          self.delta,
                                          peak_threshold = peak_threshold,
                                          verbose = self.verbose)
        
        print_comment('Max min dictionary calculated', verbose)
                                          

    __error_during_analysis = False #hacky way of doing this. TODO: fix

    @property
    def analysable_data(self):
        if self.max_min_dictionary['maxima_number'] < 3:
            analysable = False
            print_comment_v("Cannot analyse data: too few maxima (%i) in data: %s"%(self.max_min_dictionary['maxima_number'], self.max_min_dictionary))
        elif max(self.v) > 100.0:
            analysable = False
            print_comment_v("Cannot analyse data: max of v (%f) >100"%max(self.v))
        elif min(self.v) > -5.0:
            analysable = False
            print_comment_v("Cannot analyse data: min of v (%f) > -5"%min(self.v))
        elif max(self.v) < 10.0:
            analysable = False
            print_comment_v("Cannot analyse data: max of v (%f) < 10"%max(self.v))
        elif self.__error_during_analysis:
            analysable = False
            print_comment_v("Cannot analyse data: error during analysis...")
        else:
            analysable = True

        return analysable

    @analysable_data.setter
    def analysable_data(self, val):
        self.__error_during_analysis = True

    def plot_results(self):
        """
        Method represents the results visually.
        """

        import matplotlib.pyplot as plt

        minima_times = self.max_min_dictionary['minima_times']
        maxima_times = self.max_min_dictionary['maxima_times']

        for time in minima_times:
            plt.axvline(x=time)
        for time in maxima_times:
            plt.axvline(x=time,color='r')

        plt.xlabel('Time (ms)')
        plt.ylabel('Voltage (mV)')

        plt.plot(self.t,self.v)
        plt.show()

    def analyse(self):
        """If data is analysable analyses and puts all results into a dict"""    

        if self.analysable_data:
            analysis_results = {}
            max_min_dictionary=self.max_min_dictionary

            analysis_results['average_minimum'] = np.average(max_min_dictionary['minima_values'])
            analysis_results['average_maximum'] = np.average(max_min_dictionary['maxima_values'])
            analysis_results['min_peak_no'] = max_min_dictionary['minima_number']
            analysis_results['max_peak_no'] = max_min_dictionary['maxima_number']
            analysis_results['mean_spike_frequency'] = mean_spike_frequency(max_min_dictionary['maxima_times'])
            analysis_results['interspike_time_covar'] = spike_covar(max_min_dictionary['maxima_times'])
            analysis_results['first_spike_time'] = max_min_dictionary['maxima_times'][0]
            
            max_min_isi = max_min_interspike_time(max_min_dictionary['maxima_times'])
            analysis_results['max_interspike_time'] = max_min_isi[0]
            analysis_results['min_interspike_time'] = max_min_isi[1]
            
            trough_phases=minima_phases(max_min_dictionary)

            try:
                analysis_results['trough_phase_adaptation'] = exp_fit(trough_phases[0],trough_phases[1])
            except:
                logging.warning('trough_phase_adaptation raising an error')

            spike_width_list = spike_widths(self.v,self.t,max_min_dictionary,self.baseline,self.delta)

            try:
                analysis_results['spike_width_adaptation'] = exp_fit(spike_width_list[0],spike_width_list[1])
            except:
                logging.warning('spike_width_adaptation raising an exception, exp_fit looks problematic')

            spike_frequency_list = spike_frequencies(max_min_dictionary['maxima_times'])
            analysis_results['peak_decay_exponent'] = three_spike_adaptation(max_min_dictionary['maxima_times'],max_min_dictionary['maxima_values'])

            analysis_results['trough_decay_exponent'] = three_spike_adaptation(max_min_dictionary['minima_times'],max_min_dictionary['minima_values'])

            analysis_results['spike_frequency_adaptation'] = exp_fit(spike_frequency_list[0],spike_frequency_list[1])
            analysis_results['spike_broadening'] = spike_broadening(spike_width_list[1])
            analysis_results['peak_linear_gradient'] = linear_fit(max_min_dictionary["maxima_times"],max_min_dictionary["maxima_values"])

            #analysis_results['broadening_index'] = broadening_index(self.v,self.t)

            #this line here is because PPTD needs to be compared directly with experimental data:
            if self.target_data_path!=None and len(self.target_data_path)>0:
                t_experimental,v_experimental=load_csv_data(self.target_data_path)
                try:
                    analysis_results['pptd_error']=pptd_error(self.t,self.v,
                                              t_experimental,v_experimental,
                                              dvdt_threshold=self.dvdt_threshold)
                except:
                    print_comment_v('WARNING PPTD failure')
                    analysis_results['pptd_error'] = 1

            self.analysis_results=analysis_results

        else: 
            self.analysis_results = None
            print_comment_v('Data not suitable for analysis',True)
            
        print_comment('Analysis complete',self.verbose)

        return self.analysis_results


class NetworkAnalysis(object):
    """Analysis class for networks of spiking cells, mainly simulation data

    :param v: time-dependent variable (usually voltage)
    :param t: time-vector
    :param analysis_var: dictionary containing parameters to be used
        in analysis such as delta for peak detection
    :param start_analysis: time t where analysis is to start
    :param end_analysis: time in t where analysis is to end

    """

    def __init__(self,
                 volts,
                 t,
                 analysis_var,
                 start_analysis=0,
                 end_analysis=None,
                 smooth_data=False,
                 show_smoothed_data=False,
                 smoothing_window_len=11,
                 verbose=False):

        self.volts = volts
        
        if not isinstance(self.volts, dict):
            raise ValueError("NetworkAnalysis requires a dict of y values with reference vs. voltage trace")
        
        for ref in self.volts.keys():
            if not len(t)==len(self.volts[ref]):
                raise ValueError("One of the voltage traces (%s) has a different length to the time trace (%s != %s)!"%(ref, len(self.volts[ref]), len(t)))
                
        self.t = t
        
        self.verbose=verbose
        
        if smooth_data == True:
            for ref in volts.keys():
                # TODO improve this craziness
                self.volts[ref] = smooth(np.array(self.volts[ref]), window_len=smoothing_window_len).tolist()
                
            if show_smoothed_data == True:
                from matplotlib import pyplot as plt
                
                for ref in volts.keys():
                    plt.plot(self.t, self.volts[ref], label=ref)
                plt.legend()
                plt.show()


        start_index=self.__nearest_index(self.t,start_analysis)
        
        if end_analysis is None:
            end_analysis = t[-1]
            end_index=len(self.t)-1
        else:
            end_index=self.__nearest_index(self.t,end_analysis)
    

        if end_analysis!=None or start_analysis!=0:  
            self.t=t[start_index:end_index+1]
            for ref in volts.keys():
                self.volts[ref] =volts[ref][start_index:end_index+1]
                
        self.delta = analysis_var['peak_delta']
        self.baseline = analysis_var['baseline']
        self.dvdt_threshold = analysis_var['dvdt_threshold']


        if "peak_threshold" in analysis_var.keys():
            peak_threshold = analysis_var["peak_threshold"]
        else:
            peak_threshold = None

        
        self.max_min_dictionaries = {}
        for ref in self.volts.keys():
            max_min_dict = max_min_simple(self.volts[ref],
                                   self.t,
                                   self.delta,
                                   peak_threshold = peak_threshold,
                                   verbose=self.verbose)
                                   
            self.max_min_dictionaries[ref] = max_min_dict

    
    def __nearest_index(self,
        array,
        target_value):

        """Finds index of first nearest value to target_value in array"""
        nparray=np.array(array)
        differences=np.abs(nparray-target_value)
        min_difference=differences.min()
        index=np.nonzero(differences==min_difference)[0][0]
        return index

        
    '''
    targets: the standard targets to evaluate (min_peak_no, minimum, spike_broadening, etc). If None, evaluate all 
    extra_targets: used if targets==None for specifying additional targets, e.g. cell0:value_100
    '''
    def analyse(self, targets=None, extra_targets=None):
        """ Analyses and puts all results into a dict"""    

        analysis_results = {}
        
        for ref in self.volts.keys():
            max_min_dictionary=self.max_min_dictionaries[ref]
            
            print_comment('Analysing data with %i maxima, %i minima %s'%(max_min_dictionary['maxima_number'], 
                                                              max_min_dictionary['minima_number'],
                                                              '(targets: %s)'%targets if targets else ''), self.verbose)
            
            v = self.volts[ref]
            
            pre = '%s:'%(ref)
            
            max = -1 * sys.float_info.max
            min = sys.float_info.max
            for val in v:
                if val > max: max = val
                if val < min: min = val
            if targets==None or pre+'maximum' in targets:
                analysis_results[pre+'maximum'] = max
            if targets==None or pre+'minimum' in targets:
                analysis_results[pre+'minimum'] = min
            print_comment('Max: %s, min %s'%(max, min), self.verbose)    
            
            if targets==None or pre+'min_peak_no' in targets:
                analysis_results[pre+'min_peak_no'] = max_min_dictionary['minima_number']
                
            if targets==None or pre+'max_peak_no' in targets:
                analysis_results[pre+'max_peak_no'] = max_min_dictionary['maxima_number']
                
            if max_min_dictionary['maxima_number'] >= 1:
                
                if targets==None or pre+'average_maximum' in targets:
                    analysis_results[pre+'average_maximum'] = np.average(max_min_dictionary['maxima_values'])
                if targets==None or pre+'first_spike_time' in targets:
                    analysis_results[pre+'first_spike_time'] = max_min_dictionary['maxima_times'][0]
                
            if max_min_dictionary['minima_number'] >= 1:
                
                if targets==None or pre+'average_minimum' in targets:
                    analysis_results[pre+'average_minimum'] = np.average(max_min_dictionary['minima_values'])
            

            if targets==None or pre+'mean_spike_frequency' in targets:
                    
                if max_min_dictionary['maxima_number'] >= 3:
                    analysis_results[pre+'mean_spike_frequency'] = mean_spike_frequency(max_min_dictionary['maxima_times'])
                else:
                    analysis_results[pre+'mean_spike_frequency'] = 0
                
            if max_min_dictionary['maxima_number'] >= 3:
                
                if targets==None or pre+'interspike_time_covar' in targets:
                    analysis_results[pre+'interspike_time_covar'] = spike_covar(max_min_dictionary['maxima_times'])


                if targets==None or pre+'trough_phase_adaptation' in targets:
                    trough_phases=minima_phases(max_min_dictionary)

                    try:
                        analysis_results[pre+'trough_phase_adaptation'] = exp_fit(trough_phases[0],trough_phases[1])
                    except:
                        logging.warning('trough_phase_adaptation raising an error')

                if targets==None or pre+'spike_broadening' in targets or pre+'spike_width_adaptation' in targets:

                    spike_width_list = spike_widths(v,self.t,max_min_dictionary,self.baseline,self.delta)

                    if len(spike_width_list)>=2 and len(spike_width_list[0])>0:
                        if targets==None or pre+'spike_broadening' in targets:
                            analysis_results[pre+'spike_broadening'] = spike_broadening(spike_width_list[1])

                        if targets==None or pre+'spike_width_adaptation' in targets:
                            try:
                                analysis_results[pre+'spike_width_adaptation'] = exp_fit(spike_width_list[0],spike_width_list[1])
                            except:
                                logging.warning('spike_width_adaptation raising an exception, exp_fit looks problematic')
                    else:
                        logging.warning('spike_width_list does not have enough points for calculating spike_width_adaptation or spike_broadening: %s'%spike_width_list)
                        

                max_min_isi = max_min_interspike_time(max_min_dictionary['maxima_times'])
                
                if targets==None or pre+'max_interspike_time' in targets:
                    analysis_results[pre+'max_interspike_time'] = max_min_isi[0]
                if targets==None or pre+'min_interspike_time' in targets:
                    analysis_results[pre+'min_interspike_time'] = max_min_isi[1]

                if targets==None or pre+'peak_decay_exponent' in targets or pre+'spike_frequency_adaptation' in targets:
                    spike_frequency_list = spike_frequencies(max_min_dictionary['maxima_times'])
                    
                    if targets==None or pre+'peak_decay_exponent' in targets:
                        analysis_results[pre+'peak_decay_exponent'] = three_spike_adaptation(max_min_dictionary['maxima_times'],max_min_dictionary['maxima_values'])

                    if targets==None or pre+'spike_frequency_adaptation' in targets:
                        analysis_results[pre+'spike_frequency_adaptation'] = exp_fit(spike_frequency_list[0],spike_frequency_list[1])

                if targets==None or pre+'trough_decay_exponent' in targets:
                    analysis_results[pre+'trough_decay_exponent'] = three_spike_adaptation(max_min_dictionary['minima_times'],max_min_dictionary['minima_values'])

                if targets==None or pre+'peak_linear_gradient' in targets:
                    analysis_results[pre+'peak_linear_gradient'] = linear_fit(max_min_dictionary["maxima_times"],max_min_dictionary["maxima_values"])

            if targets==None or pre+'average_last_1percent' in targets:
                num_points_to_ave = int(len(v)/100.0)
                last_vs = v[len(v)-num_points_to_ave:]
                ave = 0
                for vv in last_vs: 
                    ave+=vv 
                ave = ave/len(last_vs)
                print_comment("Getting average of last %i points (%s->%s) of all %i (%s->%s): %s"%(len(last_vs),last_vs[0],last_vs[-1],len(v),v[0],v[-1], ave), self.verbose)
                analysis_results[pre+'average_last_1percent'] = ave
                
            
            other_targets = []
            
            if targets!=None: 
                other_targets.extend(targets)
            if extra_targets!=None: 
                other_targets.extend(extra_targets)
                
            for target in other_targets:
                
                # e.g. cell0:value_100 => value at 100ms
                if target.startswith(pre+"value_"):
                    target_time = float(target.split(':')[1].split('_')[1])
                    i=0
                    while self.t[i] < target_time:
                        value = v[i]
                        i+=1
                    analysis_results[target] = value
                    
                # e.g. cell0:average_100_200 => average value between 100ms & 200ms
                if target.startswith(pre+"average_"):
                    try:
                        start_time = float(target.split(':')[1].split('_')[1])
                        end_time = float(target.split(':')[1].split('_')[2])

                        average = 0
                        num = 0
                        for i in range(len(self.t)):
                            if self.t[i] >= start_time and self.t[i] <= end_time:
                                average += v[i]
                                num+=1
                        if num>0:        
                            average = average/num
                            analysis_results[target] = average
                    except ValueError:
                        # Ignoring as it could be average_last_1percent etc.
                        pass
                
                
        self.analysis_results=analysis_results

        return self.analysis_results

