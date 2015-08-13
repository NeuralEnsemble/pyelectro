
import matplotlib.pyplot as pylab

import pprint
pp = pprint.PrettyPrinter(indent=4)

from pyelectro import analysis


def _add_horizontal_line(y, times):
    
    ys = [y,y]
    xs = [times[0], times[-1]]
    pylab.plot(xs, ys, 'k--')

def simple_iclamp_analysis(volts, 
                           times, 
                           analysis_var=None, 
                           start_analysis = 0,
                           end_analysis = None,
                           plot=False, 
                           show_plot_already = True):
    
    if analysis_var == None:
        analysis_var={'peak_delta':0,
                      'baseline':0,
                      'dvdt_threshold':0,
                      'peak_threshold':0}

    analysed=analysis.IClampAnalysis(volts,
                                     times,
                                     analysis_var,
                                     start_analysis=start_analysis,
                                     end_analysis= end_analysis if end_analysis is not None else times[-1],
                                     smooth_data=False,
                                     show_smoothed_data=False,
                                     max_min_method=analysis.max_min_simple)
                                     

    analysed.analyse()

    pp.pprint(analysed.analysis_results)
    maxmin = analysed.max_min_dictionary
    #pp.pprint(maxmin)
    
    if plot:

        fig = pylab.figure()
        fig.canvas.set_window_title("Data analysed (%i traces at %i time points)"%(len(volts),len(times)))

        pylab.xlabel('Time (ms)')
        pylab.ylabel('Voltage (mV)')
        pylab.grid('on')
    
        _add_horizontal_line(analysed.analysis_results['average_maximum'], times)
        _add_horizontal_line(analysed.analysis_results['average_minimum'], times)

        for i in range(len(maxmin['maxima_times'])):
            pylab.plot(maxmin['maxima_times'][i],maxmin['maxima_values'][i],'ro')

        for i in range(len(maxmin['minima_times'])):
            pylab.plot(maxmin['minima_times'][i],maxmin['minima_values'][i],'go')


        pylab.plot(times, volts)

        if show_plot_already:
            pylab.show()
            
    return analysed.analysis_results