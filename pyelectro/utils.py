
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

    analysis.print_comment_v(pp.pformat(analysed.analysis_results))
    maxmin = analysed.max_min_dictionary
    
    if plot:

        fig = pylab.figure()
        fig.canvas.set_window_title("Data analysed (%i traces at %i time points)"%(len(volts),len(times)))

        pylab.xlabel('Time (ms)')
        pylab.ylabel('Voltage (mV)')
        pylab.grid('on')
    
        if analysed.analysis_results:
            if analysed.analysis_results.has_key('average_maximum'):
                _add_horizontal_line(analysed.analysis_results['average_maximum'], times)

            if analysed.analysis_results.has_key('average_minimum'):
                _add_horizontal_line(analysed.analysis_results['average_minimum'], times)

        if maxmin:
            for i in range(len(maxmin['maxima_times'])):
                pylab.plot(maxmin['maxima_times'][i],maxmin['maxima_values'][i],'ro')

            for i in range(len(maxmin['minima_times'])):
                pylab.plot(maxmin['minima_times'][i],maxmin['minima_values'][i],'go')


        pylab.plot(times, volts)

        if show_plot_already:
            pylab.show()
            
    return analysed.analysis_results


def simple_network_analysis(volts, 
                           times, 
                           analysis_var=None, 
                           start_analysis = 0,
                           end_analysis = None,
                           plot=False, 
                           show_plot_already = True,
                           targets=None,
                           extra_targets=None,
                           verbose=False):
    
    if analysis_var == None:
        analysis_var={'peak_delta':0,
                      'baseline':0,
                      'dvdt_threshold':0,
                      'peak_threshold':0}
                     

    analysed=analysis.NetworkAnalysis(volts,
                                     times,
                                     analysis_var,
                                     start_analysis=start_analysis,
                                     end_analysis= end_analysis if end_analysis is not None else times[-1],
                                     smooth_data=False,
                                     show_smoothed_data=False,
                                     verbose=verbose)
                                     

    analysed.analyse(targets=targets, extra_targets=extra_targets)

    analysis.print_comment_v(pp.pformat(analysed.analysis_results))
    
    if plot:
        fig = pylab.figure()
        fig.canvas.set_window_title("Data analysed (%i traces at %i time points): %s"%(len(volts.keys()),len(times), volts.keys()))

        pylab.xlabel('Time (ms)')
        pylab.ylabel('Voltage (mV)')
        pylab.grid('on')
            
        for vk in volts.keys():
            vs = volts[vk]
            maxmin = analysed.max_min_dictionaries[vk]
            pre = '%s:'%vk
            

            if analysed.analysis_results:
                if analysed.analysis_results.has_key(pre+'average_maximum'):
                    _add_horizontal_line(analysed.analysis_results[pre+'average_maximum'], times)
                    
                if analysed.analysis_results.has_key(pre+'maximum'):
                    _add_horizontal_line(analysed.analysis_results[pre+'maximum'], times)

                if analysed.analysis_results.has_key(pre+'average_minimum'):
                    _add_horizontal_line(analysed.analysis_results[pre+'average_minimum'], times)
                    
                if analysed.analysis_results.has_key(pre+'minimum'):
                    _add_horizontal_line(analysed.analysis_results[pre+'minimum'], times)

            if maxmin:
                for i in range(len(maxmin['maxima_times'])):
                    pylab.plot(maxmin['maxima_times'][i],maxmin['maxima_values'][i],'ro')

                for i in range(len(maxmin['minima_times'])):
                    pylab.plot(maxmin['minima_times'][i],maxmin['minima_values'][i],'go')


            pylab.plot(times, vs)

        if show_plot_already:
            pylab.show()
            
    return analysed.analysis_results