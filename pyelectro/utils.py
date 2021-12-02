"""Pyelectro utility functions.

This module includes a set of utility functions to simplify common tasks.
"""
import matplotlib.pyplot as pylab
import pprint
from pyelectro import analysis


pp = pprint.PrettyPrinter(indent=4)


def _add_horizontal_line(y, times):

    ys = [y, y]
    xs = [times[0], times[-1]]
    pylab.plot(xs, ys, "k--")


def simple_iclamp_analysis(
    volts,
    times,
    analysis_var=None,
    start_analysis=0,
    end_analysis=None,
    plot=False,
    show_plot_already=True,
):
    """
    A utility function to quickly carry out a simple current clamp analysis
    (IClampAnalysis).

    :param v: time-dependent variable (usually voltage)
    :type v: iterable
    :param t: time-array (1-to-1 correspondence with v_array)
    :type t: iterable
    :param start_analysis: time in v,t where analysis is to start
    :type start_analysis: float
    :param end_analysis: time in v,t where analysis is to end
    :type end_analysis: float

    :returns: dictionary of analysis results
    """
    if analysis_var is None:
        analysis_var = {
            "peak_delta": 0,
            "baseline": 0,
            "dvdt_threshold": 0,
            "peak_threshold": 0,
        }

    analysed = analysis.IClampAnalysis(
        volts,
        times,
        analysis_var,
        start_analysis=start_analysis,
        end_analysis=end_analysis if end_analysis is not None else times[-1],
        smooth_data=False,
        show_smoothed_data=False,
        max_min_method=analysis.max_min_simple,
    )

    analysed.analyse()

    analysis.print_comment_v(pp.pformat(analysed.analysis_results))
    maxmin = analysed.max_min_dictionary

    if plot:

        fig = pylab.figure()
        pylab.get_current_fig_manager().set_window_title(
            "Data analysed (%i traces at %i time points)" % (len(volts), len(times))
        )

        pylab.xlabel("Time (ms)")
        pylab.ylabel("Voltage (mV)")
        pylab.grid("on")

        if analysed.analysis_results:
            if "average_maximum" in analysed.analysis_results:
                _add_horizontal_line(
                    analysed.analysis_results["average_maximum"], times
                )

            if "average_minimum" in analysed.analysis_results:
                _add_horizontal_line(
                    analysed.analysis_results["average_minimum"], times
                )

        if maxmin:
            for i in range(len(maxmin["maxima_times"])):
                pylab.plot(maxmin["maxima_times"][i], maxmin["maxima_values"][i], "ro")

            for i in range(len(maxmin["minima_times"])):
                pylab.plot(maxmin["minima_times"][i], maxmin["minima_values"][i], "go")

        pylab.plot(times, volts)

        if show_plot_already:
            pylab.show()

    return analysed.analysis_results


def simple_network_analysis(
    volts,
    times,
    analysis_var=None,
    start_analysis=0,
    end_analysis=None,
    plot=False,
    show_plot_already=True,
    targets=None,
    extra_targets=None,
    verbose=False,
):
    """
    A utility function to quickly carry out a simple network analysis
    (IClampAnalysis).

    :param v: time-dependent variable (usually voltage)
    :type v: iterable
    :param t: time-vector
    :type t: iterable
    :param analysis_var: dictionary containing parameters to be used
        in analysis such as delta for peak detection
    :type analysis_var: dict
    :param start_analysis: time t where analysis is to start
    :type start_analysis: float
    :param end_analysis: time in t where analysis is to end
    :type end_analysis: float

    :returns: dictionary of analysis results

    """

    if analysis_var is None:
        analysis_var = {
            "peak_delta": 0,
            "baseline": 0,
            "dvdt_threshold": 0,
            "peak_threshold": 0,
        }

    analysed = analysis.NetworkAnalysis(
        volts,
        times,
        analysis_var,
        start_analysis=start_analysis,
        end_analysis=end_analysis if end_analysis is not None else times[-1],
        smooth_data=False,
        show_smoothed_data=False,
        verbose=verbose,
    )

    analysed.analyse(targets=targets, extra_targets=extra_targets)

    analysis.print_comment_v(pp.pformat(analysed.analysis_results))

    if plot:
        fig = pylab.figure()
        fig.canvas.set_window_title(
            "Data analysed (%i traces at %i time points): %s"
            % (len(volts.keys()), len(times), volts.keys())
        )

        pylab.xlabel("Time (ms)")
        pylab.ylabel("Voltage (mV)")
        pylab.grid("on")

        for vk in volts.keys():
            vs = volts[vk]
            maxmin = analysed.max_min_dictionaries[vk]
            pre = "%s:" % vk

            if analysed.analysis_results:
                if pre + "average_maximum" in analysed.analysis_results:
                    _add_horizontal_line(
                        analysed.analysis_results[pre + "average_maximum"], times
                    )

                if pre + "maximum" in analysed.analysis_results:
                    _add_horizontal_line(
                        analysed.analysis_results[pre + "maximum"], times
                    )

                if pre + "average_minimum" in analysed.analysis_results:
                    _add_horizontal_line(
                        analysed.analysis_results[pre + "average_minimum"], times
                    )

                if pre + "minimum" in analysed.analysis_results:
                    _add_horizontal_line(
                        analysed.analysis_results[pre + "minimum"], times
                    )

            if maxmin:
                for i in range(len(maxmin["maxima_times"])):
                    pylab.plot(
                        maxmin["maxima_times"][i], maxmin["maxima_values"][i], "ro"
                    )

                for i in range(len(maxmin["minima_times"])):
                    pylab.plot(
                        maxmin["minima_times"][i], maxmin["minima_values"][i], "go"
                    )

            pylab.plot(times, vs)

        if show_plot_already:
            pylab.show()

    return analysed.analysis_results
