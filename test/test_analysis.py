import os.path

from pyelectro import analysis
from pyelectro import io

import os
import pprint

pp = pprint.PrettyPrinter(indent=4)

try:
    import unittest2 as unittest
except ImportError:
    import unittest

class TestAnalysis(unittest.TestCase):


    def get_data(self, tmax, dt, vmin, vmax, spikes):

        times = []
        data = []
        t = 0
        while (t<=tmax):
            if len(spikes)>0 and t>=spikes[0]:
                v = vmax
                spikes = spikes[1:]
            else:
                v = vmin
            data.append(v)
            times.append(t)
            t+=dt

        return times, data

    def get_real_data(self):

        data_file = 'Gran_0.dat'
        delimiter = '\t'
        if os.path.isfile(data_file):
            times, data = analysis.load_csv_data(data_file, delimiter=delimiter)
        elif os.path.isfile('test/'+data_file):
            times, data = analysis.load_csv_data('test/'+data_file, delimiter=delimiter)

        print("Loaded data with %i times & %i datapoints from %s"%(len(times),len(data),data_file))
        return times, data



    def test_spike_frequencies(self):

        print("- test_spike_frequencies()")

        data = [10,20,40,140]

        times, freqs = analysis.spike_frequencies(data)

        print("Times: %s"%times)
        print("Freqs: %s"%freqs)


        self.assertEqual(freqs[0],100)
        self.assertEqual(freqs[1],50)
        self.assertEqual(freqs[2],10)


    def add_horizontal_line(self, y, times):

        import matplotlib.pyplot as pylab
        ys = [y,y]
        xs = [times[0], times[-1]]
        pylab.plot(xs, ys, 'k--') 


    def test_iclamp_analysis_data(self, show=False):

        print("- test_iclamp_analysis_data()")

        analysis_var={'peak_delta':0,'baseline':0,'dvdt_threshold':0, 'peak_threshold':0}

        times, volts = self.get_real_data()

        max_min_methods = [analysis.max_min, analysis.max_min_simple]

        for max_min_method in max_min_methods:

            print('--------------  Analysing with: %s'%max_min_method)
            analysis_data=analysis.IClampAnalysis(volts,
                                               times,
                                               analysis_var,
                                               start_analysis=0,
                                               end_analysis=1000,
                                               smooth_data=False,
                                               show_smoothed_data=False,
                                               max_min_method=max_min_method)

            analysed = analysis_data.analyse()

            pp.pprint(analysed)

            test_data = \
                {'average_maximum': 20.332122777777784,
                'average_minimum': -78.491198000000011,
                'first_spike_time': 108.44,
                'interspike_time_covar': 0.019741062134352557,
                'max_peak_no': 18,
                'mean_spike_frequency': 34.545824019508231,
                'min_peak_no': 17,
                'peak_decay_exponent': -0.064912249086890028,
                'peak_linear_gradient': -0.0020092762353974025,
                'spike_broadening': 1.0495985656104889,
                'spike_frequency_adaptation': 0.015301587514290844,
                'spike_width_adaptation': 0.0078514736435321177,
                'trough_decay_exponent': 0.0043242589967645087,
                'trough_phase_adaptation': 0.01048418950808087,
                'max_interspike_time': 31.080000000000013,
                'min_interspike_time': 28.70999999999998}

            for key in analysed.keys():
                self.assertAlmostEqual(analysed[key],test_data[key])


            if show:
                import matplotlib.pyplot as pylab
                fig = pylab.figure()
                fig.canvas.set_window_title("Data loaded (%i traces at %i time points)"%(len(volts),len(times)))

                pylab.xlabel('Time (ms)')
                pylab.ylabel('Voltage (mV)')
                pylab.grid('on')


                maxmin = analysis_data.max_min_dictionary
                pp.pprint(maxmin)

                self.add_horizontal_line(analysis_data.analysis_results['average_maximum'], times)
                self.add_horizontal_line(analysis_data.analysis_results['average_minimum'], times)

                for i in range(len(maxmin['maxima_times'])):
                    pylab.plot(maxmin['maxima_times'][i],maxmin['maxima_values'][i],'go')

                for i in range(len(maxmin['minima_times'])):
                    pylab.plot(maxmin['minima_times'][i],maxmin['minima_values'][i],'ro')

                pylab.plot(times, volts)


        if show:
            pylab.show()

    def test_net_analysis_data(self):

        print("- test_net_analysis_data()")

        analysis_var={'peak_delta':0,'baseline':0,'dvdt_threshold':0, 'peak_threshold':0}

        times, volts0 = self.get_real_data()
        volts = {'v': volts0}

        max_min_methods = [analysis.max_min_simple]

        for max_min_method in max_min_methods:

            print('--------------  Analysing with: %s in NetworkAnalysis'%max_min_method)
            analysis_data=analysis.NetworkAnalysis(volts,
                                               times,
                                               analysis_var,
                                               start_analysis=0,
                                               end_analysis=1000,
                                               smooth_data=False,
                                               show_smoothed_data=False)

            analysed = analysis_data.analyse(extra_targets=['v:value_50','v:value_200','v:average_200_201','v:average_100_200'])

            pp.pprint(analysed)

            test_data = \
                {'v:average_100_200': -55.34134969170771,
                'v:average_200_201': -73.48867227722772,
                'v:average_last_1percent': -63.559740088571395,
                'v:average_maximum': 20.332122777777784,
                'v:average_minimum': -78.491198000000011,
                'v:first_spike_time': 108.44,
                'v:interspike_time_covar': 0.019741062134352557,
                'v:max_peak_no': 18,
                'v:maximum': 23.789272,
                'v:mean_spike_frequency': 34.545824019508231,
                'v:min_peak_no': 17,
                'v:minimum': -79.47941999999999,
                'v:peak_decay_exponent': -0.064912249086890028,
                'v:peak_linear_gradient': -0.0020092762353973994,
                'v:spike_broadening': 1.0495985656104889,
                'v:spike_frequency_adaptation': 0.015301587514290829,
                'v:spike_width_adaptation': 0.0078514736435321194,
                'v:trough_decay_exponent': 0.0043242589967644922,
                'v:trough_phase_adaptation': 0.010484189508080874,
                'v:value_200': -75.04318,
                'v:value_50': -62.48935,
                'v:max_interspike_time': 31.080000000000013,
                'v:min_interspike_time': 28.70999999999998}

            for key in analysed.keys():
                self.assertAlmostEqual(analysed[key],test_data[key])



    def test_max_min(self):

        print("- test_max_min()")

        times, data = self.get_data(10,1,-80,30,[2,6])

        res = analysis.max_min(data, times)

        assert(res['minima_values'] == [-80])
        assert(res['maxima_values'][0] == 30)

        times, data = self.get_data(600,0.1,-80,40,[23,55.5,120,333.88,555.999])

        res = analysis.max_min(data, times, verbose=True)


        times, data = self.get_real_data()
        print(io.summary(data, "Real data set"))

        pts = [0, None, -1000]
        for pt in pts:
            res = analysis.max_min(data, times, peak_threshold=pt, verbose=True)

            print('Found %i maxima: %s'%(len(res['maxima_times']),res['maxima_times']))

            assert(res['maxima_times'] == [108.44, 139.52, 169.08, 198.15, 227.02, 255.80000000000004, 284.55, 313.26, 341.98, 370.7, 399.40999999999997, 428.13, 456.84999999999997, 485.58, 514.3100000000001, 543.04, 571.77, 600.54])



    def runTest(self):
        print("Running tests in TestAnalysis")

if __name__ == '__main__':

    import logging
    logging.basicConfig(level=logging.WARNING)
    ta = TestAnalysis()
    #ta.test_iclamp_analysis_data(True)
    ta.test_net_analysis_data()
