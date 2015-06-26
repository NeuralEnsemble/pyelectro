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
        
        if os.path.isfile('100pA_1.csv'):
            times, data = analysis.load_csv_data('100pA_1.csv')
        elif os.path.isfile('test/100pA_1.csv'):
            times, data = analysis.load_csv_data('test/100pA_1.csv')
        
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
        
        

    def test_iclamp_analysis_data(self):
        
        print("- test_iclamp_analysis_data()")
        
        analysis_var={'peak_delta':0,'baseline':0,'dvdt_threshold':0, 'peak_threshold':0}
        
        times, data = self.get_real_data()

        analysis_data=analysis.IClampAnalysis(data,
                                           times,
                                           analysis_var,
                                           start_analysis=0,
                                           end_analysis=1000,
                                           smooth_data=False,
                                           show_smoothed_data=False)
                                          
        analysed = analysis_data.analyse()
                                           
        pp.pprint(analysed)
        



    def test_max_min(self):
        
        print("- test_max_min()")

        times, data = self.get_data(10,1,-80,30,[2,6])
        
        res = analysis.max_min(data, times)
        
        assert(res['minima_values'] == [-80])
        assert(res['maxima_values'][0] == 30)
        
        times, data = self.get_data(600,0.1,-80,40,[23,55.5,120,333.88,555.999])
        
        res = analysis.max_min(data, times)
        
        
        times, data = self.get_real_data()
        print(io.summary(data, "Real data set"))
        
        res = analysis.max_min(data, times)
        
        assert(res['maxima_times'] == [164, 187, 210, 233, 255, 278, 299, 320, 341, 362, 383, 405, 426, 447, 467, 487, 508, 528, 549, 570, 590, 612, 633, 654, 675, 695, 716, 736, 757, 779, 820, 841, 861, 882])
        
        
        
    def runTest(self):
        print("Running tests in TestAnalysis")

if __name__ == '__main__':
    
    ta = TestAnalysis()
    ta.test_iclamp_analysis_data()