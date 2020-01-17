"""
This program reads the events in CSV format.
"""
import sys
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import iminuit
import probfit
import uproot
from tqdm import tqdm
from pkgs.fit.fit import sigmoid


def smooth(x, window_len=11, window="hanning"):

    if not window in ["flat", "hanning", "hamming", "bartlett", "blackman"]:
        raise ValueError(
            "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        )

    s = np.r_[x[window_len - 1 : 0 : -1], x, x[-2 : -window_len - 1 : -1]]
    # print(len(s))
    if window == "flat":  # moving average
        w = np.ones(window_len, "d")
    else:
        w = eval("np." + window + "(window_len)")

    y = np.convolve(w / w.sum(), s, mode="same")
    return y



# this function will take an event with two peaks
def analysis(raw_data):
    debug = False
    n_events = raw_data.size
    n_per_event = raw_data[0].size
    timesteps = np.arange(0, n_per_event, 1)
    raw_data = np.negative(raw_data)
    risetimes = []
    t1s = []
    t2s = []
    evt = []
    print(f'Total number of events: {n_events}')
    # for i in range(15, 16):
    for i in range(n_events):
        smoothed_data = smooth(raw_data[i], window_len=51, window='bartlett')
        data = np.concatenate((timesteps[:, np.newaxis],
                               smoothed_data[:n_per_event, np.newaxis]), axis=1)
        peak_range = [0.2, data[:, 1].max()]  # cut for two peaks
        peaks, properties = signal.find_peaks(data[:, 1],
                                              height=peak_range, distance=200)
        results_w = signal.peak_widths(data[:, 1], peaks, rel_height=0.95)
        if len(peaks) >= 2 and data[peaks[0], 1] > data[peaks[1], 1] and data[:, 1].max() < 1.5 and peaks[1] < 9000:
            peaks_start = [int(p) for p in (results_w[2])]
            ranged_sm_p1 = data[peaks_start[0]:peaks[0]]
            ranged_sm_p2 = data[peaks_start[1]:peaks[1]]
            temp_times = []
            """
            Fitting first peak:
            """
            t1 = np.arange(0, ranged_sm_p1[:, 0].size)
            chi2_fit = probfit.Chi2Regression(
                sigmoid,
                t1,
                ranged_sm_p1[:, 1],
            )

            minuit = iminuit.Minuit(
                chi2_fit,
                p0=np.max(ranged_sm_p1[:, 1]),
                # p1=peaks[0],
                p1=t1.max(),
                p3=np.min(ranged_sm_p1[:, 1]),
                #limit_p3=(np.min(data[a : peaks[0], 1]), None),
                pedantic=False,
                print_level=-1,
            )
            minuit.migrad()
            try:
                minuit.hesse()
            except Exception as e:
                print(str(e))
            p = minuit.values
            tz1 = p[1] - np.log(p[0] / (0.2 * ranged_sm_p1[:, 1].max() - p[3]) - 1.0) / p[2]
            temp_times.append(tz1 + ranged_sm_p1[0, 0])
            """
            Fitting second peak
            """
            t2 = np.arange(0, ranged_sm_p2[:, 0].size)
            chi2_fit_p2 = probfit.Chi2Regression(
                sigmoid,
                t2,
                ranged_sm_p2[:, 1],
            )
            
            minuit_p2 = iminuit.Minuit(
                chi2_fit_p2,
                #                 p0=np.max(data[b : peaks[1], 1]),
                #                 p1=peaks[1],
                limit_p1=(t2.max()/2, t2.max()),
                #                 p3=np.min(data[b : peaks[1], 1]),
                #                 limit_p3=(np.min(data[b : peaks[1], 1]), None),
                #                 fix_p3=np.min(data[b : peaks[1], 1]),
                pedantic=False,
                print_level=-1,
            )
            
            minuit_p2.migrad()
            try:
                minuit_p2.hesse()
            except Exception as e:
                print(str(e))
                
            pp = minuit_p2.values
            tz2 = pp[1] - np.log(pp[0] / (0.2 * ranged_sm_p2[:, 1].max() - pp[3]) - 1.0) / pp[2]
            temp_times.append(tz2 + ranged_sm_p2[0, 0])
            risetimes.append([temp_times[0], temp_times[1]])
            evt.append(i)
            t1s.append(temp_times[0])
            t2s.append(temp_times[1])
            # if i == 18:
            #     # print(tz1, temp_times[1])
            #     print([s for s in t2s])

            if debug == True:
                print(f'The value of the variable p0 is {minuit.values["p0"]}')
                print(f'The value of the variable p1 is {minuit.values["p1"]}')
                print(f'The value of the variable p2 is {minuit.values["p2"]}')
                print(f'The value of the variable p3 is {minuit.values["p3"]}')
                print(tz1, tz2)
                print(risetimes)
                plt.figure(figsize=(25, 7))
                chi2_fit_p2.draw(minuit_p2)
                plt.show()
    # print(risetimes)
    # return risetimes
    return (evt, t1s, t2s)




def main(filename):

    df_np = uproot.open(filename)["t1"].array("channel1")
    ev, test1, test2 = analysis(df_np)
    test1 = np.asarray(test1)
    test2 = np.asarray(test2)
    times = np.concatenate((test1[:, np.newaxis], test2[:, np.newaxis]), axis=1)
    for i in range(test1.size):
        print(f'{ev[i]}, {test1[i]}, {test2[i]}')
    # print(times)
    # diffs = [np.abs(d[1]-d[0]) for d in testing]
    # tp = np.all(np.isnan(diffs), axis=1)
    # print(tp)
    # print(diffs)
    # plt.hist(diffs[~tp])



if __name__ == "__main__":
    main(sys.argv[1])
