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
    n_events = raw_data.size
    n_per_event = raw_data[0].size
    timesteps = np.arange(0, n_per_event, 1)
    raw_data = np.negative(raw_data)
    risetimes = []
    for i in tqdm(range(n_events)):
        smoothed_data = smooth(raw_data[i], window_len=51, window='bartlett')
        data = np.concatenate((timesteps[:, np.newaxis],
                               smoothed_data[:n_per_event, np.newaxis]), axis=1)
        peak_range = [0.1, data[:, 1].max()]  # cut for two peaks
        peaks, properties = signal.find_peaks(data[:, 1],
                                              height=peak_range, distance=10)
        results_w = signal.peak_widths(data[:, 1], peaks, rel_height=0.95)
        if len(peaks) >= 2 and data[peaks[0], 1] > data[peaks[1], 1]:
            peaks_start = [int(p) for p in (results_w[2])]
            ranged_sm_p1 = data[peaks_start[0]:peaks[0]]
            ranged_sm_p2 = data[peaks_start[1]:peaks[1]]
            temp_times = []
            for peak_range, peak in zip([ranged_sm_p1, ranged_sm_p2], peaks):
                chi2_fit = probfit.Chi2Regression(sigmoid, peak_range[:, 0], peak_range[:, 1])
                minuit = iminuit.Minuit(
                    chi2_fit,
                    p0=np.max(peak_range[:, 1]),
                    p1=peak,
                    p3=np.min(peak_range[:, 1]),
                    # limit_p3=(0.01, 1.0),
                    pedantic=False,
                    print_level=-1,
                )
                minuit.migrad()
                try:
                    minuit.hesse()
                except Exception as e:
                    print(f'Error at event {i} :{str(e)}')
                p = minuit.values
                tz = p[1] - np.log(p[0] / (0.2 * peak_range[:, 1].max() - p[3]) - 1) / p[2]
                temp_times.append(tz)
                if tz == np.nan:
                    print(i)
            risetimes.append([temp_times[0], temp_times[1]])
    return risetimes




def main(filename):

    df_np = uproot.open(filename)["t1"].array("channel1")
    testing = analysis(df_np)
    for i, t in enumerate(testing):
        print(i, t)
    diffs = [np.abs(d[1]-d[0]) for d in testing]
    # tp = np.all(np.isnan(diffs), axis=1)
    # print(tp)
     # print(diffs)
    # plt.hist(diffs[~tp])



if __name__ == "__main__":
    main(sys.argv[1])
