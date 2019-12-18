import sys
import numpy as np
import scipy.signal as signal
from tqdm import tqdm
from ROOT import RDataFrame

RDF = RDataFrame


def smooth(x, window_len=11, window="hanning"):
    if not window in ["flat", "hanning", "hamming", "bartlett", "blackman"]:
        raise ValueError(
            "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        )

    s = np.r_[x[window_len - 1 : 0 : -1], x, x[-2 : -window_len - 1 : -1]]
    # print(len(s))
    if window == "flat":  # moving average
        w = np.ones(window_len)
    else:
        w = eval("np." + window + "(window_len)")

    y = signal.convolve(s, w / w.sum(), mode="same")

    return y


def peak_counter(filename):
    d_peaks = []
    df = RDF("t1", filename)
    df_np = df.AsNumpy()
    n_events = df_np["channel1"].size
    n_per_event = df_np["channel1"][0].size()
    timesteps = np.arange(0, n_per_event, 1)

    df = None
    for i in tqdm(range(n_events)):
        sm_ch1 = df_np["channel1"][i]
        sm_ch1 *= -1.0

        smoothed_data_ch1 = smooth(sm_ch1, window_len=51, window="bartlett")
        data_ch1 = np.concatenate(
            (smoothed_data_ch1[:n_per_event, np.newaxis], timesteps[:, np.newaxis]),
            axis=1,
        )
        peaks_ch1, _ = signal.find_peaks(
            data_ch1[:, 0], height=0.2 * data_ch1[:, 0].max()
        )
        # results_w_ch1 = signal.peak_widths(data_ch1[:, 0], peaks_ch1, rel_height=0.98)
        if len(peaks_ch1) > 1:
            d_peaks.append(i)
    return d_peaks


if __name__ == "__main__":
    double_peaks = peak_counter(sys.argv[1])
    for index in double_peaks:
        print(index)
