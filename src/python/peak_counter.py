import sys
import numpy as np
import scipy.signal as signal
from tqdm import tqdm
from smooth import smooth


def smooth_loc(x, window_len=11, window="hanning"):
    if not window in ["flat", "hanning", "hamming", "bartlett", "blackman"]:
        raise ValueError(
            "Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        )

    s = np.r_[x[window_len - 1 : 0 : -1], x, x[-2 : -window_len - 1 : -1]]
    if window == "flat":  # moving average
        w = np.ones(window_len)
    else:
        w = eval("np." + window + "(window_len)")

    y = signal.convolve(s, w / w.sum(), mode="same")

    return y


def peak_counter(filename):
    root_usage = False
    if root_usage:
        from ROOT import RDataFrame

        RDF = RDataFrame

        df = RDF("tree", filename)
        df_np = df.AsNumpy()
        n_events = df_np["channel1"].size
        n_per_event = df_np["channel1"][0].size()
        df = None
    else:
        import uproot

        df_np = uproot.open(filename)["tree"].array("channel1")
        n_events = df_np.size
        n_per_event = df_np[0].size
    df_np = np.negative(df_np)
    timesteps = np.arange(0, n_per_event, 1)
    amp_threshold = 0.1
    dist = 1000
    print(f"The threshold for a peak is: {amp_threshold}")
    print(f"Minimum distance between 2 peaks in timesteps: {dist} ")
    d_peaks = []
    for i in tqdm(range(n_events)):
        sm_ch1 = df_np[i] - np.mean(df_np[i][0:2000])
        sm_ch1[np.where(sm_ch1 < 0.0)] = 0.0

        smoothed_data_ch1 = smooth(sm_ch1, window_len=51, window="bartlett")
        data_ch1 = np.concatenate(
            (smoothed_data_ch1[:n_per_event, np.newaxis], timesteps[:, np.newaxis]),
            axis=1,
        )
        peaks_ch1, _ = signal.find_peaks(
            data_ch1[:, 0], height=[amp_threshold, None], distance=dist
        )
        # results_w_ch1 = signal.peak_widths(data_ch1[:, 0], peaks_ch1, rel_height=0.98)
        if len(peaks_ch1) > 1:
            d_peaks.append(i)
    return d_peaks


if __name__ == "__main__":
    double_peaks = peak_counter(sys.argv[1])
    for index in double_peaks:
        print(index)
