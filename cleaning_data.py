import pandas as pd
import numpy as np
import scipy.signal as signal
from ROOT import TFile, TTree, std


def smooth(x, window_len=11, window="hanning"):
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

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

    y = np.convolve(w / w.sum(), s, mode="valid")
    return y


def sub_range(dataframe):
    temp = dataframe.values
    myContainer = []
    myTimes = []

    for i in range(temp.shape[0]):
        x = temp[i].flatten()
        x *= -1.0
        smoothed_data = smooth(x, window_len=51, window="bartlett")
        timesteps = np.arange(0, 1000, 1)
        data = np.concatenate(
            (smoothed_data[0:1000, np.newaxis], timesteps[:, np.newaxis]), axis=1
        )

        peaks, _ = signal.find_peaks(data[:, 0], height=0.2 * data[:, 0].max())
        results_w = signal.peak_widths(data[:, 0], peaks, rel_height=0.95)
        min_t = int(results_w[2:][0][0])
        ranged_sm_data = data[min_t : peaks[0]]

        # peaks, _ = signal.find_peaks(smoothed_data, height=min_h)
        # ranged_sm_data = smoothed_data[
        #     smoothed_data >= 0.2 * smoothed_data[peaks].max()
        # ]
        # dataa = np.concatenate(
        #     (smoothed_data[0:1000][:, np.newaxis], timesteps_in_s[:, np.newaxis]),
        #     axis=1,
        # )
        # temp_indices = []
        # for m in ranged_sm_data:
        #     temp_indices.append(np.where(m == dataa))
        # indices = []
        # for i in range(len(temp_indices)):
        #     indices.append(temp_indices[i][0][0])

        myContainer.append(ranged_sm_data[:, 0])
        myTimes.append(ranged_sm_data[:, 1])
    return np.asarray(myContainer), np.asarray(myTimes)


def smoothed_df(dataframe):
    temp = dataframe.values
    myContainer = []
    for i in range(temp.shape[0]):
        x = temp[i].flatten()
        x *= -1.0
        smoothed_data = smooth(x, window_len=51, window="bartlett")
        myContainer.append(smoothed_data)
    return np.asarray(myContainer)


def save_smoothed(df1, df2, ranged=False):

    if ranged:
        container1, times1 = sub_range(df1)
        container2, times2 = sub_range(df2)
        f = TFile("clean_data.root", "recreate")
        subch1 = std.vector("double")()
        subch2 = std.vector("double")()
        timesteps1 = std.vector("double")()
        timesteps2 = std.vector("double")()
        tree = TTree("subrange", "")
        tree.Branch("ch1_sub", subch1)
        tree.Branch("ch2_sub", subch2)
        tree.Branch("ch1_time", timesteps1)
        tree.Branch("ch2_time", timesteps2)
        for i in range(len(container1)):
            for j in range(container1[i].size):
                subch1.push_back(container1[i][j])
                timesteps1.push_back(times1[i][j])
            for j in range(container2[i].size):
                subch2.push_back(container2[i][j])
                timesteps2.push_back(times2[i][j])
            tree.Fill()
            subch1.clear()
            subch1.shrink_to_fit()
            subch2.clear()
            subch2.shrink_to_fit()
            timesteps1.clear()
            timesteps1.shrink_to_fit()
            timesteps2.clear()
            timesteps2.shrink_to_fit()
        tree.Write()
        f.Close()
    else:
        container1 = smoothed_df(df1)
        container2 = smoothed_df(df2)
        f = TFile("smooth_data.root", "recreate")
        sm_ch1 = std.vector("double")()
        sm_ch2 = std.vector("double")()
        tree = TTree("channels", "")
        tree.Branch("sm_ch1", sm_ch1)
        tree.Branch("sm_ch2", sm_ch2)
        for i in range(len(container1)):
            for j in range(container1[i].size):
                sm_ch1.push_back(container1[i][j])
                sm_ch2.push_back(container2[i][j])
            tree.Fill()
            sm_ch1.clear()
            sm_ch1.shrink_to_fit()
            sm_ch2.clear()
            sm_ch2.shrink_to_fit()
        tree.Write()
        f.Close()


if __name__ == "__main__":
    df1 = pd.read_csv("./Labs/muon_271119/Run0/ch1_values.csv", header=None)
    df2 = pd.read_csv("./Labs/muon_271119/Run0/ch2_values.csv", header=None)
    save_smoothed(df1, df2, ranged=True)
