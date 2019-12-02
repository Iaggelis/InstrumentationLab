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


def sub_range(dataframe, channel=2):
    temp = dataframe.values
    if channel == 1:
        min_h = 0.02
    elif channel == 2:
        min_h = 0.1
    myContainer = []

    for i in range(temp.shape[0]):
        x = temp[i].flatten()
        x *= -1.0
        smoothed_data = smooth(x, window_len=51, window="bartlett")
        peaks, _ = signal.find_peaks(smoothed_data, height=min_h)
        ranged_sm_data = smoothed_data[
            smoothed_data >= 0.2 * smoothed_data[peaks].max()
        ]
        # indices = []
        # for m in ranged_sm_data:
        #   indices.append(np.where(m == smoothed_data))

        # indices = np.asarray(indices)
        myContainer.append(ranged_sm_data)
    return np.asarray(myContainer)


if __name__ == "__main__":
    df1 = pd.read_csv("./Labs/muon_271119/Run0/ch1_values.csv", header=None)
    df2 = pd.read_csv("./Labs/muon_271119/Run0/ch2_values.csv", header=None)

    container1 = sub_range(df1, channel=1)
    container2 = sub_range(df2)

    subch1 = std.vector("float")()
    subch2 = std.vector("float")()

    f = TFile("testing1.root", "recreate")
    tree = TTree("subrange", "")
    tree.Branch("ch1_sub", subch1)
    tree.Branch("ch2_sub", subch2)
    for i in range(len(container1)):
        for point in container1[i]:
            subch1.push_back(point)
        for point in container2[i]:
            subch2.push_back(point)
        tree.Fill()
        subch1.clear()
        subch1.shrink_to_fit()
        subch2.clear()
        subch2.shrink_to_fit()
    f.Write()
    f.Close()
