"""
This program reads the events in CSV format.
"""
import sys
import pandas as pd
import numpy as np
import scipy.signal as signal
from clize import run
from sigtools import modifiers


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

    y = np.convolve(w / w.sum(), s, mode="same")
    return y


def smoothed_df(dataframe):
    points_per_event = dataframe.values[0].size
    myContainer = []
    for temp_df in dataframe.values:
        smoothed_data = smooth(temp_df.flatten(), window_len=51, window="bartlett")
        myContainer.append(smoothed_data[0:points_per_event])
    return np.asarray(myContainer)


def sub_range(dataframe):
    myContainer = []
    myTimes = []
    for temp_df in dataframe.values:
        temp_df *= -1.0
        points_per_event = temp_df.flatten().size
        smoothed_data = smooth(temp_df.flatten(), window_len=51, window="bartlett")
        timesteps = np.arange(0, points_per_event, 1)
        data = np.concatenate(
            (smoothed_data[0:points_per_event, np.newaxis], timesteps[:, np.newaxis]),
            axis=1,
        )

        peaks, _ = signal.find_peaks(data[:, 0], height=0.2 * data[:, 0].max())
        results_w = signal.peak_widths(data[:, 0], peaks, rel_height=0.95)
        min_t = int(results_w[2:][0][0])
        # ranged_sm_data = data[min_t : peaks[0]]
        myContainer.append(data[min_t : peaks[0] + 1, 0])
        myTimes.append(data[min_t : peaks[0] + 1, 1])
    return np.asarray(myContainer), np.asarray(myTimes)


def save_smoothed(df1, df2, ranged=False, whole=False):
    from ROOT import TTree, std

    if ranged:
        container1, times1 = sub_range(df1)
        container2, times2 = sub_range(df2)
        subch1 = std.vector("double")()
        subch2 = std.vector("double")()
        timesteps1 = std.vector("double")()
        timesteps2 = std.vector("double")()
        tree = TTree("subrange", "")
        tree.Branch("ch1_sub", subch1)
        tree.Branch("ch2_sub", subch2)
        tree.Branch("ch1_time", timesteps1)
        tree.Branch("ch2_time", timesteps2)
        n_events = len(container1)
        for i in range(n_events):
            for ch1, t1 in zip(container1[i], times1[i]):
                subch1.push_back(ch1)
                timesteps1.push_back(t1)
            for ch2, t2 in zip(container2[i], times2[i]):
                subch2.push_back(ch2)
                timesteps2.push_back(t2)
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
    if whole:
        container1 = smoothed_df(df1)
        container2 = smoothed_df(df2)
        events = len(container1)
        points_per_event = container1[0].size
        sm_ch1 = std.vector("double")()
        sm_ch2 = std.vector("double")()
        timesteps = std.vector("double")()
        times = np.arange(0, points_per_event, 1)
        tree = TTree("channels", "")
        tree.Branch("ch1", sm_ch1)
        tree.Branch("ch2", sm_ch2)
        tree.Branch("time", timesteps)
        for i in range(events):
            for ch1, ch2, t in zip(container1[i], container2[i], times):
                sm_ch1.push_back(ch1)
                sm_ch2.push_back(ch2)
                timesteps.push_back(t)
            tree.Fill()
            sm_ch1.clear()
            sm_ch1.shrink_to_fit()
            sm_ch2.clear()
            sm_ch2.shrink_to_fit()
            timesteps.clear()
            timesteps.shrink_to_fit()
        tree.Write()


def make_csv(file1, file2):
    def merge_files(file_ch1, file_ch2):
        df1 = pd.read_csv(file_ch1, header=None, dtype="unicode")
        df2 = pd.read_csv(file_ch2, header=None, dtype="unicode")
        # res = pd.concat([df1, df2], axis=1)
        return np.asarray(pd.concat([df1, df2], axis=1))

    data_channels = merge_files(file1, file2)
    timestep_in_s = 10 * float(data_channels[7, 0])
    # timestep_in_ns = timestep_in_s * 1e9
    ayy = float(data_channels[9, 0])
    n_events = int(data_channels[10, 0])
    points_per_event = int(round(timestep_in_s / ayy))
    with open("event_information.log", "w") as f:
        f.write(f"Number of events           : {n_events} \n")
        f.write(f"Number of points per events: {points_per_event} \n")
        f.write(f"Timestep in ns             : {timestep_in_s * 1e9} \n")

    data_channels = data_channels[13:]

    x_ch1 = np.empty((n_events, points_per_event))
    x_ch2 = np.empty((n_events, points_per_event))
    for i in range(n_events):
        tempch1 = []
        tempch2 = []
        for data in data_channels[
            int(points_per_event * i) : int(points_per_event * (i + 1))
        ]:
            tempch1.append(data[0])
            tempch2.append(data[1])
        x_ch1[i] = np.transpose(np.asarray(tempch1))
        x_ch2[i] = np.transpose(np.asarray(tempch2))
    np.savetxt("ch1_test.csv", x_ch1, fmt="%1.4f", delimiter=",")
    np.savetxt("ch2_test.csv", x_ch2, fmt="%1.4f", delimiter=",")


@modifiers.kwoargs("format_csv", "format_root", "format_both")
def main(f1, f2, format_csv=False, format_root=False, format_both=False):
    """ Converts .dat files from two oscilloscope channes    
    to csv and/or .root format. In order to convert to root
    file, the input must be of csv format.

    :param f1: first input file.
    :param f2: Second input file.
    
    :param format_csv: Choose if csv output will be used.
    :param format_root: Choose if root output will be used.
    """
    if format_csv:
        make_csv(sys.argv[1], sys.argv[2])

    if format_root:
        from ROOT import TFile

        df1 = pd.read_csv(f1, header=None)
        df2 = pd.read_csv(f2, header=None)
        f = TFile("whole_data.root", "recreate")
        save_smoothed(df1, df2, ranged=True, whole=True)
        f.Write()
        f.Close()

    if format_both:
        make_csv(sys.argv[1], sys.argv[2])
        from ROOT import TFile

        df1 = pd.read_csv("ch1_test.csv", header=None)
        df2 = pd.read_csv("ch2_test.csv", header=None)
        f = TFile("whole_data.root", "recreate")
        save_smoothed(df1, df2, ranged=True, whole=True)
        f.Write()
        f.Close()


if __name__ == "__main__":
    run(main)
