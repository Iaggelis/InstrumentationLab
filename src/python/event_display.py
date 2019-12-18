import sys
import pandas as pd
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
from ROOT import TH1F, TCanvas, RDataFrame

RDF = RDataFrame


def sigmoid(x, p0, p1, p2, p3):
    return p0 / (1.0 + np.exp(-1.0 * p2 * (x - p1))) + p3


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


def display_event(filename):
    df = RDF("t1", filename)
    df_np = df.AsNumpy()
    n_events = df_np["channel1"].size
    df = None
    i = int(input("Look at event: "))
    while i != -1 and i < n_events:
        print(50 * "#")
        print(f"Looking at event {i}:")
        print(50 * "#")

        """ Plotting first channel """
        sm_ch1 = df_np["channel1"][i]
        n_per_event = sm_ch1.size()
        timesteps = np.arange(0, n_per_event, 1)
        sm_ch1 *= -1.0

        smoothed_data_ch1 = smooth(sm_ch1, window_len=51, window="bartlett")
        data_ch1 = np.concatenate(
            (smoothed_data_ch1[:n_per_event, np.newaxis], timesteps[:, np.newaxis]),
            axis=1,
        )
        peaks_ch1, _ = signal.find_peaks(
            data_ch1[:, 0], height=0.2 * data_ch1[:, 0].max()
        )
        results_w_ch1 = signal.peak_widths(data_ch1[:, 0], peaks_ch1, rel_height=0.98)
        min_t = int(results_w_ch1[2:][0][0])
        max_t = int(results_w_ch1[2:][1][0])

        print("#" * 10 + " Channel 1 " + "#" * 10)
        print(f"Peaks width: {results_w_ch1[0]}")
        for j in range(results_w_ch1[1:][0].size):
            print(f"Found peaks in range: {results_w_ch1[2][j]}, {results_w_ch1[3][j]}")

        print(
            f"Found peak at index {peaks_ch1} with amplitude: {data_ch1[peaks_ch1, 0]}"
        )

        """ Plotting second channel """
        sm_ch2 = df_np["channel2"][i]
        sm_ch2 *= -1.0
        smoothed_data_ch2 = smooth(sm_ch2, window_len=51, window="bartlett")
        data_ch2 = np.concatenate(
            (smoothed_data_ch2[:n_per_event, np.newaxis], timesteps[:, np.newaxis]),
            axis=1,
        )
        peaks_ch2, _ = signal.find_peaks(
            data_ch2[:, 0], height=0.2 * data_ch2[:, 0].max()
        )
        results_w_ch2 = signal.peak_widths(data_ch2[:, 0], peaks_ch2, rel_height=0.98)

        print("#" * 10 + " Channel 2 " + "#" * 10)
        print(f"Peaks width: {results_w_ch2[0]}")
        for j in range(results_w_ch2[1:][0].size):
            print(f"Found peaks in range: {results_w_ch2[2][j]}, {results_w_ch2[3][j]}")

        print(
            f"Found peak at index {peaks_ch2} with amplitude: {data_ch2[peaks_ch2, 0]}"
        )

        ###### Plotting ############
        fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(15, 7))

        ax[0, 0].plot(timesteps, sm_ch1, "-", lw=1)
        ax[0, 0].set_title("Raw Data, Channel 1")
        ax[0, 1].plot(data_ch1[:, 1], data_ch1[:, 0], "-", lw=1)
        ax[0, 1].plot(peaks_ch1, data_ch1[:, 0][peaks_ch1], "x")
        ax[0, 1].hlines(*results_w_ch1[1:], color="C2")
        ax[0, 1].set_title("Smoothed Data, Channel 1")

        ax[1, 0].plot(timesteps, sm_ch2, "-", lw=1)
        ax[1, 0].set_title("Raw Data, Channel 2")
        ax[1, 1].plot(data_ch2[:, 1], data_ch2[:, 0], "-", lw=1)
        ax[1, 1].plot(peaks_ch2, data_ch2[:, 0][peaks_ch2], "x")
        ax[1, 1].hlines(*results_w_ch2[1:], color="C2")
        ax[1, 1].set_title("Smoothed Data, Channel 2")

        fig.tight_layout()
        plt.show()
        i = int(input("Look at event: "))


if __name__ == "__main__":
    display_event(sys.argv[1])
