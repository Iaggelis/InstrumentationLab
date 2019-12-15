import sys
import pandas as pd
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt


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


def plot_event(df1, df2):
    with open("event_information.log", "r") as f:
        info = f.readlines()
        n_events = int(info[0].split()[-1])
        n_per_event = int(info[1].split()[-1])
    i = int(input("Look at event: "))
    while i != -1 and i < n_events:
        timesteps = np.arange(0, n_per_event, 1)
        """ Plotting first channel """
        sm_ch1 = df1.values[i].flatten()
        sm_ch1 *= -1.0
        smoothed_data_ch1 = smooth(sm_ch1, window_len=51, window="bartlett")
        data_ch1 = np.concatenate(
            (smoothed_data_ch1[0:n_per_event, np.newaxis], timesteps[:, np.newaxis]),
            axis=1,
        )
        peaks_ch1, _ = signal.find_peaks(
            data_ch1[:, 0], height=0.2 * data_ch1[:, 0].max()
        )
        results_w_ch1 = signal.peak_widths(data_ch1[:, 0], peaks_ch1, rel_height=0.95)
        print("#" * 10 + " Channel 1 " + "#" * 10)
        print(f"Peaks width: {results_w_ch1[1:][0]}")
        for j in range(results_w_ch1[1:][0].size):
            print(
                f"Found peaks in range: {results_w_ch1[1:][j + 1][0]}, {results_w_ch1[1:][j + 2]}"
            )
        print(
            f"Found peak at index {peaks_ch1} with amplitude: {data_ch1[:, 0][peaks_ch1]}"
        )
        if df2:
            """ Plotting second channel """
            sm_ch2 = df2.values[i].flatten()
            sm_ch2 *= -1.0
            smoothed_data_ch2 = smooth(sm_ch2, window_len=51, window="bartlett")
            data_ch2 = np.concatenate(
                (
                    smoothed_data_ch2[0:n_per_event, np.newaxis],
                    timesteps[:, np.newaxis],
                ),
                axis=1,
            )
            peaks_ch2, _ = signal.find_peaks(
                data_ch2[:, 0], height=0.2 * data_ch2[:, 0].max()
            )
            results_w_ch2 = signal.peak_widths(
                data_ch2[:, 0], peaks_ch2, rel_height=0.95
            )

            print("#" * 10 + " Channel 2 " + "#" * 10)
            print(f"Peaks width: {results_w_ch2[1:][0]}")
            for j in range(results_w_ch2[1:][0].size):
                print(
                    f"Found peaks in range: {results_w_ch2[1:][j + 1]}, {results_w_ch2[1:][j + 2]}"
                )

            print(
                f"Found peak at index {peaks_ch2} with amplitude: {data_ch2[:, 0][peaks_ch2]}"
            )

            fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(20, 7))

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

        else:
            _, ax = plt.subplots(nrows=1, ncols=2, figsize=(20, 7))

            ax[0].plot(timesteps, sm_ch1, "-", lw=1)
            ax[0].set_title("Raw Data, Channel 1")
            ax[1].plot(data_ch1[:, 1], data_ch1[:, 0], "-", lw=1)
            ax[1].plot(peaks_ch1, data_ch1[:, 0][peaks_ch1], "x")
            ax[1].hlines(*results_w_ch1[1:], color="C2")
            ax[1].set_title("Smoothed Data, Channel 1")

        plt.tight_layout()
        plt.show()
        i = int(input("Look at event: "))


if __name__ == "__main__":
    if len(sys.argv) == 3:
        df1 = pd.read_csv(sys.argv[1], header=None)
        df2 = pd.read_csv(sys.argv[2], header=None)
        plot_event(df1, df2)
    elif len(sys.argv) == 2:
        df1 = pd.read_csv(sys.argv[1], header=None)
        plot_event(df1, None)

