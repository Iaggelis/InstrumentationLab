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
    temp1, temp2 = df1.values, df2.values
    i = int(input("Look at event: "))
    while i != -1 and i < 1000:
        fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(20, 7))
        timesteps = np.arange(0, 1000, 1)
        """ Plotting first channel """
        sm_ch1 = temp1[i].flatten()
        sm_ch1 *= -1.0
        smoothed_data_ch1 = smooth(sm_ch1, window_len=51, window="bartlett")
        data_ch1 = np.concatenate(
            (smoothed_data_ch1[0:1000, np.newaxis], timesteps[:, np.newaxis]), axis=1
        )
        peaks_ch1, _ = signal.find_peaks(
            data_ch1[:, 0], height=0.2 * data_ch1[:, 0].max()
        )
        results_w_ch1 = signal.peak_widths(data_ch1[:, 0], peaks_ch1, rel_height=0.95)
        ax[0, 0].plot(timesteps, sm_ch1, "-", lw=1)
        ax[0, 0].set_title("Raw Data, Channel 1")
        ax[0, 1].plot(data_ch1[:, 1], data_ch1[:, 0], "-", lw=1)
        ax[0, 1].plot(peaks_ch1, data_ch1[:, 0][peaks_ch1], "x")
        ax[0, 1].hlines(*results_w_ch1[1:], color="C2")
        ax[0, 1].set_title("Smoothed Data, Channel 1")
        print("#" * 10 + " Channel 1 " + "#" * 10)
        print("Peaks width: {wd}".format(wd=results_w_ch1[1:][0]))
        for i in range(results_w_ch1[1:][0].size):
            print(
                "Found peaks in range: {rr1}, {rr2}".format(
                    rr1=results_w_ch1[1:][i + 1], rr2=results_w_ch1[1:][i + 2]
                )
            )
        print(
            "Found peak at index {t} with amplitude: {amp}".format(
                t=peaks_ch1, amp=data_ch1[:, 0][peaks_ch1]
            )
        )

        """ Plotting second channel """
        sm_ch2 = temp2[i].flatten()
        sm_ch2 *= -1.0
        smoothed_data_ch2 = smooth(sm_ch2, window_len=51, window="bartlett")
        data_ch2 = np.concatenate(
            (smoothed_data_ch2[0:1000, np.newaxis], timesteps[:, np.newaxis]), axis=1
        )
        peaks_ch2, _ = signal.find_peaks(
            data_ch2[:, 0], height=0.2 * data_ch2[:, 0].max()
        )
        results_w_ch2 = signal.peak_widths(data_ch2[:, 0], peaks_ch2, rel_height=0.95)
        ax[1, 0].plot(timesteps, sm_ch2, "-", lw=1)
        ax[1, 0].set_title("Raw Data, Channel 2")
        ax[1, 1].plot(data_ch2[:, 1], data_ch2[:, 0], "-", lw=1)
        ax[1, 1].plot(peaks_ch2, data_ch2[:, 0][peaks_ch2], "x")
        ax[1, 1].hlines(*results_w_ch2[1:], color="C2")
        ax[1, 1].set_title("Smoothed Data, Channel 2")
        print("#" * 10 + " Channel 2 " + "#" * 10)
        print("Peaks width: {wd}".format(wd=results_w_ch2[1:][0]))
        for i in range(results_w_ch2[1:][0].size):
            print(
                "Found peaks in range: {rr1}, {rr2}".format(
                    rr1=results_w_ch2[1:][i + 1], rr2=results_w_ch2[1:][i + 2]
                )
            )
        print(
            "Found peak at index {t} with amplitude: {amp}".format(
                t=peaks_ch2, amp=data_ch2[:, 0][peaks_ch2]
            )
        )

        plt.tight_layout()
        plt.show()
        i = int(input("Look at event: "))


if __name__ == "__main__":
    df1 = pd.read_csv("./Labs/muon_271119/Run0/ch1_values.csv", header=None)
    df2 = pd.read_csv("./Labs/muon_271119/Run0/ch2_values.csv", header=None)
    plot_event(df1, df2)
