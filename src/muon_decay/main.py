"""
This program is for analysis of
muon decays events
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
from scipy import integrate


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
    maxs = []
    for m in(raw_data):
        maxs.append(np.max(m))

    risetimes = []
    t1s = []
    t2s = []
    evt = []
    integrals = []
    print(f'Total number of events: {n_events}')

    for i in tqdm(range(n_events)):
        # smoothed_data = smooth(raw_data[i], window_len=51, window='bartlett')
        smoothed_data = signal.savgol_filter(raw_data[i], 51, 3)
        data = np.concatenate((timesteps[:, np.newaxis],
                               smoothed_data[:n_per_event, np.newaxis]), axis=1)
        peak_range = [0.2, data[:, 1].max()]  # cut for two peaks
        peaks, properties = signal.find_peaks(data[:, 1],
                                              height=peak_range, distance=300)
        results_w = signal.peak_widths(data[:, 1], peaks, rel_height=0.95)
        if len(peaks) >= 2 and data[peaks[0], 1] > data[peaks[1], 1] and data[:, 1].max() < np.max(maxs)-0.01 and peaks[1] < 9000:
            peaks_start = [int(p) for p in (results_w[2])]
            ranged_sm_p1 = data[peaks_start[0]:peaks[0]]
            ranged_sm_p2 = data[peaks_start[1]:peaks[1]]
            temp_times = []
            """
            Calculating Integrals
            """
            i1 = integrate.simps(ranged_sm_p1[:, 1], ranged_sm_p1[:, 0])
            i2 = integrate.simps(ranged_sm_p2[:, 1], ranged_sm_p2[:, 0])
            integrals.append([i1, i2])
            """
            Fitting first peak:
            """
            t1 = np.arange(0, ranged_sm_p1[:, 0].size)
            chi2_fit = probfit.Chi2Regression(
                sigmoid,
                # t1,
                ranged_sm_p1[:, 0],
                ranged_sm_p1[:, 1],
            )

            minuit = iminuit.Minuit(
                chi2_fit,
                p0=np.max(ranged_sm_p1[:, 1]),
                p1=np.mean(ranged_sm_p1[:, 0]),
                # p1=peaks[0],
                # p1=t1.max(),
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
            temp_times.append(tz1)
            """
            Fitting second peak
            """
            t2 = np.arange(0, ranged_sm_p2[:, 0].size)
            chi2_fit_p2 = probfit.Chi2Regression(
                sigmoid,
                # t2,
                ranged_sm_p2[:, 0],
                ranged_sm_p2[:, 1],
            )
            
            minuit_p2 = iminuit.Minuit(
                chi2_fit_p2,
                #                 p0=np.max(data[b : peaks[1], 1]),
                p1=np.mean(ranged_sm_p2[:, 0]),
                # limit_p1=(t2.max()/2, t2.max()),
                p3=np.min(ranged_sm_p2[:, 1]),
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
            temp_times.append(tz2)
            risetimes.append([temp_times[0], temp_times[1]])
            evt.append(i)
            t1s.append(temp_times[0])
            t2s.append(temp_times[1])

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

    return (evt, t1s, t2s, integrals)




def main(filename):

    tree_names = ["t1", "tree"]
    for t_n in tree_names:
        try:
            df_np = uproot.open(filename)[t_n].array("channel1")
            break
        except Exception as ex:
            print(ex)

    ev, tz1, tz2, ints = analysis(df_np)
    tz1 = np.asarray(tz1)
    tz2 = np.asarray(tz2)
    times = np.concatenate((tz1[:, np.newaxis], tz2[:, np.newaxis]), axis=1)
    np.savetxt(f'{filename[:-5]}_events.log', ev, fmt='%1i')

    diffs = tz2 - tz1
    t_step_per_point = 4 * 1e-01 # needs more details, in nanosecond
    diffs *= t_step_per_point
    plotting = 2
    if plotting == 1:
        plt.subplot(131)
        plt.hist(times[:, 0])
        plt.title("First Peak")
        plt.subplot(132)
        plt.hist(times[:, 1])
        plt.title("Second Peak")
        plt.subplot(133)
        plt.hist(diffs)
        plt.title("Difference")
        plt.show()
    elif plotting == 2:
        import ROOT
        from ROOT import TH1F, TCanvas, TF1

        expo = """
        double expon(double *x, double *par)
        {
        return par[0] * TMath::Exp(-1.0 * x[0]/ par[1]);
        };
        """
        ROOT.gInterpreter.Declare(expo)
        fit_exp = TF1("fit_exp", ROOT.expon, 0, diffs.max(), 2)
        c1 = TCanvas("Histograms", "test")
        c1.Divide(2, 2)
        p1 = TH1F("p1", " First Peak", 20, 2100, 3000)
        p2 = TH1F("p2", " Second Peak", 20, 2500, 10000)
        diff = TH1F("diff", "Difference", 30, 0, diffs.max())
        # diff = TH1F("diff", "Difference", diffs.size, 0, diffs.max())
        # for i in range(diffs.size):
        #     diff.SetBinContent(i+1, diffs[i])
        for i in range(tz1.size):
            p1.Fill(tz1[i])
            p2.Fill(tz2[i])
            diff.Fill(diffs[i])
        c1.cd(1)
        p1.Draw()
        c1.cd(2)
        p2.Draw()
        c1.cd(3)
        diff.Draw()
        c1.cd(4)
        diff.Fit(fit_exp)
        c1.Draw()
        c1.Update()

        input("Finished")
    else:
        pass

   
if __name__ == "__main__":
    main(sys.argv[1])
