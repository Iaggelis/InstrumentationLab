// g++ analysis.cpp -o analysis.exe `root-config --cflags --ldflags --glibs` -lMinuit -lRooFit -lRooFitCore
#include <ROOT/RDataFrame.hxx>
#include <TApplication.h>
#include <TMath.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <vector>
#include <iostream>
#include <memory>
#include <algorithm> // For std::minmax_element
#include <iterator>  // For global begin() and end()

using namespace std;
using RDF = ROOT::RDataFrame;

double sigmoid(double *x, double *par)
{
    return par[0] / (1.0 + TMath::Exp(-par[2] * (x[0] - par[1]))) + par[3];
}

int main(int argc, char *argv[])
{

    TApplication theApp("App", &argc, argv);
    gStyle->SetOptStat(111111);

    // int nworkers = 4;
    // if (nworkers != 1)
    // {
    // ROOT::EnableImplicitMT(nworkers);
    // }
    // RDF df("channels", "smooth_data.root");
    // RDF df("channels", "whole_data.root");

    RDF df("subrange", theApp.Argv(1));
    RDF df_whole("channels", theApp.Argv(1));

    int mode = -1;
    cout << "Choose mode: (0 for per event, 1 for total analysis)\n";
    cin >> mode;
    cout << '\n';

    auto points = 0;
    auto t = df_whole.Range(0, 1);
    t.ForeachSlot([&points](unsigned int slot, const vector<float> &ch1_data) { points = ch1_data.size(); }, {"ch1"});
    auto inverter = [](vector<float> ch1_data) {
        for (auto &point : ch1_data)
        {
            point *= (-1);
        }
        return ch1_data;
    };

    auto plothist = [&](unsigned int slot, const vector<float> &ch1_data, const vector<float> &timesteps1,
                        const vector<float> &ch2_data, const vector<float> &timesteps2, const int entry) {
        auto canvas = new TCanvas("canv", "Oscillator Channels", 900, 700);
        canvas->Divide(2, 2);
        /// Fitting 1st channel:
        canvas->cd(1);

        auto minmax_time_ch1 = minmax_element(timesteps1.begin(), timesteps1.end());
        auto minmax_amp_ch1 = minmax_element(ch1_data.begin(), ch1_data.end());
        auto min_t_ch1 = static_cast<int>(*minmax_time_ch1.first);
        auto max_t_ch1 = static_cast<int>(*minmax_time_ch1.second);
        auto max_v_ch1 = static_cast<float>(*minmax_amp_ch1.second);

        auto hist_ch1 = new TH1F("h1", "1st Channel", ch1_data.size(), min_t_ch1, max_t_ch1);
        for (int i = 0; i < ch1_data.size(); ++i)
        {
            hist_ch1->SetBinContent(hist_ch1->GetBin(i + 1), ch1_data[i]);
        }
        // hist_ch1->DrawClone("");
        auto g1 = new TF1("g1", sigmoid, min_t_ch1, max_t_ch1, 4);
        // g1->SetParameter(0, *minmax_amp.second);
        // g1->SetParameter(3, *minmax_amp.first);
        g1->SetParameter(1, max_t_ch1);
        hist_ch1->Fit(g1, "R");
        auto fit_result_ch1 = hist_ch1->GetFunction("g1");
        const auto pars_ch1 = fit_result_ch1->GetParameters();
        const auto tz_ch1 = pars_ch1[1] - 1.0 / pars_ch1[2] * TMath::Log(pars_ch1[0] / (0.2 * max_v_ch1 - pars_ch1[3]) - 1);
        // cout << pars[0] << " , " << pars[1] << " , " << pars[2] << " , " << pars[3] << " , " << max_v << '\n';
        // cout << tz << '\n';

        /// Fitting 2nd channel:
        canvas->cd(2);

        auto minmax_time_ch2 = minmax_element(timesteps2.begin(), timesteps2.end());
        auto minmax_amp_ch2 = minmax_element(ch2_data.begin(), ch2_data.end());
        auto min_t_ch2 = static_cast<int>(*minmax_time_ch2.first);
        auto max_t_ch2 = static_cast<int>(*minmax_time_ch2.second);
        auto max_v_ch2 = static_cast<float>(*minmax_amp_ch2.second);

        auto hist_ch2 = new TH1F("h2", "2nd Channel", ch2_data.size(), min_t_ch2, max_t_ch2);
        for (int i = 0; i < ch2_data.size(); ++i)
        {
            hist_ch2->SetBinContent(hist_ch2->GetBin(i + 1), ch2_data[i]);
        }
        // hist_ch2->DrawClone("");
        auto g2 = new TF1("g2", sigmoid, min_t_ch2, max_t_ch2, 4);
        // g2->SetParameter(0, *minmax_amp_ch2.second);
        // g2->SetParameter(3, *minmax_amp_ch2.first);
        g2->SetParameter(1, max_t_ch2);
        hist_ch2->Fit(g2, "R");
        auto fit_result_ch2 = hist_ch2->GetFunction("g2");
        const auto pars_ch2 = fit_result_ch2->GetParameters();
        const auto tz_ch2 = pars_ch2[1] - 1.0 / pars_ch2[2] * TMath::Log(pars_ch2[0] / (0.2 * max_v_ch2 - pars_ch2[3]) - 1);
        // cout << pars[0] << " , " << pars[1] << " , " << pars[2] << " , " << pars[3] << " , " << max_v << '\n';
        // cout << tz << '\n';
        auto temp = df_whole.Range(entry, entry + 1);
        auto h3 = new TH1F("h3", "1st Channe", points, 0, points);
        auto h4 = new TH1F("h4", "2nd Channl", points, 0, points);
        auto fill_hist = [&h3, &h4](unsigned int slot, const vector<float> &ch1_data, const vector<float> &ch2_data) {
            for (int i = 0; i < ch1_data.size(); ++i)
            {
                h3->SetBinContent(h3->GetBin(i + 1), ch1_data[i]);
                h4->SetBinContent(h4->GetBin(i + 1), ch2_data[i]);
            }
        };

        temp.ForeachSlot(fill_hist, {"ch1", "ch2"});
        canvas->cd(3);
        h3->Draw();
        fit_result_ch1->Draw("same");
        canvas->cd(4);
        h4->Draw();
        fit_result_ch2->Draw("same");
        cout << "######## Event Done. #############\n";
        cout << "Quit ROOT to continue.\n";
        theApp.Run(true);
    };

    auto fit = [](unsigned int slot, const vector<float> &ch1_data) {
        auto minmax_amp = minmax_element(ch1_data.begin(), ch1_data.end());
        auto max_v = static_cast<float>(*minmax_amp.second);

        // Hist for the whole signal
        TH1F hist_ch1("h1", "test histo", ch1_data.size(), min_t, max_t);
        for (int i = 0; i < ch1_data.size(); ++i)
        {
            hist_ch1.SetBinContent(hist_ch1.GetBin(i + 1), ch1_data[i]);
        }
    //     auto g1 = new TF1("g1", sigmoid, min_t, max_t, 4);
    //     // g1->SetParameter(0, *minmax_amp.second);
    //     // g1->SetParameter(3, *minmax_amp.first);
    //     g1->SetParameter(1, max_t);
    //     hist_ch1.Fit(g1, "RQ0");

    //     auto fit_result = hist_ch1.GetFunction("g1");
    //     const auto pars = fit_result->GetParameters();
    //     const auto tz = pars[1] - 1.0 / pars[2] * TMath::Log(pars[0] / (0.2 * max_v - pars[3]) - 1);
    //     return tz;
        return 0;
    };

    auto usage = mode == 0 ? 0 : 1;

    if (usage == 0)
    {
        int event = -2;
        while (event != -1)
        {
            cout << "Enter event number: ";
            cin >> event;
            cout << event << '\n';
            auto df_02 = df.Range(event, event + 1);
            auto df_03 = df_02.DefineSlot("entry", [&event](unsigned int slot) { return event; });
            df_03.ForeachSlot(plothist, {"ch1_sub", "ch1_time", "ch2_sub", "ch2_time", "entry"});
            theApp.Terminate();
        }
    }

    if (usage == 1)
    {
        auto aug_df = df.DefineSlot("mean_t1", fit, {"ch1_sub", "ch1_time"}).DefineSlot("mean_t2", fit, {"ch2_sub", "ch2_time"});
        auto aug_df2 = aug_df.Define("diff", "abs(mean_t2 - mean_t1)");
        auto hist = aug_df2.Histo1D({"hist", "", 50, 0, 400}, "diff");
        hist->DrawClone();
        // aug_df2.Snapshot("analysis", "analysis_output.root", {"mean_t1", "mean_t2", "diff"});
    }

    // auto polya_fit = new TF1("fit", polya, 50, 300, 3);
    // polya_fit->SetParNames("theta", "mean_time", "constant");
    // // polya_fit->SetParameter("mean_time", mean_t);

    // hist2->Fit(polya_fit);
    cout << "Program done!" << '\n';
    theApp.Run(true);
    return 0;
}
