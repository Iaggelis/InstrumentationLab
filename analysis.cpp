// g++ analysis.cpp -o analysis.exe `root-config --cflags --ldflags --glibs` -lMinuit -lRooFit -lRooFitCore -lgsl -lgslcblas -lm
#include <ROOT/RDataFrame.hxx>
#include <TApplication.h>
#include <TMath.h>

#include <TF1.h>

#include <vector>
#include <iostream>
#include <memory>
#include <algorithm> // For std::minmax_element
#include <iterator>  // For global begin() and end()

using namespace std;
using RDF = ROOT::RDataFrame;

template <typename Input1, typename Input2, typename BinaryOperation>
void zip(Input1 b1, Input1 e1, Input2 b2, BinaryOperation binOp)
{
    while (b1 != e1)
        binOp(*b1++, *b2++);
}
double sigmoid(double *x, double *par)
{
    return par[0] / (1.0 + TMath::Exp(-par[2] * (x[0] - par[1]))) + par[3];
}
double polya(double *x, double *par)
{
    /*
    par[0] is the theta parameter
    can set par[1] as mean_t and make it fixed from tf1 method
    */
    auto frac = (x[0] / par[1]);
    auto pdf = par[2] * (pow(1 + par[0], 1 + par[0]) / TMath::Gamma(1 + par[0])) * pow(frac, par[0]) * TMath::Exp(-1.0 * (1 + par[0]) * frac);
    return pdf;
}

int main(int argc, char **argv)
{
    TApplication theApp("App", &argc, argv);
    // int nworkers = 4;
    // if (nworkers != 1)
    // {
    // ROOT::EnableImplicitMT(nworkers);
    // }
    // RDF df("channels", "smooth_data.root");
    // RDF df("subrange", "clean_data.root");
    RDF df("subrange", "whole_data.root");
    // RDF df("channels", "whole_data.root");

    vector<double> mean_ch1, mean_ch2;
    auto inverter = [](vector<double> ch1_data) {
        for (auto &point : ch1_data)
        {
            point *= (-1);
        }
        return ch1_data;
    };

    auto plothist = [](unsigned int slot, const vector<double> &ch1_data, const vector<double> &timesteps) {
        auto minmax_time = minmax_element(timesteps.begin(), timesteps.end());
        auto minmax_amp = minmax_element(ch1_data.begin(), ch1_data.end());
        auto min_t = static_cast<int>(*minmax_time.first);
        auto max_t = static_cast<int>(*minmax_time.second);
        auto hist_ch1 = new TH1D("h1", "test histo", ch1_data.size(), min_t, max_t);
        for (int i = 0; i < ch1_data.size(); ++i)
        {
            hist_ch1->SetBinContent(hist_ch1->GetBin(i), ch1_data[i]);
        }
        hist_ch1->Draw("");
        auto g1 = new TF1("g1", sigmoid, min_t, max_t, 4);
        // g1->SetParameter(0, *minmax_amp.second);
        // g1->SetParameter(3, *minmax_amp.first);
        g1->SetParameter(1, max_t);
        hist_ch1->Fit(g1, "R");
    };

    auto fit = [](unsigned int slot, const vector<double> &ch1_data, const vector<double> &timesteps) {
        auto minmax_time = minmax_element(timesteps.begin(), timesteps.end());
        auto minmax_amp = minmax_element(ch1_data.begin(), ch1_data.end());
        auto min_t = static_cast<int>(*minmax_time.first);
        auto max_t = static_cast<int>(*minmax_time.second);
        auto max_v = static_cast<int>(*minmax_amp.second);
        TH1D hist_ch1("h1", "test histo", ch1_data.size(), min_t, max_t);
        for (int i = 0; i < ch1_data.size(); ++i)
        {
            hist_ch1.SetBinContent(hist_ch1.GetBin(i + 1), ch1_data[i]);
        }
        auto g1 = new TF1("g1", sigmoid, min_t, max_t, 4);
        // g1->SetParameter(0, *minmax_amp.second);
        // g1->SetParameter(3, *minmax_amp.first);
        g1->SetParameter(1, max_t);
        hist_ch1.Fit(g1, "RQ0");

        auto fit_result = hist_ch1.GetFunction("g1");
        const auto pars = fit_result->GetParameters();
        const auto tz = pars[1] - 1.0 / pars[2] * TMath::Log(pars[0] / (0.2 * max_v - pars[3]) - 1);
        return tz;
    };

    auto df_02 = df.Range(0, 1000);
    // df_02.ForeachSlot(plothist, {"ch1_sub", "ch1_time"});
    // df_02.ForeachSlot(plothist, {"ch2_sub", "ch2_time"});
    auto aug_df = df_02.DefineSlot("mean_t1", fit, {"ch1_sub", "ch1_time"}).DefineSlot("mean_t2", fit, {"ch2_sub", "ch2_time"});
    auto aug_df2 = aug_df.Define("diff", "abs(mean_t2 - mean_t1)");
    auto hist = aug_df2.Histo1D({"hist", "", 50, 0, 400}, "diff");
    hist->Draw();

    // df_02.ForeachSlot(plothist, {"ch1", "time"});

    // auto polya_fit = new TF1("fit", polya, 50, 300, 3);
    // polya_fit->SetParNames("theta", "mean_time", "constant");
    // // polya_fit->SetParameter("mean_time", mean_t);

    // hist2->Fit(polya_fit);
    cout << "Program done!" << '\n';
    theApp.Run(true);
    return 0;
}
