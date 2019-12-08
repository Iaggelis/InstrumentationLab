// g++ analysis.cpp -o analysis.exe `root-config --cflags --ldflags --glibs` -lMinuit -lRooFit -lRooFitCore -lgsl -lgslcblas -lm

#include <ROOT/RDataFrame.hxx>
#include <TApplication.h>

#include <TF1.h>

#include <vector>
#include <iostream>
#include <memory>
#include <algorithm> // For std::minmax_element
// #include <tuple>     // For std::tie
#include <iterator> // For global begin() and end()

using namespace std;
using RDF = ROOT::RDataFrame;

template <typename Input1, typename Input2, typename BinaryOperation>
void zip(Input1 b1, Input1 e1, Input2 b2, BinaryOperation binOp)
{
    while (b1 != e1)
        binOp(*b1++, *b2++);
}

int main(int argc, char **argv)
{
    TApplication theApp("App", &argc, argv);
    // int nworkers = 4;
    // if (nworkers != 1)
    // {
    //     ROOT::EnableImplicitMT(nworkers);
    // }
    // RDF df("channels", "smooth_data.root");
    RDF df("subrange", "clean_data.root");

    vector<double> mean_ch1, mean_ch2;
    auto inverter = [](vector<double> ch1_data) {
        for (auto &point : ch1_data)
        {
            point *= (-1);
        }
        return ch1_data;
    };

    auto plothist = [](const vector<double> &ch1_data, const vector<double> &timesteps) {
        auto hist_ch1 = new TH1D("h1", "test histo", ch1_data.size(), 0, 0.2);
        for (int i = 0; i < ch1_data.size(); ++i)
        {
            hist_ch1->SetBinContent(hist_ch1->GetBin(i), ch1_data[i]);
        }
        hist_ch1->Draw("L");
    };

    auto fitting1 = [&mean_ch1](unsigned int slot, const vector<double> &ch1_data, const vector<double> &timesteps) {
        auto minmax_time = minmax_element(timesteps.begin(), timesteps.end());
        auto minmax_amp = minmax_element(ch1_data.begin(), ch1_data.end());
        auto min_t = static_cast<int>(*minmax_time.first);
        auto max_t = static_cast<int>(*minmax_time.second);
        TH1D hist_ch1("h1", "test histo", ch1_data.size(), min_t, max_t);
        for (int i = 0; i < ch1_data.size(); ++i)
        {
            hist_ch1.SetBinContent(hist_ch1.GetBin(i + 1), ch1_data[i]);
        }
        auto g1 = new TF1("g1", "gaus", min_t + 5, max_t - 5);

        hist_ch1.Fit(g1, "RQ0");

        auto fit_result = hist_ch1.GetFunction("g1");
        const auto norm = fit_result->GetParameter(0);
        const auto mean = fit_result->GetParameter(1);
        const auto sigma = fit_result->GetParameter(2);
        mean_ch1.push_back(mean);
    };

    auto fitting2 = [&mean_ch2](unsigned int slot, const vector<double> &ch1_data, const vector<double> &timesteps) {
        auto minmax_time = minmax_element(timesteps.begin(), timesteps.end());
        auto minmax_amp = minmax_element(ch1_data.begin(), ch1_data.end());
        auto min_t = static_cast<int>(*minmax_time.first);
        auto max_t = static_cast<int>(*minmax_time.second);
        TH1D hist_ch1("h2", "test histo", ch1_data.size(), min_t, max_t);
        for (int i = 0; i < ch1_data.size(); ++i)
        {
            hist_ch1.SetBinContent(hist_ch1.GetBin(i + 1), ch1_data[i]);
        }
        auto g1 = new TF1("g1", "gaus", min_t + 5, max_t - 5);

        hist_ch1.Fit(g1, "RQ)");
        // hist_ch1.DrawClone();

        auto fit_result = hist_ch1.GetFunction("g1");
        const auto norm = fit_result->GetParameter(0);
        const auto mean = fit_result->GetParameter(1);
        const auto sigma = fit_result->GetParameter(2);
        mean_ch2.push_back(mean);
    };

    auto df_02 = df.Range(0, 1);
    df_02.ForeachSlot(fitting1, {"ch1_sub", "ch1_time"});
    df_02.ForeachSlot(fitting2, {"ch2_sub", "ch2_time"});
    auto hist2 = new TH1D("hist2", "times", 20, 0, 400);
    for (int i = 0; i < mean_ch1.size(); i++)
    {
        hist2->Fill(fabs(mean_ch1[i] - mean_ch2[i]));
    }
    hist2->Draw("");

    



    cout << "Program done!" << '\n';
    theApp.Run(true);
    return 0;
}
