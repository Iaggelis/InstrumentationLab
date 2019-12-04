// g++ analysis.cpp -o analysis.exe `root-config --cflags --ldflags --glibs` -lMinuit -lRooFit -lRooFitCore -lgsl -lgslcblas -lm

#include <ROOT/RDataFrame.hxx>
#include <TApplication.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooLandau.h>
#include <RooPlot.h>
#include <RooClassFactory.h>
#include <RooDataSet.h>
#include <RooAddPdf.h>
#include <Math/ProbFuncMathCore.h>
#include <RooMsgService.h>
#include <TSpectrum.h>
#include <TPolyMarker.h>

#include <vector>
#include <iostream>
#include <memory>
#include <algorithm> // For std::minmax_element
#include <tuple>     // For std::tie
#include <iterator>  // For global begin() and end()
// #include <tbb/concurrent_vector.h>

// #include <gsl/gsl_errno.h>
// #include <gsl/gsl_spline.h>

using namespace std;
using RDF = ROOT::RDataFrame;
using namespace RooFit;

shared_ptr<TTree> makeTTree(const vector<double> &old_data, const vector<double> &timesteps);

int main(int argc, char **argv)
{
    TApplication theApp("App", &argc, argv);
    int nworkers = 1;
    // if (nworkers != 1)
    // {
    //     ROOT::EnableImplicitMT(nworkers);
    // }
    RDF df("channels", "smooth_data.root");

    vector<double> measurements;
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

    auto fitting = [](unsigned int slot, const vector<double> &ch1_data, const vector<double> &timesteps) {

    };

    auto df_02 = df.Range(0, 1);
    // df_02.Foreach(plothist, {"ch1_sub", "ch1_time"});
    // df_02.ForeachSlot(fit_tree, {"ch2_sub", "ch2_time"});
    df_02.ForeachSlot(fitting, {"sm_ch1"});

    // auto hist2 = new TH1D("hist2", "times", 50, 96, 100);
    // for (const auto mes : measurements)
    // {
    //     hist2->Fill(mes);
    // }
    // hist2->Draw("");

    cout << "Program done!" << '\n';
    theApp.Run(true);
    return 0;
}
