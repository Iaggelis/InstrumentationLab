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

using namespace std;
using RDF = ROOT::RDataFrame;
using namespace RooFit;

shared_ptr<TTree> makeTTree(const vector<double> &old_data, const vector<double> &timesteps);

int main(int argc, char **argv)
{
    TApplication theApp("App", &argc, argv);
    // RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    int nworkers = 1;
    // if (nworkers != 1)
    // {
    //     ROOT::EnableImplicitMT(nworkers);
    // }
    RDF df("subrange", "clean_data.root");

    vector<double> measurements;

    auto plothist = [](const vector<double> &ch1_data, const vector<double> &timesteps) {
        auto hist_ch1 = new TH1D("h1", "test histo", ch1_data.size(), 0, 0.2);
        for (int i = 0; i < ch1_data.size(); ++i)
        {
            hist_ch1->SetBinContent(hist_ch1->GetBin(i), ch1_data[i]);
        }
        // hist_ch1->Draw("L");

        auto spec = new TSpectrum();
        const auto points = timesteps.size();
        auto d1 = new TH1D("d1", "test spectrum", points, 0, *max_element(timesteps.begin(), timesteps.end()));
        double source[points], dest[points];
        for (int i = 0; i < points; i++)
        {
            source[i] = hist_ch1->GetBinContent(i + 1);
        }
        // For smoothing
        spec->Background(source, points, 10, TSpectrum::kBackDecreasingWindow,
                         TSpectrum::kBackOrder8, kTRUE,
                         TSpectrum::kBackSmoothing5, kTRUE);
        for (int i = 0; i < points; i++)
        {
            d1->SetBinContent(i + 1, source[i]);
        }
        //  d1->SetLineColor(kRed);
        d1->Draw("L");

        auto nfound = spec->SearchHighRes(source, dest, points, 8, 2, kTRUE, 3, kTRUE, 3);
        auto xpeaks = spec->GetPositionX();
        const auto xmin = 0;
        const auto xmax = points;

        auto d = new TH1D("d", "", points, xmin, xmax);
        double a;
        int bin;
        double fPositionX[nfound], fPositionY[nfound];

        for (int i = 0; i < nfound; i++)
        {
            const auto a = xpeaks[i];
            bin = 1 + Int_t(a + 0.5);
            fPositionX[i] = d1->GetBinCenter(bin);
            fPositionY[i] = d1->GetBinContent(bin);
        }

        auto pm = new TPolyMarker(nfound, fPositionX, fPositionY);
        d1->GetListOfFunctions()->Add(pm);
        pm->SetMarkerStyle(23);
        pm->SetMarkerColor(kRed);
        pm->SetMarkerSize(1.3);

        for (int i = 0; i < points; i++)
        {
            d->SetBinContent(i + 1, dest[i]);
        }
        d->SetLineColor(kRed);
        d->Draw("SAME");
        cout << "Found " << nfound << " peaks\n";

        // hist_ch1->Draw("SAME");
    };

    auto fit_tree = [&](unsigned int slot, const vector<double> &ch1_data, const vector<double> &timesteps) {
        auto dataTree = makeTTree(ch1_data, timesteps);
        auto minmax_time = minmax_element(timesteps.begin(), timesteps.end());
        auto minmax_amp = minmax_element(ch1_data.begin(), ch1_data.end());
        auto min_t = static_cast<int>(*minmax_time.first);
        auto max_t = static_cast<int>(*minmax_time.second);
        RooRealVar v_amp("v_amp", "Amplitude on Volts", *minmax_amp.first, *minmax_amp.second);
        RooRealVar t("t", "time(ns)", min_t, max_t);

        auto mid = (minmax_amp.second - ch1_data.begin());
        RooRealVar mean("mean", "mean ", min_t, max_t);
        t.setRange("signal", min_t + 1, max_t - 2);

        RooRealVar sigma("sigma", "mass ", 25.0, 10.0, 30.0);
        RooRealVar yield("yield", "norm ", -10.0, 150.0);
        RooGaussian pdf("pdf", "Gaussian PDF", t, mean, sigma);
        RooDataSet ds("ds", "ds", RooArgSet(t, v_amp), Import(*dataTree));
        RooAddPdf totalPDF("totalPDF", "", pdf, yield);
        totalPDF.fixAddCoefRange("signal");
        totalPDF.setNormRange("signal");

        totalPDF.fitTo(ds, PrintLevel(1), Range("signal"));
        measurements.push_back(mean.getVal());
        auto xframe = t.frame();
        ds.plotOnXY(xframe, YVar(v_amp));
        totalPDF.plotOn(xframe, ProjWData(v_amp, ds));
        // totalPDF.plotOn(xframe, Range("Full"));
        xframe->Draw("SAME");
    };

    auto fitting = [&](unsigned int slot, const vector<double> &ch1_data, const vector<double> &timesteps) {
        const int points = ch1_data.size();

        auto minmax_time = minmax_element(timesteps.begin(), timesteps.end());
        auto minmax_amp = minmax_element(ch1_data.begin(), ch1_data.end());
        auto min_t = static_cast<int>(*minmax_time.first);
        auto max_t = static_cast<int>(*minmax_time.second);
        auto hist_ch1 = new TH1D("hist_ch1", "ch1 Gauss Fit", points, min_t, max_t);
        for (int i = 0; i < ch1_data.size(); ++i)
        {
            hist_ch1->SetBinContent(hist_ch1->GetBin(i + 1), ch1_data[i]);
        }
        RooRealVar t("t", "time(ns)", min_t, max_t);
        t.setRange("signal", min_t + 5, max_t - 5);
        RooRealVar mean("mean", "mean ", min_t + 5, max_t - 5);
        RooRealVar sigma("sigma", "mass ", 25.0, 10.0, 30.0);
        RooRealVar yield("yield", "norm ", 1.2, -1.0, 2.0);
        RooGaussian pdf("pdf", "Gaussian PDF", t, mean, sigma);
        RooDataHist dh("dh", "data", RooArgList(t), Import(*hist_ch1));
        RooAddPdf totalPDF("totalPDF", "", pdf, yield);
        totalPDF.fixAddCoefRange("signal");
        totalPDF.setNormRange("signal");

        auto result_b5 = totalPDF.fitTo(dh, Range("signal"), PrintLevel(1));
        // hist_ch1->Draw("L");
        auto xframe = t.frame();
        dh.plotOn(xframe);
        totalPDF.plotOn(xframe);
        xframe->Draw("SAME");
    };
    auto df_02 = df.Range(0, 1);
    df_02.ForeachSlot(fitting, {"ch2_sub", "ch2_time"});

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

shared_ptr<TTree> makeTTree(const vector<double> &old_data, const vector<double> &timesteps)
{
    const size_t points = old_data.size();

    shared_ptr<TTree> tree;
    tree = make_shared<TTree>("tree", "tree");
    double v_amp, t;
    tree->Branch("v_amp", &v_amp);
    tree->Branch("t", &t);
    for (int i = 0; i < points; ++i)
    {
        v_amp = old_data[i];
        t = timesteps[i];
        tree->Fill();
    }
    return tree;
};