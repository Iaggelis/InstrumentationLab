// g++ temp_analysis.cpp -o temp_analysis.exe `root-config --cflags --ldflags --glibs` -lMinuit -lRooFit -lRooFitCore -lSpectrum
#include <ROOT/RDataFrame.hxx>
#include <TApplication.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooPlot.h>
#include <RooClassFactory.h>
#include <RooDataSet.h>
#include <RooAbsReal.h>
#include <RooAbsPdf.h>
#include <RooFunctorBinding.h>
#include <TSpectrum.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>

#include <vector>
#include <iostream>

using namespace std;
using namespace RooFit;
using RDF = ROOT::RDataFrame;

TTree *makeTTree(const vector<float> &data);

int main(int argc, char **argv)
{
    TApplication theApp("App", &argc, argv);
    RDF df("channels", "testing1.root");

    auto inverter = [](vector<float> ch1_data) {
        for (auto &point : ch1_data)
        {
            point *= (-1);
        }
        return ch1_data;
    };

    auto plothist = [](const vector<float> &ch1_data) {
        auto h1 = new TH1F("h1", "test histo", ch1_data.size(), 0, 0.2);
        for (int i = 0; i < ch1_data.size(); ++i)
        {
            h1->SetBinContent(h1->GetBin(i), ch1_data[i]);
        }
        h1->Draw("L");
    };

    auto fit_tree = [&](const vector<float> &ch1_data) {
        auto dataTree = makeTTree(ch1_data);

        RooRealVar x("x", " ", 0.002, 0.3);

        x.setRange("signal", 0.09, 0.2);
        RooRealVar mean("mean", "mean ", 0.1, 0.05, 0.2);
        RooRealVar sigma("sigma", "mass ", 0.01, -0.1, 0.1);
        RooGaussian gauss("gauss", "Gaussian PDF", x, mean, sigma);
        RooDataSet ds("ds", "ds", RooArgSet(x), Import(*dataTree));

        auto result_b5 = gauss.fitTo(ds, Range("signal"), SumW2Error(false), PrintLevel(0));
        auto xframe = x.frame();
        ds.plotOn(xframe);
        gauss.plotOn(xframe, DrawOption("SAME"));
        xframe->Draw("SAME");
    };

    auto fitting = [&](const vector<float> &ch1_data) {
        const int points = ch1_data.size();

        auto hist_ch1 = new TH1D("hist_ch1", "ch1 Gauss Fit", points, 0.001, 0.3);
        RooRealVar x("x", " ", 0.002, 0.3);

        for (int i = 0; i < points; ++i)
        {
            hist_ch1->SetBinContent(hist_ch1->GetBin(i), ch1_data[i]);
        }

        // auto spec = new TSpectrum();
        // auto d1 = new TH1D("d1", "", points, 0.001, 0.3);
        // double source[points];

        // for (int i = 0; i < points; i++)
        // {
        //     source[i] = hist_ch1->GetBinContent(i + 1);
        // }
        // spec->Background(source, points, 10, TSpectrum::kBackDecreasingWindow,
        //                  TSpectrum::kBackOrder8, kTRUE,
        //                  TSpectrum::kBackSmoothing5, kTRUE);
        // for (int i = 0; i < points; i++)
        // {
        //     d1->SetBinContent(i + 1, source[i]);
        // }
        // d1->SetLineColor(kRed);
        // d1->Draw("L");
        // hist_ch1->Draw("SAME");

        x.setRange("signal", 0.09, 0.2);
        RooRealVar mean("mean", "mean ", 0.1, 0.05, 0.2);
        RooRealVar sigma("sigma", "mass ", 0.01, -0.1, 0.1);
        RooGaussian gauss("gauss", "Gaussian PDF", x, mean, sigma);
        RooDataHist dh("dh", "data", RooArgList(x), d1);

        auto result_b5 = gauss.fitTo(dh, Range("signal"), SumW2Error(false), PrintLevel(0));
        d1->Draw("L");
        auto xframe = x.frame();
        dh.plotOn(xframe);
        gauss.plotOn(xframe, DrawOption("SAME"));
        xframe->Draw("SAME");
    };

    // auto gcdf =  ROOT::Math::gaussian_cdf(double x,double sigma,double x0=0);
    auto df_02 = df.Define("ch1_inv", inverter, {"ch1"}).Range(1, 2);

    // df_02.Foreach(plothist, {"ch1_inv"});
    // df_02.Foreach(fitting, {"ch1_inv"});
    df_02.Foreach(fit_tree, {"ch1_inv"});

    // df.Foreach(fitting, {"ch1"});

    cout << "Program done!" << '\n';
    theApp.Run(true);
    return 0;
}

TTree *makeTTree(const vector<float> &old_data)
{
    TTree *tree = new TTree("tree", "tree");
    vector<float> data;
    tree->Branch("x", &data);
    data = old_data;
    tree->Fill();
    return tree;
};