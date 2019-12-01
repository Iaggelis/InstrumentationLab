// g++ analysis.cpp -o analysis.exe `root-config --cflags --ldflags --glibs` -lMinuit -lRooFit -lRooFitCore -ltbb
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
// #include <TStyle.h>

#include <vector>
#include <iostream>
#include <memory>
// #include <tbb/concurrent_vector.h>

using namespace std;
using RDF = ROOT::RDataFrame;
using namespace RooFit;

shared_ptr<TTree> makeTTree(const vector<float> &old_data, const vector<float> &timesteps);

int main(int argc, char **argv)
{
    TApplication theApp("App", &argc, argv);
    // RooMsgService::instance().setStreamStatus(1, false);
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    // RDF df("channels", "Labs/play_data/testing1.root");
    RDF df("channels", "testing1.root");

    // for (const auto name : df.GetColumnNames())
    // {
    //     cout << name << '\n';
    // }
    // static int i = 0;
    // auto printing = [](vector<float> ch1) {
    //     cout << ch1.size() << " at " << i++ << '\n';
    //     return true; };
    // df.Foreach(printing, {"ch1"});

    // tbb::concurrent_vector<RooRealVar> measurements; // In question!!!
    vector<double> measurements;
    auto inverter = [](vector<float> ch1_data) {
        for (auto &point : ch1_data)
        {
            point *= (-1);
        }
        return ch1_data;
    };

    auto plothist = [](const vector<float> &ch1_data) {
        auto h1 = new TH1F("h1", "test histo", ch1_data.size(), 0, 0.2);
        // for (auto &point : ch1_data)
        for (int i = 0; i < ch1_data.size(); ++i)
        {

            // h1->GetBin(i);
            h1->SetBinContent(h1->GetBin(i), ch1_data[i]);
            // h1->Fill(point);
        }
        h1->Draw("L");
    };

    auto fitting = [&](const vector<float> &ch1_data) {
        const int points = ch1_data.size();

        auto hist_ch1 = new TH1D("hist_ch1", "ch1 Fit", points, 0.001, 0.3);
        // RooDataSet data_un("data_un", "Unbinned ML", RooArgSet(x));
        for (int i = 0; i < points; ++i)
        {
            hist_ch1->SetBinContent(hist_ch1->GetBin(i), ch1_data[i]);
            // hist_ch1->Fill(ch1_data[i]);
        }
        RooRealVar x("x", "m_{Z} (GeV)", 0.002, 0.3);
        RooRealVar mean("mean", "mean Z mass", 0.15, 0.05, 0.25);
        RooRealVar sigma("sigma", "mass sigma", 0.01, -0.1, 0.1);
        // RooGaussian gauss("gauss", "Gaussian PDF", x, mean, sigma);
        RooLandau lan("lan", "Landau PDF", x, mean, sigma);
        RooDataHist dh("dh", "data", x, hist_ch1);
        // RooAddPdf pdf("pdf", "ayy", lan);
        // x.setRange("signal", 0.09, 0.4);

        // auto result_b5 = gauss.fitTo(data_un, Range("signal"),SumW2Error(false), PrintLevel(0));
        auto result_b5 = lan.fitTo(dh, Range(0.09, 0.4, false), SumW2Error(false), PrintLevel(0));
        // auto result_b5 = pdf.fitTo(dh, Range("signal"), SumW2Error(false), PrintLevel(0));

        // hist_ch1->Draw("");

        auto xframe = x.frame();
        dh.plotOn(xframe, DataError(RooAbsData::None));
        lan.plotOn(xframe);
        xframe->Draw("SAME");
    };

    auto fit_tree = [&](const vector<float> &ch1_data, const vector<float> &timesteps) {
        auto dataTree = makeTTree(ch1_data, timesteps);

        RooRealVar x("x", "points", 0.001, 1.0);
        RooRealVar t("t", "time(ns)", 75, 120);

        RooRealVar mean("mean", "mean ", 90, 75, 120);
        RooRealVar sigma("sigma", "mass ", 3.5, 2.0, 4.0);
        RooGaussian pdf("pdf", "Gaussian PDF", t, mean, sigma);
        // RooLandau pdf("pdf", "Landau PDF", t, mean, sigma);

        // t.setRange("signal", 80, 110);
        RooDataSet ds("ds", "ds", RooArgSet(t, x), Import(*dataTree));

        pdf.fitTo(ds, Range(75, 120, false), PrintLevel(-1000));
        // auto intFit = pdf.createIntegral(RooArgSet(t));
        measurements.push_back(mean.getVal());
        // cout << mean.getVal() << '\n';
        // auto xframe = t.frame();
        // ds.plotOnXY(xframe, YVar(x));
        // pdf.plotOn(xframe, DrawOption(""));
        // xframe->Draw("SAME");
    };
    auto df_02 = df.Define("ch1_inv", inverter, {"ch1"});

    // df_02.Foreach(plothist, {"ch1_inv"});
    // df_02.Foreach(fitting, {"ch1_inv"});
    df_02.Foreach(fit_tree, {"ch1_inv", "timesteps"});

    auto hist2 = new TH1D("hist2", "times", 50, 96, 100);
    for (const auto mes : measurements)
    {
        hist2->Fill(mes);
    }
    hist2->Draw("");

    // df.Foreach(fitting, {"ch1"});
    cout << "Program done!" << '\n';
    theApp.Run(true);
    return 0;
}
shared_ptr<TTree> makeTTree(const vector<float> &old_data, const vector<float> &timesteps)
{
    // TTree *tree = new TTree("tree", "tree");
    shared_ptr<TTree> tree;
    tree = make_shared<TTree>("tree", "tree");
    double *x = new double;
    double *t = new double;

    tree->Branch("x", x);
    tree->Branch("t", t);
    const int points = old_data.size();
    for (int i = 0; i < points; ++i)
    {
        *x = old_data[i];
        *t = timesteps[i];
        tree->Fill();
    }
    return tree;
};