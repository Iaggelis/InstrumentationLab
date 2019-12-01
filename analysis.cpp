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
#include <TSpectrum.h>
#include <TPolyMarker.h>

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
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    int nworkers = 1;
    // if (nworkers != 1)
    // {
    //     ROOT::EnableImplicitMT(nworkers);
    // }
    RDF df("channels", "testing1.root");

    // for (const auto name : df.GetColumnNames())
    // {
    //     cout << name << '\n';
    // }

    vector<double> measurements;
    auto inverter = [](vector<float> ch1_data) {
        for (auto &point : ch1_data)
        {
            point *= (-1);
        }
        return ch1_data;
    };

    auto plothist = [](const vector<float> &ch1_data, const vector<float> &timesteps) {
        auto hist_ch1 = new TH1F("h1", "test histo", ch1_data.size(), 0, 0.2);
        for (int i = 0; i < ch1_data.size(); ++i)
        {
            hist_ch1->SetBinContent(hist_ch1->GetBin(i), ch1_data[i]);
        }
        // hist_ch1->Draw("L");

        auto spec = new TSpectrum();
        const auto points = timesteps.size();
        auto d1 = new TH1D("d1", "test spectrum", points, 0, 0.2);
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

    auto fit_tree = [&](unsigned int slot, const vector<float> &ch1_data, const vector<float> &timesteps) {
        auto dataTree = makeTTree(ch1_data, timesteps);

        RooRealVar v_amp("v_amp", "Amplitude on Volts", 0.001, 1.0);
        RooRealVar t("t", "time(ns)", 75, 120);

        RooRealVar mean("mean", "mean ", 90, 75, 120);
        RooRealVar sigma("sigma", "mass ", 3.5, 2.0, 4.0);
        RooGaussian pdf("pdf", "Gaussian PDF", t, mean, sigma);
        // RooLandau pdf("pdf", "Landau PDF", t, mean, sigma);

        // t.setRange("signal", 80, 110);
        RooDataSet ds("ds", "ds", RooArgSet(t, v_amp), Import(*dataTree));

        pdf.fitTo(ds, Range(75, 120, false), PrintLevel(-1000));
        measurements.push_back(mean.getVal());
        // cout << mean.getVal() << '\n';
        // auto xframe = t.frame();
        // ds.plotOnXY(xframe, YVar(x));
        // pdf.plotOn(xframe, DrawOption(""));
        // xframe->Draw("SAME");
    };
    auto df_02 = df.Define("ch1_inv", inverter, {"ch1"}).Range(0, 1);
    df_02.Foreach(plothist, {"ch1_inv", "timesteps"});

    // auto df_02 = df.DefineSlot("ch1_inv", inverter, {"ch1"}).Range(0, 1);
    // df_02.ForeachSlot(fit_tree, {"ch1_inv", "timesteps"});

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

shared_ptr<TTree> makeTTree(const vector<float> &old_data, const vector<float> &timesteps)
{
    shared_ptr<TTree> tree;
    tree = make_shared<TTree>("tree", "tree");
    double v_amp, t;
    tree->Branch("v_amp", &v_amp);
    tree->Branch("t", &t);
    const int points = old_data.size();
    for (int i = 0; i < points; ++i)
    {
        v_amp = old_data[i];
        t = timesteps[i];
        tree->Fill();
    }
    return tree;
};