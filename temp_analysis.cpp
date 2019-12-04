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
    // RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    int nworkers = 1;
    // if (nworkers != 1)
    // {
    //     ROOT::EnableImplicitMT(nworkers);
    // }
    // RDF df("channels", "testing1.root");
    RDF df("subrange", "clean_data.root");

    // for (const auto name : df.GetColumnNames())
    // {
    //     cout << name << '\n';
    // }

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

        RooRealVar v_amp("v_amp", "Amplitude on Volts", *minmax_amp.first, *minmax_amp.second);
        RooRealVar t("t", "time(ns)", *minmax_time.first, *minmax_time.second);

        auto mid = (minmax_amp.second - ch1_data.begin());
        // RooRealVar mean("mean", "mean ", timesteps[mid], *minmax_time.first, *minmax_time.second);
        RooRealVar mean("mean", "mean ", *minmax_time.first, timesteps[mid]);

        RooRealVar sigma("sigma", "mass ", 2.5, 7.0);
        RooGaussian pdf("pdf", "Gaussian PDF", t, mean, sigma);
        // RooLandau pdf("pdf", "Landau PDF", t, mean, sigma);

        // t.setRange("signal", 80, 110);

        RooDataSet ds("ds", "ds", RooArgSet(v_amp, t), Import(*dataTree));

        pdf.fitTo(ds, PrintLevel(2), Range(*minmax_time.first, timesteps[mid]));
        measurements.push_back(mean.getVal());
        // cout << mean.getVal() << '\n';
        auto xframe = t.frame();
        ds.plotOnXY(xframe, YVar(v_amp));
        pdf.plotOn(xframe, ProjWData(v_amp, ds));
        xframe->Draw("SAME");
    };

    auto fitting = [&](unsigned int slot, const vector<double> &ch1_data, const vector<double> &timesteps) {
        const int points = ch1_data.size();

        auto hist_ch1 = new TH1D("hist_ch1", "ch1 Gauss Fit", points, 0.001, 0.3);
        // auto minmax_time = minmax_element(timesteps.begin(), timesteps.end());
        auto minmax_amp = minmax_element(ch1_data.begin(), ch1_data.end());
        auto maxElementIndex = max_element(ch1_data.begin(), ch1_data.end()) - ch1_data.begin();
        cout << points << '\n';
        RooRealVar t("t", "time(ns)", 0.001, 3);

        for (int i = 0; i < points; i++)
        {
            hist_ch1->SetBinContent(hist_ch1->GetBin(i), ch1_data[i]);
        }
        // t.setRange("signal", 0.09, 0.2);
        RooRealVar mean("mean", "mean ", 0.1, 0.2);
        RooRealVar sigma("sigma", "mass ", 0.05, 0.2);
        RooGaussian pdf("pdf", "Gaussian PDF", t, mean, sigma);
        RooDataHist dh("dh", "data", RooArgList(t), Import(*hist_ch1));

        auto result_b5 = pdf.fitTo(dh, SumW2Error(false), Range(0.01, 0.08), PrintLevel(3));
        // hist_ch1->Draw("L");
        auto xframe = t.frame();
        dh.plotOn(xframe);
        pdf.plotOn(xframe);
        xframe->Draw("SAME");
    };
    // auto df_02 = df.Define("ch1_inv", inverter, {"ch1"}).Range(0, 1);
    auto df_02 = df.Range(0, 1);
    // df_02.Foreach(plothist, {"ch1_sub", "ch1_time"});
    // df_02.ForeachSlot(fit_tree, {"ch2_sub", "ch2_time"});
    df_02.ForeachSlot(fitting, {"ch2_sub", "ch2_time"});

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

shared_ptr<TTree> makeTTree(const vector<double> &old_data, const vector<double> &timesteps)
{
    const size_t points = old_data.size();

    // // Interpolation process:
    // vector<double> interpolated_data, interpolated_time;
    // double xi, yi;
    // auto y = old_data.data();
    // auto x = timesteps.data();

    // gsl_interp_accel *acc = gsl_interp_accel_alloc();
    // gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, 10);

    // gsl_spline_init(spline, x, y, 10);
    // for (xi = x[0]; xi < x[points - 1]; xi += 0.01)
    // {
    //     // yi = gsl_spline_eval(spline, xi, acc);
    //     interpolated_data.push_back(gsl_spline_eval(spline, xi, acc));
    //     interpolated_time.push_back(xi);
    // }

    // gsl_spline_free(spline);
    // gsl_interp_accel_free(acc);
    // Done with interpolation

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