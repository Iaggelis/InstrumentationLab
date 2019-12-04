#include <ROOT/RDataFrame.hxx>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooLandau.h>
#include <RooPlot.h>
#include <RooClassFactory.h>
#include <RooDataSet.h>
#include <RooAddPdf.h>
#include <RooMsgService.h>
#include <ROOT/RVec.hxx>
#include <TH1F.h>

#include <vector>
#include <iostream>
#include <memory>
#include <algorithm> // For std::minmax_element
#include <iterator>  // For global begin() and end()
#include <TTree.h>
// Use O3 optimization level for just-in-time compilation
#pragma cling optimize(3)

// using namespace RooFit;
using Vec_t = const ROOT::VecOps::RVec<float> &;

void plothist(Vec_t ch1_data, Vec_t timesteps)
{
    auto hist_ch1 = new TH1F("h1", "test histo", ch1_data.size(), 0, 0.2);
    for (int i = 0; i < ch1_data.size(); ++i)
    {
        hist_ch1->SetBinContent(hist_ch1->GetBin(i), ch1_data[i]);
    }
    hist_ch1->Draw("L");
}

// void fit_tree(Vec_t ch1_data, Vec_t timesteps)
// {
//     std::cout << "ayyyyyyyy\n";
//     auto makeTTree = [&](Vec_t old_data, Vec_t timesteps) {
//         std::shared_ptr<TTree> tree;

//         tree = std::make_shared<TTree>("tree", "tree");
//         double v_amp, t;
//         tree->Branch("v_amp", &v_amp);
//         tree->Branch("t", &t);
//         const int points = old_data.size();
//         for (int i = 0; i < points; ++i)
//         {
//             v_amp = old_data[i];
//             t = timesteps[i];
//             tree->Fill();
//         }
//         return tree;
//     };
//     auto dataTree = makeTTree(ch1_data, timesteps);
//     auto minmax_time = std::minmax_element(timesteps.begin(), timesteps.end());
//     auto minmax_amp = std::minmax_element(ch1_data.begin(), ch1_data.end());
//     RooRealVar v_amp("v_amp", "Amplitude on Volts", *minmax_amp.first, *minmax_amp.second);
//     RooRealVar t("t", "time(ns)", *minmax_time.first, *minmax_time.second);

//     auto mid = (minmax_amp.second - ch1_data.begin());
//     RooRealVar mean("mean", "mean ", timesteps[mid], *minmax_time.first, *minmax_time.second);
//     RooRealVar sigma("sigma", "mass ", 0.5);
//     RooGaussian pdf("pdf", "Gaussian PDF", t, mean, sigma);
//     // RooLandau pdf("pdf", "Landau PDF", t, mean, sigma);

//     // t.setRange("signal", 80, 110);
//     RooDataSet ds("ds", "ds", RooArgSet(t, v_amp), Import(*dataTree));

//     pdf.fitTo(ds, PrintLevel(0));
//     //        measurements.push_back(mean.getVal());
//     // cout << mean.getVal() << '\n';
//     auto xframe = t.frame();
//     ds.plotOnXY(xframe, YVar(v_amp));
//     pdf.plotOn(xframe, DrawOption(""));
//     xframe->Draw("SAME");
// }
