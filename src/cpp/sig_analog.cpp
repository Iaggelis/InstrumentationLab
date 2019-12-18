#include <ROOT/RDataFrame.hxx>
#include <TApplication.h>
#include <TMath.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1.h>

#include <vector>
#include <iostream>
#include <memory>
#include <algorithm> // For std::minmax_element
#include <iterator>  // For global begin() and end()

using namespace std;
using RDF = ROOT::RDataFrame;
float sigmoid(float *x, float *par)
{
    return par[0] / (1.0 + TMath::Exp(-par[2] * (x[0] - par[1]))) + par[3];
}
int main(int argc, char *argv[])
{

    TApplication theApp("App", &argc, argv);

    // int nworkers = 4;
    // if (nworkers != 1)
    // {
    // ROOT::EnableImplicitMT(nworkers);
    // }

    RDF df("t1", theApp.Argv(1));

    int mode = -1;
    cout << "Choose mode: (0 for per event, 1 for total analysis)\n";
    cin >> mode;
    cout << '\n';

    auto plothist = [&](unsigned int slot, const vector<float> &ch1_data, const vector<float> &ch2_data) {
        auto canvas = new TCanvas("canv", "Oscillator Channels", 900, 700);
        canvas->Divide(2, 1);
        /// Fitting 1st channel:
        canvas->cd(1);

        auto minmax_amp_ch1 = minmax_element(ch1_data.begin(), ch1_data.end());
        auto max_v_ch1 = static_cast<float>(*minmax_amp_ch1.second);

        auto hist_ch1 = new TH1F("h1", "1st Channel", ch1_data.size(), 0, ch1_data.size());
        hist_ch1->Reset();
        for (int i = 0; i < ch1_data.size(); ++i)
        {
            hist_ch1->SetBinContent(hist_ch1->GetBin(i + 1), -ch1_data[i]);
        }
        hist_ch1->DrawClone();
        /// Fitting 2nd channel:
        canvas->cd(2);

        auto minmax_amp_ch2 = minmax_element(ch2_data.begin(), ch2_data.end());
        auto max_v_ch2 = static_cast<float>(*minmax_amp_ch2.second);

        auto hist_ch2 = new TH1F("h2", "2nd Channel", ch2_data.size(), 0, ch2_data.size());
        hist_ch2->Reset();
        for (int i = 0; i < ch2_data.size(); ++i)
        {
            hist_ch2->SetBinContent(hist_ch2->GetBin(i + 1), -ch2_data[i]);
        }
        hist_ch2->DrawClone("");

        cout << "######## Event Done. #############\n";
        cout << "Quit ROOT to continue.\n";
        theApp.Run(true);
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
            df_02.ForeachSlot(plothist, {"channel1", "channel2"});
            theApp.Terminate();
        }
    }

    cout << "Program done!" << '\n';
    theApp.Run(true);
    return 0;
}
