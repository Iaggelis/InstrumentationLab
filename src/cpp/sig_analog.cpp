#include <ROOT/RDataFrame.hxx>
#include <TApplication.h>
#include <TMath.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <vector>
#include <iostream>
#include <memory>
#include <algorithm> // For std::minmax_element
#include <iterator>  // For global begin() and end()

using namespace std;
using RDF = ROOT::RDataFrame;
double sigmoid(double *x, double *par)
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
    // RDF df("channels", "smooth_data.root");
    // RDF df("channels", "whole_data.root");

    RDF df("subrange", theApp.Argv(1));
    RDF df_whole("channels", theApp.Argv(1));

    cout << "Program done!" << '\n';
    theApp.Run(true);
    return 0;
}
