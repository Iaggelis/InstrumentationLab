#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <memory>
#include <math.h>

#include <TTree.h>
#include <TFile.h>
using namespace std;

struct Info
{
    string temp1;
    string temp2;
    string temp3;
    string temp4;
    string temp5; // useless text
    string temp6;
    string temp7;         // uselees again
    double timestep_in_s; // time info
    string temp9;         // also time
    double temp10;        // timestep info
    int n_events;
    string temp12;
    string temp13;
    int points;
};

shared_ptr<TTree> makeTTree(const vector<vector<double>> &data_ch1, const vector<vector<double>> &data_ch2, const size_t n_events)
{

    shared_ptr<TTree> tree;
    tree = make_shared<TTree>("tree", "tree");
    vector<double> ch1, ch2;
    tree->Branch("ch1", &ch1);
    tree->Branch("ch2", &ch2);
    const size_t points_per_event = data_ch1[0].size();

    for (size_t i = 0; i < n_events; i++)
    {
        for (size_t j = i * points_per_event; j < (i + 1) * points_per_event; j++)
        {
            ch1.push_back(data_ch1[i][j]);
            ch2.push_back(data_ch2[i][j]);
        }
        tree->Fill();
    }
    return tree;
};

int main(int argc, char const *argv[])
{
    TFile hfile("test_file.root", "RECREATE");
    vector<string> filenames;
    filenames.push_back("./Labs/muon_271119/Run1/3052B-01_1_1_113001_2010_191127-173434.dat");
    filenames.push_back("./Labs/muon_271119/Run1/3052B-01_2_2_113002_1410_191127-173434.dat");
    Info info;
    vector<vector<double>> data1, data2;
    vector<vector<vector<double>>> events;
    events.push_back(data1);
    events.push_back(data2);

    for (int s = 0; s < filenames.size(); s++)
    {
        ifstream inputf(filenames[s]);

        vector<string> info_vec;
        info_vec.reserve(13);
        if (inputf.is_open())
        {
            // loop for reading the information
            for (int i = 0; i < 13; i++)
            {
                string line;
                if (getline(inputf, line))
                {
                    info_vec.push_back(line);
                }
            }
            info.timestep_in_s = 10 * stod(info_vec[7]);
            info.temp10 = stod(info_vec[9]);
            info.n_events = stoi(info_vec[10]);
            info.points = static_cast<int>(round(info.timestep_in_s / info.temp10));
            const size_t points_per_event = info.points; // need to be know from the info part
            const size_t n_events = info.n_events;       // also from the info
            events[s].resize(n_events);
            // loop for reading the data into a 2D vector
            for (size_t i = 0; i < n_events; i++)
            {
                string line;
                for (size_t j = i * points_per_event; j < (i + 1) * points_per_event; j++)
                {
                    if (getline(inputf, line))
                    {
                        events[s][i].push_back(stod(line));
                    }
                }
            }
        }
        inputf.close();
    }

    auto myTree = makeTTree(events[0], events[1], info.n_events);
    myTree->Write();
    hfile.Close();
    return 0;
}