#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <exception>
#include <iomanip>

void fit2(){
    //	gStyle->SetOptStat(0);
    gStyle->SetOptFit(1111);

	Float_t new_v;
	int nlin=0;
    char line[200];
    int ev1=0; 
    float n1,n2,n3;
    int npoints, ppevt;
	int nlines=0;
	float intgrl1=0, intgrl2=0;
	float cut1, cut2;
	float res[4];
	
	//	vector<vector<float>> *point1, *point2;
	vector <float> evtpnts;
	vector <float> *evtpnts1, *evtpnts2;
	vector <float> *hlp1, *hlp2;
	vector <float> int1, int2;

	//TF1 *f0 = new TF1("f0","landau(0)",3100, 4000);   // for the landau fit
	//	TF1 *f0 = new TF1("f0","gaus(0)+gaus(3)",-1, 1);   // for the double gaussian fit
    TF1 *f0 = new TF1("f0","[0]+[1]/(1+exp(-1*(x-[2])/[3]))" ,170, 220);   // for the FD fit
    TF1 *FD = new TF1("FD","[0]+[1]/(1+exp(-1*(x-[2])/[3]))" ,100, 300);   // for the FD Draw 
    f0->SetParName(0, "qbase");
    f0->SetParName(1, "qmax");
    f0->SetParName(2, "tFD");
    f0->SetParName(3, "tau");


	
	TH1F *signal1 = new TH1F("signal1", ";time(.)",1000,0, 1000);
	TH1F *signal2 = new TH1F("signal2", ";time(.)",10000,0, 10000);

//	TH1F *start_pnt = new TH1F("start_pnt", ";time(ponts)",200,3200, 3400 );

	TFile *f1 = new TFile("mm_oscilN.root", "READ");

	// create tree
	TTree *t1 = (TTree*)f1->Get("t1");
	t1->SetBranchAddress("points1" ,&evtpnts1);
	t1->SetBranchAddress("points2" ,&evtpnts2);

    int events1 = (int) t1->GetEntries();

	nlines=events1;
	cout << "events= "<<nlines << endl;

////////////////////        Select event	
	TCanvas *c1 = new TCanvas("c1", "CHANELS 1+2", 100, 0, 600, 500);
	int evt_dspl=45;

//	for ( int j=0; j<1 ; j++)
    {
            signal1->Reset();signal2->Reset();
       		t1->GetEntry(evt_dspl);
	    	hlp1=evtpnts1;
            hlp2=evtpnts2;
            npoints= evtpnts1->size();
			cout <<"npoints = "<< npoints<< endl;
		
			cut1=0; for (int i=0; i<100 ; i++) cut1+=hlp1->at(i); cut1=cut1/100;
			cut2=0; for (int i=0; i<100 ; i++) cut2+=hlp2->at(i); cut2=cut2/100;

            cout <<"Cut 1=  "<< cut1 << endl;
			
	for (int i=0; i<hlp1->size() ; i++) 
	{
        signal1->SetBinContent(i+1, -1*hlp1->at(i) -(-1* cut1) );
	}


	///////////////  Display	

    f0->SetParameters(0,1.9, 140, 5);
    signal1->SetAxisRange(0,npoints,"X");
    signal1->SetAxisRange( -0.1, 1.1,"Y");
    signal1->Draw();
    signal1->Fit(f0, "R");

     float a=signal1->GetMaximum();
        cout <<"Time="<<f0->GetX(a/10)  << endl;  // this is the time at the 10% of the histo max

        
        for(int j=0; j<=3; j++){
            res[j]= f0->GetParameter(j);
        //    cout<< res[j] <<endl;
        }

        FD->SetParameters(res[0], res[1], res[2],res[3]);
        FD->SetLineColor(3);
        FD->Draw("Same");
        //	  c1 -> Update();
        cout <<"Event # "<< evt_dspl << endl;
        evt_dspl++;
 
	}

    
    
}

