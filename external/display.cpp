#include <TTree.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <exception>
#include <iomanip>

void display(){
	
	Float_t new_v;
	int nlin=0;
	char line[200];
    int ev1=0; 
    float n1=0,n2=0,n3=0;
    float npoints, ppevt;
	int nlines=0;
	float intgrl1=0;
	float cut1, cut2;
	
	vector<vector<float>> point, point2;
	vector <float> evtpnts;
	vector <float> hlp1, hlp2;

		
	TH1F *signal1 = new TH1F("signal1", ";time(.)",10000,0, 10000);
	TH1F *signal2 = new TH1F("signal2", ";time(.)",10000,0, 10000);
	
    
    FILE *fp  = fopen ("3052B-01_1_1_113001_1101_191211-195023.dat","r");
    FILE *fp2 = fopen ("3052B-01_2_2_113002_1102_191211-195023.dat","r");
//	FILE *fp  = fopen  ("3052B-1.dat","r");
//  FILE *fp2 = fopen  ("3052B-2.dat","r");

    while (fgets(&line,50,fp)) { ev1++;}
	cout<<"all #lines ="<< ev1<<endl;
    fclose(fp);
    fclose(fp2);
    
    ev1=40013;
    
    FILE *fp  = fopen ("3052B-01_1_1_113001_1101_191211-195023.dat","r");
    FILE *fp2 = fopen ("3052B-01_2_2_113002_1102_191211-195023.dat","r");
    //	FILE *fp  = fopen  ("3052B-1.dat","r");
    //    FILE *fp2 = fopen  ("3052B-2.dat","r");

	//  Read  the header  
	for (int i=0; i<7; i++) {  fgets(&line,50,fp); cout << line << endl; }
	fgets(&line,50,fp); sscanf(line,"%f",&n1);  cout << n1<<"  "<< line << endl;
	fgets(&line,50,fp); sscanf(line,"%f",&n2);  cout << n2<<"  "<< line <<endl;
	fgets(&line,50,fp); sscanf(line,"%f",&n3);  cout << n3<<"  "<< line << endl;
	for (int i=0; i<3; i++) {  fgets(&line,50,fp); }
	
	for (int i=0; i<13; i++) {  fgets(&line,50,fp2); }
	
	npoints= (n1*10/n3);
	nlines=(int)npoints;
	ppevt= (n1*10/n3);	
	cout << "ponits/signals= "<<npoints << endl;
	cout << "  # events    = "<<(ev1-13)/npoints << endl;
//////////////////////*******************************************

	for (int j=0; j<(ev1-13)/npoints ;j++) 
	{
		for (int i=0; i<nlines; i++) 
		{
			fgets(&line,50,fp);  sscanf(line,"%E",&n1); evtpnts.push_back(-1*n1);
		}
		point.push_back(evtpnts);
		evtpnts=vector<float>();
	}

	for (int j=0; j<(ev1-13)/npoints ;j++) 
	{
		for (int i=0; i<nlines; i++) 
		{
			fgets(&line,50,fp2); sscanf(line,"%E",&n1); evtpnts.push_back(-1*n1);
		}
		point2.push_back(evtpnts);
		evtpnts=vector<float>();
	}
	
	fclose(fp);
    fclose(fp2);

////////////////////        Select event	
	TCanvas *c1 = new TCanvas("c1", "CHANELS 1+2", 100, 0, 600, 500);
	int evt_dspl=0;
	int w=0;	

	while(w==0){	
	
    signal1->Reset();signal2->Reset();
	hlp1=point [evt_dspl];
	hlp2=point2[evt_dspl];

	cut1=0;	for (int i=0; i<100 ; i++) cut1+=hlp1.at(i); cut1=cut1/100;
	cut2=0; for (int i=0; i<100 ; i++) cut2+=hlp2.at(i); cut2=cut2/100;
	
	for (int i=0; i<hlp1.size() ; i++) 
	{
		signal1->SetBinContent(i+1, hlp1.at(i)-cut1 );
		signal2->SetBinContent(i+1, hlp2.at(i)-cut2 );
	}
///////////////  Display	
	c1->cd();
    signal1->SetAxisRange(0,npoints,"X");
    signal1->SetAxisRange( -0.5, 2,"Y");
//    signal1->SetAxisRange( -0.03,1,"Y");
    signal1->Draw();
	signal2->SetLineColor(2);
	signal2->Draw("SAME");
	c1 -> Update();		

	cout <<"Event # "<< evt_dspl << endl;


			cout << "Type Event number to display or negative to quit \n>" ;
			//cin.get( line );
			scanf("%d",&evt_dspl);
			if(evt_dspl<0) break;	
			w=0;
	}	// wait loop 
	
}
