/**********************************************************************
 * Course: Computational Particle Physics  
 * Author: Sidiropoulou Ourania
 * A.E.M.: 4214
 *  
 * Description of the scipt:
 *      Coincidence between 2 scintillators
 *********************************************************************/ 
void sim_frontview(void)
{ 
	gROOT -> Reset();
	
	/**************** Variables *************************************/	
	const Int_t N = 10000;                               // number of muons
	const Double_t pi = 4*atan(1);                       // pi = 3.14	
	Int_t count = 0, cntmm1=0, cntmm2=0, cntmm3=0, cntmm4=0 ;  // number of muons which pass through the surface of the rectangle
    Int_t cnt_atl = 0;
	Int_t trigger = 0;                               // flag for coincidence
    int linecol=0;
    Double_t angle_rad;                                  // the angle of muons in radians
	Double_t muons_x_position;                           // the starting x position of the muons
    int coin1=1;
    int yes1 = 0, yes2 = 0, yes3 = 0, yes4 = 0;
    Double_t x1, y1, x2, y2, x, y, a;
    float eventsTmm1 = 0, eventsTmm2 = 0, eventsTmm3 = 0, eventsTmm4 = 0;
    int   clm1=0, clm2=0,clm3=0, clm4=0;
    int   event_type=0;
    
    float dfactor =0.02;   // =2/100
    
    // Y coordinates in cm
    Double_t    Sci1_Y2=38.5 * dfactor;
    Double_t    Sci1_Y1=36.5 * dfactor;
    Double_t    MM4_Y2=23.35 * dfactor;
    Double_t    MM4_Y1=22.85 * dfactor;
    Double_t    MM3_Y2=19.6 * dfactor;
    Double_t    MM3_Y1=19.1 * dfactor;
    Double_t    MM2_Y2=16.05 * dfactor;
    Double_t    MM2_Y1=15.55 * dfactor;
    Double_t    MM1_Y2=12.5 * dfactor;
    Double_t    MM1_Y1=12 * dfactor;
    Double_t    Sci2_Y2=10 * dfactor;
    Double_t    Sci2_Y1=0 * dfactor;

    // X coordinates in cm
    Double_t    Sci1_X2=29.5 * dfactor;
    Double_t    Sci1_X1=13.5* dfactor;
    Double_t    MM4_X2=24.5* dfactor;
    Double_t    MM4_X1=15.5* dfactor;
    Double_t    MM3_X2=24.5* dfactor;
    Double_t    MM3_X1=15.5* dfactor;
    Double_t    MM2_X2=24.5* dfactor;
    Double_t    MM2_X1=15.5* dfactor;
    Double_t    MM1_X2=24.5* dfactor;
    Double_t    MM1_X1=15.5* dfactor;
    Double_t    Sci2_X2=34* dfactor;
    Double_t    Sci2_X1=2* dfactor;
    
    



    Double_t Hshift= 0*dfactor;
    Double_t S2Hshift=0*dfactor;
    
//////////////////////////////////////////////////////////////////////////////////??????????????????
    MM1_X1+=Hshift;  MM2_X1+=Hshift;  MM3_X1+=Hshift;  MM4_X1+=Hshift;
    MM1_X2+=Hshift;  MM2_X2+=Hshift;  MM3_X2+=Hshift;  MM4_X2+=Hshift;
    MM2_X1=MM1_X1; MM2_X2=MM1_X2;
    
    Sci2_X1+=S2Hshift;
    Sci2_X2+=S2Hshift;

	TCanvas *c1 = new TCanvas("c1", "front view ",0,0,1000,1000); // (x_position, y_position, width, height)
	
	TPad *pad1 = new TPad("pad1", " ",0, 0, 1.0, 1.0);  // display tracks
	pad1->Draw();
	
//	TPad pad2 ("pad2", " ",0.5, 0, 1, 0.5);   // display plot 1
//	pad2.Draw();

//    TPad pad3 ("pad3", " ",0.5, 0.5, 1.0, 1.0);   // display plot 2
//	pad3.Draw();
//    TH1F *angular_distribution = new TH1F("phi", "phi: ang. distr. of muons (Generated muons)", 100, 0, 180);
		
	TH1F *angular_distribution = new TH1F("phi", "phi: ang. distr. of muons (triggers and mm_signal)", 200, 0, 180);
    TH1F *angular_distribution2 = new TH1F("phi2", "phi: ang. distr. of muons (signal in MM)", 200, 0, 180);
	TH1F *cos_ang_distribution = new TH1F("cos(phi)", "cos(phi)", 100, -1, 1);
//    TF1 flux("flux", "0.204*cos(x)*cos(x)-0.045", -1.5,1.5.);
    TH1F *h_event_type= new TH1F("h_event_type", ";event_type ",20,0,20);
	TF1 *flux = new TF1("flux", "cos(x)*cos(x)-0.0005", -1.14, 1.14);////////////////////////////////////flux???????????????????
    
	pad1->cd();

	// Muons starting point

// 1st scintillator (width = 0.5 (10 cm), height = 0.01 (0.2 cm))
    TBox *rect1 = new TBox(Sci1_X1, Sci1_Y1 , Sci1_X2, Sci1_Y2);    // (x_init, y_init) = (0, 0.99) and (x_final, y_final) = (1.0, 1.0)
	rect1 -> SetFillColor(1);
	rect1->Draw();
// 2nd scintillator (width = 0.5 (10 cm), height = 0.01 (0.2 cm))
    TBox *rect6 = new TBox(Sci2_X1, Sci2_Y1, Sci2_X2, Sci2_Y2);    // (x_init, y_init) = (0, 0.99) and (x_final, y_final) = (1.0, 1.0)
    rect6 -> SetFillColor(1);
    rect6->Draw();

	TBox *rect2 = new TBox(MM1_X1, MM1_Y1, MM1_X2, MM1_Y2);
	rect2 -> SetFillColor(5);
	TBox *rect3 = new TBox(MM2_X1, MM2_Y1, MM2_X2, MM2_Y2);
	rect3 -> SetFillColor(5);
	TBox *rect4 = new TBox(MM3_X1, MM3_Y1, MM3_X2, MM3_Y2);
	rect4 -> SetFillColor(5);
	TBox *rect5 = new TBox(MM4_X1, MM4_Y1, MM4_X2, MM4_Y2);
	rect5 -> SetFillColor(5);

    // Fill in the arrays of starting points and angular distribution

//	for (Int_t i = 0; i<N; i++)
    while(count<N)
    {
	  yes1=0;yes2=0;yes3=0;yes4=0;
// Generate random x_positions for the muons, the y_position is always set to 0.99
		muons_x_position = gRandom -> Uniform(Sci1_X1,Sci1_X2);
//      muons_x_position = 0.5;
//      Generate random numbers which represent the angular distribution of muons
//      angle_rad        = gRandom -> Uniform(0.2, 2.9);//(0.5, 2.8);  //(30 .... -30)
        angle_rad     =  flux->GetRandom()+1.57;   //(0.5, 2.8);  //(30 .... -30)          ////////////giati????
		x1 = muons_x_position;
		y1 = Sci1_Y2;
		a =tan(angle_rad);
		x2=(a*x1-y1)/a;
		y2 = 0; 
				
		TLine *line = new TLine( x1, y1, x2, y2); // line generated
//  TRIGGER
        for ( y=Sci2_Y1; y<=Sci2_Y2; y+=0.001){
            x = ( (x2 - x1) / (y2 - y1) ) * (y - y1) + x1 ;
            if ( x>=Sci2_X1 && x<=Sci2_X2){
					angular_distribution->Fill(angle_rad*180/pi);
					cos_ang_distribution->Fill(cos(angle_rad));
                    linecol=2;
                    line->SetLineColor(linecol);
                    line->SetLineWidth(1);
                    line->SetLineStyle(3);
                    trigger=1;
					count++;////////////autes pou kanoun trigger
					break;
            }
        }
/////////////////////  track an xtuph8oun oloi
//       coin1=1;  // track if at least one mm chamber                 ///////////////// dior8wsh ok

        if(trigger==1)        {   //  Tmm1
            y=MM1_Y1;
            x = ( (x2 - x1) / (y2 - y1) ) * (y - y1) + x1 ;
            if ( x>=MM1_X1 && x<=MM1_X2){
                yes1 = 1; 		eventsTmm1++;     /////////xtuph8hke o Tmm1
            }

            y=MM2_Y1;
            x = ( (x2 - x1) / (y2 - y1) ) * (y - y1) + x1 ;
            if ( x>=MM2_X1 && x<=MM2_X2){
                yes2 = 1;       eventsTmm2++;   //////////xtuph8hke o Tmm2
            }
    
    	y=MM3_Y1;
        x = ( (x2 - x1) / (y2 - y1) ) * (y - y1) + x1 ;
        if ( x>=MM3_X1 && x<=MM3_X2){
                yes3 = 1;		eventsTmm3++;
        }

        y=MM4_Y1;
        x = ( (x2 - x1) / (y2 - y1) ) * (y - y1) + x1 ;
        if ( x>=MM4_X1 && x<=MM4_X2){
            yes4 = 1;		eventsTmm4++;
        }

    if(yes1==1) clm1=1; else clm1=0;
    if(yes2==1) clm2=2; else clm2=0;
    if(yes3==1) clm3=4; else clm3=0;
    if(yes4==1) clm4=8; else clm4=0;
    event_type=clm1+clm2+clm3+clm4;
    h_event_type->Fill(event_type);
    
    
    if (yes1*yes2*yes3*yes4==1)
    {    cntmm1++;
        angular_distribution2->Fill(angle_rad*180/pi);
    }
    if (yes1==1||yes2==1||yes3==1||yes4==1)
    {
        linecol=3;
        line->SetLineStyle(3);
        line->SetLineColor(linecol);
        line->SetLineWidth(2);
        cnt_atl++;
    }

 }
		trigger = 0;
        if(linecol!=0)		line->Draw();
        linecol=0;
        if(count%1000 == 0) cout << count << endl;

    }  //while
    
    
    
    
    rect2->Draw();
    rect3->Draw();
    rect4->Draw();
    rect5->Draw();
    rect1->Draw();
    rect6->Draw();


    TCanvas *c2 = new TCanvas("c2", " ",300,0,600,550); // (x_position, y_position, width, height)
//    c2->Divide(1,2);
//    c2->cd(1);
//152
    angular_distribution->Draw();
//	c2->cd(2);
    angular_distribution2->SetLineColor(2);
    angular_distribution2->Draw("Same");
    cout << " Total muons generated = \t" << N <<endl;
    cout << " \n Triggers = \t"  << count<< endl;
    cout << " \n signal at 4 MMs=\t"<< cntmm1;
    cout << " \n events/Triggers =\t" << 1.*cntmm1/(count*1.)<<endl;
    cout << " Events/chamber:"<<endl;
    cout << "Tmm1:\t"<<eventsTmm1 <<endl;
    cout << "Tmm2:\t"<<eventsTmm2 <<endl;
    cout << "Tmm3:\t"<<eventsTmm3 <<endl;
    cout << "Tmm4:\t"<<eventsTmm4 <<endl;

    TCanvas *c3 = new TCanvas("c3", " Event type ",300,0,600,550);
    h_event_type->Draw();

}

