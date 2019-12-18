/**********************************************************************
 * Course: Computational Particle Physics  
 * Author: Sidiropoulou Ourania
 * A.E.M.: 4214
 *  
 * Description of the scipt:
 *      Coincidence between 2 scintillators
 *********************************************************************/ 
void sim_sideview(void)
{ 
	gROOT -> Reset();
	
	/**************** Variables *************************************/	
	const Int_t N = 1000;                               // number of muons
	const Double_t pi = 4*atan(1);                       // pi = 3.14	
	Int_t count = 0, cntmm1=0;                                     // number of muons which pass through the surface of the rectangle
	Int_t trigger = 0;                               // flag for coincidence
    int linecol=0;
    Double_t angle_rad;                                  // the angle of muons in radians
	Double_t muons_x_position;                           // the starting x position of the muons
    int coin1=1;
    int yes1 = 0, yes2 = 0, yes3 = 0, yes4 = 0;
    Double_t x1, y1, x2, y2, x, y, a;
    float eventsTmm1 = 0, eventsTmm2 = 0, eventsTmm3 = 0, eventsTmm4 = 0;

//    Double_t Sc1x1=0.0, Sc1x2=0.220, Sc1y1=0.820, Sc1y2=0.860;   //X
//    Double_t Sc2x1=0.700, Sc2x2=0.910, Sc2y1=0.000, Sc2y2=0.200;   //X
    Double_t Sc1x1=0.5, Sc1x2=0.7, Sc1y1=0.71, Sc1y2=0.75;   //X                
    Double_t Sc2x1=0.5, Sc2x2=0.7, Sc2y1=0.000, Sc2y2=0.200;   //X
//     Double_t Sc1x1=0.350, Sc1x2=0.650, Sc1y1=0.422, Sc1y2=0.462;    //Y
//     Double_t Sc2x1=0.350, Sc2x2=0.950, Sc2y1=0.000, Sc2y2=0.200;    //Y


    Double_t MM1x1=0.5, MM1x2=0.7, MM1y1=0.24, MM1y2=0.25;
    Double_t MM2x1=0.5, MM2x2=0.7, MM2y1=0.31, MM2y2=0.32;
    Double_t MM3x1=0.5, MM3x2=0.7, MM3y1=0.38, MM3y2=0.39;
    Double_t MM4x1=0.5, MM4x2=0.7, MM4y1=0.45, MM4y2=0.46;



    Double_t Hshift=0.00;
    Double_t S2Hshift=0.0;
    
//////////////////////////////////////////////////////////////////////////////////??????????????????
    MM1x1+=Hshift;
    MM1x2+=Hshift;
    MM2x1=MM1x1; MM2x2=MM1x2;
    
    Sc2x1+=S2Hshift;
    Sc2x2+=S2Hshift;

	TCanvas *c1 = new TCanvas("c1", "side view ",0,0,1000,1000); // (x_position, y_position, width, height)
	
	TPad *pad1 = new TPad("pad1", " ",0, 0, 1.0, 1.0);  // display tracks
	pad1->Draw();
	
//	TPad pad2 ("pad2", " ",0.5, 0, 1, 0.5);   // display plot 1
//	pad2.Draw();

//    TPad pad3 ("pad3", " ",0.5, 0.5, 1.0, 1.0);   // display plot 2
//	pad3.Draw();
//    TH1F *angular_distribution = new TH1F("phi", "phi: ang. distr. of muons (Generated muons)", 100, 0, 180);
		
	TH1F *angular_distribution = new TH1F("phi", "phi: ang. distr. of muons (triggers and mm_signal)", 100, 0, 180);
    TH1F *angular_distribution2 = new TH1F("phi2", "phi: ang. distr. of muons (signal in MM)", 100, 0, 180);
	TH1F *cos_ang_distribution = new TH1F("cos(phi)", "cos(phi)", 100, -1, 1);
//    TF1 flux("flux", "0.204*cos(x)*cos(x)-0.045", -1.5,1.5.);
//    TF1 flux("flux", "0.204*cos(x)*cos(x)-0.045", -1.14,1.14);
    
	TF1 *flux = new TF1("flux", "cos(x)*cos(x)-0.0005", -1.14, 1.14);////////////////////////////////////flux???????????????????
    
	pad1->cd();

	// Muons starting point
	TBox *rect1 = new TBox(Sc1x1, Sc1y1, Sc1x2, Sc1y2);    // (x_init, y_init) = (0, 0.99) and (x_final, y_final) = (1.0, 1.0)
	rect1 -> SetFillColor(1);
	rect1->Draw();

    TBox *rect6 = new TBox(Sc2x1, Sc2y1, Sc2x2, Sc2y2);    // (x_init, y_init) = (0, 0.99) and (x_final, y_final) = (1.0, 1.0)
    rect6 -> SetFillColor(1);
    rect6->Draw();

    // 1st scintillator (width = 0.5 (10 cm), height = 0.01 (0.2 cm))
	TBox *rect2 = new TBox(MM1x1, MM1y1, MM1x2, MM1y2);
	rect2 -> SetFillColor(5);
//	rect2.Draw();

	// 2nd scintillator (width = 0.5 (10 cm), height = 0.01 (0.2 cm))
	TBox *rect3 = new TBox(MM2x1, MM2y1, MM2x2, MM2y2);
	rect3 -> SetFillColor(5);
//	rect3.Draw()
	TBox *rect4 = new TBox(MM3x1, MM3y1, MM3x2, MM3y2);
	rect4 -> SetFillColor(5);
	TBox *rect5 = new TBox(MM4x1, MM4y1, MM4x2, MM4y2);
	rect5 -> SetFillColor(5);

    // Fill in the arrays of starting points and angular distribution

	for (Int_t i = 0; i<N; i++){

	  yes1=0;yes2=0;yes3=0;yes4=0;
// Generate random x_positions for the muons, the y_position is always set to 0.99
		muons_x_position = gRandom -> Uniform(Sc1x1,Sc1x2);
//        muons_x_position = 0.5;

// Generate random numbers which represent the angular distribution of muons
//        angle_rad        = gRandom -> Uniform(0.2, 2.9);//(0.5, 2.8);  //(30 .... -30)
        angle_rad        =  flux->GetRandom()+1.57;//(0.5, 2.8);  //(30 .... -30)          ////////////giati????
		
		x1 = muons_x_position;
		y1 = Sc1y2;
		a =tan(angle_rad);
		x2=(a*x1-y1)/a;
		y2 = 0; 
				
		TLine *line = new TLine( x1, y1, x2, y2); // line generated
//  TRIGGER
        for ( y=Sc2y1; y<=Sc2y2; y+=0.001){
            x = ( (x2 - x1) / (y2 - y1) ) * (y - y1) + x1 ;
            if ( x>=Sc2x1 && x<=Sc2x2){
					angular_distribution->Fill(angle_rad*180/pi);
					cos_ang_distribution->Fill(cos(angle_rad));
                    linecol=2;
					line->SetLineColor(linecol);
                    line->SetLineStyle(3);
                    trigger=1;
					count++;////////////autes pou kanoun trigger
                 //   line.Draw();
					break;
            }
        }


/////////////////////LILY:track an xtuph8oun oloi
  
 //       coin1=1;  // track if at least one mm chamber                 /////////////////dior8wsh ok
if(trigger==1)        {   //  Tmm1
            y=MM1y1;

            x = ( (x2 - x1) / (y2 - y1) ) * (y - y1) + x1 ;
            if ( x>=MM1x1 && x<=MM1x2){
                //linecol=4;////////////////trigger KAI hit
                //line->SetLineColor(linecol);
                //angular_distribution2->Fill(angle_rad*180/pi);
                //cntmm1++;
                //coin1=0;
		yes1 = 1;
		eventsTmm1++;
 /////////xtuph8hke o Tmm1
            }
    //if(coin1==1){
    y=MM2y1;
    x = ( (x2 - x1) / (y2 - y1) ) * (y - y1) + x1 ;
    if ( x>=MM2x1 && x<=MM2x2){
        //linecol=4;
        //line->SetLineColor(linecol);
        //angular_distribution2->Fill(angle_rad*180/pi);
        //cntmm1++;
      yes2 = 1;
      eventsTmm2++;
 //////////xtuph8hke o Tmm2
    }
    //}
    

    	y=MM3y1;
    x = ( (x2 - x1) / (y2 - y1) ) * (y - y1) + x1 ;
    if ( x>=MM3x1 && x<=MM3x2){
      yes3 = 1;
      eventsTmm3++;}
	y=MM4y1;
    x = ( (x2 - x1) / (y2 - y1) ) * (y - y1) + x1 ;
    if ( x>=MM4x1 && x<=MM4x2){
      yes4 = 1;
      eventsTmm4++;
}

	if (yes1*yes2*yes3*yes4==1)
{

linecol=3;
line->SetLineStyle(3);
line->SetLineColor(linecol);
//line->SetLineWidth(2);
angular_distribution2->Fill(angle_rad*180/pi);
cntmm1++;



 }

 }
        
        
		trigger = 0;
        if(linecol!=0)		line->Draw();
        linecol=0;

		
	}
    rect2->Draw();
    rect3->Draw();
    rect4->Draw();
    rect5->Draw();

    TCanvas *c2 = new TCanvas("c2", " ",300,0,600,550); // (x_position, y_position, width, height)
//    c2->Divide(1,2);
//    c2->cd(1);
//152
    angular_distribution->Draw();
//	c2->cd(2);
    angular_distribution2->SetLineColor(2);
    angular_distribution2->Draw("Same");
		cout << " Total muons generated = " << N << " \n Triggers = " <<
 count<< " \n signal at least to one MM="<< cntmm1;
		cout << " \n events/Triggers = " << 1.*cntmm1/(count*1.)<<endl<<"Events/chamber:"<<endl<<"Tmm1:"<<eventsTmm1<<endl<<"Tmm2:"<<eventsTmm2<<endl<<"Tmm3:"<<eventsTmm3<<endl<<"Tmm4:"<<eventsTmm4<<endl;


}

