#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <Riostream.h>
#include <vector>
#include <TString.h>
#include <TH1F.h>
#include "hiss.h"
#include <TCanvas.h>
#include <TF1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TApplication.h>
#include "PMTStyle.h"
#include "PMType.h"
#include "Pedestal.h"
#include "SPEResponse.h"
#include "PMT.h"
#include "DFTmethod.h"
#include "SPEFitter.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "TROOT.h"
#include "TEnv.h"
#include "TBrowser.h"
#include "TMultiGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TPolyLine3D.h"
#include "TVirtualPad.h"
#include "Riostream.h"
#include "TVirtualFitter.h"
#include "TPluginManager.h"
#include "TClass.h"
#include "TMath.h"
#include "TSystem.h"
#include <stdlib.h>
#include "TCanvas.h"

#include "HFitInterface.h"
#include "Fit/DataRange.h"
#include "Math/MinimizerOptions.h"

#include <ctype.h>
using namespace std;
//main function
void runfile_working_first(TString path_to_M1)
{
  gROOT->Reset();
  
  PMTStyle::SetDefaultStyle();
  
//define Canvas event	
    TCanvas *c1 = new TCanvas( "c1", "" );
	//c1->SetLogy(); //set log Y
	Int_t HV[5] = {1250,1300,1350,1401,1450};
	int i=0;
	int j=0;
	int PED_LED_Count =2;
	int HV_Count = 5;
	int num =1;
//run through each folder and file	
TString PED = "/HVSCAN/%d/PED/F1--Trace--00000.txt";
TString LED = "/HVSCAN/%d/LED/F1--Trace--00000.txt";
//new TString definitions
TString HV_Value_PED;
TString HV_Value_LED;
TString Filename;
double Q ;
double sigma;
double amp;
int index;
double gainerror;
double sig_reduced;
double sig_reduced_err;
double xbar;
double xbarErr;
int dat;
//ofstream ff ("gains.txt"); // write the respective voltages and gains to a file in directory
cout << "Input 0 for Juno file analysis, 1 for HV folder analysis"<<endl;
cin>>index;
TString filename;
TString histname;
				
TString histname_LED;
TString histname_PED;
if (index == 1)
	{	
			ofstream labdata ("lab_data.txt");
				labdata <<"Fit data for voltages"<< endl;
	 			labdata <<"Mu Mu_err w w_err alpha alpha_err lambda lambda_err Theta Theta_err sig_reduced sig_reduced_err Gain Gain_err"<< endl;
	 			labdata <<" "<< endl;
			for (i=0; i<5; ++i)
				{	
				
				//define 2 strings to specify to hiss whether we are in PED or LED
				HV_Value_PED = TString ("Pedestal Run, V_supply = ") + Form("%d",HV[i]) + TString("V");
				HV_Value_LED = TString ("LED Run, V_supply = ") + Form("%d",HV[i]) + TString("V");
				//define the 2 histograms for PED and LED
				TH1F *histo_PED = hiss(path_to_M1 + Form(PED, HV[i]) , HV_Value_PED, index, 0);
				TH1F *histo_LED = hiss(path_to_M1 + Form(LED, HV[i]) , HV_Value_LED, index, 0);
				
				//cout << Form(PED, HV[i]) << endl; getchar();
				Q = histo_PED->GetMean(); //get Q initially
				sigma = histo_PED->GetRMS(); //get sigma initially
				
				histo_PED->GetXaxis()->SetTitle("charge (nVs)"); //set Xaxis title
				histo_PED->GetYaxis()->SetTitle("No. of Entries"); //set Yaxis title
				
				histo_LED->GetXaxis()->SetTitle("charge (nVs)"); //set Xaxis title
				histo_LED->GetYaxis()->SetTitle("No. of Entries)"); //set Yaxis title
				
				//define fit parameters
				TF1  *Fit_Gauss = new TF1("Fit_Gauss","gaus", (Q - 10*sigma), (Q + 10*sigma));
				Fit_Gauss->SetParameters(amp*histo_PED->GetBinWidth(1)*(1/(sqrt(2*M_PI)*sigma)),Q,sigma);
				Fit_Gauss->SetNpx(10000); 
				histo_PED->Draw("PE");
				histo_PED->Fit("Fit_Gauss","","", Q-3.0*sigma,Q+3*sigma);
				Q 		= histo_PED->GetFunction("Fit_Gauss")->GetParameter(1); //get Q from fit
				sigma 	= histo_PED->GetFunction("Fit_Gauss")->GetParameter(2); //get sigma from fit
				cout <<"Q = "<< Q <<" sigma = "<< sigma<< endl;
				//Fit_Gauss->Draw("same");
				
				c1->Update();
				c1->WaitPrimitive(); //ROOT waits until you hit ENTER
				
				Fit_Gauss->SetParameters(amp*histo_LED->GetBinWidth(1)*(1/(sqrt(2*M_PI)*sigma)),Q,sigma);
				Fit_Gauss->SetParLimits(1, Q-2.0*sigma, Q+2.0*sigma); //[1] is for Q, predefined by "gaus"	
				Fit_Gauss->SetParLimits(2, 0.5*sigma,  1.5*sigma); // [2] is for sigma, predefined by "gaus"
					
				histo_LED->Draw();
				histo_LED->Fit("Fit_Gauss","","", Q-5.0*sigma,Q+3*sigma);
				//Fit_Gauss -> Draw("same");
					
				Q 		= Fit_Gauss->GetParameter(1); //histo_LED->GetFunction("Fit_Gauss")->GetParameter(1); //get Q from new fit
				sigma 	= Fit_Gauss->GetParameter(2); //histo_LED->GetFunction("Fit_Gauss")->GetParameter(2); //get sigma from new fit
				cout << "Q_New = " << Q << " sigma_New = " << sigma << endl;
				
				double N_Tot	= histo_LED->Integral(); // Get N_Tot from LED
				cout <<"N_Tot = "<< N_Tot << endl;
				double N0 = Fit_Gauss->GetParameter(0);
				N0 *= sqrt(2*M_PI)*sigma/histo_LED->GetBinWidth(1);
				
				cout <<"N0 = "<< N0 << endl;
				double MU = -log(N0/N_Tot); //calculate Mu
				cout <<"Mu = " << MU << endl;
				
				
				c1->Update();
				c1->WaitPrimitive(); //ROOT waits until you hit ENTER
				
				/* ... */
				
				cout << "" << endl;
				cout << "" << endl;
				
				histo_LED->GetListOfFunctions()->Remove( histo_LED->GetFunction( "Fit_Gauss") );
				histo_LED->SetMarkerStyle( 20 );
				histo_LED->SetMarkerSize( 0.4 );
				histo_LED->SetLineColor( kBlack );
				histo_LED->SetMarkerColor( kBlack );
				histo_LED->SetStats(0);
				histo_LED->Draw( "" );
				
				
				Double_t _G = ( histo_LED->GetMean() - Q )/(MU); //calculated in nVs
				cout << " Esimated G : " << _G << endl;
				
				SPEFitter fit;
				Double_t p_test[4] = { 1.0/_G, 7.0, 1.0/(0.5*_G), 0.2 };
				SPEResponse gamma_test( PMType::GAMMA, p_test );
				
				Int_t nbins = histo_LED->GetNbinsX();
				Double_t xmin = histo_LED->GetXaxis()->GetBinLowEdge(1);
				Double_t xmax = histo_LED->GetXaxis()->GetBinUpEdge(nbins);
				cout << " No. of bins : " << nbins << endl;
				cout << " ( " << xmin << ", " << xmax << " ) " << endl;
				
				Double_t rms = histo_LED->GetRMS();
				Double_t rmse = histo_LED->GetRMSError();
				cout << " RMS : " << rms <<" +/- "<< rmse << endl;
				
				
				DFTmethod dft( 2.0*nbins, xmin, xmax, gamma_test );
				dft.wbin = histo_LED->GetBinWidth(1);
				dft.Norm = histo_LED->Integral();
				dft.Q0 = Q;
				dft.s0 = sigma;
				dft.mu = MU;
				
				fit.SetDFTmethod( dft );
				fit.FitwDFTmethod( histo_LED );
				
					
				dft.Norm = fit.vals[0];
				dft.Q0 = fit.vals[1];
				dft.s0 = fit.vals[2];
				dft.mu = fit.vals[3]; 
				Double_t p_fit[4] = { fit.vals[4], fit.vals[5], fit.vals[6], fit.vals[7] };
				dft.spef.SetParams( p_fit );
				TGraph *grBF = dft.GetGraph();
				grBF->Draw( "SAME,L" );
				
				Double_t Gfit = ( fit.vals[7]/fit.vals[6]+(1.0-fit.vals[7])/fit.vals[4] ); 
				cout << " Gain : " << Gfit/(50*1.60217662e-10) << endl;
				cout << " Gain (no. of PEs) : " << Gfit << endl;
				//ff <<HV[i]<<" "<<  Gfit/(50*1.60217662e-10) <<" "<< Gfit<<endl;  // write the respective voltages and gains to a file in directory
				cout << "" << endl;
				cout << "" << endl;
				sig_reduced = 1/sqrt(1 + fit.vals[5]);
				sig_reduced_err = 0.5*pow( (1+fit.vals[5]), -1.5 );
				gainerror = (fit.vals[7]/fit.vals[6])* ( sqrt( pow( (fit.errs[7]/fit.vals[7]),2 ) + pow( (fit.errs[6]/fit.vals[6]),2 ))   +   sqrt( pow( (fit.errs[7]/fit.vals[7]),2 ) + pow( (fit.errs[4]/fit.vals[4]),2 )) );
				//gaindata <<"angle Mu Mu_err w w_err alpha alpha_err lambda lambda_err Theta Theta_err sig_reduced sig_reduced_err Gain Gain_err"<< endl;
				labdata <<fit.vals[3]<<" "<<fit.errs[3]<<" "<<fit.vals[7]<<" "<<fit.errs[7]<<" "<<fit.vals[6]<<" "<<fit.errs[6]<<" "<<fit.vals[4]<<" "<<fit.errs[4]<<" "<<fit.vals[5]<<" "<<fit.errs[5]<<" "<< sig_reduced<<" "<<sig_reduced_err<<" "<<Gfit/(50*1.60217662e-10) <<" "<< gainerror/(50*1.60217662e-10)<< endl;
				
				//c1->Update();
				//c1->WaitPrimitive();
				}
	}
 else if (index == 0)
	 {
	 	
	 	cout <<"Input the data set you wish to analyze:"<<endl;
	 	cout <<"Input 248 for scan248"<<endl;
	 	cout <<"Input 846 for scan846"<<endl;
	 	cout <<"Input 3712 for scan3712"<<endl;
	 	cout <<"Input 3737 for scan3737"<<endl;
	 	cout <<"Input 3899 for scan3899"<<endl;
	 	cout <<"Input 4050 for scan4050"<<endl;
	 	cout <<"Input 4230 for scan4230"<<endl;
	 	cout <<"Input 4232 for scan4232"<<endl;
	 	
	 	cin >> dat;
	 	filename = TString("gain_data_scan") + Form("%d",dat);
	 	ofstream gaindata (filename);
	 	Int_t n = 24;
		Float_t X_ERR[n];
	 	Float_t ANGLES[n];
		/*Float_t POSITION_1[n];
	 	Float_t POSITION_1_ERR[n];
	 	Float_t POSITION_2[n];
	 	Float_t POSITION_2_ERR[n];
	 	Float_t POSITION_3[n];
	 	Float_t POSITION_3_ERR[n];
	 	Float_t POSITION_4[n];
	 	Float_t POSITION_4_ERR[n];
	 	Float_t POSITION_5[n];
	 	Float_t POSITION_5_ERR[n];
	 	Float_t POSITION_6[n];
	 	Float_t POSITION_6_ERR[n];
	 	Float_t POSITION_7[n];
	 	Float_t POSITION_7_ERR[n];*/
	 
	 	Float_t PMT_DATA[14][24];
	 
	 	TMultiGraph  *mg  = new TMultiGraph();
	 	for (int p = 1; p<8; ++p)
	 	{	
	 			gaindata <<"Fit data for position "<< p <<": "<< endl;
	 			gaindata <<"angle Mu Mu_err w w_err alpha alpha_err lambda lambda_err Theta Theta_err sig_reduced sig_reduced_err Gain Gain_err"<< endl;
	 			//gaindata <<"angle Theta Gain"<< endl;
	 			gaindata <<" "<< endl;
				for (int f=0; f<24; ++f) {ANGLES[f] = f*15;	X_ERR[f] = 0;}
				//gr_name = TString("GR_") + Form("%d", p);
			
	 			for (int a=0; a<24; ++a)
	 			    {
			 			//define 2 strings to specify to hiss whether we are in PED or LED
						histname = TString("position = ") + Form("%d",p) + TString(" angle = ")+ Form("%d", 15*a);
						histname_LED = histname + TString(" (LED)");
						histname_PED = histname + TString(" (PED)");
						HV_Value_PED = path_to_M1 +   TString("/scan") +Form ("%d",dat) + TString("_position") + Form("%d", p) + TString("_angle") + Form("%d", 15*a) + TString("_PED.txt");
						HV_Value_LED = path_to_M1 +   TString("/scan") +Form ("%d",dat) + TString("_position") + Form("%d", p) + TString("_angle") + Form("%d", 15*a) + TString("_LED.txt");
						cout<<HV_Value_LED<<endl;
						cout<<HV_Value_PED<<endl;
						//define the 2 histograms for PED and LED
						TH1F *histo_PED = hiss(HV_Value_PED, histname_PED, index, dat);
						TH1F *histo_LED = hiss(HV_Value_LED, histname_LED, index, dat);
						
						//cout << Form(PED, HV[i]) << endl; getchar();
						Q = histo_PED->GetMean(); //get Q initially
						sigma = histo_PED->GetRMS(); //get sigma initially
						
						histo_PED->GetXaxis()->SetTitle("charge (DUQ)"); //set Xaxis title
						histo_PED->GetYaxis()->SetTitle("No. of Entries"); //set Yaxis title
						
						histo_LED->GetXaxis()->SetTitle("charge (DUQ)"); //set Xaxis title
						histo_LED->GetYaxis()->SetTitle("No. of Entries)"); //set Yaxis title
						
						//define fit parameters
						TF1  *Fit_Gauss = new TF1("Fit_Gauss","gaus", (Q - 5*sigma), (Q + 3*sigma));
						Fit_Gauss->SetParameters(amp*histo_PED->GetBinWidth(1)*(1/(sqrt(2*M_PI)*sigma)),Q,sigma);
						Fit_Gauss->SetNpx(10000); 
						histo_PED->Draw("PE");
						histo_PED->Fit("Fit_Gauss","","", Q-5.0*sigma,Q+3*sigma);
						Q 		= histo_PED->GetFunction("Fit_Gauss")->GetParameter(1); //get Q from fit
						sigma 	= histo_PED->GetFunction("Fit_Gauss")->GetParameter(2); //get sigma from fit
						cout <<"Q = "<< Q <<" sigma = "<< sigma<< endl;
						Fit_Gauss->Draw("same");
						
						c1->Update();
						c1->WaitPrimitive(); //ROOT waits until you hit ENTER
					
						/*TString PdfName;
						PdfName = TString("scan") + Form("%d", dat) + TString("_Results.pdf");*/
						c1->Print("test.pdf(","pdf");
						
						Fit_Gauss->SetParameters(amp*histo_LED->GetBinWidth(1)*(1/(sqrt(2*M_PI)*sigma)),Q,sigma);
						Fit_Gauss->SetParLimits(1, Q-2.0*sigma, Q+2.0*sigma); //[1] is for Q, predefined by "gaus"	
						Fit_Gauss->SetParLimits(2, 0.5*sigma,  1.5*sigma); // [2] is for sigma, predefined by "gaus"
							
						histo_LED->Draw();
						histo_LED->Fit("Fit_Gauss","","", Q-5.0*sigma,Q+3*sigma);
						Fit_Gauss -> Draw("same");
							
						Q 		= Fit_Gauss->GetParameter(1); //histo_LED->GetFunction("Fit_Gauss")->GetParameter(1); //get Q from new fit
						sigma 	= Fit_Gauss->GetParameter(2); //histo_LED->GetFunction("Fit_Gauss")->GetParameter(2); //get sigma from new fit
						cout << "Q_New = " << Q << " sigma_New = " << sigma << endl;
						
						double N_Tot	= histo_LED->Integral(); // Get N_Tot from LED
						cout <<"N_Tot = "<< N_Tot << endl;
						double N0 = Fit_Gauss->GetParameter(0);
						N0 *= sqrt(2*M_PI)*sigma/histo_LED->GetBinWidth(1);
						
						cout <<"N0 = "<< N0 << endl;
						double MU = -log(N0/N_Tot); //calculate Mu
						cout <<"Mu = " << MU << endl;
						
						
						c1->Update();
						c1->WaitPrimitive(); //ROOT waits until you hit ENTER
						c1->Print("test.pdf","pdf");
						
						/* ... */
						
						cout << "" << endl;
						cout << "" << endl;
						
						histo_LED->GetListOfFunctions()->Remove( histo_LED->GetFunction( "Fit_Gauss") );
						histo_LED->SetMarkerStyle( 20 );
						histo_LED->SetMarkerSize( 0.4 );
						histo_LED->SetLineColor( kBlack );
						histo_LED->SetMarkerColor( kBlack );
						histo_LED->SetStats(0);
						histo_LED->Draw( "" );
						
						
						Double_t _G = ( histo_LED->GetMean() - Q )/(MU); //calculated in nVs
						cout << " Esimated G : " << _G << endl;
						
						SPEFitter fit;
						Double_t p_test[4] = { 1.0/_G, 7.0, 1.0/(0.5*_G), 0.2 };
						SPEResponse gamma_test( PMType::GAMMA, p_test );
						
						Int_t nbins = histo_LED->GetNbinsX();
						Double_t xmin = histo_LED->GetXaxis()->GetBinLowEdge(1);
						Double_t xmax = histo_LED->GetXaxis()->GetBinUpEdge(nbins);
						cout << " No. of bins : " << nbins << endl;
						cout << " ( " << xmin << ", " << xmax << " ) " << endl;
						
						Double_t rms = histo_LED->GetRMS();
						Double_t rmse = histo_LED->GetRMSError();
						cout << " RMS : " << rms <<" +/- "<< rmse << endl;
						
						
						DFTmethod dft( 2.0*nbins, xmin, xmax, gamma_test );
						dft.wbin = histo_LED->GetBinWidth(1);
						dft.Norm = histo_LED->Integral();
						dft.Q0 = Q;
						dft.s0 = sigma;
						dft.mu = MU;
						
						fit.SetDFTmethod( dft );
						fit.FitwDFTmethod( histo_LED );
						
							
						dft.Norm = fit.vals[0];
						dft.Q0 = fit.vals[1];
						dft.s0 = fit.vals[2];
						dft.mu = fit.vals[3]; 
						Double_t p_fit[4] = { fit.vals[4], fit.vals[5], fit.vals[6], fit.vals[7] };
						dft.spef.SetParams( p_fit );
						TGraph *grBF = dft.GetGraph();
						grBF->PaintStats(0);
						grBF->Draw( "SAME,L" );
						
						Double_t Gfit = ( fit.vals[7]/fit.vals[6]+(1.0-fit.vals[7])/fit.vals[4] ); 
						
						//ff <<HV[i]<<" "<<  Gfit/(50*1.60217662e-10) <<" "<< Gfit<<endl;  // write the respective voltages and gains to a file in directory
						cout << "" << endl;
						cout << "" << endl;
						//gaindata <<"angle xbar Q Mu w Theta Gain"<< endl;
						//gaindata << a*15 <<" "<< histo_LED->GetMean() <<" "<< Q <<" "<< MU <<" "<<fit.vals[7]<<" "<<fit.vals[5]<<" "<< Gfit << endl;
						// xbar = histo_LED->GetMean();
						// xbarErr = histo_LED->GetMeanError();
						// gainerror = Gfit*
						sig_reduced = 1/sqrt(1 + fit.vals[5]);
						sig_reduced_err = 0.5*pow( (1+fit.vals[5]), -1.5 );
						gainerror = (fit.vals[7]/fit.vals[6])* ( sqrt( pow( (fit.errs[7]/fit.vals[7]),2 ) + pow( (fit.errs[6]/fit.vals[6]),2 ))   +   sqrt( pow( (fit.errs[7]/fit.vals[7]),2 ) + pow( (fit.errs[4]/fit.vals[4]),2 )) );
						//gaindata <<"angle Mu Mu_err w w_err alpha alpha_err lambda lambda_err Theta Theta_err sig_reduced sig_reduced_err Gain Gain_err"<< endl;
						gaindata << a*15 <<" "<<fit.vals[3]<<" "<<fit.errs[3]<<" "<<fit.vals[7]<<" "<<fit.errs[7]<<" "<<fit.vals[6]<<" "<<fit.errs[6]<<" "<<fit.vals[4]<<" "<<fit.errs[4]<<" "<<fit.vals[5]<<" "<<fit.errs[5]<<" "<< sig_reduced<<" "<<sig_reduced_err<<" "<<Gfit <<" "<< gainerror<< endl;
						cout << " Gain (no. of PEs) : " << Gfit <<" +/- "<< gainerror << endl;
						/*if (p==1) {POSITION_1[a] = Gfit;	POSITION_1_ERR = gainerror;	}
						if (p==2) {POSITION_2[a] = Gfit;	POSITION_2_ERR = gainerror;}
						if (p==3) {POSITION_3[a] = Gfit;	POSITION_3_ERR = gainerror;}
						if (p==4) {POSITION_4[a] = Gfit;	POSITION_4_ERR = gainerror;}
						if (p==5) {POSITION_5[a] = Gfit;	POSITION_5_ERR = gainerror;}
						if (p==6) {POSITION_6[a] = Gfit;	POSITION_6_ERR = gainerror;}
						if (p==7) {POSITION_7[a] = Gfit;	POSITION_7_ERR = gainerror;}*/
						
						if ((fit.chi2r <= 3) && (fit.fit_status == 0))	{	PMT_DATA[2*p-2][a] = Gfit;	PMT_DATA[2*p-1][a] = gainerror;}
						cout << PMT_DATA[2*p-2][a] << "    " << PMT_DATA[2*p-1][a] << endl;
						//gaindata <<" "<< endl;
						c1->Update();
						c1->WaitPrimitive();
						c1->Print("test.pdf","pdf");
					}	
			
			
		}
	    		auto gr_1 = new TGraphErrors(n,ANGLES,PMT_DATA[0],X_ERR,PMT_DATA[1]);	gr_1->SetMarkerColor(1); gr_1->SetMarkerStyle(8); gr_1->SetName("Position 1"); gr_1->SetTitle("Position 1");
			auto gr_2 = new TGraphErrors(n,ANGLES,PMT_DATA[2],X_ERR,PMT_DATA[3]);	gr_2->SetMarkerColor(2); gr_2->SetMarkerStyle(8); gr_2->SetName("Position 2"); gr_2->SetTitle("Position 2");
			auto gr_3 = new TGraphErrors(n,ANGLES,PMT_DATA[4],X_ERR,PMT_DATA[5]);	gr_3->SetMarkerColor(3); gr_3->SetMarkerStyle(8); gr_3->SetName("Position 3"); gr_3->SetTitle("Position 3");
			auto gr_4 = new TGraphErrors(n,ANGLES,PMT_DATA[6],X_ERR,PMT_DATA[7]);	gr_4->SetMarkerColor(4); gr_4->SetMarkerStyle(8); gr_4->SetName("Position 4"); gr_4->SetTitle("Position 4");
			auto gr_5 = new TGraphErrors(n,ANGLES,PMT_DATA[8],X_ERR,PMT_DATA[9]);	gr_5->SetMarkerColor(5); gr_5->SetMarkerStyle(8); gr_5->SetName("Position 5"); gr_5->SetTitle("Position 5");
			auto gr_6 = new TGraphErrors(n,ANGLES,PMT_DATA[10],X_ERR,PMT_DATA[11]);	gr_6->SetMarkerColor(6); gr_6->SetMarkerStyle(8); gr_6->SetName("Position 6"); gr_6->SetTitle("Position 6");
			auto gr_7 = new TGraphErrors(n,ANGLES,PMT_DATA[12],X_ERR,PMT_DATA[13]);	gr_7->SetMarkerColor(7); gr_7->SetMarkerStyle(8); gr_7->SetName("Position 7"); gr_7->SetTitle("Position 7");
			mg->Add(gr_1);
			mg->Add(gr_2);
			mg->Add(gr_3);
			mg->Add(gr_4);
			mg->Add(gr_5);
			mg->Add(gr_6);
			mg->Add(gr_7);
	 		mg->GetXaxis()->SetTitle("Azimuthal Angle (Degrees)"); //set Xaxis title
			mg->GetYaxis()->SetTitle("Gain (No. of PEs))"); //set Yaxis title
	 		mg->Draw("AP");
			c1->BuildLegend();
			c1->Print("test.pdf)","pdf");
	 
	 
	 
	 }	
 else
 	 {
 	 	cout<<"Invalid value for index."<<endl;
 	 } 


		

	
	
 return;	
}
int main(int argc, char ** argv){
	if(argc < 2) return 1;
	TString path = argv[1];
	TApplication app("app", &argc, argv);
	runfile_working_first(path);
	app.Run();
	return 0;
}
