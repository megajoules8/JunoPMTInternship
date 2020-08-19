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
#include "TFile.h"
#include "DFTmethod.h"
#include "SPEFitter.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "THStack.h"

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
#include "TMatrixD.h"
#include "TMatrixDSym.h"
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
TString ECDF;
TString FCDF;
TString STACK;
TString STACK2;
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
Float_t theta;
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
	 			labdata <<"Mu  Mu_err  w  w_err  alpha  alpha_err  lambda  lambda_err  Theta  Theta_err  sig_reduced  sig_reduced_err  Gain  Gain_err  Chi2r  Fit converged?"<< endl;
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
	 	//cout <<"Input theta :"<<endl;
	 	//cin >> theta;
	 	filename = TString("gain_data_scan") + Form("%d",dat);
	 	TFile *out_file = new TFile("my_rootfile.root","RECREATE");
	 	ofstream gaindata (filename);
	 	ofstream ff ("matrices.txt");
	 	Int_t n = 24;
		Float_t X_ERR[n];
	 	Float_t ANGLES[n];
	 	Float_t BW;
	 Float_t ANGLES_1[n];
	 Float_t ANGLES_2[n];
	 Float_t ANGLES_3[n];
	 Float_t ANGLES_4[n];
	 Float_t ANGLES_5[n];
	 Float_t ANGLES_6[n];
	 Float_t ANGLES_7[n];
	 Int_t counts[7];
	 Int_t REM = 0;
	 
	 	Float_t PMT_DATA[14][24];
	 	Float_t PMT_DATA_NORM[14][24];
	 	TMultiGraph  *mg  = new TMultiGraph();
	 	TMultiGraph  *mg2  = new TMultiGraph();
	 	Float_t max_w =0;
	 	Float_t max_alpha=0;
	 	Float_t max_lambda=0;
	 	Float_t max_theta=0;
	 	Float_t max_mu=0;
	 	Float_t min_w =100;
	 	Float_t min_alpha=100;
	 	Float_t min_lambda=100;
	 	Float_t min_theta=100;
	 	Float_t min_mu=100;
	 	Float_t max_g =0;
	 	Float_t	min_g =100;
	 	Float_t max_chi =0;
	 	Float_t	min_chi =4;
	 	Float_t max_KS =0;
	 	Float_t	min_KS =1;
	 	Float_t chi_red;
	 	Float_t NDF;
	 	Float_t chi;
	 	Float_t KS;
	 	Float_t TEMP;
	 	Float_t TEMP1;
	 	Float_t TOT;
	 	Float_t TOT1;
	 	TString PdfName_end;
	 	Float_t D; 
	 	//PdfName_end = TString("scan") + Form("%d", dat) + TString("_BW_") + Form("%d", BW)+ TString("_Results.pdf)");
	 	TString PdfName_mid;
		//PdfName_mid = TString("scan") + Form("%d", dat) + TString("_BW_") + Form("%d", BW)+ TString("_Results.pdf");
	 	Float_t g_pos_1;
		Float_t g_sum_1;
	 		 	
	 	TH1F *rel_err_w 	= new TH1F("dw", "Relative error of w", 2000 , 0, 100);
	 	rel_err_w->GetXaxis()->SetTitle("Relative error of w");
		rel_err_w->GetYaxis()->SetTitle("Counts");
	 
	 	TH1F *rel_err_alpha 	= new TH1F("dalpha", "Relative error of alpha", 2000 , 0, 100);
	 	rel_err_alpha->GetXaxis()->SetTitle("Relative error of alpha");
		rel_err_alpha->GetYaxis()->SetTitle("Counts");
	 
	 	TH1F *rel_err_lambda 	= new TH1F("dlambda", "Relative error of lambda", 2000 , 0, 100);
	 	rel_err_lambda->GetXaxis()->SetTitle("Relative error of lambda");
		rel_err_lambda->GetYaxis()->SetTitle("Counts");
	 
	 	TH1F *rel_err_theta 	= new TH1F("dtheta", "Relative error of theta", 2000 , 0, 100);
	 	rel_err_theta->GetXaxis()->SetTitle("Relative error of thet");
		rel_err_theta->GetYaxis()->SetTitle("Counts");
	 
	 	TH1F *rel_err_mu 	= new TH1F("dmu", "Relative error of mu", 2000 , 0, 100);
	 	rel_err_mu->GetXaxis()->SetTitle("Relative error of mu");
		rel_err_mu->GetYaxis()->SetTitle("Counts");
	 
	 	TH1F *rel_err_g 	= new TH1F("dg", "Relative error of G", 2000 , 0, 100);
	 	rel_err_g->GetXaxis()->SetTitle("Relative error of G");
		rel_err_g->GetYaxis()->SetTitle("Counts");
	 
	 	TH1F *chisqr 	= new TH1F("chisqr", "Histogram of Chi-Square/Nbins for the Bin range (5,25)", 400 , 0, 20);
	 	chisqr->GetXaxis()->SetTitle("Chi-Square/Nbins");
		chisqr->GetYaxis()->SetTitle("Counts");
	 	
	 	TH1F *KST 	= new TH1F("KST", "Histogram of KS Statistics", 2000 , 0, 1);
	 	KST->GetXaxis()->SetTitle("KS Statistic Value");
		KST->GetYaxis()->SetTitle("Counts");
	 
	 	
	 
	 	for (int p = 1; p<8; ++p)
	 	{	
	 			gaindata <<"Fit data for position "<< p <<": "<< endl;
	 			gaindata <<"angle Mu Mu_err w w_err alpha alpha_err lambda lambda_err Theta Theta_err sig_reduced sig_reduced_err Gain Gain_err chi_sq fit.status D  is Point Good? "<< endl;
	 			//gaindata <<"angle Theta Gain"<< endl;
	 			gaindata <<" "<< endl;
				for (int f=0; f<24; ++f) {X_ERR[f] = 0;}
				//gr_name = TString("GR_") + Form("%d", p);
				Int_t count = 0;
				
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
						Q 	= histo_PED->GetFunction("Fit_Gauss")->GetParameter(1); //get Q from fit
						sigma 	= histo_PED->GetFunction("Fit_Gauss")->GetParameter(2); //get sigma from fit
						cout <<"Q = "<< Q <<" sigma = "<< sigma<< endl;
						Fit_Gauss->Draw("same");
						
						c1->Update();
						c1->WaitPrimitive(); //ROOT waits until you hit ENTER
						
						BW = ceil(histo_LED->GetBinWidth(2));
						//PdfName_end = TString("scan") + Form("%d", dat) + TString("_BW_") + Form("%.0f", BW)+ TString("_theta_") + Form ("%.2f", theta) + TString("_Results.pdf)");
						//PdfName_mid = TString("scan") + Form("%d", dat) + TString("_BW_") + Form("%.0f", BW)+ TString("_theta_") + Form ("%.2f", theta) + TString("_Results.pdf");
						PdfName_end = TString("scan") + Form("%d", dat) + TString("_BW_") + Form("%.0f", BW) + TString("_Results.pdf)");
						PdfName_mid = TString("scan") + Form("%d", dat) + TString("_BW_") + Form("%.0f", BW) + TString("_Results.pdf");
						TString PdfName_start;
						//PdfName_start = TString("scan") + Form("%d", dat) + TString("_BW_") + Form("%.0f", BW)+ TString("_theta_") + Form ("%.2f", theta) + TString("_Results.pdf(");
						PdfName_start = TString("scan") + Form("%d", dat) + TString("_BW_") + Form("%.0f", BW) + TString("_Results.pdf(");
						c1->Print(PdfName_start,"pdf");
						
						Fit_Gauss->SetParameters(amp*histo_LED->GetBinWidth(1)*(1/(sqrt(2*M_PI)*sigma)),Q,sigma);
						Fit_Gauss->SetParLimits(1, Q-2.0*sigma, Q+2.0*sigma); //[1] is for Q, predefined by "gaus"	
						Fit_Gauss->SetParLimits(2, 0.5*sigma,  1.5*sigma); // [2] is for sigma, predefined by "gaus"
							
						histo_LED->Draw();
					
						//if ((p==1)&&(a==1)) {histo_LED->Write();}
						//if ((p==5)&&(a==23)) {histo_LED->Write();}
					
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
						c1->Print(PdfName_mid ,"pdf");
						
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
						Double_t p_test[4] = { 1.0/_G, 10.0, 1/(0.1*_G), 0.2 };
						SPEResponse gamma_test( PMType::GAMMA, p_test );
						
						Int_t nbins = histo_LED->GetNbinsX();
						Double_t xmin = histo_LED->GetXaxis()->GetBinLowEdge(1);
						Double_t xmax = histo_LED->GetXaxis()->GetBinUpEdge(nbins);
						cout << " No. of bins : " << nbins << endl;
						cout << " ( " << xmin << ", " << xmax << " ) " << endl;
						
						Double_t rms = histo_LED->GetRMS();
						Double_t rmse = histo_LED->GetRMSError();
						cout << " RMS : " << rms <<" +/- "<< rmse << endl;
						
						
						DFTmethod dft( 4.0*nbins, xmin, xmax, gamma_test );
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
						c1->Update(); c1->WaitPrimitive(); c1->Print(PdfName_mid ,"pdf");
						TString STATUS;
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
						Float_t pderiv_w = (1/fit.vals[6]) - (1/fit.vals[4]);
						Float_t pderiv_alpha = -fit.vals[7]/pow(fit.vals[6],2);
						Float_t pderiv_lambda = -(1-fit.vals[7])/pow(fit.vals[4],2);
						//gainerror = (fit.vals[7]/fit.vals[6])* ( sqrt( pow( (fit.errs[7]/fit.vals[7]),2 ) + pow( (fit.errs[6]/fit.vals[6]),2 )))   +   ((1-fit.vals[7])/fit.vals[4])*(sqrt( pow( (fit.errs[7]/fit.vals[7]),2 ) + pow( (fit.errs[4]/fit.vals[4]),2 )) );
						double cov_al_lam;
						double cov_al_w;
						double cov_w_lam;
	 					cov_al_lam = fit.mFFT->CovMatrix(4,6);
						cov_al_w = fit.mFFT->CovMatrix(7,6);
						cov_w_lam = fit.mFFT->CovMatrix(4,7);
						//gainerror = sqrt ( pow(pderiv_w*fit.errs[7],2) + pow(pderiv_alpha*fit.errs[6],2) + pow(pderiv_lambda*fit.errs[4],2) );
						gainerror = sqrt ( pow(pderiv_w*fit.errs[7],2) + pow(pderiv_alpha*fit.errs[6],2) + pow(pderiv_lambda*fit.errs[4],2) + 2*cov_al_lam*pderiv_lambda*pderiv_alpha + 2*cov_al_w*pderiv_w*pderiv_alpha + 2*cov_w_lam*pderiv_lambda*pderiv_w );
						BW = histo_LED->GetBinWidth(2);
						//cout<< "8th Bin content of fit = "<<grBF->Eval(histo_LED->GetXaxis()->GetBinCenter(8))<<" +/- "<<grBF->GetErrorY(8)<<endl;
						//cout<< "8th Bin content of LED = "<<histo_LED->GetBinContent(8)<<" +/- "<<histo_LED->GetBinError(8)<<endl;
						chi = 0;
						NDF = 0;
						for (int z=5; z<26; z++) {chi += pow( ( grBF->Eval(histo_LED->GetXaxis()->GetBinCenter(z)) - histo_LED->GetBinContent(z) )/histo_LED->GetBinError(z) , 2 );	if(histo_LED->GetBinContent(z)>0) {++NDF;} }
						
						//NDF = NDF-dft.spef.nparams-4;
						chi_red = chi/(NDF);
						cout<<"NDF = "<<NDF<<endl;
						cout<< "chi_sq/Nbins for the Bin range (5,25) = "<<chi_red<<endl;
						chisqr-> Fill(chi_red); 		if (chi_red > max_chi) {max_chi = chi_red;}	if (chi_red < min_chi) {min_chi = chi_red;}
						//cout<<fit.ndof<<endl;
						
						STACK = TString("Empirical Relative Cumulative Frequency and Relative Cumulative Frequency from Fit for the LED for position = ") + Form("%d",p) + TString(" angle = ")+ Form("%d", 15*a);
						STACK2 = TString("Empirical Relative Cumulative Frequency and Relative Cumulative Frequency from Fit for the LED for position = ") + Form("%d",p) + TString(" angle = ")+ Form("%d", 15*a) + TString(";Charge (DUQ); Relative Cumulative Frequency");
						ECDF = TString("Empirical Cumulative Frequency for the LED for position = ") + Form("%d",p) + TString(" angle = ")+ Form("%d", 15*a);
						FCDF = TString("Cumulative Frequency from Fit for the LED for position = ") + Form("%d",p) + TString(" angle = ")+ Form("%d", 15*a);
						
						THStack *hs = new THStack(STACK,STACK2);
						
						TH1F *LED_CDF 	= new TH1F(ECDF, ECDF, nbins , xmin-BW, xmax);
	 					TH1F *FIT_CDF 	= new TH1F(FCDF, FCDF, nbins , xmin-BW, xmax);
					
						
	 					//hs->GetXaxis()->SetTitle("Charge (DUQ)");
						//hs->GetYaxis()->SetTitle("Relative Cumulative Frequency");
						//cout<<"Good"<<endl;
						TEMP = 0;	TEMP1 = 0;	D = 0;	TOT = 0;	TOT1 = 0;
						LED_CDF-> Fill(histo_LED->GetBinContent(0));
						FIT_CDF-> Fill(grBF->Eval(histo_LED->GetXaxis()->GetBinCenter(0)));
						for (int t=0; t<nbins; ++t){TOT+= histo_LED->GetBinContent(t);	TOT1 += grBF->Eval(histo_LED->GetXaxis()->GetBinCenter(t));}
						for (int t=0; t<nbins; ++t){TEMP += histo_LED->GetBinContent(t);	LED_CDF-> SetBinContent((t+1), TEMP/TOT);	 TEMP1 += grBF->Eval(histo_LED->GetXaxis()->GetBinCenter(t));	FIT_CDF-> SetBinContent(t+1, TEMP1/TOT1);	if ( abs((TEMP/TOT)-(TEMP1/TOT1)) > D ){D = abs((TEMP/TOT)-(TEMP1/TOT1));}	}
						
						LED_CDF->SetMarkerStyle( 20 ); LED_CDF->SetMarkerSize( 0.4 ); LED_CDF->SetLineColor( kRed ); LED_CDF->SetMarkerColor( kRed ); LED_CDF->SetStats(0);  hs->Add( LED_CDF );
	 					FIT_CDF->SetMarkerStyle( 20 ); FIT_CDF->SetMarkerSize( 0.4 ); FIT_CDF->SetLineColor( kBlue ); FIT_CDF->SetMarkerColor( kBlue ); FIT_CDF->SetStats(0);  hs->Add( FIT_CDF );
	 					hs->Draw("nostack");	c1->Update(); c1->WaitPrimitive(); c1->Print(PdfName_mid ,"pdf");
						
						//Float_t temp = 0;
						//KS = 0;
						//for (int z=5; z<26; z++) {temp = abs (grBF->Eval(histo_LED->GetXaxis()->GetBinCenter(z)) - histo_LED->GetBinContent(z)); 	if(temp > KS) {KS = temp; temp = 0;}	}
						KST-> Fill(D);		if (D > max_KS) {max_KS = D;}	if (D < min_KS) {min_KS = D;}
						cout<< "KS statistic value (D) for the position "<<p<<", angle "<<a*15<<" = "<< D <<endl;
						
						ff <<"Correlation matrix for Position = "<<p<<" , Angle = "<<a*15<<endl;
						ff <<" "<<endl;
						ff <<"            |      0     |      1      |      2     |      3      |      4     |      5      |      6     |      7      |"<<endl;
						ff <<"_________________________________________________________________________________________________________________________"<<endl;
						for (int y =0; y<8; ++y)
							{
								ff<<"     "<<y<<"      |";
								for (int x=0; x<8; ++x)
									{
										ff<<"     "<<fit.mFFT->Correlation(y,x);
									}
								ff<<"     "<<endl;
							}
						ff <<" "<<endl;
					
						ff <<"Covariance matrix for Position = "<<p<<" , Angle = "<<a*15<<endl;
						ff <<" "<<endl;
						ff <<"            |      0     |      1      |      2     |      3      |      4     |      5      |      6     |      7      |"<<endl;
						ff <<"_________________________________________________________________________________________________________________________"<<endl;
						for (int y =0; y<8; ++y)
							{
								ff<<"     "<<y<<"      |";
								for (int x=0; x<8; ++x)
									{
										ff<<"     "<<fit.mFFT->CovMatrix(y,x);
									}
								ff<<"     "<<endl;
							}
						ff <<" "<<endl;	
					
						if ((fit.chi2r <= 3) && (fit.fit_status == 0)  && (gainerror/Gfit<0.1) )	{ANGLES[count] = 15*a;  PMT_DATA[2*p-2][count] = Gfit;	PMT_DATA[2*p-1][count] = gainerror;   ++count;	STATUS = "Yes";}
						else {STATUS = "No"; ++REM;}
						gaindata << a*15 <<"  "<<fit.vals[3]<<"  "<<fit.errs[3]<<"  "<<fit.vals[7]<<"  "<<fit.errs[7]<<"  "<<fit.vals[6]<<"  "<<fit.errs[6]<<"  "<<fit.vals[4]<<"  "<<fit.errs[4]<<"  "<<fit.vals[5]<<"	"<<fit.errs[5]<<"  "<< sig_reduced<<"  "<<sig_reduced_err<<"  "\
						<<Gfit <<"  "<< gainerror<<"  "<< fit.chi2r<<"  "<<"  "<<fit.fit_status<<"  "<<D<<"  "<<STATUS<<"  "<<endl;
						cout << " Gain (DUQ) : " << Gfit <<" +/- "<< gainerror << endl;
						cout << " Bin Width : " << BW << endl;
						//gaindata <<" "<< endl;
						
						if ((fit.chi2r <= 3) && (fit.fit_status == 0))
							
						{	rel_err_w-> Fill(fit.errs[7]*100/fit.vals[7]); 		if (fit.errs[7]*100/fit.vals[7] > max_w) {max_w = fit.errs[7]*100/fit.vals[7];} 		if (fit.errs[7]*100/fit.vals[7] < min_w) {min_w = fit.errs[7]*100/fit.vals[7];}
							rel_err_alpha-> Fill(fit.errs[6]*100/fit.vals[6]); 	if (fit.errs[6]*100/fit.vals[6] > max_alpha) {max_alpha = fit.errs[6]*100/fit.vals[6];}		if (fit.errs[6]*100/fit.vals[6] < min_alpha) {min_alpha = fit.errs[6]*100/fit.vals[6];}
							rel_err_lambda-> Fill(fit.errs[4]*100/fit.vals[4]); 	if (fit.errs[4]*100/fit.vals[4] > max_lambda) {max_lambda = fit.errs[4]*100/fit.vals[4];}	if (fit.errs[7]*100/fit.vals[7] < min_lambda) {min_lambda = fit.errs[4]*100/fit.vals[4];}
							rel_err_theta-> Fill(fit.errs[5]*100/fit.vals[5]); 	if (fit.errs[5]*100/fit.vals[5] > max_theta) {max_theta = fit.errs[5]*100/fit.vals[5];}		if (fit.errs[7]*100/fit.vals[7] < min_theta) {min_theta = fit.errs[5]*100/fit.vals[5];}
							rel_err_mu-> Fill(fit.errs[3]*100/fit.vals[3]); 	if (fit.errs[3]*100/fit.vals[3] > max_mu) {max_mu = fit.errs[3]*100/fit.vals[3];}		if (fit.errs[7]*100/fit.vals[7] < min_mu) {min_mu = fit.errs[3]*100/fit.vals[3];}
						 	rel_err_g-> Fill(gainerror*100/Gfit); 			if (gainerror*100/Gfit > max_g) {max_g = gainerror*100/Gfit;}					if (gainerror*100/Gfit < min_g) {min_g = gainerror*100/Gfit;}
						}
					
						
		 				//cout<<cov<<endl;
						
						
					}
			
			if (p == 1) {	for(Int_t r=0; r<=count; ++r) {ANGLES_1[r] = ANGLES[r]; g_sum_1 += PMT_DATA[2*p-2][r];}	counts[p-1] = count; g_pos_1 = g_sum_1/count;	for(Int_t s=0; s<=count; ++s) {PMT_DATA_NORM[2*p-2][s] = PMT_DATA[2*p-2][s]/g_pos_1;	PMT_DATA_NORM[2*p-1][s] = PMT_DATA[2*p-1][s]/g_pos_1;} }
			if (p == 2) {	for(Int_t r=0; r<=count; ++r) {ANGLES_2[r] = ANGLES[r]; }	counts[p-1] = count;	for(Int_t s=0; s<=count; ++s) {PMT_DATA_NORM[2*p-2][s] = PMT_DATA[2*p-2][s]/g_pos_1;	PMT_DATA_NORM[2*p-1][s] = PMT_DATA[2*p-1][s]/g_pos_1;}}
			if (p == 3) {	for(Int_t r=0; r<=count; ++r) {ANGLES_3[r] = ANGLES[r]; }	counts[p-1] = count;	for(Int_t s=0; s<=count; ++s) {PMT_DATA_NORM[2*p-2][s] = PMT_DATA[2*p-2][s]/g_pos_1;	PMT_DATA_NORM[2*p-1][s] = PMT_DATA[2*p-1][s]/g_pos_1;}}
			if (p == 4) {	for(Int_t r=0; r<=count; ++r) {ANGLES_4[r] = ANGLES[r]; }	counts[p-1] = count;	for(Int_t s=0; s<=count; ++s) {PMT_DATA_NORM[2*p-2][s] = PMT_DATA[2*p-2][s]/g_pos_1;	PMT_DATA_NORM[2*p-1][s] = PMT_DATA[2*p-1][s]/g_pos_1;}}
			if (p == 5) {	for(Int_t r=0; r<=count; ++r) {ANGLES_5[r] = ANGLES[r]; }	counts[p-1] = count;	for(Int_t s=0; s<=count; ++s) {PMT_DATA_NORM[2*p-2][s] = PMT_DATA[2*p-2][s]/g_pos_1;	PMT_DATA_NORM[2*p-1][s] = PMT_DATA[2*p-1][s]/g_pos_1;}}
			if (p == 6) {	for(Int_t r=0; r<=count; ++r) {ANGLES_6[r] = ANGLES[r]; }	counts[p-1] = count;	for(Int_t s=0; s<=count; ++s) {PMT_DATA_NORM[2*p-2][s] = PMT_DATA[2*p-2][s]/g_pos_1;	PMT_DATA_NORM[2*p-1][s] = PMT_DATA[2*p-1][s]/g_pos_1;}}
			if (p == 7) {	for(Int_t r=0; r<=count; ++r) {ANGLES_7[r] = ANGLES[r]; }	counts[p-1] = count;	for(Int_t s=0; s<=count; ++s) {PMT_DATA_NORM[2*p-2][s] = PMT_DATA[2*p-2][s]/g_pos_1;	PMT_DATA_NORM[2*p-1][s] = PMT_DATA[2*p-1][s]/g_pos_1;}}
			
			
		}
	    		chisqr->SetMarkerStyle( 20 ); chisqr->SetMarkerSize( 0.4 ); chisqr->SetLineColor( kBlack ); chisqr->SetMarkerColor( kBlack ); chisqr->SetStats(0); chisqr->GetXaxis()->SetRangeUser(min_chi,max_chi); chisqr->Draw( "" );
	 		c1->Update(); c1->WaitPrimitive(); c1->Print(PdfName_mid ,"pdf");
	 		KST->SetMarkerStyle( 20 ); KST->SetMarkerSize( 0.4 ); KST->SetLineColor( kBlack ); KST->SetMarkerColor( kBlue ); KST->SetStats(0); KST->GetXaxis()->SetRangeUser(min_KS,max_KS); KST->Draw( "" );
			c1->Update(); c1->WaitPrimitive(); c1->Print(PdfName_mid ,"pdf");
	 		rel_err_w->SetMarkerStyle( 20 ); rel_err_w->SetMarkerSize( 0.4 ); rel_err_w->SetLineColor( kBlack ); rel_err_w->SetMarkerColor( kBlack ); rel_err_w->SetStats(0); rel_err_w->GetXaxis()->SetRangeUser(min_w,max_w); rel_err_w->Draw( "" );
			c1->Update(); c1->WaitPrimitive(); c1->Print(PdfName_mid ,"pdf");
			rel_err_alpha->SetMarkerStyle( 20 ); rel_err_alpha->SetMarkerSize( 0.4 ); rel_err_alpha->SetLineColor( kBlack ); rel_err_alpha->SetMarkerColor( kBlack ); rel_err_alpha->SetStats(0); rel_err_alpha->GetXaxis()->SetRangeUser(min_alpha,max_alpha); rel_err_alpha->Draw( "" );
			c1->Update(); c1->WaitPrimitive(); c1->Print(PdfName_mid ,"pdf");
			rel_err_lambda->SetMarkerStyle( 20 ); rel_err_lambda->SetMarkerSize( 0.4 ); rel_err_lambda->SetLineColor( kBlack ); rel_err_lambda->SetMarkerColor( kBlack ); rel_err_lambda->SetStats(0); rel_err_lambda->GetXaxis()->SetRangeUser(min_lambda,max_lambda); rel_err_lambda->Draw( "" );
			c1->Update(); c1->WaitPrimitive(); c1->Print(PdfName_mid ,"pdf");
			rel_err_theta->SetMarkerStyle( 20 ); rel_err_theta->SetMarkerSize( 0.4 ); rel_err_theta->SetLineColor( kBlack ); rel_err_theta->SetMarkerColor( kBlack ); rel_err_theta->SetStats(0); rel_err_theta->GetXaxis()->SetRangeUser(min_theta,max_theta); rel_err_theta->Draw( "" );
			c1->Update(); c1->WaitPrimitive(); c1->Print(PdfName_mid ,"pdf");
			rel_err_mu->SetMarkerStyle( 20 ); rel_err_mu->SetMarkerSize( 0.4 ); rel_err_mu->SetLineColor( kBlack ); rel_err_mu->SetMarkerColor( kBlack ); rel_err_mu->SetStats(0); rel_err_mu->GetXaxis()->SetRangeUser(min_mu,max_mu); rel_err_mu->Draw( "" );
			c1->Update(); c1->WaitPrimitive(); c1->Print(PdfName_mid ,"pdf");
	 		rel_err_g->SetMarkerStyle( 20 ); rel_err_g->SetMarkerSize( 0.4 ); rel_err_g->SetLineColor( kBlack ); rel_err_g->SetMarkerColor( kBlack ); rel_err_g->SetStats(0); rel_err_g->GetXaxis()->SetRangeUser(min_g,max_g); rel_err_g->Draw( "" );
			c1->Update(); c1->WaitPrimitive(); c1->Print(PdfName_mid ,"pdf");
	 
	 		auto gr_1 = new TGraphErrors(counts[0],ANGLES_1,PMT_DATA[0],X_ERR,PMT_DATA[1]);	gr_1->SetMarkerColor(1); gr_1->SetLineColor(1); gr_1->SetMarkerStyle(7); gr_1->SetName("Position 1"); gr_1->SetTitle("Position 1");
			auto gr_2 = new TGraphErrors(counts[1],ANGLES_2,PMT_DATA[2],X_ERR,PMT_DATA[3]);	gr_2->SetMarkerColor(2); gr_2->SetLineColor(2); gr_2->SetMarkerStyle(7); gr_2->SetName("Position 2"); gr_2->SetTitle("Position 2");
			auto gr_3 = new TGraphErrors(counts[2],ANGLES_3,PMT_DATA[4],X_ERR,PMT_DATA[5]);	gr_3->SetMarkerColor(3); gr_3->SetLineColor(3); gr_3->SetMarkerStyle(7); gr_3->SetName("Position 3"); gr_3->SetTitle("Position 3");
			auto gr_4 = new TGraphErrors(counts[3],ANGLES_4,PMT_DATA[6],X_ERR,PMT_DATA[7]);	gr_4->SetMarkerColor(4); gr_4->SetLineColor(4); gr_4->SetMarkerStyle(7); gr_4->SetName("Position 4"); gr_4->SetTitle("Position 4");
			auto gr_5 = new TGraphErrors(counts[4],ANGLES_5,PMT_DATA[8],X_ERR,PMT_DATA[9]);	gr_5->SetMarkerColor(5); gr_5->SetLineColor(5); gr_5->SetMarkerStyle(7); gr_5->SetName("Position 5"); gr_5->SetTitle("Position 5");
			auto gr_6 = new TGraphErrors(counts[5],ANGLES_6,PMT_DATA[10],X_ERR,PMT_DATA[11]); gr_6->SetMarkerColor(6); gr_6->SetLineColor(6); gr_6->SetMarkerStyle(7); gr_6->SetName("Position 6"); gr_6->SetTitle("Position 6");
			auto gr_7 = new TGraphErrors(counts[6],ANGLES_7,PMT_DATA[12],X_ERR,PMT_DATA[13]); gr_7->SetMarkerColor(7); gr_7->SetLineColor(7); gr_7->SetMarkerStyle(7); gr_7->SetName("Position 7"); gr_7->SetTitle("Position 7");
			mg->Add(gr_1);
			mg->Add(gr_2);
			mg->Add(gr_3);
			mg->Add(gr_4);
			mg->Add(gr_5);
			mg->Add(gr_6);
			mg->Add(gr_7);
	 		TString gtitle;
	 		gtitle = TString("Graph of Gain (DUQ) vs. Azimuthal angle for scan") + Form("%d",dat) + TString(", ") + Form("%d", REM) + TString(" points removed");
	 		mg->SetTitle (gtitle);
	  		mg->GetXaxis()->SetTitle("Azimuthal Angle (Degrees)"); //set Xaxis title
			mg->GetYaxis()->SetTitle("Gain (DUQ))"); //set Yaxis title
	 		mg->GetYaxis()->SetRangeUser(0,2000);
	 		mg->Draw("APL");
			c1->BuildLegend();
	 		c1->Update(); c1->WaitPrimitive();
			c1->Print( PdfName_mid ,"pdf");
		 	
	 		
	 		auto gr2_1 = new TGraphErrors(counts[0],ANGLES_1,PMT_DATA_NORM[0],X_ERR,PMT_DATA_NORM[1]);	gr2_1->SetMarkerColor(1); gr2_1->SetLineColor(1); gr2_1->SetMarkerStyle(7); gr2_1->SetName("Position 1"); gr2_1->SetTitle("Position 1");
			auto gr2_2 = new TGraphErrors(counts[1],ANGLES_2,PMT_DATA_NORM[2],X_ERR,PMT_DATA_NORM[3]);	gr2_2->SetMarkerColor(2); gr2_2->SetLineColor(2);  gr2_2->SetMarkerStyle(7); gr2_2->SetName("Position 2"); gr2_2->SetTitle("Position 2");
			auto gr2_3 = new TGraphErrors(counts[2],ANGLES_3,PMT_DATA_NORM[4],X_ERR,PMT_DATA_NORM[5]);	gr2_3->SetMarkerColor(3); gr2_3->SetLineColor(3); gr2_3->SetMarkerStyle(7); gr2_3->SetName("Position 3"); gr2_3->SetTitle("Position 3");
			auto gr2_4 = new TGraphErrors(counts[3],ANGLES_4,PMT_DATA_NORM[6],X_ERR,PMT_DATA_NORM[7]);	gr2_4->SetMarkerColor(4); gr2_4->SetLineColor(4); gr2_4->SetMarkerStyle(7); gr2_4->SetName("Position 4"); gr2_4->SetTitle("Position 4");
			auto gr2_5 = new TGraphErrors(counts[4],ANGLES_5,PMT_DATA_NORM[8],X_ERR,PMT_DATA_NORM[9]);	gr2_5->SetMarkerColor(5); gr2_5->SetLineColor(5); gr2_5->SetMarkerStyle(7); gr2_5->SetName("Position 5"); gr2_5->SetTitle("Position 5");
			auto gr2_6 = new TGraphErrors(counts[5],ANGLES_6,PMT_DATA_NORM[10],X_ERR,PMT_DATA_NORM[11]); 	gr2_6->SetMarkerColor(6); gr2_6->SetLineColor(6); gr2_6->SetMarkerStyle(7); gr2_6->SetName("Position 6"); gr2_6->SetTitle("Position 6");
			auto gr2_7 = new TGraphErrors(counts[6],ANGLES_7,PMT_DATA_NORM[12],X_ERR,PMT_DATA_NORM[13]); 	gr2_7->SetMarkerColor(7); gr2_7->SetLineColor(7); gr2_7->SetMarkerStyle(7); gr2_7->SetName("Position 7"); gr2_7->SetTitle("Position 7");
	 		mg2->Add(gr2_1);
			mg2->Add(gr2_2);
			mg2->Add(gr2_3);
			mg2->Add(gr2_4);
			mg2->Add(gr2_5);
			mg2->Add(gr2_6);
			mg2->Add(gr2_7);
	 		TString gtitlenorm;
	 		gtitlenorm = TString("Graph of Normallized Gain vs. Azimuthal angle for scan") + Form("%d",dat) + TString(", ") + Form("%d", REM) + TString(" points removed");;
	 		mg2->SetTitle (gtitlenorm);
	 		mg2->GetXaxis()->SetTitle("Azimuthal Angle (Degrees)"); //set Xaxis title
			mg2->GetYaxis()->SetTitle("Normalized Gain"); //set Yaxis title
	 		mg2->GetYaxis()->SetRangeUser(0,2);
	 		mg2->Draw("APL");
			c1->BuildLegend();
	 		c1->Update(); c1->WaitPrimitive();
			c1->Print( PdfName_end ,"pdf");
	 		cout <<"Mean relative error of w , alpha, lambda, theta, mu, G: "<<endl;
	 		cout <<rel_err_w->GetMean()<<" "<<rel_err_alpha->GetMean()<<" "<<rel_err_lambda->GetMean()<<" "<<rel_err_theta->GetMean()<<" "<<rel_err_mu->GetMean()<<" "<<rel_err_g->GetMean()<<endl;
	 		out_file->Close();

	 
	 
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
