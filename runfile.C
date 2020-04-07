#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <Riostream.h>
#include <vector>
#include <TString.h>
#include <TH1F.h>
#include "hiss_s.h"
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


using namespace std;
//main function
void runfile(TString path_to_M1)
{

  gROOT->Reset();
  
  PMTStyle::SetDefaultStyle();
  
//define Canvas event	
    TCanvas *c1 = new TCanvas( "c1", "" );
	//c1->SetLogy(); //set log Y
	//Int_t HV[5] = {1250,1300,1350,1401,1450};
	int i=0;
	int j=0;
	//int PED_LED_Count =2;
	//int HV_Count = 5;
	int num =1;

//run through each folder and file	
// TString PED = "/HVSCAN/%d/PED/F1--Trace--00000.txt";
// TString LED = "/HVSCAN/%d/LED/F1--Trace--00000.txt";

//new TString definitions
TString HV_Value_PED;
TString HV_Value_LED;
double Q ;
double sigma;
double amp;
//ofstream ff ("gains.txt"); // write the respective voltages and gains to a file in directory
	
for (i=0; i<5; ++i)
	{	
	
	//define 2 strings to specify to hiss whether we are in PED or LED
	HV_Value_PED = TString ("Pedestal Run");
	HV_Value_LED = TString ("LED Run");
	//define the 2 histograms for PED and LED
	TH1F *histo_PED = hiss_s(path_to_M1 , HV_Value_PED, -100, 200, 0);
	
	TH1F *histo_LED = hiss_s(path_to_M1 , HV_Value_LED, -100, 200, 1);
	
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
	
	//double N = Fit_Gauss->Integral(Q-5*sigma, Q+5*sigma);
	//double N = histo_LED->GetXaxis()->FindBin(50.0) - h->GetXaxis()->FindBin(0.)
	//double N = histo_LED->GetXaxis()->FindBin(Q-3*sigma) - histo_LED->GetXaxis()->FindBin(Q+3*sigma);
	//double N = histo_LED -> Integral(histo_LED->FindFixBin(Q - 5*sigma), histo_LED->FindFixBin(Q + 3*sigma), ""); //get N from LED

	//Leo told me not to use this integral, and to use a method to find the number of "Entries" within the range of the fit.
	// No. of "Entries" (N) = Sum of all amplitudes within the range of interest
	// so we need to figure out a way to get the  sum of the amplitudes within the range of interest WITHOUT using the integral (or if you do, we need to divide by the correct factors) 
	//So line 71 is what we used yesterday, but you may get something from line 70 (in comments).
	
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
	
	c1->Update();
	c1->WaitPrimitive();

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

