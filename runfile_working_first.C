#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include "Riostream.h"
#include <vector>
#include <TString.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include "hiss.h"
using namespace std;
//main function
void runfile_working_first(TString path_to_M1)
	{ 
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
TString PED = "HVSCAN/%d/PED/F1--Trace--00000.txt";
TString LED = "HVSCAN/%d/LED/F1--Trace--00000.txt";
double Q ;
double sigma;
double amp;

for (i=0; i<5; ++i)
	{	
	//define the 2 histograms for PED and LED
	TH1F *histo_PED = hiss(path_to_M1 + Form(PED, HV[i]));
	TH1F *histo_LED = hiss(path_to_M1 + Form(LED, HV[i]));
	//cout << Form(PED, HV[i]) << endl; getchar();
	Q = histo_PED->GetMean(); //get Q initially
	sigma = histo_PED->GetRMS(); //get sigma initially
	
	histo_PED->GetXaxis()->SetTitle("charge(nVs)"); //set Xaxis title
	histo_PED->GetYaxis()->SetTitle("amplitude"); //set Yaxis title
	
	histo_LED->GetXaxis()->SetTitle("charge(nVs)"); //set Xaxis title
	histo_LED->GetYaxis()->SetTitle("amplitude"); //set Yaxis title
	
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
	cout <<"Q_New = "<< Q <<" sigma_New = "<< sigma<< endl;
	
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
	}

return;	
}



//vectors storing pices of filenames	
	// vector< vector <TString> > PED_LED;
	// vector<TString> a;
	// a.clear();
	// a.push_back("/PED/");
	// a.push_back("/LED/");
	// PED_LED.push_back(a);
	// a.clear();
	
	// vector< vector<TString> > HV;
	// vector<TString> b;
	// b.clear();
	// b.push_back("1250");
	// b.push_back("1300");
	// b.push_back("1350");
	// b.push_back("1401");
	// b.push_back("1450");
	// HV.push_back(b);



// while (i<HV_Count)
	// {
		// while (j<PED_LED_Count)
		// {
			// //histogram event definition
			// TH1F *hiss_gram = hiss(path + HV[0][i] + PED_LED[0][j] , filename, j); 
			// //j=0 corresponds to PED and j=1 corresponds to LED
			// hiss_gram->GetXaxis()->SetTitle("charge(nVs)"); //set Xaxis title
			// hiss_gram->GetYaxis()->SetTitle("amplitude"); //set Yaxis title
			// hiss_gram->Draw();
			// //Mean and the SD
			// double Q ;
			// double sigma;
			// double amp;
			// //define choices for PED and LED
			
			// if (j==0) //PED
				// {
				// ofstream ff ("Q_sigma.txt");
				// Q 		= hiss_gram->GetMean();
				// sigma 	= hiss_gram->GetRMS();
				// amp		= hiss_gram->Integral();
				
				
				// // TF1  *Fit_Gauss = new TF1("Fit_Gauss","[0]* exp(-0.5*pow(x-[1])/[2]) ,2 )) ", (Q - 3*sigma), (Q + 3*sigma));
				// // Fit_Gauss-> SetParNames("amp_gauss", "Mean_Gauss", "sigma_gauss", "min", "max");
				// // Fit_Gauss->SetParameters(amp, Q, sigma);
				
				
				
			// //define the fit function instance	
				// //TF1  *Fit_Gauss = new TF1("Fit_Gauss","(1/(2*[pi*[0])) * exp(-0.5* pow( ((x-[1])/[0]) ,2 )) ", (Q - 3*sigma), (Q + 3*sigma));
				// TF1  *Fit_Gauss = new TF1("Fit_Gauss","gaus", (Q - 3*sigma), (Q + 3*sigma));
				
				// Fit_Gauss->SetParameters(amp*hiss_gram->GetBinWidth(1)*(1/(sqrt(2*M_PI)*sigma)),Q,sigma);
				// Fit_Gauss->SetNpx(10000);
				// Fit_Gauss->Draw("same");
				
				// //hiss_gram->Fit("Fit_Gauss", "RI");
				// //Fit_Gauss->Draw();
				
				// cout << "This is PED"	<<endl;
				// cout <<"Q = "	<< Q 	<< " sigma = "<< sigma << " Integral = "<< amp << endl;
				// ff   <<Q 		<< endl;
				// ff   << sigma 	<<endl;
				// ff   <<  amp		<< endl;
				
				// }
			
			// if (j==1) //LED
				// {
				// cout<< "This is LED"<<endl; 
				// ifstream scan;
				// scan.open("Q_sigma.txt");
				
				// while(scan >> Q >> sigma)
					// {
					
					// cout << "Q from previous= " << Q << " " << "sigma from previous= " << sigma << " Integral from previous= " << amp <<endl;
      
					// }
				// scan.close();
				// } 

			
			//++j;
			//++num; 
			//ROOT waits untill you hit enter

		//++i;
		//j=0;

	
	
	
	//https://root-forum.cern.ch/t/error-in-range-inf-nan-propagated-to-the-pad/29988 saving as pdf
