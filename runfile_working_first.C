

#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include "Riostream.h"
#include <vector>
#include <TString.h>
#include "hiss.C"
using namespace std;
//main function
void runfile_working_first(TString path, TString filename)
	{ 
//define Canvas event	
    TCanvas *c1 = new TCanvas( "c1", "" );
	c1->SetLogy(); //set log Y
//vectors storing pices of filenames	
	vector< vector <TString> > PED_LED;
	vector<TString> a;
	a.clear();
	a.push_back("/PED/");
	a.push_back("/LED/");
	PED_LED.push_back(a);
	a.clear();
	
	vector< vector<TString> > HV;
	vector<TString> b;
	b.clear();
	b.push_back("1250");
	b.push_back("1300");
	b.push_back("1350");
	b.push_back("1401");
	b.push_back("1450");
	HV.push_back(b);
	
	int i=0;
	int j=0;
	int PED_LED_Count =2;
	int HV_Count = 5;
	int num =1;
//run through each folder and file	
	while (i<HV_Count)
	{
		while (j<PED_LED_Count)
		{
			//histogram event definition
			TH1F *hiss_gram = hiss(path + HV[0][i] + PED_LED[0][j] , filename, j); 
			//j=0 corresponds to PED and j=1 corresponds to LED
			hiss_gram->GetXaxis()->SetTitle("charge(nVs)"); //set Xaxis title
			hiss_gram->GetYaxis()->SetTitle("amplitude"); //set Yaxis title
			hiss_gram->Draw();
			//Mean and the SD
			double Q ;
			double sigma;
			//define choices for PED and LED
			
			if (j==0) //PED
				{
				ofstream ff ("Q_sigma.txt");
				Q 		= hiss_gram->GetMean();
				sigma 	= hiss_gram->GetRMS();
			//define the fit function instance	
				TF1  *Fit_Gauss = new TF1("Fit_Gauss","(1/(2*[pi*[0])) * exp(-0.5* pow( ((x-[1])/[0]) ,2 )) ", Q - 3*sigma, Q + 3*sigma);
				//Fit_Gauss->SetParameters(Q,sigma);
				//Fit_Gauss->SetNpx(10000);
				//Fit_Gauss->Draw();
				//hiss_gram->Fit_Gauss("Fit_Gauss");

				cout<< "This is PED"<<endl;
				cout <<"Q = "<< Q << " sigma = "<< sigma <<endl;
				ff   <<Q << endl;
				ff   << sigma <<endl;
				}
			
			if (j==1) //LED
				{
				cout<< "This is LED"<<endl; 
				ifstream scan;
				scan.open("Q_sigma.txt");
				
				while(scan >> Q >> sigma)
					{
					
					cout << "Q from previous= " << Q << " " << "sigma from previous= " << sigma << endl;
      
					}
				scan.close();
				}

			c1->Update();
			c1->WaitPrimitive(); //ROOT waits untill you hit enter
			++j;
			++num;
		}
		++i;
		j=0;
	}

	
	}

//resources: https://root-forum.cern.ch/t/open-files-in-a-directory-with-a-for-loop/12471
//https://stackoverflow.com/questions/50139639/how-to-run-a-c-program-multiple-times-with-different-input-files/50139735
