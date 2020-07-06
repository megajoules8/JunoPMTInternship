#include <vector>
#include <iostream>
#include <fstream>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TVectorT.h>
#include <TApplication.h>
#include "TVectorT.h"
#include "TClass.h"
using namespace std;
int main(int argc, char ** argv)
{
	if(argc != 2)
	{
		std::cout << "Usage: " << argv[0] << " file.root" << std::endl;
		return 1;
	}
	
	TString PED;
	TString LED;
	TString Nom;
	int dat;
	
	cout <<"Input the data set you wish to analyze:"<<endl;
	cout <<"Input 248 for scan248"<<endl;
	cout <<"Input 3737 for scan3737"<<endl;
	cout <<"Input 3899 for scan3899"<<endl;
	cout <<"Input 846 for scan846"<<endl;
	cin >> dat;
		
	for (int p = 1; p<8; ++p)
	{
				
		//run through each folder and file	
		Nom = argv[1] + TString("/scan") + Form("%d",dat)+ TString("/scan") +Form("%d", dat) + TString("_position") + Form("%d", p) + TString("_angle");
		
		for (int a = 0; a<24; ++a)
			{				
				TString filename;
				TString histname;
				
				TString histname_LED;
				TString histname_PED;
				TString filename_LED;
				TString filename_PED;
				filename = Nom + Form("%d", 15*a) + TString(".root");
				histname = TString("position = ") + Form("%d",p) + TString(" angle = ")+ Form("%d", 15*a);
				
				histname_LED = histname + TString(" (LED)");
				histname_PED = histname + TString(" (PED)");
				
				filename_LED = TString("scan") +Form("%d", dat) + TString("_position") + Form("%d", p) + TString("_angle")+ Form("%d", 15*a) + TString("_LED.txt");
				filename_PED = TString("scan") +Form("%d", dat) + TString("_position") + Form("%d", p) + TString("_angle")+ Form("%d", 15*a) + TString("_PED.txt");
				cout << "accessing file "<<filename<<endl;
				TFile * f = TFile::Open(filename);
				TTree * pmt_tree = (TTree*) f->Get("pmt_tree");
				TVectorT<float> * time_vector = (TVectorT<float> *) f->Get("time_vector");
				std::vector<float> * wave_vector = 0;
				pmt_tree->SetBranchAddress("wave_vector", &wave_vector);
				//TApplication TApp("TApp", &argc, argv);
				//TCanvas * c = new TCanvas();
				int nbins = 2001;
				float t_min;
				float t_max;
				float h_max;
				float h_min;
				if (dat == 3737)
					{
						t_min = 400;
						t_max = 470;
						h_max = 1000;
						h_min = -6000;
					
					}
				else if (dat == 3899)
					{
						t_min = 220;
						t_max = 280;
						h_max = 1000;
						h_min = -6000;
					}
				else if (dat == 248)
					{
						t_min = 300;
						t_max = 350;
						h_max = 200;
						h_min = -3000;
					} 
				float t_min_PED = 0;
				float t_max_PED = t_max-t_min;
				float Integral = 0;
				float Integral_PED = 0 ;
				float bin_width = 0;
				//definition of the histogram
				TH1F *Juno = new TH1F("Juno", histname_LED, nbins , h_min, h_max);
				TH1F *JunoPED = new TH1F("JunoPED", histname_PED, nbins , h_min, h_max);
				for(int i=0; i < pmt_tree->GetEntries(); ++i)
					{
					pmt_tree->GetEntry(i);
					TVectorT<float> wave_vector_root(wave_vector->size());
					for(int iWV=0; iWV < wave_vector->size(); ++iWV)
						{
							wave_vector_root[iWV] = wave_vector->at(iWV); 
							if ( ( (*time_vector)(iWV) >= t_min ) && ( (*time_vector)(iWV) < t_max) )
								{
									Integral += wave_vector_root[iWV]; 
								}
							if ( ( (*time_vector)(iWV) >= t_min_PED ) && ( (*time_vector)(iWV) < t_max_PED) )
								{
									Integral_PED += wave_vector_root[iWV]; 
								}	
						}
					Juno-> Fill(Integral);
					JunoPED-> Fill(Integral_PED);	
					Juno->GetXaxis()->SetTitle("Integral");
					Juno->GetYaxis()->SetTitle("Counts");
					JunoPED->GetXaxis()->SetTitle("Integral");
					JunoPED->GetYaxis()->SetTitle("Counts");
					Integral = 0;
					Integral_PED = 0;
				}
					//JunoPED->Draw();
					//c->Update();
					//c->WaitPrimitive();
					//Juno->Draw();
				 ofstream ff (filename_LED);
				 //ff <<"Juno PMT data"<<endl;
				 ff <<"Histogram of Integral vs. Counts"<<endl;
				 ff <<"No. of bins = "<<nbins<<endl;
				 ff <<"position = "<< p <<" , angle = "<<15*a<<endl;
				 ff <<"Integral"<<" "<<"counts"<<endl;
				 ofstream ffP (filename_PED);
				 //ffP <<"Juno PMT data - Pedestal"<<endl;
				 ffP <<"Histogram of Integral vs. Counts - Pedestal"<<endl;
				 ffP <<"No. of bins = "<<nbins<<endl;
				 ffP <<"position = "<< p <<" , angle = "<<15*a<<endl;
				 ffP <<"Integral"<<" "<<"counts"<<endl;
				 
				 for (int i=0; i <Juno->GetNbinsX(); i++)
					{
				       		if (Juno->GetBinCenter(i)<-10000) 
							{ff << Juno->GetBinCenter(i) << "	" << 0 << endl;}
					 	else
							{ff << Juno->GetBinCenter(i) << "	" << Juno->GetBinContent(i) << endl;} //write to file
				  	}
				  	ff.close();
				 for (int i=0; i <JunoPED->GetNbinsX(); i++)
					{
				        ffP << JunoPED->GetBinCenter(i) << "	" << JunoPED->GetBinContent(i) << endl; //write to file
				  	}
				  	ffP.close();
			}		
	}			
				//TApp.Run();
				return 0;
}
