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

	TFile * f = TFile::Open(argv[1]);

	TTree * pmt_tree = (TTree*) f->Get("pmt_tree");
	TVectorT<float> * time_vector = (TVectorT<float> *) f->Get("time_vector");
	std::vector<float> * wave_vector = 0;
	pmt_tree->SetBranchAddress("wave_vector", &wave_vector);
	TApplication TApp("TApp", &argc, argv);

	TCanvas * c = new TCanvas();

	int nbins = 2000;
	float t_min = 400;
	float t_max = 500;
	float t_min_PED = 0;
	float t_max_PED = 100;
	float Integral = 0;
	float Integral_PED = 0 ;
	float bin_width = 0;

//definition of the histogram
TH1F *Juno = new TH1F("Juno", "Juno", nbins , -8000, 1000);
TH1F *JunoPED = new TH1F("JunoPED", "JunoPED", nbins , -1000, 1000);

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
		//cout<<"Integral at Entry no: "<< i << " = "<< Integral <<endl;
		//TGraph * g = new TGraph(*time_vector, wave_vector_root);
		//g->GetXaxis()->SetTitle("Time (ns)");
		//g->GetYaxis()->SetTitle("Amplitude");
		Juno->GetXaxis()->SetTitle("Integral");
		Juno->GetYaxis()->SetTitle("Counts");
		JunoPED->GetXaxis()->SetTitle("Integral");
		JunoPED->GetYaxis()->SetTitle("Counts");
		//g->SetTitle(TString::Format("Event %d", i));
		//g->SetMarkerSize(.5);
		//g->SetMarkerStyle(24);
		//g->SetMarkerColor(kBlue);
		//g->Draw("ALP");

		//delete g;
		Integral = 0;
		Integral_PED = 0;
	}

		JunoPED->Draw();
		c->Update();
		c->WaitPrimitive();

		Juno->Draw();

	 ofstream ff ("Juno_data.txt");
	 ff <<"Juno PMT data"<<endl;
	 ff <<"Histogram of Integral vs. Counts"<<endl;
	 ff <<"No. of bins = "<<nbins<<endl;
	 ff <<"************************************"<<endl;
	 ff <<"Integral"<<" "<<"counts"<<endl;

	 ofstream ffP ("Juno_data_PED.txt");
	 ffP <<"Juno PMT data - Pedestal"<<endl;
	 ffP <<"Histogram of Integral vs. Counts - Pedestal"<<endl;
	 ffP <<"No. of bins = "<<nbins<<endl;
	 ffP <<"************************************"<<endl;
	 ffP <<"Integral"<<" "<<"counts"<<endl;
	 
	 for (int i=0; i <Juno->GetNbinsX(); i++)
		{
	        ff << Juno->GetBinCenter(i) << "	" << Juno->GetBinContent(i) << endl; //write to file
	  	}
	  	ff.close();

	 for (int i=0; i <JunoPED->GetNbinsX(); i++)
		{
	        ffP << JunoPED->GetBinCenter(i) << "	" << JunoPED->GetBinContent(i) << endl; //write to file
	  	}
	  	ffP.close();


	TApp.Run();
	return 0;
}
