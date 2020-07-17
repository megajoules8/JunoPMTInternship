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
	/*if(argc != 2)
	{
		std::cout << "Usage: " << argv[0] << " file.root" << std::endl;
		return 1;
	}*/
	int pos;
	int ang;
	int sc;
	cout << "Input the dataset you wish to analyze:"<< endl;
	cout << "248 for scan248" << endl;
	cout << "846 for scan846" << endl;
	cout << "3712 for scan3712" << endl;
	cout << "248 for scan248" << endl;
	cout << "3737 for scan3737" << endl;
	cout << "3899 for scan3899" << endl;
	cout << "4050 for scan4050" << endl;
	cout << "4230 for scan4230" << endl;
	cout << "4232 for scan4232" << endl;
	cin >> sc ;
	cout << " " << endl;
	cout << "Input the position" << endl;
	cin >> pos ;
	cout << " " << endl;
	cout << "Input the angle" << endl;
	cin >> ang ;
	cout << " " << endl;
	
	TString fname = argv[1] + TString("/scan") + Form("%d",sc) + TString ("_position") + Form("%d",pos) + TString("_angle") + Form ("%d", ang) + TString(".root");
	TFile * f = TFile::Open(fname);

	TTree * pmt_tree = (TTree*) f->Get("pmt_tree");
	TVectorT<float> * time_vector = (TVectorT<float> *) f->Get("time_vector");
	std::vector<float> * wave_vector = 0;
	pmt_tree->SetBranchAddress("wave_vector", &wave_vector);
	TApplication TApp("TApp", &argc, argv);

	TCanvas * c = new TCanvas();

	int nbins = 2000;
	float t_min = 0;
	float t_max = 1000;
	//float Integral = 0;
	float bin_width = 0;

//definition of the histogram
TH1F *Juno = new TH1F("Juno", "Event Integration: Position = 1, angle = 0", nbins , t_min, t_max);

	for(int i=0; i < pmt_tree->GetEntries(); ++i)
		{

		pmt_tree->GetEntry(i);

		TVectorT<float> wave_vector_root(wave_vector->size());
		for(int iWV=0; iWV < wave_vector->size(); ++iWV)
			{
				wave_vector_root[iWV] = wave_vector->at(iWV);
				Juno-> Fill(iWV, wave_vector_root[iWV]); 

			}
			
		//cout<<"Integral at Entry no: "<< i << " = "<< Integral <<endl;
		//TGraph * g = new TGraph(*time_vector, wave_vector_root);
		//g->GetXaxis()->SetTitle("Time (ns)");
		//g->GetYaxis()->SetTitle("Amplitude");
		Juno->GetXaxis()->SetTitle("Time(ns)");
		Juno->GetYaxis()->SetTitle("Total Counts");
		//g->SetTitle(TString::Format("Event %d", i));

		//g->Draw("ALP");
		
		//c->Update();
		//c->WaitPrimitive();
		//delete g;
		//Integral = 0;
	}
		Juno->SetMarkerSize(1.0);
		Juno->SetMarkerStyle(24);
		Juno->SetMarkerColor(kBlue);
	Juno->Draw();
	 // ofstream ff ("Juno_data.txt");
	 // ff <<"Juno PMT data"<<endl;
	 // ff <<"Histogram of Integral vs. Counts"<<endl;
	 // ff<<"No. of bins = "<<nbins<<endl;
	 // ff<<"************************************"<<endl;
	 // ff<<"Integral"<<" "<<"counts"<<endl;
	 
	 // for (int i=0; i <Juno->GetNbinsX(); i++)
		// {
	 //        ff << Juno->GetBinCenter(i) << "	" << Juno->GetBinContent(i) << endl; //write to file
	 //  	}
	 //  	ff.close();


	TApp.Run();
	return 0;
}
