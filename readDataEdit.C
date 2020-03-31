#include <vector>
#include <iostream>
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

	// for (int k =0; k<10; ++k)
	// {
	// 	cout<< time_vector[k]<<endl;
	// }
	std::vector<float> * wave_vector = 0;

	pmt_tree->SetBranchAddress("wave_vector", &wave_vector);

	TApplication TApp("TApp", &argc, argv);

	TCanvas * c = new TCanvas();

	int count = 0;
	float t_min = 0;
	float t_max = 0;
	float Integral = 0;
	float bin_width = 0;

//determining the bin width, the times (x value) corresponding to 400ns and 550ns, 
//and the no. of bins (count) within the range 

	for (int j=0; j < wave_vector->size(); ++j)
			{
				if (( time_vector[j] > 400 ) && ( time_vector[j] < 550))
					{
						if (count ==0)
						{
							cout<<"hall!"<<endl;
							//t_min =  time_vector[j];
						}
						if (count ==1)
						{
							//bin_width =  time_vector[j]- t_min;
							cout<<"hall! hall!"<<endl;
						}

						++count;
						//t_max = time_vector[j];
					}
				
			}


//definition of the histogram using previous findings
	TH1F* Juno = new TH1F("Juno", "Juno", count , t_min - (bin_width/2), t_max +(bin_width/2));

	for(int i=0; i < pmt_tree->GetEntries(); ++i)
		{

		pmt_tree->GetEntry(i);

		TVectorT<float> wave_vector_root(wave_vector->size());
		for(int iWV=0; iWV < wave_vector->size(); ++iWV)
			{
				wave_vector_root[iWV] = wave_vector->at(iWV);
				if (( time_vector[iWV] >= t_min ) && ( time_vector[iWV]<= t_max ))
					{
						Integral += wave_vector->at(iWV);
					}
			}

		cout<<"Integral at Entry no: "<< i << " = "<< Integral <<endl;
		TGraph * g = new TGraph(*time_vector, wave_vector_root);
		g->GetXaxis()->SetTitle("Time (ns)");
		g->GetYaxis()->SetTitle("Amplitude");
		g->SetTitle(TString::Format("Event %d", i));
		g->SetMarkerSize(.5);
		g->SetMarkerStyle(24);
		g->SetMarkerColor(kBlue);
		g->Draw("ALP");

		c->Update();
		c->WaitPrimitive();
		delete g;
	}

	TApp.Run();
	return 0;
}
