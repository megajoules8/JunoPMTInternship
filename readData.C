#include <vector>
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TVectorT.h>
#include <TApplication.h>

int main(int argc, char ** argv){
	if(argc != 2){
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

	for(int i=0; i < pmt_tree->GetEntries(); ++i){
		pmt_tree->GetEntry(i);

		TVectorT<float> wave_vector_root(wave_vector->size());
		for(int iWV=0; iWV < wave_vector->size(); ++iWV){
			wave_vector_root[iWV] = wave_vector->at(iWV);
		}

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
