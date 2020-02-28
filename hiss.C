#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
using namespace std;

//main function
TH1F* hiss(TString path, TString filename){
  
  TString histname = filename;
  histname.ReplaceAll(".",""); //remove all "."s from histname
  histname.ReplaceAll("/",""); //remove all "/"s from histname
  ifstream read;
    read.open((path+filename).Data());
    //read.open("test.txt");
    double charge;
    double amplitude;
    double charge_min;
    double charge_max;
    double bin_width;
    int count=0;
//ignore first 5 lines
    for (int a=0; a<5; a++){
        read.ignore(1000000000000, '\n');
        }	
//read to find charge_max and charge_min
    while(read >> charge >> amplitude)
	{
     	//cout << "charge = " << charge << " " << "amplitude= " << amplitude << endl;
		
        if (count == 0)
		{
            charge_min = charge;
			//cout << "charge min = " << charge_min << endl;
        }
        if (count == 1)
		{
            bin_width = charge - charge_min;
		//	cout << "bin_width = " << bin_width << endl;
        }

        charge_max = charge;
        ++count;
    }
	
    //cout << charge_min << " " << charge_max <<" "<<bin_width<< endl;

//start calculations
	double charge_min_temp;
	charge_min_temp = charge_min;
    charge_min = -1*charge_max*1.0e+9;
    charge_max = -1*charge_min_temp*1.0e+9;
    bin_width  = bin_width*1.0e+9;
 
	//cout << charge_min << " " << charge_max <<" "<<bin_width<< "  " << count << endl;
//create and allocate histogram

   read.close();
   
   
//start histogram event		
  TH1F* hist = new TH1F(path+histname, histname, count , charge_min-(bin_width/2), charge_max+(bin_width/2));
	
	ifstream scan;
    scan.open((path+filename).Data());
//ignore first 5 lines	
	for (int a=0; a<5; a++){
        scan.ignore(1000000000000, '\n');
        }
//fill	histogram while calculating each element	
    while(scan >> charge >> amplitude){
      
      charge = -1*charge*1.0e+9;
	 // cout << "charge = " << charge << " " << "amplitude= " << amplitude << endl;
      hist -> Fill(charge, amplitude/2);
  }
//std::string name ="result_" + std::string (filename) + ".root";
//scan2.close();
  TFile* file = new TFile(path + "root_file.root","RECREATE"); //new ROOT file event
  hist->GetXaxis()->SetTitle("charge(nC)"); //set Xaxis title
  hist->GetYaxis()->SetTitle("amplitude"); //set Yaxis title
  hist->Draw();
  file->Write(); //write to ROOT file
  file->Close();
//modify and update the canvas  
  gPad->Modified(); 
  gPad->Update(); 
  gSystem->ProcessEvents();
  return hist;
  return (0);
  
}
