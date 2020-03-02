#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
using namespace std;

//main function
TH1F* hiss(TString Full_path){
  
  TString histname = Full_path;
  histname.ReplaceAll(".",""); //remove all "."s from histname
  histname.ReplaceAll("/",""); //remove all "/"s from histname
  ifstream read;
    read.open(Full_path);
    //read.open("test.txt");
	//cout << Full_path;
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
	//int i=0;
 
	//cout << charge_min << " " << charge_max <<" "<<bin_width<< "  " << count << endl;
//create and allocate histogram

   read.close();
    
//start histogram event		
  TH1F* hist = new TH1F(Full_path, histname, count , charge_min-(bin_width/2), charge_max+(bin_width/2));
	ifstream scan;
    scan.open(Full_path);
//ignore first 5 lines	
	for (int a=0; a<5; a++)
		{
        scan.ignore(1000000000000, '\n');
        }
		
//fill	histogram while calculating each element	



   while(scan >> charge >> amplitude){
      
      charge = -1*charge*1.0e+9;
	 // cout << "charge = " << charge << " " << "amplitude= " << amplitude << endl;
      hist -> Fill(charge, amplitude);
  }
  
  	for ( int l=1; l <= hist->GetXaxis()->GetNbins(); l++ ) 
		{
		
		
		// hist -> Fill(charge, amplitude);
		hist -> SetBinError(l, sqrt( hist->GetBinContent( l ) ));
		
		//cout << hist->GetBinContent( l ) << endl;
		//if ( hist->GetBinContent( l )!=0 )getchar();
		
		// Double_t BC = hist -> GetBinContent(l);
        // Double_t BE = hist -> GetBinError(l);		
		// cout << "charge = " << charge << " " << "amplitude= " << amplitude << " error = "<< BE <<" Bin content = "<< BC<<endl;	
		// if ( BC!=0 )getchar();
		//hist->Sumw2();
		}
		
		
  //hist->Sumw2();
//std::string name ="result_" + std::string (filename) + ".root";
//scan2.close();
  //TFile* file = new TFile(path + "root_file.root","RECREATE"); //new ROOT file event
  //hist->GetXaxis()->SetTitle("charge(nVs)"); //set Xaxis title
  //hist->GetYaxis()->SetTitle("amplitude"); //set Yaxis title
  //hist->Draw();
  
 // if (num ==0)
  // {
	// ofstream ff ("Q_sigma.txt");
	// double Q = hist->GetMean();
	// double sigma = hist->GetRMS();
	// cout<< "This is PED"<<endl;
	// cout <<"Q = "<< Q << " sigma = "<< sigma <<endl;
	// ff   <<Q << endl;
	// ff   << sigma <<endl;
  // }
 
 // if (num==1)
 // {
	// cout<< "This is LED"<<endl; 
 // }
 
  // file->Write(); //write to ROOT file
  // file->Close();
// //modify and update the canvas  
  // gPad->Modified(); 
  // gPad->Update(); 
  // gSystem->ProcessEvents();
  return hist;
  //return (0);
  
}


// for ( int l=1; l<=hSG->GetXaxis()->GetNbins(); l++ )
// hSG->SetBinError( l, sqrt( hSG->GetBinContent( l ) ) );
