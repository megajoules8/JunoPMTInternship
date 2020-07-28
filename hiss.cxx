#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TStyle.h>
#include "hiss.h"
using namespace std;
//main function
TH1F* hiss(TString Full_path, TString HV_Value, int index, int sc )
{
  
  TString histname = Full_path;
//  histname.ReplaceAll(".",""); //remove all "."s from histname
//  histname.ReplaceAll("/",""); //remove all "/"s from histname
  ifstream read;
    read.open(Full_path);
    if (!read){
      cout << "the file does not exist" << endl;
    } 
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
    if (index == 0)
    {
      charge_min_temp = charge_min;
      charge_min = -1*charge_max;
      charge_max = -1*charge_min_temp;
      //int i=0;
    }
    else if (index == 1)
    {
    charge_min_temp = charge_min;
    charge_min = -1*charge_max*1.0e+9;
    charge_max = -1*charge_min_temp*1.0e+9;
    bin_width  = bin_width*1.0e+9;
    }
	//cout << charge_min << " " << charge_max <<" "<<bin_width<< "  " << count << endl;
//create and allocate histogram
   read.close();
    
//start histogram event	//The new argument has been set as the name of the histogram	
	//new edit: changing the integration window to [-0.1,0.55] (observing previous fit data)
  TH1F* hist = new TH1F(HV_Value, HV_Value, count , charge_min - (bin_width/2), charge_max +(bin_width/2));
   gStyle->SetOptFit(1111); //used to display the stat window on the canvas
	ifstream scan;
    scan.open(Full_path);
//ignore first 5 lines	
	for (int a=0; a<5; a++)
		{
        scan.ignore(1000000000000, '\n');
        }
		
//fill	histogram while calculating each element	
   while(scan >> charge >> amplitude)
   {
      
      if (index == 1)
      {
        charge = charge*1.0e+9;
      }
      // cout << "charge = " << charge << " " << "amplitude= " << amplitude << endl;
      hist -> Fill(-1*charge, amplitude);
  }
  int Bin_Size;
 //rebinning 
    if (index == 1)	{Bin_Size = 8;}
    if (index == 0)	{Bin_Size = 4;}
    if (sc == 3737)	{Bin_Size = 1;}
    if (sc == 3899)	{Bin_Size = 4;}
    if (sc == 248) 	{Bin_Size = 25;}
    if (sc == 4050)	{Bin_Size = 10;}
    if (sc == 4232)	{Bin_Size = 2;}
    hist-> Rebin(Bin_Size); 
	
  	for ( int l=1; l <= hist->GetXaxis()->GetNbins(); l++ ) 
		{
		
		
		// hist -> Fill(charge, amplitude);
		hist -> SetBinError(l, sqrt( hist->GetBinContent( l ) ));
		}
		
		
  return hist;
  //return (0);
  
}
// for ( int l=1; l<=hSG->GetXaxis()->GetNbins(); l++ )
// hSG->SetBinError( l, sqrt( hSG->GetBinContent( l ) ) );
