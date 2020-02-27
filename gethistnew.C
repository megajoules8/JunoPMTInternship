#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
using namespace std;


TH1F* hiss(TString path, TString filename){
	
  TString histname = filename;
  histname.ReplaceAll(".","");
  histname.ReplaceAll("/","");
  ifstream read;
    read.open((path+filename).Data());
    //read.open("test.txt");
    double charge;
    double amplitude;
    double charge_min;
    double charge_max;
    double bin_width;
    int count=0;

    for (int a=0; a<5; a++){
        read.ignore(1000000000000, '\n');
        }	

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
    cout << charge_min << " " << charge_max <<" "<<bin_width<< endl;

	double charge_min_temp;
	charge_min_temp = charge_min;
    charge_min = -1*charge_max*1000000000 / 50;
    charge_max = -1*charge_min_temp*1000000000 / 50;
    bin_width  = bin_width*1000000000 / 50;
 
	cout << charge_min << " " << charge_max <<" "<<bin_width<< endl;
//create and allocate histogram
  TH1F* hist = new TH1F(histname, histname, count , charge_min-bin_width/2, charge_max-bin_width/2);


	ifstream scan;
    scan.open((path+filename).Data());
  while(scan >> charge >> amplitude){
      
      charge = -1*charge*1000000000  / 50;
	  cout << "charge = " << charge << " " << "amplitude= " << amplitude << endl;
      hist -> Fill(charge, amplitude);
  }
  //binwidth histogram and draw
  //hist->Rebin(binwidth);
  hist->Draw();
  
  return hist;
  return (0);
}
