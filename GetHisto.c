/*
This macro creates a root histogram from an ascii data file of format
---------------------------
Channel0 CountsPerChannel0
Channel1 CountsPerChannel1
Channel2 CountsPerChannel2
etc...
---------------------------

Execution example:
.x GetHistogram.C("../path/to/the/file/","name-of-the-file.asc",nb-of-bins-ADC, nb-of-merged-bins )

Practical info:
nb-of-bins-ADC = pow(2,15) = 32768,         /!\  FIXED by the ADC
nb-of-merged-bins = binwidth = pow(2,n)     /!\  n = User defined integer

*/

#include <iostream>
#include <bits/stdc++.h>
using namespace std;
//Global parameter
int bins = 50000;


//Main function
TH1F* GetHisto(TString path, TString filename, int binwidth){

  //build histogram titles
  TString histname = filename;
  histname.ReplaceAll(".","");
  histname.ReplaceAll("/","");
  
  
  
ifstream f( (path+filename).Data() );

double number, max = INT_MIN, min = INT_MAX; 

	int number_of_lines;
    std::string line;
    std::ifstream myfile((path+filename).Data());

    while (std::getline(myfile, line))
	{
        ++number_of_lines;
	}

    std::cout << "Number of lines in text file: " << number_of_lines<< endl;
   
double time[number_of_lines];
double amplitude[number_of_lines];
double sorted_time[number_of_lines];
double sorted_amplitude[number_of_lines];
int i=0;
 
 ifstream infile;
    infile.open((path+filename).Data());
	
for (int a=0; a<5; a++)
{
 infile.ignore(1000000000000, '\n');
}	

 while (!infile.eof()) {
	  infile >> time[i] >> amplitude[i];
	  time[i] = time[i]*(-1000000000/50);
	  //cout <<"time "<< time[i] << " amplitude "<< amplitude[i] <<endl;
 	++i;
    }


for (int b=0; b<i-1; ++b)
	{
		if (time[b]>max)
		     max=time[b];
		
		if (time[b]<min)
		     min=time[b];
	}
cout <<  " Max = " << max << " min = " << min << endl; 
//sort into new arrays
 
 //define size of array
  
  int n = sizeof(time)/sizeof(time[0]); 
    sort(time, time+n); //adding a sorting function to time
   // cout << "\nArray after sorting using "
     //    "default sort is : \n"<<endl;  
    //for (int i = 0; i < n; ++i)
        //cout << time[i] << " "<<endl;
    









 //create and allocate histogram
  TH1F* hist = new TH1F(histname,histname,2*(i-1),min,max);
 
  //open data file
  ifstream read;
  read.open("LED1350.txt");
  double counts=-1;  
  
  //read data file
  int bin = -1 ;

  while (read>>counts){ //reading
   // for(int l = 0 ; l < counts; l++) //filling histogram
      //hist->Fill(bin);  
    hist->SetBinContent(bin,counts); //filling histogram, alternative method
    bin++;
  }
  
  //print number of read lines
  cout << bin << " bins were found in file " << filename << endl;
  
  //binwidth histogram and draw
  //hist->Rebin(binwidth);
  hist->Draw();
  
  return hist;

return(0);
}

