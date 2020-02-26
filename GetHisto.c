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
//Global parameter
int bins = 50000;


//Main function
TH1F* GetHisto(TString path, TString filename, int binwidth){

//COUNTING THE TOTAL NUMBER OF ENTRIES
ifstream f( (path+filename).Data() );
	int number_of_lines;
    std::string line;
    std::ifstream myfile((path+filename).Data());

    while (std::getline(myfile, line))
	{
        ++number_of_lines;
	}

    std::cout << "Number of lines in text file: " << number_of_lines<< endl;


 //IGNORING THE FIRST 5 LINES
 
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

//SAVING DATA INTO SEPERATE ARRAYS WHILE CALCULATING

 while (!infile.eof()) {
	  infile >> time[i] >> amplitude[i];
	  time[i] = time[i]*(-1000000000/50);
	 // cout <<"time "<< time[i] << " amplitude "<< amplitude[i] <<endl;
 	i++;
    }
//FINDING THE MAXIMUM AND MINIMUM

double number, max = INT_MIN, min = INT_MAX;
for (int b=0; b<i; ++b)
	{
		if (time[b]>max)
		     max=time[b];
		
		if (time[b]<min)
		     min=time[b];
	}
cout <<  " Max = " << max << " min = " << min << endl; 

//SORT INTO NEW ARRAYS
int nbins = 3999;
TString histo = filename;
 

  //create and allocate histogram
  TH1F* hist = new TH1F(histo,histo,nbins,min,max);
 
  //open data file
  infile.open((path+filename).Data());
  double counts=-1;  
  
  //read data file
  int bin = -1 ;
  while (read>>counts){ //reading
    for(int k = 0 ; k < counts; k++) //filling histogram
      hist->Fill(bin);  
    //hist->SetBinContent(bin,counts); //filling histogram, alternative method
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

