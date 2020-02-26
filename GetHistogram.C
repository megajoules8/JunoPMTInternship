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


//Global parameter
int nbins = 50000; //just a starting point, can be any value 


//Main function
TH1F* GetHistogram(TString path, TString filename, int binwidth){

  //build histogram titles
  TString histname = filename;
  histname.ReplaceAll(".","");
  histname.ReplaceAll("/","");

 //create and allocate histogram
  TH1F* hist = new TH1F(histname,histname,nbins,0,nbins); //we want to automate it and count number of lines before made into a histogram
 
	// so I think we need to move this to open the data file, read it, skip the five lines, then have the code execute creating a 
	//histogram
  //open data file 
  ifstream read;
  read.open((path+filename).Data()); 
  double counts=-1;  

  //read data file
  int bin = -1 ;
  while (read>>counts){ //reading
    for(int i = 0 ; i < counts; i++) //filling histogram hist to find the no. of bins
    // hist->Fill(bin);  
    hist->SetBinContent(bin,counts); 
    bin++;
  }
  
  
//read the maximum and minimum values of data  

double number, max = INT_MIN, min = INT_MAX; 
  ifstream read34;
  read34.open((path+filename).Data());
 
  while (read34>>number)
	{
		if ((number>max)&&(number<1))
		     max=number;
		
		if (number<min)
		     min=number;
	}
  
 cout <<  " Max = " << max << " min = " << min << endl; 
 
//print number of read lines
  cout << bin << " bins were found in file " << filename << endl;

//defining new histogram hist12 
   TH1F* hist12 = new TH1F(histname,histname,bin,min,max);
 
  //open data file
  ifstream read12;
  read12.open((path+filename).Data());
  double counts12=-1;  

  //read data file
  int bins = -1 ;
  while (read12>>counts12){ //reading
    for(int j = 0 ; j < counts12; j++) //filling histogram
    // hist->Fill(bin);  
    hist12->SetBinContent(bins,counts12); //filling histogram, alternative method
    bins++;
  }

//draw histogram
  hist12->Draw();  
  return hist12;
  
 
}
