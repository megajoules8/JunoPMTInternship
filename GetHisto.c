
#include <fstream>
#include <bits/stdc++.h>
#include <vector>
#include <stdlib.h>
using namespace std;

//Global parameter
int bins = 50000;


//Main function
TH1F* GetHisto(TString path, TString filename, int binwidth){

  //build histogram titles
  TString histname = filename;
  histname.ReplaceAll(".","");
  histname.ReplaceAll("/","");
 
  //count no. of rows 
int rows = 0, cols = 0;
   string line, item;

   ifstream file( (path+filename).Data());
   while ( getline( file, line ) )
   {
      rows++;
      if ( rows == 1 )                 // First row only: determine the number of columns
      {
         stringstream ss( line );      // Set up up a stream from this line
         while ( ss >> item ) cols++;  // Each item delineated by spaces adds one to cols
      }

       // .... //                      // Do any processing on a line here

      //cout << "Line read is " << line << endl;
   }
   file.close();

   cout << "\nFile had " << rows-5 << " rows and " << cols-1 << " columns" << endl;
 
//fill data into arrays
 
double time[rows-5];
double amplitude[rows-5];
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

//finding the max and min values
double min = INT_MAX, max = INT_MIN;
int b=0;	
  for (int b=0; b<rows-5; ++b)
	{
		if (time[b] > max)
		     max = time[b];
		
		if (time[b]<min)
		     min=time[b];
	}
  
 cout <<  " Max = " << max << " min = " << min << endl; 

//sort into new arrays
 std::string name ="result_" + std::string (filename) + ".txt";
 ofstream FF (name);
 vector< pair <double,int> > vect;
 int n = sizeof(time)/sizeof(time[0]);
 for (int c=0; c<n; c++) 
        vect.push_back( make_pair(time[c],amplitude[c]) );

sort(vect.begin(), vect.end()); 
  
     // Printing the sorted vector(after using sort()) 
    cout << "The vector after sort operation is:\n" ; 
    for (int d=0; d<n; d++) 
    { 
        // "first" and "second" are used to access 
        // 1st and 2nd element of pair respectively 
        FF << vect[d].first << " "<< vect[d].second << endl; 
			 
		 
    } 
 FF.close();	
 
 
 
 //create and allocate histogram
  TH1F* hist = new TH1F(histname,histname,2*(rows-5),min,max);
 
  //open data file
  ifstream read;
  
  read.open(name);
  double counts=-1;  
  
  //read data file
  int bin = -1 ;

  while (read>>counts){ //reading
   //for(int l = 0 ; l < counts; l++) //filling histogram
   //   hist->Fill(bin);  
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

