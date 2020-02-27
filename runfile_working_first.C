
#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include "Riostream.h"
#include <vector>
#include <TString.h>
#include "hiss.C"
using namespace std;

void runfile_working_first(TString path, TString filename)
	{ 
		
	vector< vector<TString> > PED_LED;
	vector<TString> a;
	a.clear();
	a.push_back("/PED/");
	a.push_back("/LED/");
	PED_LED.push_back(a);
	a.clear();
	
	vector< vector<TString> > HV;
	vector<TString> b;
	b.clear();
	b.push_back("1250");
	b.push_back("1300");
	b.push_back("1350");
	b.push_back("1401");
	b.push_back("1450");
	HV.push_back(b);
	
	int i=0;
	int j=0;
	int PED_LED_Count =2;
	int HV_Count = 5;
	while (i<HV_Count)
	{
		while (j<PED_LED_Count)
		{
			hiss(path + HV[0][i] + PED_LED[0][j] , filename);
			++j;
		}
		++i;
		j=0;
	}

	
	}