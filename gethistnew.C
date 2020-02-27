#include <iostream>
#include <string>
#include <fstream>
using namespace std;


int main(){

    string path;
    string filename;
    ifstream read;
//    read.open((path+filename).Data());
    read.open("test.txt");
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
     	cout << "charge = " << charge << " " << "amplitude= " << amplitude << endl;
		
        if (count == 0)
		{
            charge_min = charge;
			cout << "charge min = " << charge_min << endl;
        }
        if (count == 1)
		{
            bin_width = charge - charge_min;
			cout << "bin_width = " << bin_width << endl;
        }

        charge_max = charge;
        ++count;
    }
    cout << charge_min << " " << charge_max <<" "<<bin_width<< endl;

return (0);
}
