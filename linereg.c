#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TMatrixD.h"
#include <TFitResult.h>


void linereg()
{
Double_t x[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
Double_t y[] = {2, 5, 6, 7, 10, 11, 14, 16, 19.5, 20};
Double_t yerrs[] = {0.1, 0.15, 0.09, 0.22, 0.17, 0.05, 0.12, 0.11, 0.35, 0.2};  

TGraphErrors g(10,x,y, nullptr, yerrs);
g.SetTitle("Measurement XYZ;lenght [cm];Arb.Units");
g.SetMarkerStyle(kOpenCircle);
g.SetMarkerColor(kBlue);
g.SetLineColor(kBlue);
auto mycanvas = new TCanvas();
g.DrawClone("APE");
  
//TF1 f("Linear law","[0]+x*[1]",.5,10.5);
//f.SetLineColor(kRed); f.SetLineStyle(2);
//g.Fit(&f);
//f.DrawClone("Same");

TF1 *f1 = new TF1("f1","pol1",0.5,10.5);
g.Fit("f1","R");
f1->DrawClone("Same");
TFitResultPtr r = g.Fit("pol1","S"); 
// Double_t a = r->Value(0);
// Double_t b = r->Value(1);
// TF1 *Fn = new TF1("Fn",lin,0.5,10.5,2);
// Fn->SetNpx(500);
// Fn->SetParameter(0,a);
// Fn->SetParameter(1,b);
// Fn->DrawClone("Same");
TMatrixD cov = r->GetCorrelationMatrix();
TMatrixD cor = r->GetCovarianceMatrix();
cov.Print();
cor.Print();
return 0;
}

int main(){
linereg();
}
