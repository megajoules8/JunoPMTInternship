#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"

void linereg()
{
Double_t x[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
Double_t y[] = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20};

TGraph g(10,x,y);
graph.SetTitle("Measurement XYZ;lenght [cm];Arb.Units");
auto mycanvas = new TCanvas();
g.DrawClone("APE");
TF1 f("Linear law","[0]+x*[1]",.5,10.5);
g.Fit(&f);
f.DrawClone("Same");
TMatrixD cov = f->GetCorrelationMatrix();
TMatrixD cor = f->GetCovarianceMatrix();
cov.Print();
cor.Print();
return 0;
}

int main(){
linereg();
}
