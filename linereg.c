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
Double_t y[] = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20};

TGraph g(10,x,y);
g.SetTitle("Measurement XYZ;lenght [cm];Arb.Units");
auto mycanvas = new TCanvas();
g.DrawClone("APE");
/*TFitResultPtr r = g->Fit(“pol1”, “S”);
r.DrawClone("Same");  
TMatrixD cov = r->GetCorrelationMatrix();
TMatrixD cor = r->GetCovarianceMatrix();
cov.Print();
cor.Print();*/
return 0;
}

int main(){
linereg();
}
