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
g.SetMarkerStyle(kOpenCircle);
g.SetMarkerColor(kBlue);
g.SetLineColor(kBlue);
auto mycanvas = new TCanvas();
g.DrawClone("APE");
TF1 f("Linear law","[0]+x*[1]",.5,10.5);
f.SetLineColor(kRed); f.SetLineStyle(2);
g.Fit(&f);
f.DrawClone("Same");
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
