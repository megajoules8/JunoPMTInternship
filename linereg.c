#include “TGraph.h”
#include “TFitResult.h”
#include “TMatrixD.h”

int QA2()
{
Double_t x[] = {1, 2, 3, 4, 5};
Double_t y[] = {1.1, 1.9, 3.2, 3.9, 5.5};
TGraph g = new TGraph((sizeof(x) / sizeof(Double_t)), x, y);
g->Draw("A");
TFitResultPtr r = g->Fit(“pol1”, “S”);
TMatrixD cov = r->GetCorrelationMatrix();
TMatrixD cor = r->GetCovarianceMatrix();
cov.Print();
cor.Print();
return 0;
}
