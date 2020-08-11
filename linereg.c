
int QA2()
{
Double_t x[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
Double_t y[] = {2, 4, 6, 8, 10, 12, 14, 16, 18, 20};
TGraph g = new TGraph((sizeof(x) / sizeof(Double_t)), x, y);
g->Draw("A");
TFitResultPtr r = g->Fit(“pol1”, “S”);
TMatrixD cov = r->GetCorrelationMatrix();
TMatrixD cor = r->GetCovarianceMatrix();
cov.Print();
cor.Print();
return 0;
}
