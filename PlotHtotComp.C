
void PlotHtotComp(TString fn1="", TString fn2="", Double_t nfct=1.0, TString lgnd1="", TString lgnd2=""){

  TFile *f1 = TFile::Open(fn1.Data());  
  TFile *f2 = TFile::Open(fn2.Data());  

  TH1F *h1 = (TH1F*) f1->Get("h_total");
  TH1F *h2 = (TH1F*) f2->Get("h_total");
  h1->Draw();
  h2->SetLineColor(kRed);
  h2->Scale(nfct);
  h2->Draw("sames");

  leg = new TLegend(0.5,0.55,0.9,0.75);
  leg->AddEntry(h1,lgnd1,"l");
  leg->AddEntry(h2,lgnd2,"l");
  leg->SetFillColor(kWhite);
  leg->Draw("same");

}
