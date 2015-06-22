{


  //alluminio
  //TFile* file = new TFile("../roottople/aluminum_new.root");

  //titanio
  //TFile* file = new TFile("../roottople/titanium_new.root");

  //mylar
  // TFile* file = new TFile("../roottople/run684psi.root");
  
  //steel
  // TFile* file = new TFile("../roottople/run687psi.root");

  //pb
  // TFile* file = new TFile("../roottople/run690psi.root");

  //h2o
   TFile* file = new TFile("../roottople/run695psi.root");

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();

  Waveform->Draw("b3_tmean6-b2_tmean6>>h_time(100,-1000,3000)","b3_tmean6>0&&b2_tmean6>0&&b1_tmean7==-999&&b2_charge7<50");
  h_time->SetTitle("Time difference Sum 9 - Veto_Hole; Time [ns]; Entries / 40 ns");
  TF1 *fit0 = new TF1("fit0","pol0",-1000.,-50.);
  h_time->Fit("fit0","RBO");
  double costante = fit0->GetParameter(0);
  TF1 *fit = new TF1("fit","([0]+[1]*exp(-(x-[2])/[3]))", 50., 4000);
  
 fit->SetParName(0,"offset");
 fit->SetParName(1,"N_{0}");
 fit->SetParName(2,"#tau");
 fit->SetParName(3,"t_{0}");
 fit->SetParameter(0,200.);
 fit->SetParameter(1,500);
 fit->SetParameter(2,600);
 fit->FixParameter(0,costante);
 h_time->Fit("fit","RBO"); 
 h_time->Draw();
 fit0->Draw("same");
 
 TCanvas *c2 = new TCanvas("c2","c2");
 c2->cd();
 Waveform->Draw("b3_tmean6-b1_tmean7>>h_time2(1000,-1000,3000)","b3_tmean6>0&&b2_tmean6==-999&&b1_tmean7>0&&b2_charge7<50&&b3_charge5<-2000");
 TF1 *fit02 = new TF1("fit02","pol0",-1000.,-50.);
 h_time2->Fit("fit02","RBO");
 h_time2->SetTitle("Time difference Sum 9 - muon beam (charged particles); Time [ns]; Entries / 40 ns");
 double costante2 = fit02->GetParameter(0);
 TF1 *fit2 = new TF1("fit2","([0]+[1]*exp(-(x)/[2]))", -10., 4000);
 fit2->SetParName(0,"offset");
 fit2->SetParName(1,"N_{0}");
 fit2->SetParName(2,"#tau");
 fit2->SetParameter(0,200.); 
 fit2->FixParameter(0,costante2);
 fit2->SetParameter(1,500.);
 fit2->SetParameter(2,60);
 
 h_time2->Fit("fit2","RBO");
 h_time2->Draw();
 fit02->Draw("same");
}
