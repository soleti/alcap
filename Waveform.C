 #define Waveform_cxx
#include "Waveform.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1F.h"
#include "TProfile.h"
#include "Riostream.h"
void Waveform::Loop(Int_t typeCut, Int_t sum, TString postfix)
{
//   In a ROOT session, you can do:
//      Root > .L Waveform.C
//      Root > Waveform t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
     if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   
   double matrix_data[5][5];
   double matrix_time[5][5];
   double beam_time;

   double sum_8, sum_16, diff_8, diff_16, sum16_data;
   
   TH1F *h_cells[5][5];
   for(int JLoop=0; JLoop<5; JLoop++){
       for(int ILoop=0; ILoop<5; ILoop++){ 
           h_cells[ILoop][JLoop] = new TH1F(Form("h_%d_%d", ILoop, JLoop), Form("Cell %d,%d;Energy [MeV];Entries / 1 MeV", ILoop, JLoop), 45, -5, 40);
       }
   }
   
   TH1F *h_inner = new TH1F("h_inner", "Inner ring; Energy [MeV]; Entries / MeV", 500, -10, 490);
   TH1F *h_outer = new TH1F("h_outer", "Outer ring; Energy [MeV]; Entries / MeV", 500, -10, 490);
   TH1F *h_sum9 = new TH1F("h_sum9", "Sum 9; Energy [MeV]; Entries / MeV", 500, -10, 490);
  TH1F *h_sum16 = new TH1F("h_sum16", "Sum 16; Energy [MeV]; Entries / MeV", 500, -10, 490);
   TH1F *h_total = new TH1F("h_total",Form("Entire matrix (Run %i); Energy [MeV];Entries / 1 MeV",run_number), 500, -10, 490);
   TH1F *h_diff_inn = new TH1F("h_diff_inn", "diff inner ; Charge [pC]; Entries / 2 pC", 400, -400, 400);
   TH1F *h_diff_outer = new TH1F("h_diff_outer", "diff outer ; Charge [pC]; Entries / 2 pC", 400, -400, 400);
   TH1F *h_trg_cosmic = new TH1F("h_trg_cosmic", "h_trg_cosmic ; Charge [pC]; Entries / 10 pC", 800, -4000, 4000);

   TH1F *h_time_diff = new TH1F("h_time_diff", "h_time_diff ; Time [ns]; Entries / 1 ns", 4000, -1000, 3000);

   TH2F *h_energy_time = new TH2F("h_energy_time", "; Energy [MeV]; t_{TSc} - t_{LYSO}", 250, 0, 250, 250, -1000, 3000);
   //prova raffa
TH2F *h_energy_time2 = new TH2F("h_energy_time", "; Energy [MeV]; t_{TSc} - t_{LYSO}", 250, 0, 250, 250, -1000, 3000);

   TProfile *correlation_8 = new TProfile("correlation_8","correlation_8; Sum8 [pC] - EneSum8 [pC]",100,-300,1700,-300,1700,"");
    TProfile *correlation_16 = new TProfile("correlation_16","correlation_16; Sum16 [pC] - EneSum16 [pC]",100,-300,1700,-300,1700,"");  
   int i,j,k;
   //nentries = 10000;
   //519349 
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      
      // matrix_data[ROW][COLUMN] rear view, [0][0] is the bottom left corner, as in the logbook
      matrix_data[0][0] = b2_charge5;
      matrix_data[1][0] = b2_charge4;
      matrix_data[2][0] = b2_charge3;
      matrix_data[3][0] = b2_charge1;
      matrix_data[4][0] = b2_charge0;
  
      matrix_data[0][1] = b0_charge4;
      matrix_data[1][1] = b0_charge3;
      matrix_data[2][1] = b0_charge2;
      matrix_data[3][1] = b0_charge1;
      matrix_data[4][1] = b0_charge0;
  
      matrix_data[0][2] = b1_charge1;
      matrix_data[1][2] = b1_charge0;
      matrix_data[2][2] = b0_charge7;
      matrix_data[3][2] = b0_charge6;
      matrix_data[4][2] = b0_charge5;
  
      matrix_data[0][3] = b1_charge6;
      matrix_data[1][3] = b1_charge5;
      matrix_data[2][3] = b1_charge4;
      matrix_data[3][3] = b1_charge3;
      matrix_data[4][3] = b1_charge2;
  

      // WRONG SIGN IN MAKEROOTTOPLA
      matrix_data[0][4] = -b3_charge4;
      matrix_data[1][4] = -b3_charge3;
      matrix_data[2][4] = -b3_charge2;
      matrix_data[3][4] = -b3_charge1;
      matrix_data[4][4] = -b3_charge0;

      matrix_time[0][0] = b2_tmean5;
      matrix_time[1][0] = b2_tmean4;
      matrix_time[2][0] = b2_tmean3;
      matrix_time[3][0] = b2_tmean1;
      matrix_time[4][0] = b2_tmean0;
  
      matrix_time[0][1] = b0_tmean4;
      matrix_time[1][1] = b0_tmean3;
      matrix_time[2][1] = b0_tmean2;
      matrix_time[3][1] = b0_tmean1;
      matrix_time[4][1] = b0_tmean0;
  
      matrix_time[0][2] = b1_tmean2;
      matrix_time[1][2] = b1_tmean0;
      matrix_time[2][2] = b0_tmean7;
      matrix_time[3][2] = b0_tmean6;
      matrix_time[4][2] = b0_tmean5;
  
      matrix_time[0][3] = b1_tmean6;
      matrix_time[1][3] = b1_tmean5;
      matrix_time[2][3] = b1_tmean4;
      matrix_time[3][3] = b1_tmean3;
      matrix_time[4][3] = b1_tmean2;
  
      matrix_time[0][4] = b3_tmean4;
      matrix_time[1][4] = b3_tmean3;
      matrix_time[2][4] = b3_tmean2;
      matrix_time[3][4] = b3_tmean1;
      matrix_time[4][4] = b3_tmean0;

      beam_time = b1_tmean7;

      sum_8 = b3_charge6;
      sum_16 = b3_charge7;

      h_trg_cosmic->Fill(b2_charge7);
      bool trigger = true;
            
      if (trigger) {
	double outer = 0, inner = 0, total = 0, sum9 = 0;
	double time9 = 0;
          for (int i = 0; i < 5; i++) {
              outer += matrix_data[i][0];
              outer += matrix_data[i][4];
          }
      
          outer += matrix_data[4][1] + matrix_data[4][2] + matrix_data[4][3];
          outer += matrix_data[0][1] + matrix_data[0][2] + matrix_data[0][3];
      
          for (int i = 1; i < 4; i++) {
              inner += matrix_data[i][1];
              inner += matrix_data[i][3];
	      time9 += matrix_time[i][1]*matrix_data[i][1];
	      time9 += matrix_time[i][3]*matrix_data[i][3];
          }
      
          inner += matrix_data[1][2] + matrix_data[3][2];
          time9 += matrix_data[1][2]*matrix_time[1][2] + matrix_data[3][2]*matrix_time[3][2];

	  sum9 = inner + matrix_data[2][2];
          time9 += matrix_data[2][2]*matrix_time[2][2];

	  //sum16_data = sum9 + outer;
	  sum16_data = outer;

	  Bool_t diffCut=false;	   
	  diff_8 = 0;
	  switch(sum){
	   case 1:
	     diff_8 = inner - sum_8*4;
	     diffCut = diff_8>-20 && diff_8<20;
	     break;
	   case 2:	  
	     inner = inner + matrix_data[2][2];
	     diff_8 = inner - sum_8*4;
	     diffCut = diff_8>-20 && diff_8<20;
	     break;	  
	   }
	   diff_16 = outer - sum_16*4;
	   h_diff_inn->Fill(diff_8);
	   h_diff_outer->Fill(diff_16);
	   time9 = matrix_time[2][2];
	   
	  Bool_t addCut=false;
	   
	  switch(typeCut){
	  case 1:
	    addCut=(b3_charge5<-2000);//charged particle
	    //addCut=(b3_charge7<-10);	    
	    postfix="_palettaSottileOn";
	    break;
	  case 2:
	    addCut=(b3_charge5>-2000);//neutral particle
	    //addCut=(b3_charge7>-10);	    
	    postfix="_palettaSottileOff";
	    break;
	  case 3:
	    addCut=(b2_charge7>2000);
	    postfix="_palettaSpessaOn";
	    break;
	  case 4:
	    addCut=(b2_charge7<2000);
	    postfix="_palettaSpessaOff";
	    break;
	    //RUN VECCHI
	  case 5:
	    addCut=(b3_charge5<-2000);
	    postfix="_palettaSottileOn";
	    break;
	  case 6:
	    addCut=(b3_charge5>-2000);
	    postfix="_palettaSottileOff";
	    break;
	  case 7:
	    addCut=(b2_charge7<-2000);
	    postfix="_palettaSpessaOn";
	    break;
	  case 8:
	    addCut=(b2_charge7>-2000);
	    postfix="_palettaSpessaOff";
	    break;
	    
	  default:
	    addCut=true;
	    postfix="";
	    break;    
	  }
	  //cout<<"cut "<<addCut<<" case "<<postfix<<endl;
	  //Bool_t csm_led_cut = b2_charge6>-2000 && b2_charge7<1000;
	  Bool_t cosmics_trigger = b2_charge7 < 50;
	  //per il run 591
	  //	  if( diff_8>-20 && diff_8<20 /*&& diff_16>-100 && diff_16<50*/)

	  double calib_factor = 6.1;


	  Bool_t top_signal = false;
	  Bool_t bottom_signal = false;
	  Bool_t left_signal = false;
	  Bool_t right_signal = false;

	  double cosmic_cut = 20;

	  for (int i = 0; i < 5; i++) {
	    if (matrix_data[4][i]/calib_factor > cosmic_cut)
	      top_signal = true;
	    if (matrix_data[0][i]/calib_factor > cosmic_cut)
	      bottom_signal = true;
	    if (matrix_data[i][0]/calib_factor > cosmic_cut)
	      left_signal = true;
	    if (matrix_data[i][4]/calib_factor > cosmic_cut)
	      right_signal = true;
	  }

	  Bool_t cosmic_track = top_signal || bottom_signal || left_signal || right_signal;

	  Bool_t hole_counter = b2_tmean6 == -999;
	  Bool_t beam_counter = b1_tmean7 > 0;
	  //cosmic_track = false;
	  if(addCut && diffCut && !cosmic_track && cosmics_trigger && diff_16>-40 && diff_16<40 && hole_counter && beam_counter){
	    
	    h_time_diff->Fill(b3_tmean6-b1_tmean7);

	    h_inner->Fill(inner/calib_factor);
	    h_outer->Fill(outer/calib_factor);
	    h_sum9->Fill(sum9/calib_factor);
	    h_sum16->Fill(sum16_data/calib_factor);
	    correlation_8->Fill(sum_8*4,inner);
	    correlation_16->Fill(sum_16*4,outer);
	    
	    for (int i = 0; i < 5; i++) {
              for (int j = 0; j < 5; j++) {
		h_cells[i][j]->Fill(matrix_data[i][j]/calib_factor);
		
		total += matrix_data[i][j];
              }
	    }
	    h_total->Fill(total/calib_factor);
	    if (total/calib_factor > 150)
	      cout << jentry << " " << total/calib_factor << endl;
	    h_energy_time->Fill(total/calib_factor, jentry);
	    h_energy_time2->Fill(total/calib_factor,b3_tmean6-b1_tmean7);
	  }//end if cosmic trg
      }      
   }
   
   TCanvas *c_en_time = new TCanvas("c_en_time","c_en_time");
   c_en_time->cd();
   h_energy_time->Draw();
   TCanvas *c_en_time2 = new TCanvas("c_en_time2","c_en_time2");
   c_en_time2->cd();
   h_energy_time2->Draw();

   TCanvas *c_cells = new TCanvas("c_cells","Cells spectra");
   c_cells->Divide(5,5);
   
   TCanvas *pad = new TCanvas();
   pad->Divide(2,2);
   
   TCanvas *c_single_crys = new TCanvas("c_single_crys","c_single_crys"); 
   c_single_crys->cd();
   h_cells[2][2]->Draw("e");
   // h_cells[2][2]->SetLineColor(kGreen+2);
   h_cells[2][1]->Draw("same");
   // h_cells[2][3]->SetLineColor(kRed);
   h_cells[2][3]->Draw("same");
   k=0;
   for (j=4;j>=0;j--) {
     for (i=0; i<5;i++) {
       k++;
       c_cells->cd(k);
       //h_cells[j][i]->Draw();
       // c_single_crys->cd();
       //if(i>0 && i<4) {
       if (matrix_data[j][i]/calib_factor > 7) {
	 h_cells[j][i]->Draw();
	 h_cells[j][i]->SetLineColor(kRed+1);
       } else {
	 h_cells[j][i]->Draw();
	 h_cells[j][i]->SetLineColor(kBlack);
       }
	 //}
	 //if(i>0 && i<4){
	 //h_cells[3][i]->SetLineColor(10+i*2);
	 //h_cells[3][i]->Draw("same");
	 //}
     }
   }
   
   
   //   TCanvas *c_central = new TCanvas("c_central","Central Crystal");
   pad->cd(1);
   h_cells[2][2]->Draw();

   //TCanvas *c_inner = new TCanvas("c_inner","Inner ring");
   pad->cd(2);   
   h_inner->Draw();
   
   //TCanvas *c_outer = new TCanvas("c_outer","Outer ring");
   pad->cd(3);   
   h_outer->Draw();
   
   //TCanvas *c_total = new TCanvas("c_total","Entire matrix");
   pad->cd(4);   
   h_total->Draw();

   TCanvas *c_diff = new TCanvas("c_diff","c_diff");
   c_diff->Divide(2,2);
   c_diff->cd(1);
   h_diff_inn->Draw();
   c_diff->cd(2);
   h_diff_outer->Draw();
   c_diff->cd(3);
   correlation_8->Draw();
   c_diff->cd(4);
   correlation_16->Draw();

   TCanvas *c_trg = new TCanvas("c_trg");
   c_trg->cd();
   h_trg_cosmic->Draw();

   TCanvas *c_sum9 = new TCanvas("sum9","sum9");
   c_sum9->cd();
   h_sum9->Draw();

   TCanvas *c_time = new TCanvas("time","time");
   c_time->cd();
   h_time_diff->Draw();
   TCanvas *c_sum16 = new TCanvas("sum16","sum16");
   c_sum16->cd();
   h_sum16->Draw();

   TCanvas *c_total = new TCanvas("total","total");
   c_total->cd();
   h_total->GetYaxis()->SetTitleOffset(1.08);
   // h_total->Scale(1./h_total->Integral());
   h_total->Draw();


   h_total->SaveAs(Form("hist_%i%s.root",run_number,postfix.Data()));
   h_sum9->SaveAs(Form("h_sum9_%i%s.root",run_number,postfix.Data()));
}
