#include <vector>			
#include <string>
#include <fstream>
#include <iomanip>
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"

#include "TString.h"
#include <sys/stat.h>
#include "logn.h"
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include "TF1.h"

using namespace std;
  
void makeroottopla_psi_2015_nocorr(TString outfile, TString FileNameAll,
				int events,
				TString datadir="/home/mu2e/mainz14/daq/",
				TString roottopledir="/home/mu2e/mainz14/roottople/")
{
  /*questo programma crea un vettore onda[i], che fornisce l'intero profilo della wave senza offset e quindi 'normalizzato';
    prevede l'utilizzo di 2 cicli for di lettura del file per determinare in maniera più accurata il piedistallo. Ogni piedistallo viene registrato per singlo evento. L'integrale di carica è calcolato mediante la seguente formula:: (dt/R)*(Somma di tutti i samples del singolo evento). Ogni samples è stato moltiplicato per 0,98 (mV) e per un fattore -3,3. */  

  TObjArray* FNameSplit = outfile.Tokenize("_");

  Int_t SubRunNum = ((TObjString*) FNameSplit->At(1))->String().Atoi();
  //std::cout<<"Sub run number "<<SubRunNum<<std::endl;
  
  //  gROOT.Reset();
  
  // ########## VARIABILI INTERFACCIA UTENTE ######

  //  char FileNameAll[100];
  //  string outfile;
  //  string outfile2;
  // ostringstream nameroot(ostringstream::out);
  // ostringstream nameroot2(ostringstream::out);
  TString nameroot, nameroot2;
  TString FileName = datadir;
  FileName += FileNameAll;
  //printf("txt file name: %s \n", FileName.Data());
  // ##############################################	

  int val=0;
  float sum[4][8];
  float sum_ped[4][8];
  float dt;
  float camp;	
  //  Int_t events=0;
  Int_t TimeSamples;
 
  //###########################################################

  //cout << "Sampling frequency (GHz):" << endl;
  //cin >> camp;
  camp=1;
  camp=0.25;
  dt=(1/camp);

  TimeSamples = 1024;

  //--------------------------------------------------------------

  // cout << "Enter file name" << endl;     
  // cin >> FileNameAll;
  // cout << "Reading from: " << FileNameAll << endl;

  //--------------------------------------------------------------

  //  cout << "Enter root file name:" << endl;
  //  cin >> outfile;
  //  nameroot << outfile << ".root";
  //cout<<"Output file: " << nameroot.str() << endl;
  //nameroot2 << outfile << "full.root";
  //cout <<"Output file full: " << nameroot2.str() << endl;
  
  // 02-11-2014 Giani changes
  
  nameroot = roottopledir;
  nameroot += outfile;
  nameroot += "psi.root";
  
  nameroot2 = roottopledir;
  nameroot2 += outfile;
  nameroot2 += "psi_full.root";
  
  //--------------------------------------------------------------

  //cout << "Enter number of events" << endl;
  //  cin >> events;
    
  //#############################################################

  //  TFile *OutputFile = new TFile(nameroot.str().c_str(),"recreate");
  TFile *OutputFile = new TFile(nameroot.Data(),"recreate");
  TTree *wavetree = new TTree("Waveform","Waveform");
  //02-11-2014 Giani changes
  wavetree->SetAutoSave(1000);//auto-save the TTree each 1 kb
  //------------------------------------------------------------//
  const int max_num=TimeSamples;//10000;  // Constant: cannot be modified
  Double_t time[max_num];
  Double_t ampl,piede[4][8][max_num],onda[4][8][max_num];
  Double_t bkg[4][8],charge[4][8],bkg_charge[4][8];
  Double_t tmean[4][8], tave[4][8], tmax[4][8];
  Double_t off[4][8];

  Double_t tempo=0.;	
  Int_t ntrig=0;
  Int_t cntb=0;
  Int_t nval[4][8];
  Int_t nev;
  int jmax[4][8];
  float ondamax[4][8];
  string none;
  char String1[20],String2[20];

  int Nboards=-1;
  int BoardNum[4]={0,0,0,0};
  int Nchan[4]={0,0,0,0};
  int ChanNum[4][8];

  //================================================================
  // Read DAQ configuration
  //================================================================

  ifstream WaveAll;
  //  WaveAll.open(FileNameAll);
  WaveAll.open(FileName.Data());
  
  if(!WaveAll.is_open()) 
    printf("file %s not open correctly\n",FileName.Data());
  
  char buff[500];
  if (SubRunNum>9) {
    //WaveAll >> none >> none >> none >> none >> none >> none >> none;
    //WaveAll.ignore(500,'\n');
    WaveAll.getline(buff,500);
  } else {
    //WaveAll >> none >> none >> none >> none;
    //WaveAll >> none >> none >> none >> none;
    //WaveAll.ignore(500,'\n');
    //WaveAll.ignore(500,'\n');
    WaveAll.getline(buff,500);
    WaveAll.getline(buff,500);
  }
  WaveAll >> none >> Nboards;

  printf("Nboards = %i \n", Nboards);

  for ( int Iboard=0; Iboard<Nboards; Iboard++ )
    WaveAll >> BoardNum[Iboard];

  for ( int Iboard=0; Iboard<Nboards; Iboard++ ){
    WaveAll >> none >> Nchan[Iboard];
    for( int Ichan=0; Ichan<Nchan[Iboard]; Ichan++ ){
      WaveAll >> ChanNum[Iboard][Ichan];
      //cout << "Board/Chan " << Iboard << " " << Ichan << " " << ChanNum[Iboard][Ichan] << endl;
    }
  }

  // Create ROOT tree

  wavetree->Branch("ntrig",&ntrig,"ntrig/I");
  wavetree->Branch("evnum",&nev,"evnum/I");

  /////////
  ///////// to-be-done: inseart board#/chan# info in root tree
  /////////

  for (int Iboard=0; Iboard<Nboards; Iboard++){
    for(int Ichan=0; Ichan<Nchan[Iboard]; Ichan++){
      sprintf(String1,"b%d_bkg%d",Iboard,Ichan);
      sprintf(String2,"b%d_bkg%d/D",Iboard,Ichan);
      wavetree->Branch(String1,&bkg[Iboard][Ichan],String2);
      sprintf(String1,"b%d_charge%d",Iboard,Ichan);
      sprintf(String2,"b%d_charge%d/D",Iboard,Ichan);
      wavetree->Branch(String1,&charge[Iboard][Ichan],String2);
      sprintf(String1,"b%d_ped%d",Iboard,Ichan);
      sprintf(String2,"b%d_ped%d/D",Iboard,Ichan);
      wavetree->Branch(String1,&bkg_charge[Iboard][Ichan],String2);
      sprintf(String1,"b%d_tmean%d",Iboard,Ichan);
      sprintf(String2,"b%d_tmean%d/D",Iboard,Ichan);
      wavetree->Branch(String1,&tmean[Iboard][Ichan],String2);
      sprintf(String1,"b%d_tave%d",Iboard,Ichan);
      sprintf(String2,"b%d_tave%d/D",Iboard,Ichan);
      wavetree->Branch(String1,&tave[Iboard][Ichan],String2);
      sprintf(String1,"b%d_tmax%d",Iboard,Ichan);
      sprintf(String2,"b%d_tmax%d/D",Iboard,Ichan);
      wavetree->Branch(String1,&tmax[Iboard][Ichan],String2);
    }
  }

  //  TFile *OutputFile2 = new TFile(nameroot2.str().c_str(),"recreate");
  TFile *OutputFile2 = new TFile(nameroot2.Data(),"recreate");
  TTree *wavefull = new TTree("Wavefull","Wavefull");
  //02-11-2014 Giani changes 
  wavefull->SetAutoSave(1000);//allows to auto-save the TTree each 1 kb 
  //------------------------------------------------------------//
  wavefull->Branch("ntrig",&ntrig,"ntrig/I");
  wavefull->Branch("evnum",&nev,"evnum/I");
  wavefull->Branch("num",&TimeSamples,"num/I");
  wavefull->Branch("time",time,"time[num]/D");

  for (int Iboard=0; Iboard<Nboards; Iboard++){
    for(int Ichan=0; Ichan<Nchan[Iboard]; Ichan++){
      sprintf(String1,"b%d_onda%d",Iboard,Ichan);
      sprintf(String2,"b%d_onda%d[num]/D",Iboard,Ichan);
      wavefull->Branch(String1,onda[Iboard][Ichan],String2);
      sprintf(String1,"b%d_ped%d",Iboard,Ichan);
      sprintf(String2,"b%d_ped%d/D",Iboard,Ichan);
      wavefull->Branch(String1,&bkg_charge[Iboard][Ichan],String2);
      sprintf(String1,"b%d_charge%d",Iboard,Ichan);
      sprintf(String2,"b%d_charge%d/D",Iboard,Ichan);
      wavefull->Branch(String1,&charge[Iboard][Ichan],String2);
      sprintf(String1,"b%d_tmean%d",Iboard,Ichan);
      sprintf(String2,"b%d_tmean%d/D",Iboard,Ichan);
      wavefull->Branch(String1,&tmean[Iboard][Ichan],String2);
      sprintf(String1,"b%d_tave%d",Iboard,Ichan);
      sprintf(String2,"b%d_tave%d/D",Iboard,Ichan);
      wavefull->Branch(String1,&tave[Iboard][Ichan],String2);
    }
  }
 
  //#################### CALCOLO DEL PIEDISTALLO ########################

  int ped_events=1;

  // REGION OF INTEGRATION  ...
  // Make pedwindow as large as chargewindow in units of sampling time
  // scintillators (PLASTIC+LYSO) typical 0-100 (i.e. 0-400 ns)
  // CsI  Typical 400-800 (1.6 musec)
  // Lyso Typical 400-500
  // float Sign[4] = {1., 1., -1., -1.};//first two boards single ended for matrix, 
  //float Sign[4] = {-1., -1., 1., 1.}; //modifica del 21-5-2014 Roberto prova
  float Sign[4] = {1., 1., 1., -1.}; // SM 26/8/2014.All positive .. 2 boards full Matrix + 2 boards 5chan/eacj for Matrix
  //third and fourth boards differentials for cosmic scintillators , 
  // finger scintillators and pin diode and single crystal
  
  float CalFact = 2000./4096.; // 2V over 4096 counts

  // window for charge and pedestal integration in ns
  float pedwin_min = 0;
  float pedwin_max = 500;
  //float charge_min = 1250; // Miscetti con ritardo di 500 ns sul trigger LED+CR 
  // float charge_max = 1650; // (prima era 1800-2200) then 1300,1700, today 2/Set 1250-1650.
  //  float charge_min=1650; // Misc 2/set after run46 with Coincidence 7*Beam
  //float charge_max=2050; 
//
  // float charge_min=1550; //Misc 3/set after run 68 changed coincidence width
  // float charge_max=2050; 
//------------------------------------------------------------------------------------------
//2014-11-18 BTF group changed the treshold according to resulkt shown on run300
  //float charge_min=1750;
  //float charge_max=2100; //2200

  // 2015 11-june MISC .. per run di calorimetro .. finestra spostata ..nuova coincidenza
  float charge_min= 3000;
  float charge_max= 3500;

  //float charge_fingers_min=1820;
  //float charge_fingers_max=1880; 


  //
  //float charge_min=1650; //Misc 3/set after run 68 changed coincidence width
  //float charge_max=1800; // for testing 250 ns gate.
  // charge_min = 1700; //ivano CsI
  // charge_max = 1880; //ivano CsI
   
  //--------------------------------------------------

  for (int k=0; k<ped_events; k++){
    
    if(WaveAll.eof()) break;

    // Read event header

    WaveAll >> none >> nev;
    
    if(none!="event"){
      cout << "ERROR: wrong data decoding during PED evaluation\n";
      exit(0);
    }
    if(WaveAll.eof()) goto EndOfFilePed;
    //cout << "nev " << nev << endl;

    // Initialize variables 

    for (int Iboard=0; Iboard<Nboards; Iboard++){
      for (int Ichan=0; Ichan<Nchan[Iboard]; Ichan++){
 	nval[Iboard][Ichan] = 0;
 	off[Iboard][Ichan]  = 0.;
      }
    }

    // Loop over wave samples

    for (int i=0; i<TimeSamples ;i++){

      tempo=dt*double(i+1);

      for (int Iboard=0; Iboard<Nboards; Iboard++){
 	for (int Ichan=0; Ichan<Nchan[Iboard]; Ichan++){
 	  WaveAll >> val;
 	  //cout << "Board " << Iboard << " Chan " << Ichan << " Value " << val << endl;

 	  if ( tempo>pedwin_min &&tempo < pedwin_max ){
 	    nval[Iboard][Ichan]++;
 	    off[Iboard][Ichan] = off[Iboard][Ichan]+val;
 	  }

 	}
      }
    }

    for (int Iboard=0; Iboard<Nboards; Iboard++){
      for (int Ichan=0; Ichan<Nchan[Iboard]; Ichan++){
 	piede[Iboard][Ichan][cntb]=((off[Iboard][Ichan])/(nval[Iboard][Ichan]));
 	cout << "Nev: " << cntb << " Board: " << Iboard 
	     << " Chan: " << Ichan << " Ped: " << piede[Iboard][Ichan][cntb] 
	     << endl;
      }
    }
    cntb=cntb+1;
  }
 EndOfFilePed:

  WaveAll.close();

  cout << "PEDs evaluated with " << ped_events << " events" << endl ;
  cout << " " << endl;

  //################## CALCOLO DELL'INTEGRALE DI CARICA #########################
  
  //  WaveAll.open(FileNameAll);
  WaveAll.open(FileName.Data());

  // Skip file header
  if (SubRunNum>9) {
    //WaveAll >> none >> none >> none >> none >> none >> none >> none;
    //WaveAll.ignore(500,'\n');
    WaveAll.getline(buff,500);
  } else {
    //WaveAll >> none >> none >> none >> none;
    //WaveAll >> none >> none >> none >> none;
    //WaveAll.ignore(500,'\n');
    //WaveAll.ignore(500,'\n');
    WaveAll.getline(buff,500);
    WaveAll.getline(buff,500);
  }
  WaveAll >> none >> none;

  for ( int Iboard=0; Iboard<Nboards; Iboard++ )
    WaveAll >> none;

  for ( int Iboard=0; Iboard<Nboards; Iboard++ ){
    WaveAll >> none >> none;
    for( int Ichan=0; Ichan<Nchan[Iboard]; Ichan++ ){
      WaveAll >> none;
      //cout << "Board/Chan " << Iboard << " " << Ichan << " " << none << endl;
    }
  }

  //float time_all = dt*TimeSamples;
  //float time_ped = time_all - 800.;     

  // Loop on events

  for (int k2=0; k2<events; k2++){

    if( k2%500==0) cout << "Event Number " << k2 << endl;
    //cout << "Event Number " << k2 << endl;
    ntrig++;

    // Read event header

    WaveAll >> none >> nev;
    if(none!="event"){
      cout << "ERROR: wrong data decoding\n";
      exit(0);
    }

    if(WaveAll.eof()) goto EndOfFile;

    //cout << "nev " << nev << endl;

    // Initialize variables 

    for (int Iboard=0; Iboard<Nboards; Iboard++){
      for (int Ichan=0; Ichan<Nchan[Iboard]; Ichan++){

 	tmax[Iboard][Ichan]  = 0.;
 	tave[Iboard][Ichan]  = 0.;
 	tmean[Iboard][Ichan] = 0.;

 	sum[Iboard][Ichan]     = 0.;
 	sum_ped[Iboard][Ichan] = 0.;

 	jmax[Iboard][Ichan]    = 0; 

	ondamax[Iboard][Ichan] = -10000.;

 	bkg[Iboard][Ichan] = piede[Iboard][Ichan][0]*CalFact; 

      }
    }


    // Loop over wave samples

    for (int j=0; j<TimeSamples; j++){

      time[j] = dt*double(j+1);

      for (int Iboard=0; Iboard<Nboards; Iboard++){
 	for (int Ichan=0; Ichan<Nchan[Iboard]; Ichan++){
 	  WaveAll >> val;
 	  //cout << "Board " << Iboard << " Chan " << Ichan << " Value " << val << endl;
	  /* 
	  float ymin, ymax;
	  ymin = 210+Sign[Iboard]*bkg[Iboard][Ichan];
	  ymax = 224+Sign[Iboard]*bkg[Iboard][Ichan]; */
	  ampl=(val)*CalFact*Sign[Iboard]; 
	  
	  /*
	    if( Iboard == 0 && ampl > 13000) ampl =ampl-14000.;
	    if( Iboard == 1 && ampl > 13000) ampl =ampl-14000.;
	    if( Iboard == 2 && ampl > 13000) ampl =ampl-14000.;
	    if( Iboard == 3 && ampl > 13000) ampl =ampl-14000.;
	    if( Iboard == 2 && ampl > ymin && ampl < ymax) ampl = Sign[Iboard]*bkg[Iboard][Ichan];
	    if( Iboard == 3 && ampl > ymin && ampl < ymax) ampl = Sign[Iboard]*bkg[Iboard][Ichan];
	  */
	  if (Iboard == 3 && Ichan == 7)
	    onda[Iboard][Ichan][j]=-(ampl-Sign[Iboard]*bkg[Iboard][Ichan]); 
	  else
	    onda[Iboard][Ichan][j]=(ampl-Sign[Iboard]*bkg[Iboard][Ichan]); 
	  

// onda is the normalized sample, with pedestal subtracted and normalized in mV
 	  //	  if(time[j]>160. && time[j]<240.){   //  giani & stefano: per integrare la carica lontano da rumore di eventi fuori tempo

//2014-11-18 BTF group changed the cut for calculating the finger charge in the correct way
/*	  if( (Iboard == 3) && (Ichan == 7 || Ichan == 6) ){
	    if(time[j]>charge_fingers_min && time[j]<charge_fingers_max)
	      sum[Iboard][Ichan] = sum[Iboard][Ichan] + onda[Iboard][Ichan][j];
	      }else { */

	    if(time[j]>charge_min && time[j]<charge_max)
	      sum[Iboard][Ichan] = sum[Iboard][Ichan] + onda[Iboard][Ichan][j];
	  


 	  if( time[j]>pedwin_min && time[j]<pedwin_max ){       //ivano: ped in carica 
 	    sum_ped[Iboard][Ichan] += onda[Iboard][Ichan][j];
 	  }

 	  // 1/10/2013 (SG): timing
	  // SM put to 1 the old cut that was introducing a bias
	  // 	  if( onda[Iboard][Ichan][j]>10 ){
	  if( onda[Iboard][Ichan][j]>1 ){
	    if( onda[Iboard][Ichan][j]>ondamax[Iboard][Ichan] ){
	      jmax[Iboard][Ichan] = j;
	      ondamax[Iboard][Ichan] = onda[Iboard][Ichan][j];
	      tmax[Iboard][Ichan] = time[j];
	    }
	  }
	  //cout << "onda" << Iboard << " " << Ichan << " " << onda[Iboard][Ichan][j] << endl;
	}
	
	
      } // End loop over wave samples
    }
    for (int Iboard=0; Iboard<Nboards; Iboard++){
      for (int Ichan=0; Ichan<Nchan[Iboard]; Ichan++){

 	charge[Iboard][Ichan] = (sum[Iboard][Ichan]*dt/50.);  // Calcolo dell'interale di carica. Dove 'sum' è la somma di tutti i samples dell'evento
 	bkg_charge[Iboard][Ichan] = (sum_ped[Iboard][Ichan]*dt/50.);      
 	charge[Iboard][Ichan] = charge[Iboard][Ichan]-bkg_charge[Iboard][Ichan];
	//cout << Iboard <<  " " << Ichan << " " << charge[Iboard][Ichan] << endl;

	// Custom time for b1c7 17/6/2015 at PSI
	bool beamtime = false;
	
	if (Iboard == 1 && Ichan == 7 && beamtime) {
	  tmean[Iboard][Ichan] = -99;
	  for (int j = 0; j < TimeSamples; j++) {
	    if (onda[Iboard][Ichan][j] < -70) {
	      tmean[Iboard][Ichan] = time[j];
	      j = TimeSamples;
	    }
	  }
	} else {

	  int minval = jmax[Iboard][Ichan]-10;
	  int maxval = jmax[Iboard][Ichan]+10;
	  
	  float denom = 0.;
	  

	  if( (jmax[Iboard][Ichan]-10)<0 ) minval = 0;
	  if( (jmax[Iboard][Ichan]+10)>TimeSamples ) maxval = TimeSamples;
	  for (int j=minval; j<maxval; j++){
	    if( onda[Iboard][Ichan][j]>(ondamax[Iboard][Ichan]/3.)){
	      tmean[Iboard][Ichan] += time[j]*onda[Iboard][Ichan][j];
	      denom += onda[Iboard][Ichan][j];
	    } 
	  }
	  if(denom>0. && ondamax[Iboard][Ichan] > 20 ){
	      tmean[Iboard][Ichan] = tmean[Iboard][Ichan]/denom;
	      
	  } else {
	    tmean[Iboard][Ichan] = -999;
	    

	  }

	}
      } 
    }
    

    wavetree->Fill();
    //if( k2%100==0) wavefull->Fill();
    //if( k2%10==0) wavefull->Fill();
    wavefull->Fill();

  }
 EndOfFile:
  
  // Closing files

  cout << "Closing files" << endl;

  WaveAll.close();

  OutputFile->cd();
  wavetree->Write();
  OutputFile->Close(); 
  OutputFile2->cd();
  wavefull->Write();
  OutputFile2->Close();

  cout << "Files closed " << endl;
}
