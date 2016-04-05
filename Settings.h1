#ifndef RAAINPUTSETTINGS
#define RAAINPUTSETTINGS

#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

class Settings {
  public:
 
  static const int ntrkBins = 37;
  double xtrkbins[ntrkBins+1] = { 0.5,0.6, 0.7 , 0.8 , 0.9 , 1.0 , 1.1 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.2 , 2.4 , 3.2 , 4.0 , 4.8 , 5.6 , 6.4 , 7.2 , 9.6 , 12.0, 14.4,19.2, 24.0, 28.8, 35.2, 41.6, 48.0, 60.8,73.6,86.4,103.6,120.8,140,165,250,400};
  
  static const int nTriggers = 4;
  double triggerBins[nTriggers+1] = {0,60,80,100,1200};
  double triggerOverlapBins[nTriggers] = {0,60,80,100};
  static const int HInTriggers = 5;
  double HItriggerBins[HInTriggers+1] = {0,60,80,100,120,1200};
  double HItriggerOverlapBins[HInTriggers] = {0,60,80,100,120};
  
  static const int njetBins = 240;
  static const int maxJetBin = 1200;
  
  static const int nTriggers_trk = 6;
  double triggerBins_trk[nTriggers_trk+1] = {0,20,26,36,47,55,500};
  double triggerOverlapBins_trk[nTriggers_trk] = {0,20,26,36,47,55};
  static const int HInTriggers_trk = 5;
  double HItriggerBins_trk[HInTriggers_trk+1] = {0,14,20,35,50,500};
  double HItriggerOverlapBins_trk[HInTriggers_trk] = {0,14,20,35,50};

  static const int nTrktriggerBins = 500;
  static const int maxTrktriggerBin = 500;
  
  static const int nCentBins = 32;
  int lowCentBin[nCentBins] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,0,6,10,2,6,10,14,16,0,0,14,0};
  int highCentBin[nCentBins] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,20,20,20,6,10,14,16,18,6,10,18,2};
  double nColl[nCentBins] = {1819,1433,1127,882,685.9,526.5,399.3,297.5,217.1,155.1,107.9,73.51,48.76,31.46,19.69,12.02,7.042,3.974,2.12,1.164,392.4,98.33,30.74,805.7,267.25,65.415,15.87,5.502,1079.1,754,10.686, 1626};
  double TAA[nCentBins] = {25.98,20.46,1127./70.0,882./70.0,685.9/70.0,526.5/70.0,399.3/70.0,297.5/70.0,217.1/70.0,155.1/70.0,107.9/70.0,73.51/70.0,48.76/70.0,31.46/70.0,19.69/70.0,12.02/70.0,7.042/70.0,3.974/70.0,2.12/70.0,1.164/70.0,392.4/70.0,98.33/70.0,30.74/70.0,11.51,3.819,0.9345,15.87/70.0,5.502/70.0,1079.1/70.0,754/70.0,0.1525, 1626/70.0};//Ncoll/70.0for TAA in most places, otherwise we are using offical numbers from cent group if no 70.0is there
  double TAAuncert[nCentBins] = {1.7,1.8,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,2.3,5.0,9.5,5,5,5,5,16.,5};//assume 5% uncert for 'unofficial' values for now (5.0 is official, 5 is not in this list)
  /*static const int nCentBins = 31;
  int lowCentBin[nCentBins] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,0,6,10,2,6,10,14,16,0,0,14};
  int highCentBin[nCentBins] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,20,20,20,6,10,14,16,18,6,10,18};
  double nColl[nCentBins] = {1819,1433,1127,882,685.9,526.5,399.3,297.5,217.1,155.1,107.9,73.51,48.76,31.46,19.69,12.02,7.042,3.974,2.12,1.164,392.4,98.33,30.74,805.7,267.25,65.415,15.87,5.502,1079.1,754,10.686};
  double TAA[nCentBins] = {25.98,20.46,1127./70.0,882./70.0,685.9/70.0,526.5/70.0,399.3/70.0,297.5/70.0,217.1/70.0,155.1/70.0,107.9/70.0,73.51/70.0,48.76/70.0,31.46/70.0,19.69/70.0,12.02/70.0,7.042/70.0,3.974/70.0,2.12/70.0,1.164/70.0,392.4/70.0,98.33/70.0,30.74/70.0,11.51,3.819,0.9345,15.87/70.0,5.502/70.0,1079.1/70.0,754/70.0,0.1525};//Ncoll/70.0for TAA in most places, otherwise we are using offical numbers from cent group if no 70.0is there
  double TAAuncert[nCentBins] = {1.7,1.8,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,2.3,5.0,9.5,5,5,5,5,16.};//assume 5% uncert for 'unofficial' values for now (5.0 is official, 5 is not in this list)*/

  TH2D *spec[nTriggers],               *HIspec[HInTriggers][nCentBins];
  TH1D *evtCount[nTriggers],           *HIevtCount[HInTriggers][nCentBins];
  TH2D *evtCount_JetVars[nTriggers],   *HIevtCount_JetVars[HInTriggers][nCentBins];
  TH1D *nVtxMB,                        *HInVtxMB[nCentBins];
  TH2D *spec_trk[nTriggers_trk],       *HIspec_trk[HInTriggers_trk][nCentBins];
  TH1D *evtCount_trk[nTriggers_trk],   *HIevtCount_trk[HInTriggers_trk][nCentBins];
  TH1D *nVtxMB_trk,                    *HInVtxMB_trk[nCentBins];

  TH1D * pp,                           *HI[nCentBins];
  TH1D * ppByTrigger[nTriggers],       *HIByTrigger[HInTriggers][nCentBins];
  TH1D * ppUsedByTrigger[nTriggers],   *HIUsedByTrigger[HInTriggers][nCentBins];
  TH1D * ppJets,                       *HIJets[nCentBins];
  TH1D * ppJetsByTrigger[nTriggers],   *HIJetsByTrigger[HInTriggers][nCentBins];
  TH1D * pp_perMBTrigger,              *HI_perMBTrigger[nCentBins];               

  TH1D * pp_trk,                              *HI_trk[nCentBins];                                  
  TH1D * ppByTrigger_trk[nTriggers_trk],      *HIByTrigger_trk[HInTriggers_trk][nCentBins];            
  TH1D * ppUsedByTrigger_trk[nTriggers_trk],  *HIUsedByTrigger_trk[HInTriggers_trk][nCentBins];        
  TH1D * ppMaxtrk,                            *HIMaxtrk[nCentBins];                              
  TH1D * ppMaxtrkByTrigger[nTriggers_trk],    *HIMaxtrkByTrigger[HInTriggers_trk][nCentBins];        
  TH1D * pp_perMBTrigger_trk,                 *HI_perMBTrigger_trk[nCentBins];                     

  TH1D * RAA[nCentBins];
  TH1D * RAA_trk[nCentBins];

  TH2D *h_scale, *h_scale_trk, *h_HIscale, *h_HIscale_trk;
  TH2D *h_normErr, *h_normErr_trk, *h_HInormErr, *h_HInormErr_trk;
  TH1D *h_normSyst, *h_normSyst_trk, *h_HInormSyst[nCentBins], *h_HInormSyst_trk[nCentBins];
  TH1D *RAA_totSyst[nCentBins], *PbPb_totSyst[nCentBins], *pp_totSyst; 
 
  Settings();
};
 
Settings::Settings()
{
  std::cout << "Getting setting.." << std::endl;
  return;
}

#endif
