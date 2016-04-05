#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TBranch.h"
#include "TMath.h"
#include "TAttMarker.h"
#include "TAttLine.h"
#include "getTrkCorr.h"
#include "Settings.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void countTracks(std::vector<std::string> inputFiles, int jobNum, int isPP, bool isTest = false)
{
  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();

  float jetEtaSelection = 2;
 
  Settings s; 
  if(isPP){
    for(int i = 0; i<s.nTriggers; i++)
    {
      s.spec[i] = new TH2D(Form("spectrum_trigger%d",i),"",s.njetBins,0,s.maxJetBin,s.ntrkBins,s.xtrkbins);
      s.evtCount[i] = new TH1D(Form("evtCount%d",i),";max jet p_{T};N",s.njetBins,0,s.maxJetBin);
      s.evtCount[i]->SetMarkerColor(i);
      s.evtCount_JetVars[i] = new TH2D(Form("evtCount_JetVars%d",i),"max jet #eta;max jet p_{T};N",10,-2,2,16,40,120);
    }
    s.nVtxMB = new TH1D("nVtxMB","nVtx;N Events",12,0,12);
    for(int i = 0; i<s.nTriggers_trk; i++)
    {
      s.spec_trk[i] = new TH2D(Form("spectrum_trigger%d_trk",i),"",s.nTrktriggerBins,0,s.maxTrktriggerBin,s.ntrkBins,s.xtrkbins);
      s.evtCount_trk[i] = new TH1D(Form("evtCount%d_trk",i),";max jet p_{T};N",s.nTrktriggerBins,0,s.maxTrktriggerBin);
      s.evtCount_trk[i]->SetMarkerColor(i);
    }
    s.nVtxMB_trk = new TH1D("nVtxMB_trk","nVtx;N Events",12,0,12);
  }  //end of PbPb loop
//******************************************************************************************************************************
//******************************************************************************************************************************
//******************************************************************************************************************************
  int nTrk;
  int nVtx;
  int nTrkTimesnVtx;
  bool highPurity[50000];
  float trkPt[50000];
  float trkPtError[50000];
  float trkEta[50000];
  float trkPhi[50000];
  float trkMVA[50000];
  float trkDxy1[50000];
  float trkDxyError1[50000];
  float trkDz1[50000];
  float trkDzError1[50000];
  float trkDzOverDzError[500000];
  float trkDxyOverDxyError[500000];
  float pfEcal[50000];
  float pfHcal[50000];
  float trkChi2[50000];
  float zVtx[20];
  unsigned char trkNHit[50000];
  unsigned char trkNlayer[50000];
  unsigned char trkNdof[50000];
  unsigned char trkAlgo[50000];
  unsigned char trkOriginalAlgo[50000];

  unsigned int run=0;
  unsigned int lumi=0;

  int pVtx;
  int pBeamScrape;
  //int NoiseFilter; 
  int pclusterCompatibilityFilter; 
  int pprimaryVertexFilter;  
  int phfCoincFilter3;
  int hiBin = 1;
  float hiHF = 0;

  int nref;
  float jtpt[200];
  float jteta[200];
  float jtphi[200];
  float rawpt[200];
  float chargedSum[200];
  float ecalSum[200];
  float hcalSum[200];

  int MB[20]={0};
  int j40=0;
  int j60=0;
  int j80=0;
  int t18=0;
  int t24=0;
  int t34=0;
  int t45=0;
  int t53=0;

  TFile * inputFile;
  TTree * trkCh;
  TTree * jetCh;
  TTree * evtCh;
  TTree * hltCh;
  TTree * hiCh;

  //for documenting which PD a file comes out of to avoid overlaps between PDs
  //0 is MB, 1 is jet40/60, 2 is jet80
  int PDindx[5000];
  for(unsigned int i = 0; i<inputFiles.size(); i++)
  {
    if(isPP){
      if((inputFiles.at(i).find("MinimumBias") != std::string::npos) || (inputFiles.at(i).find("MinBias") != std::string::npos)) PDindx[i]=0;
      else if(inputFiles.at(i).find("Jet40") != std::string::npos) PDindx[i]=1;
      else if(inputFiles.at(i).find("Jet80") != std::string::npos) PDindx[i]=2;
      else PDindx[i]=-1;
    }
  }

  for(int nFile = 0; nFile<inputFiles.size(); nFile++){
    inputFile = TFile::Open(inputFiles.at(nFile).c_str(),"read");
    if(isPP) trkCh = (TTree*)inputFile->Get("ppTrack/trackTree");
    
    trkCh->SetBranchAddress("nTrk",&nTrk);
    trkCh->SetBranchAddress("nVtx",&nVtx);
    trkCh->SetBranchAddress("trkPt",&trkPt);
    trkCh->SetBranchAddress("trkEta",&trkEta);
    trkCh->SetBranchAddress("trkPhi",&trkPhi);
    trkCh->SetBranchAddress("highPurity",&highPurity);
    trkCh->SetBranchAddress("trkMVA",&trkMVA);
    trkCh->SetBranchAddress("trkNHit",&trkNHit);
    trkCh->SetBranchAddress("trkPtError",&trkPtError);
    trkCh->SetBranchAddress("pfHcal",&pfHcal);
    trkCh->SetBranchAddress("pfEcal",&pfEcal);
    trkCh->SetBranchAddress("trkDxy1",&trkDxy1);
    trkCh->SetBranchAddress("trkDxyError1",&trkDxyError1);
    trkCh->SetBranchAddress("trkDz1",&trkDz1);
    trkCh->SetBranchAddress("trkDzError1",&trkDzError1);
    trkCh->SetBranchAddress("trkChi2",&trkChi2);
    trkCh->SetBranchAddress("trkNlayer",&trkNlayer);
    trkCh->SetBranchAddress("trkNdof",&trkNdof);
    trkCh->SetBranchAddress("trkAlgo",&trkAlgo);
    trkCh->SetBranchAddress("trkOriginalAlgo",&trkOriginalAlgo);
    trkCh->SetBranchAddress("zVtx",&zVtx);
    if(isPP){
      trkCh->SetBranchAddress("nTrkTimesnVtx",&nTrkTimesnVtx);
      trkCh->SetBranchAddress("trkDzOverDzError",&trkDzOverDzError);
      trkCh->SetBranchAddress("trkDxyOverDxyError",&trkDxyOverDxyError); 
    }
  
    if(isPP) jetCh = (TTree*)inputFile->Get("ak3PFJetAnalyzer/t");
    jetCh->SetBranchAddress("nref",&nref);
    jetCh->SetBranchAddress("jtpt",&jtpt);
    jetCh->SetBranchAddress("jteta",&jteta);  
    jetCh->SetBranchAddress("jtphi",&jtphi);  
    jetCh->SetBranchAddress("rawpt",&rawpt);
    jetCh->SetBranchAddress("chargedSum",&chargedSum);  
    jetCh->SetBranchAddress("ecalSum",&ecalSum);
    jetCh->SetBranchAddress("hcalSum",&hcalSum);  
    trkCh->AddFriend(jetCh);
  
    evtCh = (TTree*)inputFile->Get("skimanalysis/HltTree");
    if(isPP){
      evtCh->SetBranchAddress("pprimaryvertexFilter",&pVtx);
      evtCh->SetBranchAddress("pBeamScrapingFilter",&pBeamScrape);
      //evtCh->SetBranchAddress("pHBHENoiseFilterResultProducer",&NoiseFilter);
    }
   
    hltCh = (TTree*)inputFile->Get("hltanalysis/HltTree");
    if(isPP){
      for(int i = 0; i<20; i++) hltCh->SetBranchAddress(Form("HLT_L1MinimumBiasHF1OR_part%d_v1",i),&(MB[i]));
      hltCh->SetBranchAddress("HLT_PAJet40_NoJetID_v1",&j40);
      hltCh->SetBranchAddress("HLT_PAJet60_NoJetID_v1",&j60);
      hltCh->SetBranchAddress("HLT_PAJet80_NoJetID_v1",&j80);
    }
    trkCh->AddFriend(hltCh);
  //***********************************************************************************
  //***********************************************************************
    std::cout << "starting event loop" << std::endl;
    std::cout << trkCh->GetEntries() << std::endl;

    int numberEvt = trkCh->GetEntries();
    if(numberEvt>10000) numberEvt = 10000;
    for(int i = 0; i<numberEvt; i++)
    {
      //if(i%1000==0) std::cout << i<<"/"<<trkCh->GetEntries()<<" "<<std::endl;
      evtCh->GetEntry(i);
      //if(!NoiseFilter) continue;
      if(isPP && (!pVtx || !pBeamScrape)) continue;

      trkCh->GetEntry(i);
      
      bool hasGoodVtx = 0; 
      for(int vtx = 0; vtx<nVtx; vtx++){
        if(TMath::Abs(zVtx[vtx])<15) hasGoodVtx = true;
      }
      if(hasGoodVtx==false) continue;

      bool MinBias = 0;
      for(int j = 0; j<21; j++) MinBias = MinBias || ((isPP)?(bool)MB[j]:(bool)HIMB[j]);
      if(isPP && !MinBias && !j40 && !j60 && !j80 && !t18 && !t24 && !t34 && !t45 && !t53) continue;
  
      //**************************************************
      //for trigger combination with jet triggers
      float maxJtPt = 0;
      float maxJtEta = -99;
      for(int j=0; j<nref; j++)
      {
        if(isPP && (chargedSum[j]/rawpt[j]<0.01 || TMath::Abs(jteta[j])>jetEtaSelection)) continue;
        if(jtpt[j]>maxJtPt){
          maxJtPt = jtpt[j];
          maxJtEta = jteta[j];
        }
      }//end maxJt
  
 
      int PD = PDindx[nFile];
      if(MinBias && PD==0)
      {
        if(isPP){
          s.evtCount[0]->Fill(maxJtPt); 
          s.evtCount_JetVars[0]->Fill(maxJtEta,maxJtPt);
          s.nVtxMB->Fill(nVtx);
          s.evtCount_trk[0]->Fill(maxTrackPt); 
          s.nVtxMB_trk->Fill(nVtx);
        }
      }
      if(j40 && PD==1){s.evtCount_JetVars[1]->Fill(maxJtEta,maxJtPt); s.evtCount[1]->Fill(maxJtPt);}  
      if(j60 && PD==1){s.evtCount_JetVars[2]->Fill(maxJtEta,maxJtPt); s.evtCount[2]->Fill(maxJtPt);}  
      if(j80 && PD==2){s.evtCount_JetVars[3]->Fill(maxJtEta,maxJtPt); s.evtCount[3]->Fill(maxJtPt);} 
      
      if(PD!=3)
      {
        for(int j = 0; j<nTrk; j++)
        { 
          if(trkPt[j]<0.5 || trkPt[j]>=400) continue;
          if(TMath::Abs(trkEta[j])>1) continue;
          if(highPurity[j]!=1) continue;
          if( trkPtError[j]/trkPt[j]>0.1) continue;       
          if(TMath::Abs(trkDz1[j]/trkDzError1[j])>3 || TMath::Abs(trkDxy1[j]/trkDxyError1[j])>3) continue;

         // if((maxJtPt>50 && trkPt[j]>maxJtPt) || (maxJtPt<=50 && trkPt[j]>50)) continue;

          //dividing by pt at bin center instead of track by track pt (just a convention)
          float binCenter;
          if(isPP) binCenter = s.spec[0]->GetYaxis()->GetBinCenter(s.spec[0]->GetYaxis()->FindBin(trkPt[j]));
          if(isPP){
            if(MinBias && PD==0) s.spec[0]->Fill(maxJtPt,trkPt[j],1.0/binCenter); 
            if(j40 && PD==1)     s.spec[1]->Fill(maxJtPt,trkPt[j],1.0/binCenter); 
            if(j60 && PD==1)     s.spec[2]->Fill(maxJtPt,trkPt[j],1.0/binCenter); 
            if(j80 && PD==2)     s.spec[3]->Fill(maxJtPt,trkPt[j],1.0/binCenter); 
          }
        } //end trk loop
      }//end if statement  
    }//end event loop
    inputFile->Close();
  }//end file loop

  //for pp
  TFile * outF;
  outF = TFile::Open(Form("276pp_output_%d.root",jobNum),"recreate");
  outF->cd();
  if(isPP){
    for(int i = 0; i<s.nTriggers; i++)
    {
      s.spec[i]->Write();
      s.evtCount[i]->Write();
      s.evtCount_JetVars[i]->Write();
    }
    s.nVtxMB->Write();
    s.nVtxMB_trk->Write();
  }
}



//*************************************************************************************
//*************************************************************************************
//*************************************************************************************
int main(int argc, const char* argv[])
{
  if(argc != 5)
  {
    std::cout << "Usage: countTracks <fileList>  <job>" << std::endl;
    return 1;
  }  

  std::string fList = argv[1];
  int job = std::atoi(argv[2]);
  int totalJobs = std::atoi(argv[3]);
  int isPP = std::atoi(argv[4]);
  std::string buffer;
  std::vector<std::string> listOfFiles;
  std::ifstream inFile(fList.data());

  if(!inFile.is_open())
  {
    std::cout << "Error opening jet file. Exiting." <<std::endl;
    return 1;
  }
  else
  {
    int line = 0;
    while(true)
    {
      inFile >> buffer;
      if(inFile.eof()) break;
      if(line%totalJobs==job) listOfFiles.push_back(buffer);
      line++;
    }
  }
   
  countTracks(listOfFiles,job,isPP);
  return 0; 
}
