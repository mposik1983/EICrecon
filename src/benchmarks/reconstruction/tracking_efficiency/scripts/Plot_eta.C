// Eta coverage of detectors
// Shyam Kumar:INFN Bari, shyam.kumar@ba.infn.it; shyam055119@gmail.com

#include "TGraphErrors.h"
#include "TF1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#define mpi 0.139  // 1.864 GeV/c^2

void Plot_eta()
{

////////////////////////////////////////////////////// style of the plot
   gStyle->SetPalette(1);
   gStyle->SetOptTitle(0);
   gStyle->SetTitleOffset(.85,"X");gStyle->SetTitleOffset(1.0,"Y");
   gStyle->SetTitleSize(.05,"X");gStyle->SetTitleSize(.05,"Y");
   gStyle->SetLabelSize(.04,"X");gStyle->SetLabelSize(.04,"Y");
   gStyle->SetHistLineWidth(2);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);


   TFile *f = TFile::Open("sim.edm4hep.root");
   TTree *tree = (TTree*)f->Get("events");

   TCanvas * c1 = new TCanvas("c1","coutput",1400,1000);
   c1->SetMargin(0.10, 0.05 ,0.1,0.05);
   c1->SetGridx();

   TH1D *h0 = new TH1D("h0","h0",200,-4.,4.);
   TH1D *h1 = new TH1D("h1","h1",200,-4.,4.); // Barr1
   TH1D *h2 = new TH1D("h2","h2",200,-4.,4.); // BMT1
   TH1D *h3 = new TH1D("h3","h3",200,-4.,4.);
   TH1D *h4 = new TH1D("h4","h4",200,-4.,4.);
   TH1D *h5 = new TH1D("h5","h5",200,-4.,4.);


   h0->GetXaxis()->SetTitle("#eta");
   h0->GetYaxis()->SetTitle("Entries (a.u.)");
   h0->GetXaxis()->CenterTitle();
   h0->GetYaxis()->CenterTitle();

   TString vtx_hits = "-TMath::Log(TMath::Tan((TMath::ATan2(sqrt(VertexBarrelHits.position.x*VertexBarrelHits.position.x+VertexBarrelHits.position.y*VertexBarrelHits.position.y),VertexBarrelHits.position.z))/2))>>h0";
   TString Barrel_hits = "-TMath::Log(TMath::Tan((TMath::ATan2(sqrt(SiBarrelHits.position.x*SiBarrelHits.position.x+SiBarrelHits.position.y*SiBarrelHits.position.y),SiBarrelHits.position.z))/2))>>h1";
   TString BMM_hits = "-TMath::Log(TMath::Tan((TMath::ATan2(sqrt(MPGDBarrelHits.position.x*MPGDBarrelHits.position.x+MPGDBarrelHits.position.y*MPGDBarrelHits.position.y),MPGDBarrelHits.position.z))/2))>>h2";
   TString ToF_hits = "-TMath::Log(TMath::Tan((TMath::ATan2(sqrt(TOFBarrelHits.position.x*TOFBarrelHits.position.x+TOFBarrelHits.position.y*TOFBarrelHits.position.y),TOFBarrelHits.position.z))/2))>>h3";
   TString ETracker_hits = "-TMath::Log(TMath::Tan((TMath::ATan2(sqrt(TrackerEndcapHits.position.x*TrackerEndcapHits.position.x+TrackerEndcapHits.position.y*TrackerEndcapHits.position.y),TrackerEndcapHits.position.z))/2))>>h4";
   TString ETOF_hits = "-TMath::Log(TMath::Tan((TMath::ATan2(sqrt(TOFEndcapHits.position.x*TOFEndcapHits.position.x+TOFEndcapHits.position.y*TOFEndcapHits.position.y),TOFEndcapHits.position.z))/2))>>h5";


   tree->Draw(vtx_hits.Data(),"VertexBarrelHits.position.y>0","goff"); // Vtx layers
   tree->Draw(Barrel_hits.Data(),"SiBarrelHits.position.y>0 ","goff"); // Barrel layers
   tree->Draw(BMM_hits.Data(),"MPGDBarrelHits.position.y>0 ","goff"); // Barrel layers
   tree->Draw(ToF_hits.Data(),"TOFBarrelHits.position.y>0","goff"); // BMT layers1
   tree->Draw(ETracker_hits.Data(),"TrackerEndcapHits.position.y>0","goff"); // BMT layers1
   tree->Draw(ETOF_hits.Data(),"TOFEndcapHits.position.y>0","goff"); // DIRC layers


   c1->cd();
   c1->SetLogy();
   h0->Scale(1./h0->Integral());
   h1->Scale(1./h1->Integral());
   h2->Scale(1./h2->Integral());
   h3->Scale(1./h3->Integral());
   h4->Scale(1./h4->Integral());
   h5->Scale(1./h5->Integral());
   h0->SetLineColor(kBlue);
   h1->SetLineColor(kMagenta);
   h2->SetLineColor(kRed);
   h3->SetLineColor(kBlack);
   h4->SetLineColor(kGreen);
   h5->SetLineColor(kCyan);


   h0->SetMaximum(0.1);
   h0->GetXaxis()->SetRangeUser(-4.0,4.0);
   h0->Draw("hist");
   h1->Draw("hist-same");
   h2->Draw("hist-same");
   h3->Draw("hist-same");
   h4->Draw("hist-same");
   h5->Draw("hist-same");

  TLegend *l1= new TLegend(0.70,0.75,0.90,0.93);
  l1->SetTextSize(0.03);
  l1->SetBorderSize(0);
  l1->AddEntry(h0,"VertexBarrelHits");
  l1->AddEntry(h1,"SiBarrelHits");
  l1->AddEntry(h2,"MPGDBarrelHits");
  l1->AddEntry(h3,"TOFBarrelHits");
  l1->AddEntry(h4,"TrackerEndcapHits");
  l1->AddEntry(h5,"TOFEndcapHits");
  c1->cd();
  l1->Draw();
  c1->SaveAs("eta_DD4HEP.png");

}
