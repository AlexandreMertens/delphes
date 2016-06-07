/*
This macro shows how to compute jet energy scale.
root -l examples/Example4.C'("delphes_output.root", "plots.root")'

The output ROOT file contains the pT(MC)/pT(Reco) distributions for
various pT(Reco) and |eta| bins. The peak value of such distribution is
interpreted as the jet energy correction to be applied for that
given pT(Reco), |eta| bin.

This can be done by modifying the "ScaleFormula" input parameter to
the JetEnergyScale module in the delphes_card_XXX.tcl

e.g  a smooth function:

  set ScaleFormula { sqrt(3.0 - 0.1*(abs(eta)))^2 / pt + 1.0) }

or a binned function:

  set ScaleFormula {(abs(eta) > 0.0 && abs(eta) <= 2.5) * (pt > 20.0 && pt <= 50.0)  * (1.10) +
                    (abs(eta) > 0.0 && abs(eta) <= 2.5) * (pt > 50.0 && pt <= 100.0) * (1.05) +
                    (abs(eta) > 0.0 && abs(eta) <= 2.5) * (pt > 100.0)               * (1.00) +
                    (abs(eta) > 2.5 && abs(eta) <= 5.0) * (pt > 20.0 && pt <= 50.0)  * (1.10) +
                    (abs(eta) > 2.5 && abs(eta) <= 5.0) * (pt > 50.0 && pt <= 100.0) * (1.05) +
                    (abs(eta) > 2.5 && abs(eta) <= 5.0) * (pt > 100.0)               * (1.00)}

Be aware that a binned jet energy scale can produce "steps" in the corrected
jet pt distribution ...
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif

#include "TSystem.h"
#include <TH1.h>
#include "TString.h"
#include "vector"
#include <TMath.h>
#include <iostream>
#include "TGraph.h"
#include "TGraphErrors.h"
#include <typeinfo>
#include "TLorentzVector.h"

//------------------------------------------------------------------------------

double ptrangemin = 5;
double ptrangemax = 1000;
static const int Nbins = 100;

struct resolPlot
{
    TH1 *cenResolHist;
    TH1 *fwdResolHist;
    double ptmin;
    double ptmax;
    TString obj;

    resolPlot();
    resolPlot(double ptdown, double ptup, TString object);
    void set(double ptdown, double ptup, TString object);
    void print(){std::cout << ptmin << std::endl;}
};


resolPlot::resolPlot()
{
}

resolPlot::resolPlot(double ptdown, double ptup, TString object)
{
    this->set(ptdown,ptup,object);
}

void resolPlot::set(double ptdown, double ptup, TString object){
    ptmin = ptdown;
    ptmax = ptup;
    obj = object;

    cenResolHist = new TH1D(obj+"_delta_pt_"+Form("%4.2f",ptmin)+"_"+Form("%4.2f",ptmax)+"_cen", obj+"_delta_pt_"+Form("%4.2f",ptmin)+"_"+Form("%4.2f",ptmax)+"_cen", 100, -2.0, 2.0);
    fwdResolHist = new TH1D(obj+"_delta_pt_"+Form("%4.2f",ptmin)+"_"+Form("%4.2f",ptmax)+"_fwd", obj+"_delta_pt_"+Form("%4.2f",ptmin)+"_"+Form("%4.2f",ptmax)+"_fwd", 100, -2.0, 2.0);

}

void HistogramsCollection(std::vector<resolPlot> *histos, double ptmin, double ptmax, TString obj)
{
    double width;
    double ptdown;
    double ptup;
    resolPlot ptemp;

    for (int i = 0; i < Nbins; i++)
    {
        width = (ptmax - ptmin) / Nbins;
        ptdown = TMath::Power(10,ptmin + i * width );
        ptup = TMath::Power(10,ptmin + (i+1) * width );
        ptemp.set(ptdown, ptup, obj);
        histos->push_back(ptemp);
    }
}

//------------------------------------------------------------------------------

class ExRootResult;
class ExRootTreeReader;

//------------------------------------------------------------------------------

void BinLogX(TH1*h)
{

   TAxis *axis = h->GetXaxis();
   int bins = axis->GetNbins();

   Axis_t from = axis->GetXmin();
   Axis_t to = axis->GetXmax();
   Axis_t width = (to - from) / bins;
   Axis_t *new_bins = new Axis_t[bins + 1];

   for (int i = 0; i <= bins; i++) {
     new_bins[i] = TMath::Power(10, from + i * width);

   }
   axis->Set(bins, new_bins);
   delete new_bins;
} 


//------------------------------------------------------------------------------

template<typename T> 
std::pair <TH1D,TH1D> GetEff(TString branch, T* recoObj, int pdgID, ExRootTreeReader *treeReader)
{
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchReco = treeReader->UseBranch(branch);

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  //Jet *jet, *genjet;
  GenParticle *particle;
  //TObject *reco;

  TLorentzVector recoMomentum, genMomentum, bestRecoMomentum;

  Float_t deltaR;
  Float_t pt, eta;
  Long64_t entry;

  Int_t i, j;

  TH1D histGenPt  = TH1D(branch+" gen spectra Pt",branch+" gen spectra Pt", Nbins, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax));
  TH1D histRecoPt = TH1D(branch+" reco spectra Pt",branch+" reco spectra Pt", Nbins, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax));
  TH1D histGenEta  = TH1D(branch+" gen spectra Eta",branch+" gen spectra Eta", 12, -3, 3);
  TH1D histRecoEta = TH1D(branch+" reco spectra Eta",branch+" reco spectra Eta", 12, -3, 3);


  BinLogX(&histGenPt);
  BinLogX(&histRecoPt);

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    if(entry%10000 == 0) cout << "Event number: "<< entry <<endl;

    // Loop over all reconstructed jets in event
    for(i = 0; i < branchParticle->GetEntriesFast(); ++i)
    {

      particle = (GenParticle*) branchParticle->At(i);
      genMomentum = particle->P4();

      deltaR = 999;
   
      if (particle->PID == pdgID && genMomentum.Pt() > 1)
      {
    
        // Loop over all reco object in event
        for(j = 0; j < branchReco->GetEntriesFast(); ++j)
        {
          recoObj = (T*)branchReco->At(j);
          recoMomentum = recoObj->P4();
          // this is simply to avoid warnings from initial state particle
          // having infite rapidity ...
     	  //if(Momentum.Px() == 0 && genMomentum.Py() == 0) continue;
     
          // take the closest parton candidate
          if(genMomentum.DeltaR(recoMomentum) < deltaR)
          {
            deltaR = genMomentum.DeltaR(recoMomentum);
            bestRecoMomentum = recoMomentum;
          }
        }

        pt  = genMomentum.Pt();
        eta = genMomentum.Eta();

        histGenPt.Fill(pt);
        histGenEta.Fill(eta);

        if(deltaR < 0.5)
        {
          histRecoPt.Fill(pt);
          histRecoEta.Fill(eta); 
        }
      }
    }
  }

  std::pair <TH1D,TH1D> histos;

  histRecoPt.Divide(&histGenPt);
  histRecoEta.Divide(&histGenEta);

  histos.first = histRecoPt;
  histos.second = histRecoEta;

  return histos;
}

template<typename T>
void GetEres(std::vector<resolPlot> *histos, TString branch, int pdgID, ExRootTreeReader *treeReader)
{
  T* recoObj;
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchReco = treeReader->UseBranch(branch);

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;
  //TObject *object;

  TLorentzVector recoMomentum, genMomentum, bestGenMomentum;

  Float_t deltaR;
  Float_t pt, eta;
  Long64_t entry;

  Int_t i, j, bin;

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    if(entry%100000 == 0) cout << "Event number: "<< entry <<endl;

    // Loop over all reconstructed jets in event
    for(i = 0; i < branchReco->GetEntriesFast(); ++i)
    {

      recoObj = (T*) branchReco->At(i);
      recoMomentum = recoObj->P4();

      deltaR = 999;

     // Loop over all hard partons in event
     for(j = 0; j < branchParticle->GetEntriesFast(); ++j)
     {
        particle = (GenParticle*) branchParticle->At(j);

        if (particle->PID == pdgID)
        {
            genMomentum = particle->P4();

            // this is simply to avoid warnings from initial state particle
            // having infite rapidity ...
    	    if(genMomentum.Px() == 0 && genMomentum.Py() == 0) continue;
    
            // take the closest parton candidate
            if(genMomentum.DeltaR(recoMomentum) < deltaR)
            {
               deltaR = genMomentum.DeltaR(recoMomentum);
               bestGenMomentum = genMomentum;
            }
        }
      }

      if(deltaR < 0.3)
      {
        pt  = bestGenMomentum.Pt();
        eta = TMath::Abs(bestGenMomentum.Eta());

        for (bin = 0; bin < Nbins; bin++)
        {
            //std::cout << "ptmin : " << std::endl;
            //std::cout << "histos " << typeid(histos).name() << '\n';
            //std::cout << "*histos " << typeid(*histos).name() << '\n';
            //histos[bin]->();
            //std::cout << histos->at(bin).ptmin << std::endl;
            if(pt > histos->at(bin).ptmin && pt < histos->at(bin).ptmax && eta > 0.0 && eta < 2.5) 
            {
                //std::cout << histos->at(bin).ptmin << std::endl;
                if (eta < 0.5) {histos->at(bin).cenResolHist->Fill((bestGenMomentum.E()-recoMomentum.E())/bestGenMomentum.E());}
                else if (eta < 2.5) {histos->at(bin).fwdResolHist->Fill((bestGenMomentum.E()-recoMomentum.E())/bestGenMomentum.E());}
            }
        }
      }
    }
  }
}

void GetJetsEres(std::vector<resolPlot> *histos, ExRootTreeReader *treeReader)
{
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  Jet *jet, *genjet;
  GenParticle *particle;
  //TObject *object;

  TLorentzVector jetMomentum, genJetMomentum, bestGenJetMomentum;

  Float_t deltaR;
  Float_t pt, eta;
  Long64_t entry;

  Int_t i, j, bin;

  // Loop over all events
  for(entry = 0; entry < 100; ++entry) //allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    if(entry%10000 == 0) cout << "Event number: "<< entry <<endl;

    // Loop over all reconstructed jets in event
    for(i = 0; i < branchJet->GetEntriesFast(); ++i)
    {

      jet = (Jet*) branchJet->At(i);
      jetMomentum = jet->P4();

      deltaR = 999;

     // Loop over all hard partons in event
     for(j = 0; j < branchGenJet->GetEntriesFast(); ++j)
     {
        genjet = (Jet*) branchGenJet->At(j);

        genJetMomentum = genjet->P4();

	// this is simply to avoid warnings from initial state particle
        // having infite rapidity ...
	if(genJetMomentum.Px() == 0 && genJetMomentum.Py() == 0) continue;

        // take the closest parton candidate
        if(genJetMomentum.DeltaR(jetMomentum) < deltaR)
        {
           deltaR = genJetMomentum.DeltaR(jetMomentum);
           bestGenJetMomentum = genJetMomentum;
        }
      }

      if(deltaR < 0.3)
      {
        pt  = genJetMomentum.Pt();
        eta = TMath::Abs(genJetMomentum.Eta());

        for (bin = 0; bin < Nbins; bin++)
        {
            //std::cout << "ptmin : " << std::endl;
            //std::cout << "histos " << typeid(histos).name() << '\n';
            //std::cout << "*histos " << typeid(*histos).name() << '\n';
            //histos[bin]->();
            //std::cout << histos->at(bin).ptmin << std::endl;
            if(pt > histos->at(bin).ptmin && pt < histos->at(bin).ptmax && eta > 0.0 && eta < 2.5) 
            {
                //std::cout << histos->at(bin).ptmin << std::endl;
                histos->at(bin).cenResolHist->Fill((bestGenJetMomentum.E()-jetMomentum.E())/bestGenJetMomentum.E());
            }
        }
      }
    }
  }
}

TGraphErrors EresGraph(std::vector<resolPlot> *histos, bool central)
{
    Int_t bin;
    TGraphErrors gr = TGraphErrors(Nbins);
    for (bin = 0; bin < Nbins; bin++)
    {
        if (central == true) 
        {
            gr.SetPoint(bin,(histos->at(bin).ptmin+histos->at(bin).ptmax)/2.0, histos->at(bin).cenResolHist->GetRMS());
            gr.SetPointError(bin,0, histos->at(bin).cenResolHist->GetRMSError());
        }
        else 
        {
            gr.SetPoint(bin,(histos->at(bin).ptmin+histos->at(bin).ptmax)/2.0, histos->at(bin).fwdResolHist->GetRMS());
            gr.SetPointError(bin,0, histos->at(bin).fwdResolHist->GetRMSError());
        }
    }
    return gr;
}

//------------------------------------------------------------------------------

void Validation(const char *inputFile, const char *outputFile)
{
  //gSystem->Load("libDelphes");

  std::cout << "input file : " << inputFile << " " << " , output file : " << outputFile << std::endl;

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  //ExRootResult *result = new ExRootResult();

  //TestPlots *plots = new TestPlots;
  //BookHistograms(result, plots);

  //jethistos = 
  //std::vector<resolPlot> plots = HistogramsCollection ( TMath::Log10(ptmin), TMath::Log10(ptmax), Nbins, "jets"); 

  std::vector<resolPlot> plots;
  HistogramsCollection(&plots, 1, 3, "jets");

  GetJetsEres( &plots, treeReader);

  TGraphErrors gr = EresGraph(&plots, true);

  TString elecs = "Electron";
  int elID = 11;

  Electron *elec;

  std::pair <TH1D,TH1D> histos = GetEff(elecs, elec, elID, treeReader);

  TString tracks = "Track";
  Track *track;

  std::pair <TH1D,TH1D> histos_track = GetEff(tracks, track, elID, treeReader);


  // Electron Energy Resolution
  std::vector<resolPlot> plots_el;
  HistogramsCollection(&plots_el, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "electrons");
  GetEres<Electron>( &plots_el, "Electron", 11, treeReader);
  TGraphErrors gr_el = EresGraph(&plots_el, true);
  gr_el.SetName("Electron");

  // Electron Track Energy Resolution
  std::vector<resolPlot> plots_eltrack;
  HistogramsCollection(&plots_eltrack, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "electronsTracks");
  GetEres<Track>( &plots_eltrack, "Track", 11, treeReader);
  TGraphErrors gr_eltrack = EresGraph(&plots_eltrack, true);
  gr_eltrack.SetName("ElectronTracks");

  // Electron Tower Energy Resolution
  std::vector<resolPlot> plots_eltower;
  HistogramsCollection(&plots_eltower, TMath::Log10(ptrangemin), TMath::Log10(ptrangemax), "electronsTower");
  GetEres<Tower>( &plots_eltower, "Tower", 11, treeReader);
  TGraphErrors gr_eltower = EresGraph(&plots_eltower, true);
  gr_eltower.SetName("ElectronTower");


  TFile *fout = new TFile(outputFile,"recreate");
  gr.Write();
  histos.first.Write();
  histos.second.Write();
  histos_track.first.Write();
  histos_track.second.Write();
  gr_el.Write();
  gr_eltrack.Write();
  gr_eltower.Write();
  fout->Write();

  
  /*
  AnalyseEvents(treeReader, plots);

  result->Write(outputFile);
  */

  cout << "** Exiting..." << endl;

  //delete plots;
  //delete result;
  delete treeReader;
  delete chain;
}

//------------------------------------------------------------------------------
