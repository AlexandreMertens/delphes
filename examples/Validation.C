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
#include <typeinfo>
#include "TLorentzVector.h"

//------------------------------------------------------------------------------

double ptrangemin = 10;
double ptrangemax = 1000;
static const int Nbins = 10;

struct resolPlot
{
    TH1 *cenResolHist;
    TH1 *fwdResolHist;
    int ptmin;
    int ptmax;
    TString obj;

    resolPlot();
    resolPlot(double ptdown, double ptup, TString object);
    void set(double ptdown, double ptup, TString object);
    print(){std::cout << ptmin << std::endl;}
};


resolPlot::resolPlot()
{
}

resolPlot::resolPlot(double ptdown, double ptup, TString object)
{
    this->set(ptdown,ptup,object);
}

void resolPlot::set(double ptdown, double ptup, TString object){
    ptmin = int(ptdown);
    ptmax = int(ptup);
    obj = object;

    cenResolHist = new TH1D(obj+"_delta_pt_"+Form("%d",ptmin)+"_"+Form("%d",ptmax)+"_cen", obj+"_delta_pt_"+Form("%d",ptmin)+"_"+Form("%d",ptmax)+"_cen", 100, -2.0, 2.0);
    fwdResolHist = new TH1D(obj+"_delta_pt_"+Form("%d",ptmin)+"_"+Form("%d",ptmax)+"_fwd", obj+"_delta_pt_"+Form("%d",ptmin)+"_"+Form("%d",ptmax)+"_fwd", 100, -2.0, 2.0);

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

    if(entry%10 == 0) cout << "Event number: "<< entry <<endl;

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
        pt  = jetMomentum.Pt();
        eta = TMath::Abs(jetMomentum.Eta());

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
                histos->at(bin).cenResolHist->Fill(bestGenJetMomentum.Pt()/jetMomentum.Pt());
            }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

void Validation(const char *inputFile, const char *outputFile)
{
  gSystem->Load("libDelphes");

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

  std::cout << "reading " << std::endl;

  for (int i = 0; i < Nbins; i++)
  {
      std::cout << plots[i].ptmin <<" " << plots[i].ptmax <<  std::endl;
  }
  
  GetJetsEres( &plots, treeReader);
  for (int i = 0; i < Nbins; i++)
  {

      std::cout << "          " << typeid(plots[i].cenResolHist).name() << '\n';
      std::cout << "entries : " << plots[i].cenResolHist->GetEntries() << std::endl;
  }

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
