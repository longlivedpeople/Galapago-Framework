/// \file
/// \ingroup tutorial_fit
/// \notebook -js
/// Example for fitting signal/background.
/// This example can be executed with:
///
/// ~~~{.cpp}
/// root > .x FittingDemo.C  (using the CINT interpreter)
/// root > .x FittingDemo.C+ (using the native complier via ACLIC)
/// ~~~
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \author Rene Brun

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooTFnBinding.h"
#include "RooCMSShape.h"
#include "RooCBShape.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TFile.h"

#include <fstream>

using namespace RooFit;

double fitMass(TH1F* h, TString name, TString title) {

  RooRealVar x("x", "Mass (GeV)", 60, 120);
  RooDataHist histo("histo", "histo", x, Import(*h));

  RooPlot *frame = x.frame(Title(title));
  histo.plotOn(frame);

  RooRealVar acmsF("acmsF", "acmsF", 60, 50, 190);
  RooRealVar betaF("betaF", "betaF", 0.05, 0.01, 0.08);
  RooRealVar gammaF("gammaF", "gammaF", 0.1, -2, 2);
  RooRealVar peakF("peakF", "peakF", 91, 60, 120);
  RooCMSShape bkgFail("bkgFail", "bkgFail", x, acmsF, betaF, gammaF, peakF);

  RooRealVar alphaCB("alphaCB", "alpha", 1, 0, 10);
  RooRealVar nCB("nCB", "n", 1.5, 1, 5);
  RooRealVar sigmaCB("sigmaCB", "sigma", 3, 1, 5);
  RooRealVar muCB("mu", "mu", 90, 88, 92);
  RooCBShape signal("signal", "CBShape", x, muCB, sigmaCB, alphaCB, nCB);

  RooRealVar bkgfrac("bkgfrac", "fraction of background", 0.5, 0.0, 1.0);
  RooAddPdf model("model", "model",  RooArgList(bkgFail, signal), bkgfrac);

  model.fitTo(histo);
  model.plotOn(frame);
  model.plotOn(frame, Components(bkgFail), LineColor(kRed));

  //std::cout << bkgfrac.getVal() << std::endl;
  std::string outdir = "Fits/";
  std::string format = ".png";
  TCanvas c1("canvas", "canvas", 600, 600);
  c1.cd();
  frame->Draw();
  c1.SaveAs(outdir + name + format);

  return 1 - bkgfrac.getVal();
}

void FitData() {

  ofstream outfile("FitResults.txt");


  std::vector<std::string> Passed_;
  std::vector<std::string> Failed_;
  Passed_.push_back("TH1F_tnp_mass_passed_DATA");
  Failed_.push_back("TH1F_tnp_mass_failed_DATA");
  //Passed_.push_back("TH1F_tnp_mass_passed_pt2");
  //Failed_.push_back("TH1F_tnp_mass_failed_pt2");
  //Passed_.push_back("TH1F_tnp_mass_passed_pt2");
  //Failed_.push_back("TH1F_tnp_mass_failed_pt2");
  //Passed_.push_back("TH1F_tnp_mass_passed_pt2");
  //Failed_.push_back("TH1F_tnp_mass_failed_pt2");


  TFile file_ = TFile("plots_condor/SingleElectronRun2016G.root");
  for (unsigned int i = 0; i < Passed_.size(); i++) {
    TH1F *hfail = (TH1F*) file_.Get(TString(Failed_[i]));
    TH1F *hpass = (TH1F*) file_.Get(TString(Passed_[i]));

    double failFraction = fitMass(hfail, TString(Failed_[i]), "Failed probes fit (2016 data)");
    double passFraction = fitMass(hpass, TString(Passed_[i]), "Passing probes fit (2016 data)");

    double nfail = failFraction*hfail->Integral();
    double npass = passFraction*hpass->Integral();

    double eff = npass / (npass + nfail);

    outfile << ">>> From failed file: " << TString(Failed_[i]) << std::endl;
    outfile << ">>> From passed file: " << TString(Passed_[i]) << std::endl;
    outfile << "    failFraction: " << failFraction << std::endl;
    outfile << "    passFraction: " << passFraction << std::endl;
    outfile << "    total fail: " << hfail->Integral() << std::endl;
    outfile << "    total pass: " << hpass->Integral() << std::endl;
    outfile << "    fitted fail: " << nfail << std::endl;
    outfile << "    fitted pass: " << npass << std::endl;
    outfile << "    Efficiency: " << eff << std::endl;
  

    TFile file_mc = TFile("plots_condor/DYJetsToLL_M50_dispEle_2016.root");
    TH1F *hfailmc = (TH1F*) file_mc.Get(TString(Failed_[i]));
    TH1F *hpassmc = (TH1F*) file_mc.Get(TString(Passed_[i]));

    double eff_mc = hpassmc->Integral() / (hpassmc->Integral() + hfailmc->Integral());

    outfile << "    mc total fail: " << hfailmc->Integral() << std::endl;
    outfile << "    mc total pass: " << hpassmc->Integral() << std::endl;
    outfile << "    Efficiency mc: " << eff_mc  << std::endl;
    outfile << "    Correction factor in mc: " << eff/eff_mc << std::endl;
  }

}
