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
#include "TColor.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TLatex.h"

#include <fstream>

using namespace RooFit;
Int_t scolor = TColor::GetColor("#f18f01");
Int_t bcolor = TColor::GetColor("#048ba8");
Int_t tcolor = TColor::GetColor("#2e4057");

//TGaxis::SetMaxDigits(3);

double fitMass(TH1F* h, TString name, TString title) {


  TCanvas c1("canvas", "canvas", 600, 600);

  // Fit to Data with preliminar values extracted from MC fit
  RooRealVar x("x", "Mass (GeV)", 60, 120);
  RooDataHist histo("histo", "histo", x, Import(*h));

  //RooPlot *frame = x.frame(Title(title));
  RooPlot *frame = x.frame(Title(title + ";Mass (GeV); Pairs"));
  histo.plotOn(frame);

  RooRealVar acmsF("acmsF", "acmsF", 60, 50, 190);
  RooRealVar betaF("betaF", "betaF", 0.05, 0.01, 1.5);
  RooRealVar gammaF("gammaF", "gammaF", 0.1, -2, 2);
  RooRealVar peakF("peakF", "peakF", 60, 91, 120);
  RooCMSShape bkg("bkgFail", "bkgFail", x, acmsF, betaF, gammaF, peakF);
  
  //RooRealVar tau("tau","tau",-0.1,-5,+5) ;
  //RooExponential bkg("exp_mc","exp_mc",x,tau) ;

  RooRealVar alphaCB("alphaCB", "alpha", 0.5, 1, 5);
  RooRealVar nCB("nCB", "n", 1.5, 0, 20);
  RooRealVar sigmaCB("sigmaCB", "sigma", 3, 0, 10);
  RooRealVar muCB("mu", "mu", 90, 85, 95);
  RooCBShape signal("signal", "CBShape", x, muCB, sigmaCB, alphaCB, nCB);


  RooRealVar bkgfrac("bkgfrac", "fraction of background", 0.1, 0.0, 1.0);
  RooAddPdf model("model", "model",  RooArgList(bkg, signal), bkgfrac);

  //signal.fitTo(histo);
  //signal.plotOn(frame, LineColor(tcolor));

  model.fitTo(histo);
  model.plotOn(frame, LineColor(tcolor));
  model.plotOn(frame, Components(bkg), LineColor(bcolor));
  model.plotOn(frame, Components(signal), LineColor(scolor));

  TH1F *hsig = new TH1F("hsig", "", 1, 60, 120);
  hsig->SetLineColor(scolor);
  hsig->SetLineWidth(2);
  TH1F *hbkg = new TH1F("hbkg", "", 1, 60, 120);
  hbkg->SetLineColor(bcolor);
  hbkg->SetLineWidth(2);
  TH1F *htot = new TH1F("htot", "", 1, 60, 120);
  htot->SetLineColor(tcolor);
  htot->SetLineWidth(2);

  // Define legend
  TLegend *leg1 = new TLegend(0.55,0.73,0.86,0.87);
  leg1->SetTextSize(0.03);
  leg1->SetFillStyle(0);
  leg1->SetLineWidth(0);
  leg1->SetBorderSize(0);
  leg1->AddEntry("histo","Data", "P");
  //leg1->AddEntry("signal","Signal","L");
  //leg1->AddEntry("bkg","Background", "L");
  //leg1->AddEntry("model","Signal + background","L");
  leg1->AddEntry(htot,"Signal + Background","L");
  leg1->AddEntry(hsig,"Signal","L");
  leg1->AddEntry(hbkg,"Background", "L");

  // Define the label
  TLatex latex; 
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);
  latex.SetTextFont(42);
  latex.SetTextAlign(11);
  latex.SetTextSize(0.05);
  TLatex latexb;
  latexb.SetNDC();
  latexb.SetTextAngle(0);
  latexb.SetTextColor(kBlack);
  latexb.SetTextFont(42);
  latexb.SetTextAlign(11);
  latexb.SetTextSize(0.033);
  

  c1.cd();
  frame->Draw();
  //hsig.Draw("HIST,SAME");
  //hbkg.Draw("HIST,SAME");
  leg1->Draw();
  latex.DrawLatex(0.17, 0.84, "#bf{Private work}");
  latexb.DrawLatex(0.17, 0.8, "#it{CMS data}");
  std::string outdir = "Fits/";
  std::string format = ".pdf";
  c1.SaveAs(outdir + name + format);

  return 1 - bkgfrac.getVal();
}

void FitData() {

  ofstream outfile("FitResults.txt");

  TGaxis::SetMaxDigits(3);

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
  TFile file_mc = TFile("plots_condor/DYJetsToLL_M50_dispEle_2016.root");
  TString datal = "_data_2016";
  TString mcl = "_mc_2016";

  for (unsigned int i = 0; i < Passed_.size(); i++) {

    TH1F *hfail = (TH1F*) file_.Get(TString(Failed_[i]));
    TH1F *hpass = (TH1F*) file_.Get(TString(Passed_[i]));
    TH1F *hfailmc = (TH1F*) file_mc.Get(TString(Failed_[i]));
    TH1F *hpassmc = (TH1F*) file_mc.Get(TString(Passed_[i]));

    double failFraction = fitMass(hfail, TString(Failed_[i]) + datal, "Failed probes fit (Data)");
    double passFraction = fitMass(hpass, TString(Passed_[i]) + datal, "Passing probes fit (Data)");
    double failFraction_mc = fitMass(hfailmc, TString(Failed_[i]) + mcl, "Failed probes fit (MC)");
    double passFraction_mc = fitMass(hpassmc, TString(Passed_[i]) + mcl, "Passing probes fit (MC)");

    double nfail = failFraction*hfail->Integral();
    double npass = passFraction*hpass->Integral();
    double nfail_mc = failFraction_mc*hfailmc->Integral();
    double npass_mc = passFraction_mc*hpassmc->Integral();

    double eff = npass / (npass + nfail);
    double eff_mc = npass_mc / (npass_mc + nfail_mc);

    double eff_cac = hpass->Integral() / (hpass->Integral() + hfail->Integral());
    double eff_mc_cac = hpassmc->Integral() / (hpassmc->Integral() + hfailmc->Integral());

    std::cout << eff_cac/eff_mc_cac << std::endl;
    std::cout << eff/eff_mc << std::endl;


    outfile << ">>> From failed file: " << TString(Failed_[i]) << std::endl;
    outfile << ">>> From passed file: " << TString(Passed_[i]) << std::endl;
    outfile << "    Fail fraction (Data): " << failFraction << std::endl;
    outfile << "    Pass fraction (Data): " << passFraction << std::endl;
    outfile << "    Fail count (Data): " << hfail->Integral() << std::endl;
    outfile << "    Pass count (Data): " << hpass->Integral() << std::endl;
    outfile << "    N fail (Data): " << nfail << std::endl;
    outfile << "    N pass (Data): " << npass << std::endl;
    outfile << "    Efficiency fit (Data): " << eff << std::endl;
    outfile << "    Efficiency count (Data): " << eff_cac << std::endl;
    outfile << "   " << std::endl;
    outfile << "    Fail fraction (MC): " << failFraction_mc << std::endl;
    outfile << "    Pass fraction (MC): " << passFraction_mc << std::endl;
    outfile << "    Fail count (MC): " << hfailmc->Integral() << std::endl;
    outfile << "    Pass count (MC): " << hpassmc->Integral() << std::endl;
    outfile << "    N fail (MC): " << nfail_mc << std::endl;
    outfile << "    N pass (MC): " << npass_mc << std::endl;
    outfile << "    Efficiency fit (MC): " << eff_mc << std::endl;
    outfile << "    Efficiency count (MC): " << eff_mc_cac << std::endl;


    TCanvas c2("canvas_mc", "canvas_mc", 600, 600);
    c2.cd();
    hfailmc->Draw("Hist");
    c2.SaveAs("MC_failed.png");
    hpassmc->Draw("Hist");
    c2.SaveAs("MC_passed.png");

    //double eff_mc = hpassmc->Integral() / (hpassmc->Integral() + hfailmc->Integral());

  }

}
