
#include <iostream>
#include <iomanip>

#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <fstream>

//#include <TGraph.h>
//#include <TGraphErrors.h>
//#include <TH1F.h>
//#include <TH1I.h>
//#include <TStyle.h>
//#include <TAxis.h>
//#include <TLegend.h>
//#include <THStack.h>
//#include <TFile.h>
//#include <TLine.h>
//#include <TSystem.h>
//#include <TDirectory.h>

#include <sash/HESSArray.hh>
#include <sash/DataSet.hh>
#include <sash/Telescope.hh>
#include <sash/TelescopeConfig.hh>
#include <sash/Time.hh>
#include <sash/TimeDiff.hh>
#include <sash/RunHeader.hh>
#include <sash/Pointer.hh>
#include <sash/PointerSet.hh>
#include <sash/PointerSetIterator.hh>
#include <sashfile/HandlerC.hh>
#include <trigger/Block.hh>
#include <trigger/Header.hh>
#include <trigger/PackedEvent.hh>
#include <trigger/Monitor.hh>
#include <trigger/TelescopeMonitor.hh>
#include <trigger/TelescopeDelayMonitor.hh>
#include <display/MonitorList.hh>

#include <utilities/Statistics.hh>
#include <utilities/TextStyle.hh>
#include <utilities/SimpleStats.hh>
#include <utilities/Flux.hh>

#include <display/Histogram.hh>
#include <display/SkyHistogram2D.hh>
#include <display/SkyHistogram3D.hh>

#include <parisanalysis/RunList.hh>
#include <parisanalysis/RunStat.hh>
#include <parisanalysis/Analysis2DResults.hh>
#include <parisanalysis/AnalysisConfig.hh>
#include <parisanalysis/FieldOfViewAcceptance.hh>
#include <parisanalysis/AcceptanceMap.hh>
#include <parisanalysis/DataStorageRun.hh>

//#include <spectrum/SpectrumTableFinderMulti.hh>
//#include <spectrum/DataStorageLinearTable.hh>
//#include <spectrum/AcceptanceTableEfficiency.hh>
//#include <spectrum/ResolutionTableEfficiency.hh>
//#include <spectrum/SpectrumBase.hh>
//#include <spectrum/SpectrumPowerLaw.hh>

//#include <morphology/MorphologyTableFinderMulti.hh>
////#include <morphology/AngularResolutionTableEfficiency.hh>

//ROOT
#include <TROOT.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TRandom3.h>
#include <TNtuple.h>
#include <TNtupleD.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TRolke.h>

double DrawONOFFTest_atPosition(const char* Resfile, const char* PerRunMapsFile,  double l, double b, bool display, const char *OutfilePrefix = "out",bool UseConfigTarget=true)
{
  //Opening the results file
  TFile *fileResults = TFile::Open(Resfile);
  gROOT->cd();
  Sash::DataSet *results = (Sash::DataSet *)fileResults->Get("results");
  results->GetEntry(0);
  
  Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();
  ParisAnalysis::AnalysisConfig *Config =  hess->Handle("", (ParisAnalysis::AnalysisConfig *) 0);
  Config->LoadAllMembers();
  
  double UserLambda,UserBeta;

  UserLambda = l;
  UserBeta = b;

  if(UseConfigTarget){
    std::cout << "Using target position from Results File" << std::endl;
    UserLambda = Config->GetTargetRa();
    UserBeta = Config->GetTargetDec();
    std::cout << "UserLambda = " << UserLambda << ", UserBeta = " << UserBeta << std::endl;
  }


  Sash::DataSet *run_results = (Sash::DataSet *)fileResults->Get("run_results");
  gROOT->cd();

  Stash::Coordinate targetpos(Stash::Lambda(UserLambda,Stash::Angle::Degrees),
			      Stash::Beta(UserBeta,Stash::Angle::Degrees),
			      Config->GetSystem());
  
  TFile *file = TFile::Open(PerRunMapsFile);


  int count = 0;  
  std::map<int,double> vAccs, vEvents, vAllowed, vDate, vRun;

  for(int i = 0; i < run_results->GetEntries(); ++i)
    {
      run_results->GetEntry(i);
      const ParisAnalysis::DataStorageRun* dsrun = hess->Handle("GammaFOVStorage", (ParisAnalysis::DataStorageRun *) 0);
      
      std::ostringstream fname_Events;
      fname_Events << "Events_" << i;
      std::ostringstream fname_Acceptance;
      fname_Acceptance << "Acceptance_" << i;
      std::ostringstream fname_Excess;
      fname_Excess << "Excess_" << i;
      std::ostringstream fname_Significance;
      fname_Significance << "Significance_" << i;

      Display::SkyHistogram2D *hEvents = (Display::SkyHistogram2D *)file->Get(fname_Events.str().c_str());
      Display::SkyHistogram2D *hAcceptance = (Display::SkyHistogram2D *)file->Get(fname_Acceptance.str().c_str());
      Display::SkyHistogram2D *hExcess = (Display::SkyHistogram2D *)file->Get(fname_Excess.str().c_str());
      Display::SkyHistogram2D *hSignificance = (Display::SkyHistogram2D *)file->Get(fname_Significance.str().c_str());

      if(hEvents->Integral() == 0){
	continue;
      }

      int sky_bin = hEvents->FindBinPosition(targetpos);      

      double RunNumber = dsrun->GetRunNumber();
      double Events = hEvents->GetCircularIntegral(targetpos,Stash::Lambda(0.1,Stash::Angle::Degrees));
      double Exposure = hAcceptance->GetCircularIntegral(targetpos,Stash::Lambda(0.1,Stash::Angle::Degrees));
      double MJD_Start = dsrun->GetFirstEventTime().GetUTC().GetModifiedJulianDate();
      double MJD_End = dsrun->GetLastEventTime().GetUTC().GetModifiedJulianDate();
      
      //std::cout << "Run " << RunNumber << std::setprecision(8) << ":  (" << MJD_Start << ", " << MJD_End << ") -> " << Exposure << " " << " +/- " << Events << std::endl; 

      if(Exposure > 0){
	vDate[count] = (MJD_Start+MJD_End)*0.5;
	vAccs[count] = Exposure;
	vEvents[count] = Events;
	vAllowed[count] = 1;
	vRun[count] = RunNumber;
	++count;
      }

    }
  
  int iterations = 0;
  int changes = 1;
  std::map<int,double> vON, vOFF, vAlphaON, vAlphaOFF, vExcess, vSig;

  while(changes != 0){
    //std::cout << "Got In" << std::endl;
    int countchanges = 0;
    vON.clear();
    vOFF.clear();
    vAlphaON.clear();
    vAlphaOFF.clear(); 
    vExcess.clear(); 
    vSig.clear();
    
    std::map<int,double>::iterator itD = vDate.begin();
    std::map<int,double>::iterator itA = vAccs.begin();
    std::map<int,double>::iterator itE = vEvents.begin();
    std::map<int,double>::iterator itX = vAllowed.begin();
    
    for(; itD != vDate.end(); ++itD, ++itA, ++itE, ++itX){
    
      //      std::cout << itD->first << " " << itD->second << std::endl;
    
      double this_entry = itD->first;
      vON[itD->first]= itE->second;
      vAlphaON[itD->first]= itA->second;
      std::map<int,double>::iterator itDIn = vDate.begin();
      std::map<int,double>::iterator itAIn = vAccs.begin();
      std::map<int,double>::iterator itEIn = vEvents.begin();
      std::map<int,double>::iterator itXIn = vAllowed.begin();
      double sumAcc = 0;
      double sumEv = 0;
    
      for(; itDIn != vDate.end(); ++itDIn, ++itAIn, ++itEIn, ++itXIn){
	if(itDIn->first != this_entry && itXIn->second != 0){
	  sumAcc += itAIn->second;
	  sumEv += itEIn->second;
	}
      }
    
      vOFF[itD->first]= sumEv;
      vAlphaOFF[itD->first]= sumAcc;
    }
  
    //
    std::map<int,double>::iterator itON = vON.begin();
    std::map<int,double>::iterator itOFF = vOFF.begin();
    std::map<int,double>::iterator itAON = vAlphaON.begin();
    std::map<int,double>::iterator itAOFF = vAlphaOFF.begin();
    std::map<int,double>::iterator itX2 = vAllowed.begin();
    std::map<int,double>::iterator itRun = vRun.begin();
    std::map<int,double>::iterator itDate = vDate.begin();
  
    for(; itON != vON.end(); ++itON, ++itOFF, ++itAON, ++itAOFF, ++itX2){
      int this_entry = itON->first;
      double non = itON->second;
      double noff = itOFF->second;
      double alpha = itAON->second*1./itAOFF->second;
      double sig = Utilities::Statistics::Significance(non,noff,alpha);

      if(sig > 4.5 && itX2->second == 1){
	vAllowed[itX2->first-2] = 0;
	vAllowed[itX2->first-1] = 0;
	vAllowed[itX2->first] = 0;
	vAllowed[itX2->first+1] = 0;
	vAllowed[itX2->first+2] = 0;
		
	++countchanges;
      }
      //    std::cout << non-alpha*noff << " " << sig << " " << countchanges << std::endl;
    }
    //std::cout << iterations << " " << countchanges << std::endl;
    changes = countchanges;
    ++iterations;
  }

  std::cout << "Size " << vON.size() << std::endl;
  std::map<int,double>::iterator itON = vON.begin();
  std::map<int,double>::iterator itOFF = vOFF.begin();
  std::map<int,double>::iterator itAON = vAlphaON.begin();
  std::map<int,double>::iterator itAOFF = vAlphaOFF.begin();
  std::map<int,double>::iterator itX2 = vAllowed.begin();
  std::map<int,double>::iterator itRun = vRun.begin();
  std::map<int,double>::iterator itDate = vDate.begin();
  
  std::map<int,double>::iterator itONBegin = vON.begin();
  std::map<int,double>::iterator itONEnd = vON.end();
  TH1F *hON = new TH1F("ON","ON",vON.size(),itONBegin->first-0.5,itONEnd->first-0.5);
  TH1F *hAOFF = new TH1F("AOFF","AOFF",vON.size(),itONBegin->first-0.5,itONEnd->first-0.5);
  TH1F *hAllowed = new TH1F("Allowed","Allowed",vON.size(),itONBegin->first-0.5,itONEnd->first-0.5);
  TH1F *hExcess = new TH1F("Excess","Excess",vON.size(),itONBegin->first-0.5,itONEnd->first-0.5);
  TH1F *hSignif = new TH1F("Signif","Signif",vON.size(),itONBegin->first-0.5,itONEnd->first-0.5);
  TH1F *hSigDist = new TH1F("SignifD","SignifD",80,-10,10);
  TH1F *hSigDist_Allowed = new TH1F("SignifDA","SignifDA",80,-10,10);
  
  double maxsig = 0;
  
  for(; itON != vON.end(); ++itON, ++itOFF, ++itAON, ++itAOFF, ++itX2, ++itDate, ++itRun){
    int this_entry = itON->first;
    double non = itON->second;
    double noff = itOFF->second;
    double alpha = itAON->second*1./itAOFF->second;
    double sig = Utilities::Statistics::Significance(non,noff,alpha);
    double allowed = itX2->second;

    std::cout << "Bin " << this_entry << ", Run = " << vRun[this_entry] << " (MJD : " << vDate[this_entry] << "), Excess = " << non-alpha*noff << "\t -> Sig = " << sig  << std::endl;

    if(sig > maxsig)
      maxsig = sig;

    hON->Fill(itON->first, non);
    hAOFF->Fill(itON->first, alpha*noff);
    hExcess->Fill(itON->first, non-alpha*noff);
    hSignif->Fill(itON->first, sig);
    hAllowed->Fill(itON->first, allowed);
    hSigDist->Fill(sig);
    if(allowed == 1)
      hSigDist_Allowed->Fill(sig);

  }

  //std::cout << itONBegin->first << " " << itONEnd->first << std::endl;
  if(display){
  TCanvas *cONOFFTest = new TCanvas("cONOFF","cONOFF",1400,600);
  cONOFFTest->Divide(3,1);
  cONOFFTest->cd(1);
  hON->Draw();
  hAOFF->SetLineColor(2);
  hAOFF->Draw("same");
  hExcess->SetLineColor(4);
  hExcess->Draw("same");
  TLegend *leg = new TLegend(0.15,0.75,0.4,0.85);
  leg->AddEntry(hON,"ON","l");
  leg->AddEntry(hAOFF,"alpha * OFF","l");
  leg->AddEntry(hExcess,"Excess","l");
  leg->SetFillColor(0);
  leg->Draw();

  cONOFFTest->cd(2);
  hSignif->Draw();
  hAllowed->SetLineColor(7);
  hAllowed->Draw("same");
  TLegend *leg2 = new TLegend(0.15,0.75,0.4,0.85);
  leg2->AddEntry(hSignif,"Significance","l");
  leg2->AddEntry(hAllowed,"Exclusions","l");
  leg2->SetFillColor(0);
  leg2->Draw();

  cONOFFTest->cd(3)->SetLogy();
  hSigDist_Allowed->SetLineColor(2);
  hSigDist_Allowed->Draw();
  hSigDist->Draw("same");
  hSigDist_Allowed->Fit("gaus","Q","SAMES");      
  TF1 *f = hSigDist_Allowed->GetFunction("gaus");
  if(f)
    f->SetLineColor(kRed+2);

    std::ostringstream canvfile;
    canvfile << OutfilePrefix << "_OnOfftestplot.png";
    cONOFFTest->Print(canvfile.str().c_str());	

  }
  fileResults->Close();


  std::cout<<"Max significance found = "<<maxsig<<std::endl;
  return maxsig;
}
