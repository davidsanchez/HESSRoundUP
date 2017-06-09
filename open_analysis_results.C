#include <parisanalysis/AnalysisResults.hh>
#include <parisanalysis/Analysis2DResults.hh>
#include <parisanalysis/AcceptanceMap.hh>
#include <parisanalysis/RunStatistics.hh>
#include <parisanalysis/AnalysisConfig.hh>
#include <parisanalysis/PlotCollection.hh>
#include <parisanalysis/FieldOfViewAcceptance.hh>
#include <parisanalysis/RunList.hh>
#include <display/SkyHistogram2D.hh>
#include <sash/HESSArray.hh>
#include <sash/Telescope.hh>
#include <sash/TelescopeConfig.hh>
#include <sash/RunHeader.hh>
#include <sash/DataSet.hh>
#include <utilities/StringTools.hh>
#include <utilities/TextStyle.hh>
#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
#include <TPaletteAxis.h>
#include <TColor.h>
#include <TEllipse.h>
#include <TLatex.h>
#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <vector>
#include <string>

void open_analysis_results(TString filename, TString SaveName="Save"){
 
//TString filename = "Results_PKS_2155_304_ModelPlus_Combined_Std.root";
TString Sourcename="";
//TString SaveName="test";

 


//Size Canvas
double yC=800.;
double xC=(0.87/0.82)*yC; //To use with SetMargin(0.15,0.03,0.1,0.03);


//contour palette
const Int_t NRGBs = 5;
const Int_t NCont = 255;
Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
//gStyle->SetNumberContours(NCont);

	
double mytheta=0.01;	
double maxL=210;	
double minL=120;	
	
	
TFile *_file0 = new TFile(filename);
gROOT->cd();
Sash::DataSet * results = (Sash::DataSet *) _file0->Get("results");
results->GetEntry(0);

Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();




const ParisAnalysis::AnalysisConfig * fConf = hess->Get((ParisAnalysis::AnalysisConfig *)0);
double RA = fConf->GetTargetRa();
double Dec = fConf->GetTargetDec();
const ParisAnalysis::RunList * fRunList = hess->Get((ParisAnalysis::RunList *)0);
const ParisAnalysis::AnalysisResults * fResults = hess->Get("MultipleOff", (ParisAnalysis::AnalysisResults *) 0);
const ParisAnalysis::Analysis2DResults * fResults2D = hess->Get((ParisAnalysis::Analysis2DResults *) 0);
const ParisAnalysis::AcceptanceMap * fAccMap = hess->Get((ParisAnalysis::AcceptanceMap *) 0);

	std::cout<<"-------------"<<std::endl;
	std::cout<<"-------------"<<std::endl;
	fRunList->Print();
	std::cout<<"-------------"<<std::endl;
	std::cout<<"-------------"<<std::endl;



if(fConf){
	RA=-1.*fConf->GetTargetRa();
	Dec=fConf->GetTargetDec();
//	mytheta=fConf->GetTheta2Max();
	mytheta=fConf->GetTargetSize().GetDegrees();
	if(Sourcename=="")Sourcename=fConf->GetTargetName();

	std::cout<<""<<std::endl;
	std::cout<<""<<std::endl;
	fConf->Print();
	std::cout<<""<<std::endl;
	std::cout<<""<<std::endl;
}




double minxP=RA+fConf->GetMapExtensionX().GetDegrees()*0.8;
double maxxP=RA+fConf->GetMapExtensionX().GetDegrees()*0.9;
double minyP=Dec-fConf->GetMapExtensionY().GetDegrees()*0.5;
double maxyP=Dec+fConf->GetMapExtensionY().GetDegrees()*0.5;

double cosbeta = cos(Dec*2.*TMath::Pi()/360.); 



 //    First plots //
 //
 //..................................................
 //
 //
TCanvas * cTheta2 = fResults->DisplaySingle();
cTheta2->SaveAs(SaveName+"_roundup_Theta2.png");
//}


std::cout<<""<<std::endl;
std::cout<<""<<std::endl;
fResults->DisplaySigOverTime();
TCanvas * cSigvsNoff = fResults->DisplaySigVsNOFF();
cSigvsNoff->SaveAs(SaveName+"_roundup_SigvsNoff.png");

fResults->DisplaySigVsTheta2();

ParisAnalysis::Analysis2DResults *sigplot = hess->Handle("",(ParisAnalysis::Analysis2DResults *)0);
TCanvas* cansig = sigplot->DisplayRingModel();
cansig->Print(SaveName+"_roundup_SigMap.png");

////////////////////////////////////////////////////////////////////////////////////////////////////////
//                    Control map
//
//.....................................................................
//

TCanvas *cGammaAcceptance =  new TCanvas("CAcceptance","Acceptance Info",600,300);
cGammaAcceptance->Clear();
cGammaAcceptance->Divide(2,1,0.001,0.001);
cGammaAcceptance->cd(1);
Display::SkyHistogram2D * GammaAccMap=fAccMap->GetSkyGammaAcceptanceMap();

GammaAccMap->SetContour(256);
GammaAccMap->Draw("COLZ");

cGammaAcceptance->SaveAs(SaveName+"_roundup_GammaAcceptanceMap.png");

cGammaAcceptance->cd(2);
Display::SkyHistogram2D * BckAccMap=fAccMap->GetSkyBackgroundAcceptanceMap();
BckAccMap->SetContour(256);
BckAccMap->Draw("COLZ");

cGammaAcceptance->SaveAs(SaveName+"_roundup_AcceptanceMap.png");






}
