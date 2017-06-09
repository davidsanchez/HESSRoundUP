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

#include <utilities/Statistics.hh>

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


double GetNRun(TString filename){
//{
//TString filename = "Results_PKS2155m304_Mono.root";
	
TFile *_file0 = new TFile(filename);
gROOT->cd();
Sash::DataSet * results = (Sash::DataSet *) _file0->Get("results");
results->GetEntry(0);

Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();

const ParisAnalysis::RunList * fRunList = hess->Get((ParisAnalysis::RunList *)0);
return fRunList->GetSize();
}

double GetSig(TString filename){
	
TFile *_file0 = new TFile(filename);
gROOT->cd();
Sash::DataSet * results = (Sash::DataSet *) _file0->Get("results");
results->GetEntry(0);

Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();

const ParisAnalysis::AnalysisResults * fResults = hess->Get("MultipleOff", (ParisAnalysis::AnalysisResults *) 0);

Float_t Norm  = (fResults->GetTotalLiveTimeOff() > 0 ? fResults->GetTotalLiveTimeOn()/fResults->GetTotalLiveTimeOff() : 1);
UInt_t on   = fResults->GetNTotalEventsOn()[0];
UInt_t off  = fResults->GetNTotalEventsOff()[0];
return   Utilities::Statistics::Significance(on,off,Norm);

}


void MakeSummary(TString filename){

std::ifstream infile(filename);
    std::cout<<"filename significance #run "<<std::endl;
std::string line;
while (std::getline(infile, line))
{
    std::istringstream iss(line);
    double sig = GetSig(line);
    double nrun =  GetNRun(line);
    std::cout<<line<<"\t"<<sig<<"\t"<<nrun<<"\t"<<std::endl;
}


}


