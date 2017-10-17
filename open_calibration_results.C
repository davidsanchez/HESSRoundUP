#include <parisanalysis/AnalysisResults.hh>
#include <parisanalysis/Analysis2DResults.hh>
#include <parisanalysis/AcceptanceMap.hh>
#include <parisanalysis/RunStatistics.hh>
#include <parisanalysis/SkyMonitor.hh>
#include <parisanalysis/AnalysisConfig.hh>
#include <parisanalysis/PlotCollection.hh>
#include <parisanalysis/FieldOfViewAcceptance.hh>
#include <parisanalysis/TelescopeCogMap.hh>
#include <parisanalysis/RunList.hh>
#include <display/SkyHistogram2D.hh>
#include <sash/HESSArray.hh>
#include <sash/Telescope.hh>
#include <sash/TelescopeConfig.hh>
#include <sash/RunHeader.hh>
#include <sash/DataSet.hh>
#include <sash/Pointer.hh>
#include <sash/PointerSet.hh>
#include <sash/PointerSetIterator.hh>
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

void open_calibration_results(TString filename, TString SaveName="Save"){
 
//TString filename = "Results_PKS_2155_304_ModelPlus_Combined_Std.root";
//TString SaveName="test";
	


TFile *_file0 = new TFile(filename);
gROOT->cd();
Sash::DataSet * results = (Sash::DataSet *) _file0->Get("results");
results->GetEntry(0);

Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();
// const ParisAnalysis::TelescopeCogMap* Get(unsigned int telId, const char* name, ParisAnalysis::TelescopeCogMap*) const
// const ParisAnalysis::TelescopeCogMap *cog =  hess->Handle(5,ParisAnalysis::TelescopeCogMap*);

const ParisAnalysis::AnalysisConfig * fConf = hess->Get((ParisAnalysis::AnalysisConfig *)0);

const Sash::PointerSet<Sash::Telescope> &tels = fConf->GetTelsInAnalysis();

const ParisAnalysis::TelescopeCogMap *cog;

  Sash::NonConstPointer<Sash::Telescope>  first_tel = *tels.begin();
  for(Sash::PointerSet<Sash::Telescope>::iterator it = tels.begin();
      it != tels.end();++it)
    {
      UInt_t TelId = (*it)->GetId();
	cog =  hess->Handle<ParisAnalysis::TelescopeCogMap>(TelId);
	cog->Display();
  }


// const ParisAnalysis::SkyMonitor* NSB = hess->Get(ParisAnalysis::SkyMonitor*) ;

// const ParisAnalysis::SkyMonitor *NSBMonitor = hess->Get("NSB",ParisAnalysis::SkyMonitor*);
const ParisAnalysis::SkyMonitor *NSBMonitor = hess->Get<ParisAnalysis::SkyMonitor>("NSB");

// NSBMonitor->OverSample(0.1);
Display::SkyHistogram2D   *NSBMap =  NSBMonitor->GetMap();
// std::cout<<&NSBMap<<std::endl;
// NSBMonitor->Display();

TCanvas *NSBcanvas =  new TCanvas("NSB","NSB all tels",600,600);
NSBMap->SetContour(256);
NSBMap->OverSample(0.1);
NSBMap->Draw("COLZ");

NSBMap->LoadStarMap(14);
NSBMap->SetStarsVisible();
}
