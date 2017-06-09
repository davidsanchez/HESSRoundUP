#include <TCanvas.h>
#include <TRandom.h>
#include <TEllipse.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TPad.h>
#include <TEllipse.h>
#include <TText.h>
#include <TMarker.h>
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
#include <TFile.h>
#include <display/Histogram.hh>
#include <display/SkyHistogram2D.hh>

#include <iostream>
#include <iomanip>

#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <fstream>


void DrawMap(const char* Prefix, const char *OutfilePrefix = "out" )
{
                std::ostringstream MapFile;
                MapFile << Prefix<< "_roundup_RingBgMaps.root";      
      		TFile *rootfile_All = TFile::Open(MapFile.str().c_str());
      		Display::SkyHistogram2D *hSig = (Display::SkyHistogram2D *)rootfile_All->Get("SignificanceMap");
            	Display::SkyHistogram2D *hSigAllowed = (Display::SkyHistogram2D *)rootfile_All->Get("AllowedSignificanceMap");
		TCanvas *Canvas = new TCanvas("Canvas","Canvas",0,0,1000,500);
		Canvas->Divide(2,1);
		Canvas->cd(1);
                hSig->SetStats(000);
		hSig->Draw("COLZ");

		Canvas->cd(2);
               hSigAllowed->SetStats(000);
		hSigAllowed->Draw("COLZ");
    std::ostringstream canvfile;
    canvfile << OutfilePrefix << "_Maps.png";
    Canvas->Print(canvfile.str().c_str());	
}
