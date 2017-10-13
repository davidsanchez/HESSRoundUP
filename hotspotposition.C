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

int FindSpot(TH2F* SigMap,float SigMin,float RegionSize, float vpx[1000], float vpy[1000], float vxmax[1000]) {

	int nbspot=0;

	int MapSizeX=SigMap->GetNbinsX();
	int MapSizeY=SigMap->GetNbinsY();
	double BinSizeX = (SigMap->GetXaxis()->GetXmax() - SigMap->GetXaxis()->GetXmin())/SigMap->GetXaxis()->GetNbins();
	double BinSizeY = (SigMap->GetYaxis()->GetXmax() - SigMap->GetYaxis()->GetXmin())/SigMap->GetYaxis()->GetNbins();
	
	int RegionSizeX = ceil(RegionSize/BinSizeX);
	int RegionSizeY = ceil(RegionSize/BinSizeY);
				
	int posx,posy,posz;

        double Max=0;

	while((SigMap->GetMaximum()>SigMin) && (nbspot<1000))
	{
		SigMap->GetMaximumBin(posx,posy,posz);
		float COGx=0;
		float COGy=0;
		float weight=0;
		for (int x=posx-RegionSizeX; x<posx+RegionSizeY+1; x++)
		{
			if((x<MapSizeX)&&(x>0))
			{
				for (int y=posy-RegionSizeY; y<posy+RegionSizeY+1; y++)
				{
					if((y<MapSizeY)&&(y>0))
					{
						if (((x-posx)*(x-posx)*BinSizeX*BinSizeX + (y-posy)*(y-posy)*BinSizeY*BinSizeY) < RegionSize*RegionSize)
						{	
							COGx+=x*SigMap->GetBinContent(x,y);
							COGy+=y*SigMap->GetBinContent(x,y);														
							weight+=SigMap->GetBinContent(x,y);							
                                                        if (Max<SigMap->GetBinContent(x,y)){	Max =SigMap->GetBinContent(x,y);	}
						}
					}
				}
			}
		}

		posx=floor(COGx/weight);
		posy=floor(COGy/weight);
		// pos = numero du bin
		// vp = coordonnées RA/Dec en degrees
		vpx[nbspot]=SigMap->GetXaxis()->GetBinCenter(posx);
		vpy[nbspot]=SigMap->GetYaxis()->GetBinCenter(posy);
		vxmax[nbspot]=Max;
                Max = 0.;
		//std::cout << "HOT SPOT " << nbspot+1 << " " << SigMap->GetXaxis()->GetBinCenter(posx) <<" "<<SigMap->GetYaxis()->GetBinCenter(posy) << std::endl;	
		
		for (int x=posx-RegionSizeX; x<posx+RegionSizeY+1; x++)
		{
			if((x<MapSizeX) && (x>0))
			{
				for (int y=posy-RegionSizeY; y<posy+RegionSizeY+1; y++)
				{
					if((y<MapSizeY) && (y>0))
						if (((x-posx)*(x-posx)*BinSizeX*BinSizeX + (y-posy)*(y-posy)*BinSizeY*BinSizeY) < RegionSize*RegionSize)
						{
							SigMap->SetBinContent(x,y,0);
						}
				}
			}
		}
		nbspot++;

		if (nbspot==999)
		{
			cerr << "Warning Maximum sources !!!!!!" << endl;
		}
	}
	return nbspot;
}

void hotspotposition(const char* Prefix ,double MinSig = 4.5)
{
                std::ostringstream MapFile;
                MapFile << Prefix<< "_roundup_RingBgMaps.root";      
      		TFile *rootfile_All = TFile::Open(MapFile.str().c_str());
      		Display::SkyHistogram2D *hSig = (Display::SkyHistogram2D *)rootfile_All->Get("SignificanceMap");


		float vpx[1000];
		float vpy[1000];
		float vmax[1000];


		TH2F *SpotMap = new TH2F(*hSig);
		
		int nbspot = FindSpot(SpotMap,MinSig,0.5,vpx,vpy,vmax);

		
                std::cout  << MapFile.str().c_str()<< std::endl;
		for (int i=0;i<nbspot;i++) 
		{
			std::cout  << " " << i+1 << " " << -vpx[i] << " " << vpy[i] <<" Max Sig "<<vmax[i]<< std::endl;
		}
                if (nbspot==0) 	std::cout  << "No Hot spot found at "<<MinSig<<" sigma, Observe More! "<<std::endl;
	
}

void hotspotpositionONOFF(const char* Prefix ,double MinSig = 4.5)
{
                std::ostringstream MapFile;
                MapFile << Prefix<< "_roundup_ONOFFTestMaps.root";      
      		TFile *rootfile_All = TFile::Open(MapFile.str().c_str());
      		Display::SkyHistogram2D *hSig = (Display::SkyHistogram2D *)rootfile_All->->Get("ONOFFTest");


		float vpx[1000];
		float vpy[1000];
		float vmax[1000];


		TH2F *SpotMap = new TH2F(*hSig);
		
		int nbspot = FindSpot(SpotMap,MinSig,0.5,vpx,vpy,vmax);

		
                std::cout  << MapFile.str().c_str()<< std::endl;
		for (int i=0;i<nbspot;i++) 
		{
			std::cout  << " " << i+1 << " " << -vpx[i] << " " << vpy[i] <<" Max Sig "<<vmax[i]<< std::endl;
		}
                if (nbspot==0) 	std::cout  << "No Hot spot found at "<<MinSig<<" sigma, Observe More! "<<std::endl;
	
}