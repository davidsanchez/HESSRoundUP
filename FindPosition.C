bool FindPosition(double l = 329.716938, double b = -30.225588	, const char* Prefix = "")
{
// double Eth

//  if(l > 180){
//    l = l-360;
//    std::cout << "automatic change : l -> " << l << std::endl;
//  }
  Stash::Coordinate targetpos(Stash::Lambda(l,Stash::Angle::Degrees),
			      Stash::Beta(b,Stash::Angle::Degrees),
			      Crash::GetRADecJ2000System());
  
  std::ostringstream fradialfile;
  fradialfile << Prefix<< "_roundup_RadialMaps.root";
   std::cout << fradialfile.str().c_str() << std::endl;
  TFile *fradial = TFile::Open(fradialfile.str().c_str());
  Display::SkyHistogram2D *hAcc = (Display::SkyHistogram2D *)fradial->Get("SkyAcceptance"); 
  Display::SkyHistogram2D *hTime = (Display::SkyHistogram2D *)fradial->Get("AveragedLiveTimeMap");
  Display::SkyHistogram2D *hZen = (Display::SkyHistogram2D *)fradial->Get("ZenithAngleMap");
  Display::SkyHistogram2D *hOff = (Display::SkyHistogram2D *)fradial->Get("OffAxisMap");
  Display::SkyHistogram2D *hEff = (Display::SkyHistogram2D *)fradial->Get("MuonEfficiencyMap");
  Display::SkyHistogram2D *hMinTh = (Display::SkyHistogram2D *)fradial->Get("MinSafeThresholdMap");
  Display::SkyHistogram2D *hMeanTh = (Display::SkyHistogram2D *)fradial->Get("AveragedSafeThresholdMap");
  
  int binx, biny;
  hTime->FindBinPosition(targetpos,binx,biny);
  if(hAcc->GetBinContent(binx,biny) == 0)
    {
      std::cout << "Null acceptance at this position" << std::endl;
      return 0;
    }
  
  std::cout << "Bins -> " << binx << " " << biny << std::endl;
  std::cout << "- Obs. Time = " << hTime->GetBinContent(binx,biny) << " hours" << std::endl;
  std::cout << "- (OffAxis, Zenith, Efficiency)  = (" << hOff->GetBinContent(binx,biny) << ", " 
	    << hZen->GetBinContent(binx,biny) << ", " << hEff->GetBinContent(binx,biny) << ")" << std::endl;
//  std::cout << "- Min Threshold = " << hMinTh->GetBinContent(binx,biny) << " TeV, Average Threshold = " << hMeanTh->GetBinContent(binx,biny) << " TeV" << std::endl;
    
  std::ostringstream fmapsfile;
  fmapsfile << Prefix<< "_roundup_RingBgMaps.root";

  TFile *fmaps = TFile::Open(fmapsfile.str().c_str());
  Display::SkyHistogram2D *hON = (Display::SkyHistogram2D *)fmaps->Get("ONCorrelated");
  Display::SkyHistogram2D *hOFF = (Display::SkyHistogram2D *)fmaps->Get("OFFCorrelated");
  Display::SkyHistogram2D *hAlpha = (Display::SkyHistogram2D *)fmaps->Get("AlphaCorrelated");
  Display::SkyHistogram2D *hExcess = (Display::SkyHistogram2D *)fmaps->Get("ExcessCorrelated");
  Display::SkyHistogram2D *hSignificance = (Display::SkyHistogram2D *)fmaps->Get("SignificanceMap");

  std::cout << std::endl;
  std::cout << "ON = " << hON->GetBinContent(binx,biny) << ", OFF = " << hOFF->GetBinContent(binx,biny) << ", Alpha = " << hAlpha->GetBinContent(binx,biny) << std::endl;
  std::cout << "Excess = " << hExcess->GetBinContent(binx,biny) << ", Significance = " << hSignificance->GetBinContent(binx,biny) << std::endl;
  


  return 1;

}
