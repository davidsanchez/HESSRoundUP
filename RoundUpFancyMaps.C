
#include "RoundUpFancyMaps.hh"

//HESS
#include <sash/HESSArray.hh>
#include <sash/Folder.hh>
#include <sash/EnvelopeEntry.hh>
#include <sash/DataSet.hh>
#include <sash/NominalPointing.hh>
#include <sash/Time.hh>
#include <sash/TimeDiff.hh>

#include <crash/RotatedSystem.hh>

#include <utilities/Statistics.hh>
#include <utilities/Parameter.hh>
#include <utilities/TextStyle.hh>
#include <utilities/TFftConv.hh>
#include <utilities/TextStyle.hh>
#include <utilities/Flux.hh>
#include <mathutils/Fourier.hh>

#include <display/Histogram.hh>
#include <display/SkyHistogram2D.hh>

#include <parisanalysis/RunList.hh>
#include <parisanalysis/RunStat.hh>
#include <parisanalysis/Analysis2DResults.hh>
#include <parisanalysis/AnalysisConfig.hh>
#include <parisanalysis/FieldOfViewAcceptance.hh>
#include <parisanalysis/AcceptanceMap.hh>
#include <parisanalysis/DataStorageRun.hh>
#include <parisanalysis/RadialAcceptance.hh>

#include <spectrum/SpectrumTableFinderMulti.hh>
#include <spectrum/DataStorageLinearTable.hh>
#include <spectrum/AcceptanceTableEfficiency.hh>
#include <spectrum/ResolutionTableEfficiency.hh>
#include <spectrum/SpectrumBase.hh>
#include <spectrum/SpectrumPowerLaw.hh>

#include <morphology/MorphologyTableFinderMulti.hh>
#include <morphology/AngularResolutionTableEfficiency.hh>

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

//STL
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>

SurveySuite::MapMaker::MapMaker():
  fEventsMap(0),
  fAcceptanceMap(0),
  fExclusionMap(0),
  fExpectedCountsMap(0),
  fExtendedExpectedCountsMap(0),
  fExcessMap(0),
  fUncorrelatedExcessMap(0),
  fSignificanceMap(0),
  fAlphaMap(0),
  fOffMap(0),
  fOffAxisMap(0),
  fZenithAngleMap(0),
  fMuonEfficiencyMap(0),
  fMinSafeThresholdMap(0),
  fAveragedSafeThresholdMap(0),
  fAveragedLiveTimeMap(0)
{
}
 
/** Destructor */
SurveySuite::MapMaker::~MapMaker()
{  
  delete fEventsMap;
  delete fAcceptanceMap;
  delete fExclusionMap;
  delete fExpectedCountsMap;
  delete fExtendedExpectedCountsMap;
  delete fExcessMap;
  delete fUncorrelatedExcessMap;
  delete fSignificanceMap;
  delete fAlphaMap;
  delete fOffMap;
  delete fOffAxisMap;
  delete fZenithAngleMap;
  delete fMuonEfficiencyMap;
  delete fMinSafeThresholdMap;
  delete fAveragedSafeThresholdMap;
  delete fAveragedLiveTimeMap;
}

void SurveySuite::MapMaker::Clear()
{  
  fEventsAndAcceptanceFromFile = 0; 
  fEventsAndAcceptanceFromRadial = 0;
  fEventsAndAcceptanceFromRadialFile = 0;
  fExclusionFromFits = 0;
  fExclusionFromRegionFile = 0;
  fUseConfigTarget = 0; 
  fUseConfigMapParams = 0; 
  fUserLambda = 0; 
  fUserBeta = 0; 
  fBinSize = 0; 
  fExtX = 0; 
  fExtY = 0; 
  fPsiCut = 0; 
  fOSRadius = 0; 
  fAdaptFFT = 0; 
  fAdaptCut_Alpha = 0; 
  fConstantArea = 0; 
  fConstantThickness = 0; 
  fConstantInnerRadius = 0; 
  fSmoothRings = 0; 
  fRingMinimalRadius = 0; 
  fInnerRingMax = 0; 
  fOuterRingMax = 0; 
  fRingStep = 0; 
  fRingParam_AreaOverPi = 0; 
  fRingParam_ExcFracMax = 0; 
  fRingParam_Thickness = 0; 
  fRingParam_AlphaMax = 0; 
  fStandardRingRadius = 0; 
  fStandardRingThickness = 0; 
  fAverageOff = 0;
  fCorrectZenith = 0;
  fEMin = 0;
  fEMax = 0;
  fFitAcceptance = 0;
  fSelectAllRunsContributingToTheMap = 0;
  fProduceFluxProducts = 0;
  fSpectralIndex = 0;
  fPointLikeFluxMaps = 0;
  fProduceSurfaceBrightnessMap = 0;
  fExposureMapsFromFits = 0;
  fSaveResults = 0;
  fVerbose = 0;
  //Configuration flags
  fConfigureOutputsFlag = 0;
  fConfigureFluxProductsFlag = 0; 
  fConfigureRingMethodFlag = 0;
  fConfigureMapsFlag = 0; 
  fConfigureAcceptanceFlag = 0; 
  fConfigureExclusionsFlag = 0;
  fStartConfigureFlag = 0;

  fEntriesVector.clear();
  fScaleFactor.clear();
  fRadialIntegrationEbin.clear();
}

/****************************************************************************************************************/
/************************************************* CONFIGURATION ************************************************/
/****************************************************************************************************************/

/**
 * Starting configuration : clear everything
 */
void SurveySuite::MapMaker::StartConfigure()
{
  Clear(); 
  GetStartConfigureFlag() = true;
}


/**
 * Configuration of the exclusion regions. If both parameters are false, the regions from the analysis result file will be taken.
 * \param ExclusionFromFits : if true, take the specified FITS mask
 * \param ExclusionFromRegionFile : if true, take circular regions from an ascii file with l, b, radius
 * 
 */
bool SurveySuite::MapMaker::ConfigureExclusions(bool ExclusionFromFits, 
						bool ExclusionFromRegionFile)
{
  GetExclusionFromFits() = ExclusionFromFits;
  GetExclusionFromRegionFile() = ExclusionFromRegionFile;
  
  if(!ExclusionFromRegionFile && !ExclusionFromFits)
    std::cout << "WARNING: Eclusion regions will be taken from Result File" << std::endl;
  
  GetConfigureExclusionsFlag() = true;
  return 1;
}


/**
 * Configuration of the Acceptance calculation.
 *
 * \param EventsAndAcceptanceFromFile : if true, take events and acceptance maps from the analysis result file. 
 * \param EventsAndAcceptanceFromRadial : if true, (events and) acceptance maps will be computed from radial lookups. 
 * \param EventsAndAcceptanceFromRadialFile : if true, events and acceptance will be taken from a file previously generated from this program (radial file or ringbg maps)
 * If the first three parameters are set to 0, the events and acceptance maps are assumed to be given in FITS format.
 * \param PsiCut : PsiCut value to apply in case of radial acceptance computation
 * \param FitAcceptance : wether the radial acceptance histogram is fitted or not
 * \param CorrectZenith : apply zenith angle correction to the acceptance maps (time consuming!)
 * \param SelectAllRunsContributingToTheMap : if custom source position or custom map size is used, do we select all runs contributing to the map or only those contributing to the source position.
 * \param ApplySafeThreshold : Applies the safe threshold cut 
 * \param SafeThresholdFromAcceptance : two methods are available threshold from acceptance or from bias. If true, take threshold from acceptance, if false threshold is taken from bias
 * \param SafeThresholdFromPsiCut : the threshold is computed at a given offset. If true, automatically takes compute the threshold at the offset value taken for the PsiCut
 * \param SafeThresholdFromPsiValue : the threshold is computed at a given offset. The value here is to chose the offset value for the safe threshold calculation. The parameter SafeThresholdFromPsiCut should be set to false
 * \param SafeThresholdRatioParam : this parameter allows to chose the percentage of the maximum acceptance or the bias value that defines the safe threshold energy cut
 *
 */
bool SurveySuite::MapMaker::ConfigureAcceptance(bool EventsAndAcceptanceFromFile, 
						bool EventsAndAcceptanceFromRadial,
						bool EventsAndAcceptanceFromRadialFile,
						double PsiCut, 
						bool FitAcceptance, 			     
						bool CorrectZenith,
						bool SelectAllRunsContributingToTheMap,
						bool ApplySafeThreshold,
						bool SafeThresholdFromAcceptance,
						bool SafeThresholdFromPsiCut,
						double SafeThresholdFromPsiValue,
						double SafeThresholdRatioParam)
{
  //if the next 3 values are == 0 -> Ev and Acc from FITS file!!!
  GetEventsAndAcceptanceFromFile() = EventsAndAcceptanceFromFile;
  GetEventsAndAcceptanceFromRadial() = EventsAndAcceptanceFromRadial;
  GetEventsAndAcceptanceFromRadialFile() = EventsAndAcceptanceFromRadialFile;
  GetPsiCut() = PsiCut;
  GetFitAcceptance() = FitAcceptance;
  GetCorrectZenith() = CorrectZenith;
  GetSelectAllRunsContributingToTheMap() = SelectAllRunsContributingToTheMap;
  GetApplySafeThreshold() = ApplySafeThreshold;
  GetSafeThresholdFromAcceptance() = SafeThresholdFromAcceptance;
  GetSafeThresholdFromPsiCut() = SafeThresholdFromPsiCut;
  GetSafeThresholdFromPsiValue() = SafeThresholdFromPsiValue;
  GetSafeThresholdRatioParam() = SafeThresholdRatioParam;
  
  if(!EventsAndAcceptanceFromRadial && !EventsAndAcceptanceFromFile)
    std::cout << "WARNING: Events and Acceptance maps will be taken from FITS Files" << std::endl;
  if(EventsAndAcceptanceFromRadial)
    if(PsiCut == 0){
      std::cout << "ERROR: Events and Acceptance Map from Radial but PsiCut = 0" << std::endl;
      GetConfigureAcceptanceFlag() = false;
      return 0;
    }
  GetConfigureAcceptanceFlag() = true;
  return 1;
}


/**
 * Configuration of the Maps calculation.
 *
 * \param UseConfigTarget : if true, take target position from the analysis result file. It defines the center of the maps. 
 * \param UseConfigMapParams : if true, take maps parameters (bin size, extensions) from the analysis result file. 
 * \param UserLambda : custom Lambda. 
 * \param UserBeta : custom Beta.
 * \param BinSize : custom BinSize.
 * \param ExtX : Map extension X.
 * \param ExtY : Map extension Y.
 * \param OSRadius : Oversampling radius.
 * \param EMin : energy min.
 * \param EMax : energy max. 
 *
 */
bool SurveySuite::MapMaker::ConfigureMaps(bool UseConfigTarget, 
					  bool UseConfigMapParams, 
					  double UserLambda, 
					  double UserBeta, 
					  double BinSize, 
					  double ExtX, 
					  double ExtY, 
					  double OSRadius, 
					  double EMin, 
					  double EMax)
{
  GetUseConfigTarget() = UseConfigTarget;
  GetUseConfigMapParams() = UseConfigMapParams;
  GetUserLambda() = UserLambda;
  GetUserBeta() = UserBeta;
  GetBinSize() = BinSize;
  GetExtX() = ExtX;
  GetExtY() = ExtY;
  GetOSRadius() = OSRadius;
  GetEMin() = EMin;
  GetEMax() = EMax;

  if(!UseConfigTarget)
    std::cout << "WARNING: Custom target position -> make sure Lambda and Beta are well filled" << std::endl;
  if(!UseConfigMapParams)
    if(ExtX == 0 || ExtY == 0){
      std::cout << "ERROR: Map parameters not taken from Config but ExtX = 0 or ExtY = 0" << std::endl;
      GetConfigureMapsFlag() = false;
      return 0;
    }
  GetConfigureMapsFlag() = true;
  return 1;
}


/**
 * Configuration of the Ring background method.
 *
 * \param AdaptFFT : do we use the adaptive ring method 
 * In case the adaptive ring method is used, the following parameters are used : 
 * \param AdaptCut_Alpha : if true, cut on alpha value to stop the increase of the ring. Otherwise, cut on the fraction of excluded area. 
 * \param ConstantArea : rings have constant area
 * \param ConstantThickness : rings have constant thickness
 * \param ConstantInnerRadius : rings have constant inner radius
 * \param SmoothRings : do we smooth the rings number to take in order to avoid sharp edges
 * \param RingMinimalRadius : minimal radius to start with
 * \param InnerRingMax : inner ring max radius
 * \param OuterRingMax : outer ring max radius
 * \param RingStep : size if the rings increase
 * \param RingParam_AreaOverPi : defines the ring size : area divided by pi. 
 * \param RingParam_ExcFracMax : maximal allowed excluded fraction
 * \param RingParam_Thickness : defines the ring size with the thickness
 * \param RingParam_AlphaMax : maximal allowed alpha value
 * In case the adaptive ring is not used : 
 * \param StandardRingRadius : ring radius (center)  
 * \param StandardRingThickness :  ring thickness (ring size is ring radius +/- thickness)
 * \param AverageOff
 *
 */
bool SurveySuite::MapMaker::ConfigureRingMethod(bool AdaptFFT, 
						bool AdaptCut_Alpha, 
						bool ConstantArea, 
						bool ConstantThickness, 
						bool ConstantInnerRadius,
						bool SmoothRings, 
						double RingMinimalRadius, 
						double InnerRingMax, 
						double OuterRingMax, 
						double RingStep, 
						double RingParam_AreaOverPi, 
						double RingParam_ExcFracMax, 
						double RingParam_Thickness, 
						double RingParam_AlphaMax, 
						double StandardRingRadius, 
						double StandardRingThickness, 
						bool AverageOff)
{
  GetAdaptFFT() = AdaptFFT;
  GetAdaptCut_Alpha() = AdaptCut_Alpha;
  GetConstantArea() = ConstantArea;
  GetConstantThickness() = ConstantThickness;
  GetConstantInnerRadius() = ConstantInnerRadius;
  GetSmoothRings() = SmoothRings;
  GetRingMinimalRadius() = RingMinimalRadius;
  GetInnerRingMax() = InnerRingMax;
  GetOuterRingMax() = OuterRingMax;
  GetRingStep() = RingStep;
  GetRingParam_AreaOverPi() = RingParam_AreaOverPi;
  GetRingParam_Thickness() = RingParam_Thickness;
  GetRingParam_ExcFracMax() = RingParam_ExcFracMax;
  GetRingParam_AlphaMax() = RingParam_AlphaMax;
  GetStandardRingRadius() = StandardRingRadius;
  GetStandardRingThickness() = StandardRingThickness;
  GetAverageOff() = AverageOff;

  if(!AdaptFFT){
    if(StandardRingRadius < 0.01 || StandardRingThickness < 0.01){
      std::cout << "WARNING: standard ring parameters are strange" << std::endl;
    }
  }
  else{
    if(!ConstantInnerRadius && !ConstantThickness && !ConstantArea){
      std::cout << "ERROR: No method chosen for Adaptive ring calculation" << std::endl;
      GetConfigureRingMethodFlag() = false;      
      return 0;
    }
  }
  GetConfigureRingMethodFlag() = true;
  return 1;
}

/**
 * Configuration of the flux products.
 *
 * \param ProduceFluxProducts : do we produce the flux maps
 * \param SpectralIndex : spectral index to be taken to compute the expected counts maps.
 * \param PointLikeFluxMaps : if true, the point-like spectral tables will be used.
 * \param ProduceSurfaceBrightnessMap : deprecated
 * \param ExposureMapsFromFits : if true, the exposure maps are expected to be in FITS format
 *
 */
bool SurveySuite::MapMaker::ConfigureFluxProducts(bool ProduceFluxProducts,
						  double SpectralIndex, 
						  bool PointLikeFluxMaps, 
						  bool ProduceSurfaceBrightnessMap,
						  bool ExposureMapsFromFits,
						  const char *AnalysisConfig)
{
  GetProduceFluxProducts() = ProduceFluxProducts;
  GetSpectralIndex() = SpectralIndex;
  GetPointLikeFluxMaps() = PointLikeFluxMaps;
  GetProduceSurfaceBrightnessMap() = ProduceSurfaceBrightnessMap;
  GetExposureMapsFromFits() = ExposureMapsFromFits;

  std::cout << ProduceFluxProducts << " " << PointLikeFluxMaps << " " << AnalysisConfig << std::endl;
  std::map<std::string, double> fMapOS; 
  fMapOS["Loose"]=TMath::Sqrt(0.0125);
  fMapOS["Std"]=TMath::Sqrt(0.01);
  fMapOS["Faint"]=TMath::Sqrt(0.005);

  if(!ProduceSurfaceBrightnessMap){
    if(ProduceFluxProducts && PointLikeFluxMaps){
      fOSRadius = fMapOS[AnalysisConfig];
      std::cout << "Flux Maps requested with PointLike Configuration (" << AnalysisConfig << ") -> Updating OSRadius to " << fOSRadius << std::endl;
    }
  }
  else{
    fOSRadius = fMapOS[AnalysisConfig];
    fPointLikeFluxMaps = 0;
  }
  GetConfigureFluxProductsFlag() = true;
  return 1;
  
}

/**
 * Configuration of the output of the program.
 *
 * \param SaveResults : if true, the results will be saved in root format
 * \param Verbose
 *
 */
bool SurveySuite::MapMaker::ConfigureOutputs(bool SaveResults, bool Verbose)
{ 
  GetSaveResults() = SaveResults;
  GetVerbose() = Verbose;
  GetConfigureOutputsFlag() = true;
  return 1;
}

/**
 * End of Configuration.
 *
 */

bool SurveySuite::MapMaker::EndConfigure(const char *RadialConfig, const char *Resfile, int RunNumberMin, int RunNumberMax, bool UseRunListToMatch, bool UseRunsToForbid, const char *RunListToMatch, const char *RunsToForbid, std::string table_path)
{
  if(GetEventsAndAcceptanceFromRadialFile()){
    return 1;
  }
  else if(GetConfigureOutputsFlag() && GetConfigureFluxProductsFlag() && GetConfigureRingMethodFlag() && GetConfigureMapsFlag() && GetConfigureAcceptanceFlag() && GetConfigureExclusionsFlag() && GetStartConfigureFlag()){  
    if(ConfigureRuns(Resfile, RunNumberMin, RunNumberMax, UseRunListToMatch, UseRunsToForbid, RunListToMatch, RunsToForbid, table_path)){
      if(GetApplySafeThreshold()){
	ComputeSafeThresholdPerRun(RadialConfig, Resfile, table_path);
      }
      return 1;
    }
    else{
      return 0;
    }
  }
  else{
    return 0;
  }
}

/**
 *
 * Configuration of the runs to used.
 * 
 * One can use a list of runs to match or a list of runs to forbid. 
 * If the production of flux products is required, an additionnal check will be made to automatically forbid runs which have suspicious spectral tables
 *
 */
bool SurveySuite::MapMaker::ConfigureRuns(const char *Resfile, int RunNumberMin, int RunNumberMax, bool UseRunListToMatch, bool UseRunsToForbid, const char *RunListToMatch, const char *RunsToForbid, std::string table_path)
{
  if(RunNumberMax != -1){
    if(RunNumberMin >= RunNumberMax){
      std::cout << "ERROR : Check the run range!!!" << std::endl;
      return 0;
    }
  }

  std::vector<int> ForbidList, MatchList;
  if(UseRunsToForbid){
    std::ifstream ff(RunsToForbid);
    int run;
    while(!ff.eof())
    {
      std::string myline;
      std::getline(ff,myline);
      if (myline.size()==0) continue;
      std::istringstream iss_myline(myline);
      
      iss_myline >> run;
      ForbidList.push_back(run);
    }
  }
  if(UseRunListToMatch){
    std::ifstream fm(RunListToMatch);
    int run;
    while(!fm.eof())
    {
      std::string myline;
      std::getline(fm,myline);
      if (myline.size()==0) continue;
      std::istringstream iss_myline(myline);
      
      iss_myline >> run;
      MatchList.push_back(run);
    }
  } 

  TFile *fileResults = TFile::Open(Resfile);
  gROOT->cd();
  Sash::DataSet *results = (Sash::DataSet *)fileResults->Get("results");
  results->GetEntry(0);
  
  Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();
  
  ParisAnalysis::AnalysisConfig *Config =  hess->Handle("", (ParisAnalysis::AnalysisConfig *) 0);
  Config->LoadAllMembers();
  
  Sash::DataSet *run_results = (Sash::DataSet *)fileResults->Get("run_results");
  gROOT->cd();
  
  for(int i = 0; i < run_results->GetEntries(); ++i)
    {
      run_results->GetEntry(i);
      const ParisAnalysis::DataStorageRun* dsrun = hess->Handle("GammaFOVStorage", (ParisAnalysis::DataStorageRun *) 0);
      int runnr = dsrun->GetRunNumber();
      std::cout << i+1 << "/" << run_results->GetEntries() << "-> Run Nr = " << runnr << " " << RunNumberMin << " " << RunNumberMax << std::endl;
      if(RunNumberMin != 0 && runnr <= RunNumberMin) 
	continue;
      if(RunNumberMax != -1 && runnr >= RunNumberMax)
	continue;

      if(UseRunsToForbid && IsRunInFile(runnr, ForbidList, 0)){
	continue;
      }
      if(UseRunListToMatch && IsRunInFile(runnr, MatchList, 1)){
	continue;
      }
      std::cout << "adding run " << runnr << " to the list of runs to use" << std::endl;
      fEntriesVector.push_back(i);
    }  

  if(fProduceFluxProducts){
    RemoveBadTablesRuns(Resfile, table_path);
  }

  std::cout << Utilities::TextStyle::Green() <<  "NUMBER OF RUNS >> END OF CONFIGURE RUNS" << Utilities::TextStyle::Reset() <<  std::endl;
  std::cout << "Final runlist has " << fEntriesVector.size() << " runs" << std::endl;  
  return 1;
}

/*
 * Function to check if the run is in the given file
 */
bool SurveySuite::MapMaker::IsRunInFile(int RunNumber, std::vector<int> File, bool ReverseAnswer)
{
  //  std::ifstream f(File);
  int run;
  bool runinfile = false;
  bool runnotinfile = false;
  //  while(!f.eof())
  std::vector<int>::iterator it = File.begin();
  for(; it != File.end(); ++it)
    {
      /*      std::string myline;
      std::getline(f,myline);
      if (myline.size()==0) continue;
      std::istringstream iss_myline(myline);
      
      iss_myline >> run;*/
      run = *it;
      std::cout << run << " " << RunNumber << std::endl;
      if(fabs(run - RunNumber) < 1){
	std::cout << " -> run in file" << std::endl;
	runinfile = true;
	if(!ReverseAnswer)
	  continue;
      }
    }
  
  if(!ReverseAnswer) return runinfile;
  else return !runinfile;
  
}

/*
 * Check the table corresponding to each run and if the masimum acceptance is suspiciously low, forbid the run. 
 */
void SurveySuite::MapMaker::RemoveBadTablesRuns(const char *Resfile, std::string table_path)
{
  std::cout << Utilities::TextStyle::Green() <<  "REMOVE BAD TABLES RUNS" << Utilities::TextStyle::Reset() <<  std::endl;
  double fSourceSize = 0.;
  if(!GetPointLikeFluxMaps()){
    fSourceSize = 1.;
    std::cout << "Extended tables will be taken!" << std::endl;
  }
  
  TFile *fileResults = TFile::Open(Resfile);
  gROOT->cd();
  Sash::DataSet *results = (Sash::DataSet *)fileResults->Get("results");
  results->GetEntry(0);
  
  Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();
  
  ParisAnalysis::AnalysisConfig *Config =  hess->Handle("", (ParisAnalysis::AnalysisConfig *) 0);
  Config->LoadAllMembers();
  
  Sash::DataSet *run_results = (Sash::DataSet *)fileResults->Get("run_results");
  gROOT->cd();

  //Retrieve tables
  char *s = gSystem->ExpandPathName(table_path.c_str());
  setenv("SPECTRUM_PATH",s,1);
  delete[] s;

  std::string fAcceptanceNameBase = "Combined";
  Bool_t fAutoAzimut=false;
  Float_t fSourceExtension=0.;
  Float_t fTheta2Cut = fSourceSize;
  Float_t fMuonEfficiency = 0.;
  Bool_t fCheckAnalysisVersion=false;
  Bool_t fVerbose=false;

  std::map<Int_t,std::string> fMapAzimuth;
  fMapAzimuth[0]   = "North";
  fMapAzimuth[180] = "South";
  
  std::map<Int_t, Spectrum::SpectrumTableFinderMulti*> fMapTableFinder;
  std::map<Int_t, std::map<Int_t, const Spectrum::AcceptanceTableEfficiency*> > fMapTableAcceptance;
  std::map<Int_t, std::map<Int_t, const Spectrum::ResolutionTableEfficiency*> > fMapTableResolution;

  for (std::map<Int_t,std::string>::const_iterator it = fMapAzimuth.begin(); it!=fMapAzimuth.end(); ++it) {
    std::string AcceptanceName = fAcceptanceNameBase + it->second;
    Float_t Azimuth = Float_t(it->first);
    
    Spectrum::SpectrumTableFinderMulti *table = new Spectrum::SpectrumTableFinderMulti(Spectrum::SpectrumTableFinderMulti::Models,AcceptanceName.c_str(),Azimuth,fAutoAzimut,fSourceExtension,fTheta2Cut,fMuonEfficiency,fCheckAnalysisVersion,fVerbose);
    table->Process(Sash::Folder::GetFolder("results"));
    // std::cout << Utilities::TextStyle::Blue() << "SpectrumTableFinderMulti[" << Azimuth << "] = " << table << Utilities::TextStyle::Reset() << std::endl;    
    fMapTableFinder[Azimuth] = table;
    std::cout << Utilities::TextStyle::Blue() << "SpectrumTableFinderMulti[" << Azimuth << "] = Done !" << Utilities::TextStyle::Reset() << std::endl;    
    std::cout << (void*)table << std::endl;
    
    for (int itel=3;itel<=4;++itel) {
      std::ostringstream oss_tabname;
      oss_tabname << AcceptanceName << "_Tels" << itel;
      fMapTableAcceptance[Azimuth][itel] = hess->Get<Spectrum::AcceptanceTableEfficiency>(oss_tabname.str().c_str());
      std::cout << Utilities::TextStyle::Blue() << "Retrieve  AcceptanceTableEfficiency[" << Azimuth << "][" << itel << "] Named : " << oss_tabname.str().c_str() << " (" << (void*)fMapTableAcceptance[Azimuth][itel] << ")" << Utilities::TextStyle::Reset() << std::endl;
      fMapTableResolution[Azimuth][itel] = hess->Get<Spectrum::ResolutionTableEfficiency>(oss_tabname.str().c_str());
      std::cout << Utilities::TextStyle::Blue() << "Retrieve  ResolutionTableEfficiency[" << Azimuth << "][" << itel << "] Named : " << oss_tabname.str().c_str() << " (" << (void*)fMapTableResolution[Azimuth][itel] << ")" << Utilities::TextStyle::Reset() << std::endl;
      oss_tabname.str("");
    }
  }
  
  std::vector<int>::iterator it = fEntriesVector.begin();
  
  //  for(int i = 0; i < run_results->GetEntries(); ++i)
  int EntryCounter = 0;
  for(; it != fEntriesVector.end(); ++it)
    {
      int i = *it;
      run_results->GetEntry(i);
      const ParisAnalysis::DataStorageRun* dsrun = hess->Handle("GammaFOVStorage", (ParisAnalysis::DataStorageRun *) 0);

      std::cout << i+1 << "/" << run_results->GetEntries() << std::endl;

      double MeanZen = dsrun->GetMeanZenith();
      double meanCosZen = cos(MeanZen*TMath::Pi()/180.);
      double RelativeEfficiency = dsrun->GetMuonEfficiency();
      
      // Retrieve Azimuth
      Double_t MeanAzimuth = dsrun->GetMeanAzimuth();
      if (MeanAzimuth > 360. || MeanAzimuth<0.) {
	std::cout << Utilities::TextStyle::Red() << "Problem MeanAzimuth = " << MeanAzimuth << " BADLY TAKEN INTO ACCOUNT (I TRY A FIX, BUT CHECK IF THIS IS OK) !!!" << Utilities::TextStyle::Reset() << std::endl;
      }
      
      if (MeanAzimuth<0) {MeanAzimuth+=360.;}
      if (MeanAzimuth>360) {MeanAzimuth=MeanAzimuth-360.;}
      Int_t AzimuthCode = ( (MeanAzimuth<=90. || MeanAzimuth>270) ? 0 : 180 );

      // Retrieve LiveTime
      Double_t LiveTime = dsrun->GetOnLiveTime();
      
      // Retrieve TelescopeNumber and Pattern !
      Int_t NTels = dsrun->GetTelsInRun().size();
    
      // Get the Acceptance/Resolution Tables for the current run
      const Spectrum::AcceptanceTableEfficiency *AccTab = fMapTableAcceptance[AzimuthCode][NTels];
      if (!AccTab) {
	std::cout << Utilities::TextStyle::Yellow() << " Can't find AcceptanceTableEfficiency for AzimuthCode : " << AzimuthCode << " NTels : " << NTels << Utilities::TextStyle::Reset() << std::endl;
	continue;
      }
      
      const Spectrum::ResolutionTableEfficiency *ResTab = fMapTableResolution[AzimuthCode][NTels];
      if (!ResTab) {
	std::cout << Utilities::TextStyle::Yellow() << " Can't find ResolutionTableEfficiency for AzimuthCode : " << AzimuthCode << " NTels : " << NTels << Utilities::TextStyle::Reset() << std::endl;
	continue;
      }
      std::cout << meanCosZen << " " << RelativeEfficiency << " " << LiveTime << " " << NTels << std::endl;

      double MuonEff = dsrun->GetMuonEfficiency();
      double offaxisAngle = fPsiCut;
      if(!fSafeThresholdFromPsiCut)
	offaxisAngle = fSafeThresholdFromPsiValue;


      Double_t lnEmin = log(0.02);
      Double_t lnEmax = log(20.);
      Double_t lnEstep = 0.01;
      Double_t MaxAcceptance = 0;
      Double_t CurrentAcceptance = 0;
      Double_t CurrentBias = 1;
      Double_t CurrentBiasFromGraph = 1;
      Double_t lnE = lnEmin;
      Double_t lnEThresh = lnEmin;
      Double_t lnEThreshBiasFromGraph = lnEmin;
      Double_t SafeThreshold = 0;
      
      //----Acceptance method      
      while(lnE < lnEmax) {
	CurrentAcceptance = AccTab->GetAcceptance(lnE,meanCosZen,offaxisAngle,MuonEff);
	if(CurrentAcceptance>MaxAcceptance)
	  MaxAcceptance = CurrentAcceptance;     
	lnE += lnEstep;
      }
      
      std::cout << MaxAcceptance << std::endl;
      
      if(MaxAcceptance < 1){	
	std::cout << "Strange table .... run will not be used !!!" << std::endl;
	fEntriesVector.erase(it);	
	--it;
      }
    }
  std::cout << fEntriesVector.size() << std::endl;
}

void SurveySuite::MapMaker::ComputeSafeThresholdPerRun(const char *RadialConfig, const char *Resfile, std::string table_path)
{
 double fSourceSize = 0.;
  if(!GetPointLikeFluxMaps()){
    fSourceSize = 1.;
    std::cout << "Extended tables will be taken!" << std::endl;
  }

  std::cout << RadialConfig << std::endl;
  std::ostringstream f;
  f << "RadialTables/RadialTablesEnergy_" << RadialConfig << ".root";  
  TFile *RadialFile = TFile::Open(f.str().c_str());


  TFile *fileResults = TFile::Open(Resfile);
  gROOT->cd();
  Sash::DataSet *results = (Sash::DataSet *)fileResults->Get("results");
  results->GetEntry(0);
  
  Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();
  
  ParisAnalysis::AnalysisConfig *Config =  hess->Handle("", (ParisAnalysis::AnalysisConfig *) 0);
  Config->LoadAllMembers();
  
  Sash::DataSet *run_results = (Sash::DataSet *)fileResults->Get("run_results");
  gROOT->cd();

  //Retrieve tables
  char *s = gSystem->ExpandPathName(table_path.c_str());
  setenv("SPECTRUM_PATH",s,1);
  delete[] s;

  std::string fAcceptanceNameBase = "Combined";
  Bool_t fAutoAzimut=false;
  Float_t fSourceExtension=0.;
  Float_t fTheta2Cut = fSourceSize;
  Float_t fMuonEfficiency = 0.;
  Bool_t fCheckAnalysisVersion=false;
  Bool_t fVerbose=false;

  std::map<Int_t,std::string> fMapAzimuth;
  fMapAzimuth[0]   = "North";
  fMapAzimuth[180] = "South";
  
  std::map<Int_t, Spectrum::SpectrumTableFinderMulti*> fMapTableFinder;
  std::map<Int_t, std::map<Int_t, const Spectrum::AcceptanceTableEfficiency*> > fMapTableAcceptance;
  std::map<Int_t, std::map<Int_t, const Spectrum::ResolutionTableEfficiency*> > fMapTableResolution;
  
  for (std::map<Int_t,std::string>::const_iterator it = fMapAzimuth.begin(); it!=fMapAzimuth.end(); ++it) {
    std::string AcceptanceName = fAcceptanceNameBase + it->second;
    Float_t Azimuth = Float_t(it->first);
    
    Spectrum::SpectrumTableFinderMulti *table = new Spectrum::SpectrumTableFinderMulti(Spectrum::SpectrumTableFinderMulti::Models,AcceptanceName.c_str(),Azimuth,fAutoAzimut,fSourceExtension,fTheta2Cut,fMuonEfficiency,fCheckAnalysisVersion,fVerbose);
    table->Process(Sash::Folder::GetFolder("results"));
    // std::cout << Utilities::TextStyle::Blue() << "SpectrumTableFinderMulti[" << Azimuth << "] = " << table << Utilities::TextStyle::Reset() << std::endl;    
    fMapTableFinder[Azimuth] = table;
    std::cout << Utilities::TextStyle::Blue() << "SpectrumTableFinderMulti[" << Azimuth << "] = Done !" << Utilities::TextStyle::Reset() << std::endl;    
    std::cout << (void*)table << std::endl;
    
    for (int itel=3;itel<=4;++itel) {
      std::ostringstream oss_tabname;
      oss_tabname << AcceptanceName << "_Tels" << itel;
      fMapTableAcceptance[Azimuth][itel] = hess->Get<Spectrum::AcceptanceTableEfficiency>(oss_tabname.str().c_str());
      std::cout << Utilities::TextStyle::Blue() << "Retrieve  AcceptanceTableEfficiency[" << Azimuth << "][" << itel << "] Named : " << oss_tabname.str().c_str() << " (" << (void*)fMapTableAcceptance[Azimuth][itel] << ")" << Utilities::TextStyle::Reset() << std::endl;
      fMapTableResolution[Azimuth][itel] = hess->Get<Spectrum::ResolutionTableEfficiency>(oss_tabname.str().c_str());
      std::cout << Utilities::TextStyle::Blue() << "Retrieve  ResolutionTableEfficiency[" << Azimuth << "][" << itel << "] Named : " << oss_tabname.str().c_str() << " (" << (void*)fMapTableResolution[Azimuth][itel] << ")" << Utilities::TextStyle::Reset() << std::endl;
      oss_tabname.str("");
    }
  }
  
  std::vector<int>::iterator it = fEntriesVector.begin();
  
  //  for(int i = 0; i < run_results->GetEntries(); ++i)
  for(; it != fEntriesVector.end(); ++it)
    {
      int i = *it;
      run_results->GetEntry(i);
      const ParisAnalysis::DataStorageRun* dsrun = hess->Handle("GammaFOVStorage", (ParisAnalysis::DataStorageRun *) 0);

      std::cout << i+1 << "/" << run_results->GetEntries() << std::endl;

      double MeanZen = dsrun->GetMeanZenith();
      double meanCosZen = cos(MeanZen*TMath::Pi()/180.);
      double RelativeEfficiency = dsrun->GetMuonEfficiency();
      
      // Retrieve Azimuth
      Double_t MeanAzimuth = dsrun->GetMeanAzimuth();
      if (MeanAzimuth > 360. || MeanAzimuth<0.) {
	std::cout << Utilities::TextStyle::Red() << "Problem MeanAzimuth = " << MeanAzimuth << " BADLY TAKEN INTO ACCOUNT (I TRY A FIX, BUT CHECK IF THIS IS OK) !!!" << Utilities::TextStyle::Reset() << std::endl;
      }
      
      if (MeanAzimuth<0) {MeanAzimuth+=360.;}
      if (MeanAzimuth>360) {MeanAzimuth=MeanAzimuth-360.;}
      Int_t AzimuthCode = ( (MeanAzimuth<=90. || MeanAzimuth>270) ? 0 : 180 );

      // Retrieve LiveTime
      Double_t LiveTime = dsrun->GetOnLiveTime();
      
      // Retrieve TelescopeNumber and Pattern !
      Int_t NTels = dsrun->GetTelsInRun().size();
    
      // Get the Acceptance/Resolution Tables for the current run
      const Spectrum::AcceptanceTableEfficiency *AccTab = fMapTableAcceptance[AzimuthCode][NTels];
      if (!AccTab) {
	std::cout << Utilities::TextStyle::Yellow() << " Can't find AcceptanceTableEfficiency for AzimuthCode : " << AzimuthCode << " NTels : " << NTels << Utilities::TextStyle::Reset() << std::endl;
	continue;
      }
      
      const Spectrum::ResolutionTableEfficiency *ResTab = fMapTableResolution[AzimuthCode][NTels];
      if (!ResTab) {
	std::cout << Utilities::TextStyle::Yellow() << " Can't find ResolutionTableEfficiency for AzimuthCode : " << AzimuthCode << " NTels : " << NTels << Utilities::TextStyle::Reset() << std::endl;
	continue;
      }
      std::cout << meanCosZen << " " << RelativeEfficiency << " " << LiveTime << " " << NTels << std::endl;

      std::cout << RadialConfig << std::endl;
      TH2F *hEnergy = (TH2F *)RadialFile->Get("hR2D_4Tels_0");
      double int_emin = fEMin;
      double int_emax = fEMax;
      int binemin = 0;
      int binemax = -1;
      if(fEMin != 0){
	double log10emin = TMath::Log10(fEMin);
	binemin = hEnergy->GetXaxis()->FindFixBin(log10emin);
	int_emin = pow(10,hEnergy->GetXaxis()->GetBinLowEdge(binemin));
      }
      else{
	int_emin = 0.01;	
      }
      if(fEMax != -1){
	double log10emax = TMath::Log10(fEMax);
	binemax = hEnergy->GetXaxis()->FindFixBin(log10emax);
	int_emax = pow(10,hEnergy->GetXaxis()->GetBinLowEdge(binemax)+hEnergy->GetXaxis()->GetBinWidth(binemax));
      }
      else{
	int_emax = 200;
      }

      unsigned int nen = 20;
      TGraph *gBias = new TGraph(nen);

      double en[] = {0.02,0.03,0.05,0.08,0.125,0.2,0.3,0.5,0.8,1.25,
		     2.,3.,5.,8.,12.5,20.,30.,50.,80.,125.};
      std::vector<double> EnergyMC; ///< mc energy
      EnergyMC.assign(en,en+nen);
      
      double MuonEff = dsrun->GetMuonEfficiency();
      double offaxisAngle = fPsiCut;
      if(!fSafeThresholdFromPsiCut)
	offaxisAngle = fSafeThresholdFromPsiValue;

      //Loop on acceptance values at current zenith angle     
      Double_t lnEmin = log(0.02);
      Double_t lnEmax = log(20.);
      Double_t lnEstep = 0.01;
      Double_t MaxAcceptance = 0;
      Double_t CurrentAcceptance = 0;
      Double_t CurrentBias = 1;
      Double_t CurrentBiasFromGraph = 1;
      Double_t lnE = lnEmin;
      Double_t lnEThresh = lnEmin;
      Double_t lnEThreshBiasFromGraph = lnEmin;
      Double_t SafeThreshold = 0;
      
	if(GetSafeThresholdFromAcceptance()){
	  //----Acceptance method      
	  while(lnE < lnEmax) {
	    CurrentAcceptance = AccTab->GetAcceptance(lnE,meanCosZen,offaxisAngle,MuonEff);
	    if(CurrentAcceptance>MaxAcceptance)
	      MaxAcceptance = CurrentAcceptance;     
	    lnE += lnEstep;
	  }
	  //Find energy where acceptance is a given fraction  of its maximum
	  CurrentAcceptance = AccTab->GetAcceptance(lnEThresh,meanCosZen,offaxisAngle,MuonEff);
	  while(lnEThresh < lnEmax && CurrentAcceptance < MaxAcceptance * GetSafeThresholdRatioParam()) {
	    lnEThresh += lnEstep;
	    CurrentAcceptance = AccTab->GetAcceptance(lnEThresh,meanCosZen,offaxisAngle,MuonEff);
	  }
	  SafeThreshold = exp(lnEThresh);
	  if(MaxAcceptance < 1)
	    {
	      std::cout << "Strange Table -> taking threshold value from parametrisation" << std::endl;
	      SafeThreshold = 0.25+0.0085*MeanZen+exp(-7.94+0.14*MeanZen);
	    }
      	}
	else{
	  Int_t mcnpointres(200); // overkill ?
	  Double_t log_ereco_over_etrue_min = -3.;
	  Double_t log_ereco_over_etrue_max = 3.;
	  Double_t log_e_step = (log_ereco_over_etrue_max-log_ereco_over_etrue_min)/Double_t(mcnpointres-1.);
	  
	  for(unsigned int i = 0; i < nen; ++i)
	    {
	      TH1F *hBias = new TH1F("Bias","Bias",200,-3,3);
	      Double_t etrue = EnergyMC[i]; 
	      for (int imcbin=0;imcbin<mcnpointres;++imcbin) {
		Double_t log_ereco_over_etrue = log_ereco_over_etrue_min + Double_t(imcbin)*log_e_step;
		Double_t mce = etrue*TMath::Exp(log_ereco_over_etrue);
		double pdf = 0;
		pdf = ResTab->GetPDFvalue( TMath::Log(etrue), TMath::Log(mce), meanCosZen,offaxisAngle,MuonEff);
		hBias->SetBinContent(imcbin+1,pdf);
	      }	    
	      gBias->SetPoint(i,log(EnergyMC[i]),hBias->GetMean());
	      delete hBias;
	    }
	  
	  //----Bias Method
	  lnE = lnEmin;
	  CurrentBiasFromGraph = gBias->Eval(lnE);
	  
	  if(gBias->Eval(lnE) < GetSafeThresholdRatioParam()){
	    while(fabs(gBias->Eval(lnE+lnEstep)-CurrentBiasFromGraph)<0.001){
	      CurrentBiasFromGraph = gBias->Eval(lnE);
	      lnE += lnEstep;
	    } 
	  }
	  else{
	    while(CurrentBiasFromGraph > GetSafeThresholdRatioParam()){
	      CurrentBiasFromGraph = gBias->Eval(lnE);
	      lnE += lnEstep;
	    }
	  }
	  lnEThreshBiasFromGraph = lnE;
	  SafeThreshold = exp(lnEThreshBiasFromGraph);
	}
	std::cout << "SafeThreshold = " << SafeThreshold << std::endl;	
	fPerRunSafeThreshold[dsrun->GetRunNumber()]=SafeThreshold;
    }  
}

void SurveySuite::MapMaker::PrintConfigSummary()
{
  //Print Config Summary//
  std::cout << "**** CONFIG ****" << std::endl;
  std::cout << "OverSampling Radius = " << fOSRadius << std::endl;
  std::cout << "EMin = " << fEMin << " - EMax = " << fEMax << std::endl;
  std::cout << "** Files:" << std::endl;
  if(fEventsAndAcceptanceFromRadial){
    std::cout << "- Events and Acceptance from RADIAL" << std::endl;
    std::cout << "\t- PsiCut = " << fPsiCut << std::endl;
    std::cout << "\t- Zenith Correct = " << fCorrectZenith << std::endl;
  }
  else if(fEventsAndAcceptanceFromFile)
    std::cout << "- Events and Acceptance from FILE" << std::endl;
  else
    std::cout << "- Events and Acceptance from FITS" << std::endl;
  
  if(fExclusionFromRegionFile)
    std::cout << "- Exclusion from REGION FILE" << std::endl;
  else if(fExclusionFromFits)
    std::cout << "- Exclusion from FITS" << std::endl;
  else
    std::cout << "- Exclusion from FILE" << std::endl;
  
  std::cout << "** Maps:" << std::endl;
  if(fUseConfigTarget)
    std::cout << "- Target taken from result file" << std::endl;
  else 
    std::cout << "- Custom Target : Lambda = " << fUserLambda << " - Beta = " << fUserBeta << std::endl;
  
  if(fUseConfigMapParams)
    std::cout << "- Map parameters taken from result file" << std::endl;
  else{
    std::cout << "- Custom Map parameters : " << std::endl;
    std::cout << "\t" << "- BinSize = " << fBinSize << std::endl;
    std::cout << "\t" << "- ExtX = " << fExtX << ", ExtY = " << fExtY << std::endl;
  }
  
  std::cout << "** Ring Parameters:" << std::endl;
  if(!fAdaptFFT){
    std::cout << "- Standard Ring Bg method" << std::endl;
    std::cout << "\t" << "- RingRadius = " << fStandardRingRadius << std::endl;
    std::cout << "\t" << "- RingThickness = " << fStandardRingThickness << std::endl;
  }
  else{
    std::cout << "- Adaptive Ring Bg method" << std::endl;
    if(fAdaptCut_Alpha)
      std::cout << "- Cut will be applied on Alpha value : " << fRingParam_AlphaMax << std::endl;
    else
      std::cout << "- Cut will be applied on Excluded fraction of area : " << fRingParam_ExcFracMax << std::endl;
    
    if(fConstantArea){
      std::cout << "> ConstantArea : " << fRingParam_AreaOverPi << std::endl;
      std::cout << "\t+ InnerRingMin : " << fRingMinimalRadius << std::endl;
      std::cout << "\t+ InnerRingMax : " << fInnerRingMax << std::endl;
      std::cout << "\t+ InnerRingStep : " << fRingStep << std::endl;
    }
    else if(fConstantThickness){
      std::cout << "> ConstantThickness : " << fRingParam_Thickness << std::endl;
      std::cout << "\t+ InnerRingMin : " << fRingMinimalRadius << std::endl;
      std::cout << "\t+ InnerRingMax : " << fInnerRingMax << std::endl;
      std::cout << "\t+ InnerRingStep : " << fRingStep << std::endl;
    }
    else{
      if(fConstantInnerRadius){
	std::cout << "> ConstantInnerRadius : " << std::endl;
	std::cout << "\t+ InnerRingMin : " << fRingMinimalRadius << std::endl;
	std::cout << "\t+ OuterRingMin : " << fRingMinimalRadius+fRingParam_Thickness << std::endl;
	std::cout << "\t+ OuterRingMax : " << fOuterRingMax << std::endl;
	std::cout << "\t+ OuterRingStep : " << fRingStep << std::endl;
      }
    }
    if(fProduceFluxProducts){
      std::cout << "Produce Flux Maps" << std::endl;
      if(fEMin != 0 || fEMax != -1)
	if(!fEventsAndAcceptanceFromRadial)
	  std::cout << "WARNING!!! Non all energy range -> Flux Maps might be wrong!!!" << std::endl; 
    }
  }
}

/****************************************************************************************************************/
/*************************************************** UTILITIES **************************************************/
/****************************************************************************************************************/


void SurveySuite::MapMaker::RetrieveRunStat(const char *Resfile)
{
  TFile *fileResults = TFile::Open(Resfile);
  gROOT->cd();
  Sash::DataSet *results = (Sash::DataSet *)fileResults->Get("results");
  results->GetEntry(0);
  
  Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();
  ParisAnalysis::AnalysisConfig *Config =  hess->Handle("", (ParisAnalysis::AnalysisConfig *)0);
  ParisAnalysis::Analysis2DResults * r = hess->Handle("", (ParisAnalysis::Analysis2DResults *)0);
  ParisAnalysis::AcceptanceMap * acc = hess->Handle("", (ParisAnalysis::AcceptanceMap *)0);

  Config->LoadAllMembers();
  r->LoadAllMembers();
  acc->LoadAllMembers();

  std::vector<ParisAnalysis::RunStat> runs = acc->GetRuns();
  std::cout <<   (*acc->FindRun(18422)).GetRunNum() << std::endl;
  std::cout << "Printing RUNNUMBERS !!!! " << std::endl;
  for(std::vector<ParisAnalysis::RunStat>::const_iterator run = runs.begin();
      run != runs.end(); ++run) {
    std::cout << run->GetRunNum() << " " << run->GetGammaAcceptanceGradientX() << " " << run->GetGammaAcceptanceGradientY() << std::endl;;
  }
}

void SurveySuite::MapMaker::DisplayEventsMap(const char *Resfile)
{
  TFile *fileResults = TFile::Open(Resfile);
  gROOT->cd();
  Sash::DataSet *results = (Sash::DataSet *)fileResults->Get("results");
  results->GetEntry(0);
  
  Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();
  ParisAnalysis::AnalysisConfig *Config =  hess->Handle("", (ParisAnalysis::AnalysisConfig *)0);
  ParisAnalysis::Analysis2DResults * r = hess->Handle("", (ParisAnalysis::Analysis2DResults *)0);
  Config->LoadAllMembers();
  r->LoadAllMembers();
  
  Display::SkyHistogram2D *hGammaMap = r->GetGammaCandidatesMap();
  hGammaMap->Draw("colz");
}

void SurveySuite::MapMaker::DrawRadialLookup(const char *Config, double ZenithAngle)
{
  std::cout << Config << std::endl;
  
  std::ostringstream f;
  f << "RadialTables/RadialTablesEnergy_" << Config << ".root";

  TFile *RadialFile = TFile::Open(f.str().c_str());
  TH1D *hR3 = new TH1D("hR3","hR3",100,0,9);    
  TH1D *hR4 = new TH1D("hR4","hR4",100,0,9);    
  TH1D *hR3Fit = new TH1D("hR3Fit","hR3Fit",100,0,9);    
  TH1D *hR4Fit = new TH1D("hR4Fit","hR4Fit",100,0,9);   
  TH1F *hRes = new TH1F("hRes","hRes",100,0,9);    
  TH1F *hRes3Fit = new TH1F("hRes3Fit","hRes3Fit",100,0,9);    
  TH1F *hRes4Fit = new TH1F("hRes4Fit","hRes4Fit",100,0,9);    

  std::ostringstream fzenrad3;
  std::ostringstream fzenrad4;
  fzenrad3 << "hR2D_3Tels_" << int(ZenithAngle*1./10)*10;
  fzenrad4 << "hR2D_4Tels_" << int(ZenithAngle*1./10)*10;
  TH2F *h3Tels = (TH2F *)RadialFile->Get(fzenrad3.str().c_str());
  TH2F *h4Tels = (TH2F *)RadialFile->Get(fzenrad4.str().c_str());
  
  double emin = fEMin;
  double emax = fEMax;
  int binemin = 0;
  int binemax = -1;
  if(fEMin != 0){
    double log10emin = TMath::Log10(fEMin);
    binemin = h3Tels->GetXaxis()->FindFixBin(log10emin);
    emin = pow(10,h3Tels->GetXaxis()->GetBinLowEdge(binemin));
  }
  if(fEMax != -1){
    double log10emax = TMath::Log10(fEMax);
    binemax = h3Tels->GetXaxis()->FindFixBin(log10emax)-1;
    emax = pow(10,h3Tels->GetXaxis()->GetBinLowEdge(binemax)+h3Tels->GetXaxis()->GetBinWidth(binemax));
  }
  std::cout << fEMin << " - " << fEMax << " -> " << emin << "(bin " << binemin << ")" << " - " << emax << "(bin " << binemax << ")" << std::endl;  

  hR3 = h3Tels->ProjectionY("p3",binemin,binemax);
  hR4 = h4Tels->ProjectionY("p4",binemin,binemax);
  std::cout << hR3->GetMaximum() << " " << hR4->GetMaximum() << std::endl;
  hR3->Scale(1./hR3->GetMaximum());
  hR4->Scale(1./hR4->GetMaximum());
  std::cout << hR3->GetMaximum() << " " << hR4->GetMaximum() << std::endl;

  //Fit 
  TF1 *fR3 = new TF1("fR3","pol6",0,9);
  hR3->Fit(fR3,"QN");
  TF1 *fR4 = new TF1("fR4","pol6",0,9);
  hR4->Fit(fR4,"QN");
  
  for(int i = 1; i <= 100; ++i)
    hR3Fit->SetBinContent(i,fR3->Eval(hR3->GetBinCenter(i)));
  for(int i = 1; i <= 100; ++i)
    hR4Fit->SetBinContent(i,fR4->Eval(hR4->GetBinCenter(i)));


  for(int k = 1; k <= hR3->GetNbinsX(); ++k)
    {
      if(hR4->GetBinContent(k) != 0)
	hRes->SetBinContent(k,(hR3->GetBinContent(k)-hR4->GetBinContent(k))*1./hR4->GetBinContent(k));
    }

  for(int k = 1; k <= hR3->GetNbinsX(); ++k)
    {
      if(hR4->GetBinContent(k) != 0)
	hRes3Fit->SetBinContent(k,(hR3Fit->GetBinContent(k)-hR3->GetBinContent(k))*1./hR3->GetBinContent(k));
    }

  for(int k = 1; k <= hR3->GetNbinsX(); ++k)
    {
      if(hR4->GetBinContent(k) != 0)
	hRes4Fit->SetBinContent(k,(hR4Fit->GetBinContent(k)-hR4->GetBinContent(k))*1./hR4->GetBinContent(k));
    }

  TCanvas *cRadial = new TCanvas("Radial","Radial",1000,800);
  cRadial->Divide(3,2);
  cRadial->cd(1);
  hR3Fit->Draw();
  hR4Fit->SetLineColor(2);
  hR4Fit->Draw("same");
  cRadial->cd(4);
  hRes->Draw();
  cRadial->cd(2);
  hR3->Draw();
  hR3Fit->Draw("same");
  cRadial->cd(5);
  hRes3Fit->Draw();
  cRadial->cd(3);
  hR4->Draw();
  hR4Fit->SetLineColor(2);
  hR4Fit->Draw("same");
  cRadial->cd(6);
  hRes4Fit->Draw();

}

void SurveySuite::MapMaker::DisplayEventsAndAcceptanceMaps()
{
  TCanvas *cEvAcc = new TCanvas("EventsAndAcc","Events And Acceptance Maps",1000,400);
  cEvAcc->Divide(2,1);
  cEvAcc->cd(1);
  GetEventsMap()->Draw("colz");  
  cEvAcc->cd(2);
  GetAcceptanceMap()->Draw("colz");  
}

void SurveySuite::MapMaker::CreateEventsMap_FromFITS(const char *EventsMapFile)
{
  Display::SkyHistogram2D *hEv = new Display::SkyHistogram2D(EventsMapFile); 
  GetEventsMap() = hEv;
}

void SurveySuite::MapMaker::CreateAcceptanceMap_FromFITS(const char *AccMapFile)
{
  Display::SkyHistogram2D *hAcc = new Display::SkyHistogram2D(AccMapFile); 
  GetAcceptanceMap() = hAcc;
}

void SurveySuite::MapMaker::CreateExclusionMap_FromFITS(const char *ExclusionMapFile)
{
  Display::SkyHistogram2D *hEx = new Display::SkyHistogram2D(ExclusionMapFile); 
  GetExclusionMap() = hEx;
}

void SurveySuite::MapMaker::CreateEventsMap_FromResultFile(const char *Resfile)
{
  TFile *fileResults = TFile::Open(Resfile);
  gROOT->cd();
  Sash::DataSet *results = (Sash::DataSet *)fileResults->Get("results");
  results->GetEntry(0);
  
  Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();
  ParisAnalysis::Analysis2DResults * r = hess->Handle("", (ParisAnalysis::Analysis2DResults *)0);
  r->LoadAllMembers();
  
  Display::SkyHistogram2D *hGammaMap = r->GetGammaCandidatesMap();
  GetEventsMap() = hGammaMap;  

}

void SurveySuite::MapMaker::CreateAcceptanceMap_FromResultFile(const char *Resfile)
{
  TFile *fileResults = TFile::Open(Resfile);
  gROOT->cd();
  Sash::DataSet *results = (Sash::DataSet *)fileResults->Get("results");
  results->GetEntry(0);
  
  Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();
  ParisAnalysis::AcceptanceMap * acc = hess->Handle("", (ParisAnalysis::AcceptanceMap *)0);
  acc->LoadAllMembers();
  
  Display::SkyHistogram2D *hAcceptance = acc->GetSkyGammaAcceptanceMap();
  GetAcceptanceMap() = hAcceptance;
}

void SurveySuite::MapMaker::CreateExposureMaps_FromFITS(const char *ExposureMapFile)
{
  Display::SkyHistogram2D *hExpo = new Display::SkyHistogram2D(ExposureMapFile); 
  GetExpectedCountsMap() = hExpo;
}

void SurveySuite::MapMaker::CreateExposureMaps_FromFile(const char *ExposureMapFile)
{
  TFile *fileResults = TFile::Open(ExposureMapFile);
  gROOT->cd();
  Display::SkyHistogram2D *hExpo = (Display::SkyHistogram2D *)fileResults->Get("AeraTimeMap");
  GetExpectedCountsMap() = hExpo;
}

void SurveySuite::MapMaker::CreateExtendedExposureMaps_FromFITS(const char *ExposureMapFile)
{
  Display::SkyHistogram2D *hExpo = new Display::SkyHistogram2D(ExposureMapFile); 
  GetExtendedExpectedCountsMap() = hExpo;
}

void SurveySuite::MapMaker::CreateExtendedExposureMaps_FromFile(const char *ExposureMapFile)
{
  TFile *fileResults = TFile::Open(ExposureMapFile);
  gROOT->cd();
  Display::SkyHistogram2D *hExpo = (Display::SkyHistogram2D *)fileResults->Get("AeraTimeMap");
  GetExtendedExpectedCountsMap() = hExpo;
}

void SurveySuite::MapMaker::CreateExcessMaps_FromFITS(const char *ExcessMapFile)
{
  Display::SkyHistogram2D *hExcess = new Display::SkyHistogram2D(ExcessMapFile); 
  GetExcessMap() = hExcess;
}

void SurveySuite::MapMaker::CreateAlphaMaps_FromFITS(const char *AlphaMapFile)
{
  Display::SkyHistogram2D *hAlpha = new Display::SkyHistogram2D(AlphaMapFile); 
  GetAlphaMap() = hAlpha;
}

void SurveySuite::MapMaker::CreateOffMaps_FromFITS(const char *OffMapFile)
{
  Display::SkyHistogram2D *hOff = new Display::SkyHistogram2D(OffMapFile); 
  GetOffMap() = hOff;
}

/*
 * Create exclusion map depending on the choice in the configuration
 */
void SurveySuite::MapMaker::CreateExclusionMask(const char *Resfile, const char *Regionsfile)
{
  TFile *fileResults = TFile::Open(Resfile);
  gROOT->cd();
  Sash::DataSet *results = (Sash::DataSet *)fileResults->Get("results");
  results->GetEntry(0);
  Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();
  
  ParisAnalysis::AnalysisConfig *Config =  hess->Handle("", (ParisAnalysis::AnalysisConfig *) 0);
  Config->LoadAllMembers();
  
  if(fUseConfigTarget){
    std::cout << "Using target position from Results File" << std::endl;
    fUserLambda = Config->GetTargetRa();
    fUserBeta = Config->GetTargetDec();
    std::cout << "fUserLambda = " << fUserLambda << ", fUserBeta = " << fUserBeta << std::endl;
  }
  Stash::Coordinate targetpos(Stash::Lambda(fUserLambda,Stash::Angle::Degrees),
			      Stash::Beta(fUserBeta,Stash::Angle::Degrees),
			      Config->GetSystem());
  if(fUseConfigMapParams){
    std::cout << "Using Maps parameters from Results File" << std::endl;
    fBinSize = Config->GetMapBinSize().GetDegrees();
    fExtX  = Config->GetMapExtensionX().GetDegrees();
    fExtY  = Config->GetMapExtensionY().GetDegrees();
    std::cout << "fBinSize = " << fBinSize << ", fExtX = " << fExtX << ", fExtY = " << fExtY << std::endl;
  }

  //Creating Exclusion regions map
  if(fExclusionFromFits){
    std::cout << "Exclusion map generated from FITS file" << std::endl;
    CreateExclusionMap_FromFITS(Regionsfile);
    if(!fUseConfigMapParams){
      Display::SkyHistogram2D *hEx = new Display::SkyHistogram2D("Exc","Exc",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
      
      for(int k = 1; k <= hEx->GetNbinsX(); ++k)
	{
	  for(int l = 1; l <= hEx->GetNbinsY(); ++l)
	    {
	      Stash::Coordinate pop = hEx->GetBinCoordinate(k,l);	
	      Int_t bin = GetExclusionMap()->FindBinPosition(pop);
	      hEx->SetBinContent(k,l,GetExclusionMap()->GetBinContent(bin));
	    }
	}
      GetExclusionMap() = hEx;
    }

  }
  else if(fExclusionFromRegionFile){
    std::cout << "Exclusion map generated from List File" << std::endl;
    std::cout << "Coordinates must be in the same system !" << std::endl;
    std::vector<Stash::Coordinate> excpos;
    std::vector<double> excradius;
    std::ifstream reg(Regionsfile);
    double l, b, rad;
      
    Stash::Coordinate postarget(Stash::Lambda(fUserLambda,Stash::Angle::Degrees),
				Stash::Beta(fUserBeta,Stash::Angle::Degrees),
				Config->GetSystem());
    
    //Get useful exclusion regions
    double Extension = 0;
    if(fExtX > fExtY)
      Extension = fExtX+5;
    else
      Extension = fExtY+5;
      
    while(!reg.eof())
      {
	std::string myline;
	std::getline(reg,myline);
	if (myline.size()==0) continue;
	std::istringstream iss_myline(myline);
	  
	iss_myline >> l >> b >> rad;
	Stash::Coordinate posreg(Stash::Lambda(l,Stash::Angle::Degrees),
				 Stash::Beta(b,Stash::Angle::Degrees),
				 Config->GetSystem());
	if(posreg.GetAngularDistance(postarget).GetDegrees() < Extension){
	  excpos.push_back(posreg);
	  excradius.push_back(rad);
	}	  
      }

    //Actual Calculation of the regions map
    Display::SkyHistogram2D *hEx = new Display::SkyHistogram2D("Exc","Exc",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
      
    for(int k = 1; k <= hEx->GetNbinsX(); ++k)
      {
	for(int l = 1; l <= hEx->GetNbinsY(); ++l)
	  {
	    Stash::Coordinate pop = hEx->GetBinCoordinate(k,l);	
	    bool exc = 0;
	    for(int m = 0; m < excpos.size(); ++m)
	      {
		if(excpos[m].GetAngularDistance(pop).GetDegrees() < excradius[m]){
		  exc = 1;
		  break;
		}
	      }
	    if(exc)
	      hEx->SetBinContent(k,l,1);
	    else
	      hEx->SetBinContent(k,l,0);
	  }
      }
    GetExclusionMap() = hEx;      
  }
  else
    {
      std::cout << "Exclusion map generated from Results File" << std::endl;
      Display::SkyHistogram2D *hEx = new Display::SkyHistogram2D("Exc","Exc",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
      
      for(int k = 1; k <= hEx->GetNbinsX(); ++k)
	{
	  for(int l = 1; l <= hEx->GetNbinsY(); ++l)
	    {
	      Stash::Coordinate pop = hEx->GetBinCoordinate(k,l);	
	      if(Config->IsExcluded(pop,0,0))
		hEx->SetBinContent(k,l,1);
	      else
		hEx->SetBinContent(k,l,0);
	    }
	}
      GetExclusionMap() = hEx;
    }
}

/*
 * Create events and acceptance maps depending on the requested configuration
 */
void SurveySuite::MapMaker::CreateEventsAndAcceptanceMaps(const char *RadialConfig, const char *Resfile, const char *Regionsfile, const char *Evtfile, const char *Accfile, std::string table_path, const char *OutfilePrefix)
{
  TFile *fileResults = TFile::Open(Resfile);
  gROOT->cd();
  Sash::DataSet *results = (Sash::DataSet *)fileResults->Get("results");
  results->GetEntry(0);
  Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();
  
  ParisAnalysis::AnalysisConfig *Config =  hess->Handle("", (ParisAnalysis::AnalysisConfig *) 0);
  Config->LoadAllMembers();
  
  if(fUseConfigTarget){
    std::cout << "Using target position from Results File" << std::endl;
    fUserLambda = Config->GetTargetRa();
    fUserBeta = Config->GetTargetDec();
    std::cout << "fUserLambda = " << fUserLambda << ", fUserBeta = " << fUserBeta << std::endl;
  }
  Stash::Coordinate targetpos(Stash::Lambda(fUserLambda,Stash::Angle::Degrees),
			      Stash::Beta(fUserBeta,Stash::Angle::Degrees),
			      Config->GetSystem());
  if(fUseConfigMapParams){
    std::cout << "Using Maps parameters from Results File" << std::endl;
    fBinSize = Config->GetMapBinSize().GetDegrees();
    fExtX  = Config->GetMapExtensionX().GetDegrees();
    fExtY  = Config->GetMapExtensionY().GetDegrees();
    std::cout << "fBinSize = " << fBinSize << ", fExtX = " << fExtX << ", fExtY = " << fExtY << std::endl;
  }

  if(fEventsAndAcceptanceFromFile){
    std::cout << "Creating Events and Acceptance Maps from Result file" << std::endl;
    CreateEventsMap_FromResultFile(Evtfile);
    CreateAcceptanceMap_FromResultFile(Accfile);
        
    Display::SkyHistogram2D *hEvTemp = new Display::SkyHistogram2D("EvTemp","EvTemp",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
    Display::SkyHistogram2D *hAccTemp = new Display::SkyHistogram2D("AccTemp","AccTemp",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));

    for(int k = 1; k <= hEvTemp->GetNbinsX(); ++k)
      {
	for(int l = 1; l <= hEvTemp->GetNbinsY(); ++l)
	  {
	    Stash::Coordinate pop = hEvTemp->GetBinCoordinate(k,l);
	    double LambdaMod = pop.GetLambda().GetDegrees();
	    double BetaMod = pop.GetBeta().GetDegrees();
	    if(LambdaMod > 180)
	      LambdaMod = LambdaMod-360;
	    int a = GetEventsMap()->GetXaxis()->FindFixBin(-LambdaMod);
	    if(a == 0 || a > GetEventsMap()->GetNbinsX())
	      a = GetEventsMap()->GetXaxis()->FindFixBin(-(LambdaMod+360));
	    if(a == 0 || a > GetEventsMap()->GetNbinsX())
	      continue;
	    
	    int b = GetEventsMap()->GetYaxis()->FindFixBin(BetaMod);
	    hEvTemp->SetBinContent(k,l,GetEventsMap()->GetBinContent(a,b));
	  }
      }
    GetEventsMap()->Reset();
    GetEventsMap() = hEvTemp;
    
    for(int k = 1; k <= hAccTemp->GetNbinsX(); ++k)
      {
	for(int l = 1; l <= hAccTemp->GetNbinsY(); ++l)
	  {
	    Stash::Coordinate pop = hAccTemp->GetBinCoordinate(k,l);
	    double LambdaMod = pop.GetLambda().GetDegrees();
	    double BetaMod = pop.GetBeta().GetDegrees();
	    if(LambdaMod > 180)
	      LambdaMod = LambdaMod-360;
	    int a = GetAcceptanceMap()->GetXaxis()->FindFixBin(-LambdaMod);
	    if(a == 0 || a > GetAcceptanceMap()->GetNbinsX())
	      a = GetAcceptanceMap()->GetXaxis()->FindFixBin(-(LambdaMod+360));
	    if(a == 0 || a > GetAcceptanceMap()->GetNbinsX())
	      continue;
	    
	    int b = GetAcceptanceMap()->GetYaxis()->FindFixBin(BetaMod);
	    hAccTemp->SetBinContent(k,l,GetAcceptanceMap()->GetBinContent(a,b));
	  }
      }
    GetAcceptanceMap()->Reset();
    GetAcceptanceMap() = hAccTemp;
  }
  else if(fEventsAndAcceptanceFromRadial){
    std::cout << "Creating Events and Acceptance Maps from Radial Lookup" << std::endl;
    CreateEventsAndAcceptanceMaps_FromRadialLookups(RadialConfig, Resfile, Regionsfile, table_path, OutfilePrefix);
  }
  else if(fEventsAndAcceptanceFromRadialFile){
    std::cout << "Creating Events and Acceptance Maps from Radial File" << std::endl;
    TFile *fileResults = TFile::Open(Evtfile);
    gROOT->cd();
    GetAcceptanceMap() = (Display::SkyHistogram2D *)fileResults->Get("SkyAcceptance");
    GetEventsMap() = (Display::SkyHistogram2D *)fileResults->Get("SkyEvents");
    GetOffAxisMap() = (Display::SkyHistogram2D *)fileResults->Get("OffAxisMap");    
    GetZenithAngleMap() = (Display::SkyHistogram2D *)fileResults->Get("ZenithAngleMap");
    GetMuonEfficiencyMap() = (Display::SkyHistogram2D *)fileResults->Get("MuonEfficiencyMap");   
    GetMinSafeThresholdMap() = (Display::SkyHistogram2D *)fileResults->Get("MinSafeThresholdMap");
    GetAveragedSafeThresholdMap() = (Display::SkyHistogram2D *)fileResults->Get("AveragedSafeThresholdMap");
    GetAveragedLiveTimeMap() = (Display::SkyHistogram2D *)fileResults->Get("AveragedLiveTimeMap");
  }
  else{
    std::cout << "Creating Events and Acceptance Maps from FITS Files" << std::endl;
    CreateEventsMap_FromFITS(Evtfile);
    CreateAcceptanceMap_FromFITS(Accfile);    

    Display::SkyHistogram2D *hEvTemp = new Display::SkyHistogram2D("EvTemp","EvTemp",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
    Display::SkyHistogram2D *hAccTemp = new Display::SkyHistogram2D("AccTemp","AccTemp",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));

    for(int k = 1; k <= hEvTemp->GetNbinsX(); ++k)
      {
	for(int l = 1; l <= hEvTemp->GetNbinsY(); ++l)
	  {
	    Stash::Coordinate pop = hEvTemp->GetBinCoordinate(k,l);
	    double LambdaMod = pop.GetLambda().GetDegrees();
	    double BetaMod = pop.GetBeta().GetDegrees();
	    if(LambdaMod > 180)
	      LambdaMod = LambdaMod-360;
	    int a = GetEventsMap()->GetXaxis()->FindFixBin(-LambdaMod);
	    if(a == 0 || a > GetEventsMap()->GetNbinsX())
	      a = GetEventsMap()->GetXaxis()->FindFixBin(-(LambdaMod+360));
	    if(a == 0 || a > GetEventsMap()->GetNbinsX())
	      continue;
	    
	    int b = GetEventsMap()->GetYaxis()->FindFixBin(BetaMod);
	    hEvTemp->SetBinContent(k,l,GetEventsMap()->GetBinContent(a,b));
	  }
      }
    GetEventsMap()->Reset();
    GetEventsMap() = hEvTemp;
    
    for(int k = 1; k <= hAccTemp->GetNbinsX(); ++k)
      {
	for(int l = 1; l <= hAccTemp->GetNbinsY(); ++l)
	  {
	    Stash::Coordinate pop = hAccTemp->GetBinCoordinate(k,l);
	    double LambdaMod = pop.GetLambda().GetDegrees();
	    double BetaMod = pop.GetBeta().GetDegrees();
	    if(LambdaMod > 180)
	      LambdaMod = LambdaMod-360;
	    int a = GetAcceptanceMap()->GetXaxis()->FindFixBin(-LambdaMod);
	    if(a == 0 || a > GetAcceptanceMap()->GetNbinsX())
	      a = GetAcceptanceMap()->GetXaxis()->FindFixBin(-(LambdaMod+360));
	    if(a == 0 || a > GetAcceptanceMap()->GetNbinsX())
	      continue;
	    
	    int b = GetAcceptanceMap()->GetYaxis()->FindFixBin(BetaMod);
	    hAccTemp->SetBinContent(k,l,GetAcceptanceMap()->GetBinContent(a,b));
	  }
      }
    GetAcceptanceMap()->Reset();
    GetAcceptanceMap() = hAccTemp;
  }
}

/*
 * Generate runlist given the input parameters
 *
 */
void SurveySuite::MapMaker::CreateRunList(const char *Resfile, const char *Regionsfile, const char *OutfilePrefix, double SourceRadius, int sourceNumber)
{
  //Open result file to get proper config and access each run's parameter
  TFile *fileResults = TFile::Open(Resfile);
  gROOT->cd();
  Sash::DataSet *results = (Sash::DataSet *)fileResults->Get("results");
  results->GetEntry(0);
  
  Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();
  
  ParisAnalysis::AnalysisConfig *Config =  hess->Handle("", (ParisAnalysis::AnalysisConfig *) 0);
  Config->LoadAllMembers();
  
  Sash::DataSet *run_results = (Sash::DataSet *)fileResults->Get("run_results");
  gROOT->cd();

  //  std::ifstream fileList("SpectralAnalysis_New/HGPS3_Catalog_targets_mod_New.txt");
  //std::ifstream fileList("SpectralAnalysis_PublishedPositions/HGPS3_Association_targets_mod.txt");
  std::ifstream fileList("SpectralAnalysis_052015/HGPS_SourceList_SpecTestIn.txt");
  //Sash::HESSArray* hess = Sash::Folder::GetFolder(0)->GetHESSArray();
  int NUMBER;
  double GLON, GLAT, RADIUS, RADIUS2;
  while(!fileList.eof())
    {
      //      std::cout << "Getting in..." << std::endl;
      std::string myline;
      std::getline(fileList,myline);
      if (myline.size()==0) continue;
      std::istringstream iss_myline(myline);
      
      iss_myline >> NUMBER >> GLON >> GLAT >> RADIUS >> RADIUS2;
      
      std::cout << NUMBER << " " << GLON << " " << GLAT << " " << RADIUS << " " << RADIUS2 << std::endl;
            
      if(fUseConfigTarget){
	std::cout << "Using target position from Results File" << std::endl;
	fUserLambda = Config->GetTargetRa();
	fUserBeta = Config->GetTargetDec();
	std::cout << "fUserLambda = " << fUserLambda << ", fUserBeta = " << fUserBeta << std::endl;
      }
      
      fUserLambda = GLON;
      fUserBeta = GLAT;
      
      Stash::Coordinate targetpos(Stash::Lambda(fUserLambda,Stash::Angle::Degrees),
				  Stash::Beta(fUserBeta,Stash::Angle::Degrees),
				  Config->GetSystem());
      
      int skippedruns = 0;
      double maxExtension = sqrt(pow(fExtX,2)+pow(fExtY,2));
      double RunSelectionCut = 0;
      if(fSelectAllRunsContributingToTheMap)
	RunSelectionCut = fPsiCut+maxExtension;
      else
	RunSelectionCut = fPsiCut-RADIUS;
      
      std::vector<int>::iterator it = fEntriesVector.begin();
      
      //  for(int i = 0; i < run_results->GetEntries(); ++i)
      std::vector<int> selectedRuns;
      for(; it != fEntriesVector.end(); ++it)
	{
	  int i = *it;
	  run_results->GetEntry(i);
	  std::cout << i+1 << "/" << run_results->GetEntries() << std::endl;
	  const ParisAnalysis::DataStorageRun* dsrun = hess->Handle("GammaFOVStorage", (ParisAnalysis::DataStorageRun *) 0);
	  Stash::Coordinate pointingpos = dsrun->GetObsPos();      
	  
	  if(pointingpos.GetAngularDistance(targetpos).GetDegrees() > RunSelectionCut){
	    //	std::cout << "skipping run" << std::endl;
	    std::cout << Utilities::TextStyle::Yellow() <<  "SKIPPING RUN!!!" << Utilities::TextStyle::Reset() <<  std::endl;
	    
	    ++skippedruns;
	    continue;
	  }
	  
	  std::cout << dsrun->GetRunNumber() << std::endl;
	  selectedRuns.push_back(dsrun->GetRunNumber());
	}
      
      if(fSaveResults){
	//-----------------------------------------------------------------------------
	//Generate File Name according to analysis parameters    
	//and write in file
	std::ostringstream foutfile;
	
	foutfile << OutfilePrefix << "_Source_SpecTest.txt";

	
	std::ofstream ofile(foutfile.str().c_str());
	
	for(std::vector<int>::iterator itout = selectedRuns.begin(); itout != selectedRuns.end(); ++itout)
	  {
	    ofile << *itout << std::endl;
	  }
      }

    }

}

/*
 * Write runs parameters in a file
 *
 */
void SurveySuite::MapMaker::CreateRunFileInfos(const char *Resfile, const char *Regionsfile, const char *OutfilePrefix)
{
  //Open result file to get proper config and access each run's parameter
  TFile *fileResults = TFile::Open(Resfile);
  gROOT->cd();
  Sash::DataSet *results = (Sash::DataSet *)fileResults->Get("results");
  results->GetEntry(0);
  
  Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();
  
  ParisAnalysis::AnalysisConfig *Config =  hess->Handle("", (ParisAnalysis::AnalysisConfig *) 0);
  Config->LoadAllMembers();
  
  Sash::DataSet *run_results = (Sash::DataSet *)fileResults->Get("run_results");
  gROOT->cd();
  
  if(fUseConfigTarget){
    std::cout << "Using target position from Results File" << std::endl;
    fUserLambda = Config->GetTargetRa();
    fUserBeta = Config->GetTargetDec();
    std::cout << "fUserLambda = " << fUserLambda << ", fUserBeta = " << fUserBeta << std::endl;
  }
  
  Stash::Coordinate targetpos(Stash::Lambda(fUserLambda,Stash::Angle::Degrees),
			      Stash::Beta(fUserBeta,Stash::Angle::Degrees),
			      Config->GetSystem());
  
  int skippedruns = 0;
  double maxExtension = sqrt(pow(fExtX,2)+pow(fExtY,2));
  double RunSelectionCut = 0;
  if(fSelectAllRunsContributingToTheMap)
    RunSelectionCut = fPsiCut+maxExtension;
  else
    RunSelectionCut = fPsiCut;//-RADIUS;
      
  std::vector<int>::iterator it = fEntriesVector.begin();

  std::ostringstream foutfile;
  foutfile << OutfilePrefix << "_RunsInfos.txt";
  std::ofstream ofile(foutfile.str().c_str());
    
  for(; it != fEntriesVector.end(); ++it)
    {
      int i = *it;
      run_results->GetEntry(i);
      std::cout << i+1 << "/" << run_results->GetEntries() << std::endl;
      const ParisAnalysis::DataStorageRun* dsrun = hess->Handle("GammaFOVStorage", (ParisAnalysis::DataStorageRun *) 0);
      Stash::Coordinate pointingpos = dsrun->GetObsPos();      
      
      if(pointingpos.GetAngularDistance(targetpos).GetDegrees() > RunSelectionCut){
	//	std::cout << "skipping run" << std::endl;
	std::cout << Utilities::TextStyle::Yellow() <<  "SKIPPING RUN!!!" << Utilities::TextStyle::Reset() <<  std::endl;
	
	++skippedruns;
	continue;
      }
      
      std::cout << dsrun->GetRunNumber() << std::endl;

      ofile << dsrun->GetRunNumber() << " " << dsrun->GetOnLiveTime() << std::endl;
      
    }

}

/****************************************************************************************************************/
/*********************************************** MAPS COMPUTATION ***********************************************/
/****************************************************************************************************************/

/*
 * Computation of area time maps
 *
 */
void SurveySuite::MapMaker::CreateExposureMaps(const char *RadialConfig, const char *Resfile, const char *Regionsfile, std::string table_path, const char *OutfilePrefix)
{
  //Choose either point-like or extended configuration
  double fSourceSize = 0.;
  if(!GetPointLikeFluxMaps()){
    fSourceSize = 1.;
    std::cout << "Extended tables will be taken!" << std::endl;
  }
  
  std::cout << RadialConfig << std::endl;
  std::ostringstream f;
  f << "RadialTables/RadialTablesEnergy_" << RadialConfig << ".root";  
  TFile *RadialFile = TFile::Open(f.str().c_str());

  //Open result file to get proper config and access each run's parameter
  TFile *fileResults = TFile::Open(Resfile);
  gROOT->cd();
  Sash::DataSet *results = (Sash::DataSet *)fileResults->Get("results");
  results->GetEntry(0);
  
  Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();
  
  ParisAnalysis::AnalysisConfig *Config =  hess->Handle("", (ParisAnalysis::AnalysisConfig *) 0);
  Config->LoadAllMembers();
  
  Sash::DataSet *run_results = (Sash::DataSet *)fileResults->Get("run_results");
  gROOT->cd();

  //Prepare maps
  if(fUseConfigTarget){
    std::cout << "Using target position from Results File" << std::endl;
    fUserLambda = Config->GetTargetRa();
    fUserBeta = Config->GetTargetDec();
    std::cout << "fUserLambda = " << fUserLambda << ", fUserBeta = " << fUserBeta << std::endl;
  }

  Stash::Coordinate targetpos(Stash::Lambda(fUserLambda,Stash::Angle::Degrees),
			      Stash::Beta(fUserBeta,Stash::Angle::Degrees),
			      Config->GetSystem());
  
  if(fUseConfigMapParams){
    std::cout << "Using Maps parameters from Results File" << std::endl;
    fBinSize = Config->GetMapBinSize().GetDegrees();
    fExtX  = Config->GetMapExtensionX().GetDegrees();
    fExtY  = Config->GetMapExtensionY().GetDegrees();
    std::cout << "fBinSize = " << fBinSize << ", fExtX = " << fExtX << ", fExtY = " << fExtY << std::endl;
  }

  Double_t ExtX_Run = fPsiCut;
  Double_t ExtY_Run = fPsiCut;

  //Global Maps ...
  Display::SkyHistogram2D *hExpectedCountsMap = new Display::SkyHistogram2D("ExpectedCountsMap","ExpectedCountsMap",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
  Display::SkyHistogram2D *hAreaTimeMap = new Display::SkyHistogram2D("AeraTimeMap","AreaTimeMap",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));

  //Retrieve tables
  char *s = gSystem->ExpandPathName(table_path.c_str());
  setenv("SPECTRUM_PATH",s,1);
  delete[] s;

  std::string fAcceptanceNameBase = "Combined";
  Bool_t fAutoAzimut=false;
  Float_t fSourceExtension=0.;
  Float_t fTheta2Cut = fSourceSize;
  Float_t fMuonEfficiency = 0.;
  Bool_t fCheckAnalysisVersion=false;
  Bool_t fVerbose=false;

  std::map<Int_t,std::string> fMapAzimuth;
  fMapAzimuth[0]   = "North";
  fMapAzimuth[180] = "South";
  
  std::map<Int_t, Spectrum::SpectrumTableFinderMulti*> fMapTableFinder;
  std::map<Int_t, std::map<Int_t, const Spectrum::AcceptanceTableEfficiency*> > fMapTableAcceptance;
  std::map<Int_t, std::map<Int_t, const Spectrum::ResolutionTableEfficiency*> > fMapTableResolution;
  
  for (std::map<Int_t,std::string>::const_iterator it = fMapAzimuth.begin(); it!=fMapAzimuth.end(); ++it) {
    std::string AcceptanceName = fAcceptanceNameBase + it->second;
    Float_t Azimuth = Float_t(it->first);
    
    Spectrum::SpectrumTableFinderMulti *table = new Spectrum::SpectrumTableFinderMulti(Spectrum::SpectrumTableFinderMulti::Models,AcceptanceName.c_str(),Azimuth,fAutoAzimut,fSourceExtension,fTheta2Cut,fMuonEfficiency,fCheckAnalysisVersion,fVerbose);
    table->Process(Sash::Folder::GetFolder("results"));
    // std::cout << Utilities::TextStyle::Blue() << "SpectrumTableFinderMulti[" << Azimuth << "] = " << table << Utilities::TextStyle::Reset() << std::endl;    
    fMapTableFinder[Azimuth] = table;
    std::cout << Utilities::TextStyle::Blue() << "SpectrumTableFinderMulti[" << Azimuth << "] = Done !" << Utilities::TextStyle::Reset() << std::endl;    
    std::cout << (void*)table << std::endl;

    for (int itel=3;itel<=4;++itel) {
      std::ostringstream oss_tabname;
      oss_tabname << AcceptanceName << "_Tels" << itel;
      fMapTableAcceptance[Azimuth][itel] = hess->Get<Spectrum::AcceptanceTableEfficiency>(oss_tabname.str().c_str());
      std::cout << Utilities::TextStyle::Blue() << "Retrieve  AcceptanceTableEfficiency[" << Azimuth << "][" << itel << "] Named : " << oss_tabname.str().c_str() << " (" << (void*)fMapTableAcceptance[Azimuth][itel] << ")" << Utilities::TextStyle::Reset() << std::endl;
      fMapTableResolution[Azimuth][itel] = hess->Get<Spectrum::ResolutionTableEfficiency>(oss_tabname.str().c_str());
      std::cout << Utilities::TextStyle::Blue() << "Retrieve  ResolutionTableEfficiency[" << Azimuth << "][" << itel << "] Named : " << oss_tabname.str().c_str() << " (" << (void*)fMapTableResolution[Azimuth][itel] << ")" << Utilities::TextStyle::Reset() << std::endl;
      oss_tabname.str("");
    }
  }

  //AreaTime Map will be computed under the assumption of a power law spectrum 
  Spectrum::SpectrumPowerLaw *spl = hess->Handle("",(Spectrum::SpectrumPowerLaw*)0);
  std::cout << spl << std::endl;
  
  spl->SetParameters(1.,GetSpectralIndex());
  spl->GetRefEnergy() = 1.0;

  int skippedruns = 0;
  double maxExtension = sqrt(pow(fExtX,2)+pow(fExtY,2));
  double RunSelectionCut = 0;
  if(fSelectAllRunsContributingToTheMap)
    RunSelectionCut = fPsiCut+maxExtension;
  else
    RunSelectionCut = fPsiCut;

  std::ostringstream foutfile;
  foutfile << OutfilePrefix << "_RunsAreaInfos.txt";
  std::ofstream ofile(foutfile.str().c_str());

  ofile << "RunNumber" << " " << "Lambda" << " " << "Beta" << " " << "MeanZen" << " " << "RelativeEfficiency" << " " << "LiveTime" << " " << "SafeThreshold" << " " << "EnergyIntegrationBound" << " " << "TheoricRateAtCenter*LiveTime*1e8" << std::endl;      

  std::vector<int>::iterator it = fEntriesVector.begin();
  for(; it != fEntriesVector.end(); ++it)
    {
      int i = *it;
      run_results->GetEntry(i);
      std::cout << i+1 << "/" << run_results->GetEntries() << std::endl;
      const ParisAnalysis::DataStorageRun* dsrun = hess->Handle("GammaFOVStorage", (ParisAnalysis::DataStorageRun *) 0);
      Stash::Coordinate pointingpos = dsrun->GetObsPos();      
  
      if(pointingpos.GetAngularDistance(targetpos).GetDegrees() > RunSelectionCut){
	//	std::cout << "skipping run" << std::endl;
	std::cout << Utilities::TextStyle::Yellow() <<  "SKIPPING RUN!!!" << Utilities::TextStyle::Reset() <<  std::endl;
  
	++skippedruns;
	continue;
      }

	      
      Display::SkyHistogram2D *hRunExpectedCountsMap = new Display::SkyHistogram2D("RunExpectedCountsMap","RunExpectedCountsMap",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
      Display::SkyHistogram2D *hRunAreaTimeMap = new Display::SkyHistogram2D("RunAreaTimeMap","RunAreaTimeMap",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));

      double MeanZen = dsrun->GetMeanZenith();
      double meanCosZen = cos(MeanZen*TMath::Pi()/180.);
      double RelativeEfficiency = dsrun->GetMuonEfficiency();
      
      // Retrieve Azimuth
      Double_t MeanAzimuth = dsrun->GetMeanAzimuth();
      if (MeanAzimuth > 360. || MeanAzimuth<0.) {
	std::cout << Utilities::TextStyle::Red() << "Problem MeanAzimuth = " << MeanAzimuth << " BADLY TAKEN INTO ACCOUNT (I TRY A FIX, BUT CHECK IF THIS IS OK) !!!" << Utilities::TextStyle::Reset() << std::endl;
      }
      
      if (MeanAzimuth<0) {MeanAzimuth+=360.;}
      if (MeanAzimuth>360) {MeanAzimuth=MeanAzimuth-360.;}
      Int_t AzimuthCode = ( (MeanAzimuth<=90. || MeanAzimuth>270) ? 0 : 180 );

      // Retrieve LiveTime
      Double_t LiveTime = dsrun->GetOnLiveTime();
      
      // Retrieve TelescopeNumber and Pattern !
      Int_t NTels = dsrun->GetTelsInRun().size();
    
      // Get the Acceptance/Resolution Tables for the current run
      const Spectrum::AcceptanceTableEfficiency *AccTab = fMapTableAcceptance[AzimuthCode][NTels];
      if (!AccTab) {
	std::cout << Utilities::TextStyle::Yellow() << " Can't find AcceptanceTableEfficiency for AzimuthCode : " << AzimuthCode << " NTels : " << NTels << Utilities::TextStyle::Reset() << std::endl;
	continue;
      }
      
      const Spectrum::ResolutionTableEfficiency *ResTab = fMapTableResolution[AzimuthCode][NTels];
      if (!ResTab) {
	std::cout << Utilities::TextStyle::Yellow() << " Can't find ResolutionTableEfficiency for AzimuthCode : " << AzimuthCode << " NTels : " << NTels << Utilities::TextStyle::Reset() << std::endl;
	continue;
      }
      std::cout << meanCosZen << " " << RelativeEfficiency << " " << LiveTime << " " << NTels << std::endl;

      //Retrieve the radial acceptance table to have consistent energy range
      std::cout << RadialConfig << std::endl;
      TH2F *hEnergy = (TH2F *)RadialFile->Get("hR2D_4Tels_0");
      double int_emin = fEMin;
      double int_emax = fEMax;
      int binemin = 0;
      int binemax = -1;
      if(fEMin != 0){
	double log10emin = TMath::Log10(fEMin);
	binemin = hEnergy->GetXaxis()->FindFixBin(log10emin);
	int_emin = pow(10,hEnergy->GetXaxis()->GetBinLowEdge(binemin));
      }
      else{
	int_emin = 0.01;	
      }
      if(fEMax != -1){
	double log10emax = TMath::Log10(fEMax);
	binemax = hEnergy->GetXaxis()->FindFixBin(log10emax);
	int_emax = pow(10,hEnergy->GetXaxis()->GetBinLowEdge(binemax)+hEnergy->GetXaxis()->GetBinWidth(binemax));
      }
      else{
	int_emax = 200;
      }
      Double_t SafeThreshold = 0;
      if(GetApplySafeThreshold()){
	SafeThreshold = fPerRunSafeThreshold[dsrun->GetRunNumber()];
	std::cout << "SafeThreshold = " << SafeThreshold << std::endl;	

	int_emin = fEMin;
	binemin = 0;
	if(fEMax != -1 && SafeThreshold > fEMax){
	  std::cout << Utilities::TextStyle::Yellow() <<  "SKIPPING RUN!!!" << Utilities::TextStyle::Reset() <<  std::endl;
		
	  std::cout << "Skipping run : safe threshold above max energy." << std::endl; 
	  ++skippedruns;
	  continue;
	}
	if(SafeThreshold > fEMin){
	  double log10eth = TMath::Log10(SafeThreshold);
	  binemin = hEnergy->GetXaxis()->FindFixBin(log10eth);
	  int_emin = pow(10,hEnergy->GetXaxis()->GetBinLowEdge(binemin));	
	}
      }
      delete hEnergy;

      std::cout << "Integration will be performed between " << int_emin << " " << int_emax << std::endl;

      int nbinsArea = fPsiCut*1./0.15 + 1;
      std::cout << "nbins to compute theoric rate = " << nbinsArea << std::endl;
      //extending histogram beyond fPsiCut in order to ensure proper interpolation at PsiCut value. 
      TH1F *hRTemp = new TH1F("hRTemp","hRTemp",nbinsArea,0,fPsiCut+0.15);
      for(int r = 1; r <= nbinsArea; ++r){
	hRTemp->SetBinContent(r,spl->TheoricRate(TMath::Log(int_emin),TMath::Log(int_emax),meanCosZen,hRTemp->GetBinCenter(r),RelativeEfficiency,AccTab,ResTab));
      }

      double TheoricRateAtCenter = spl->TheoricRate(TMath::Log(int_emin),TMath::Log(int_emax),meanCosZen,0,RelativeEfficiency,AccTab,ResTab);
          
      for(int k = 1; k <= hRunExpectedCountsMap->GetNbinsX(); ++k)
	{
       	  for(int l = 1; l <= hRunExpectedCountsMap->GetNbinsY(); ++l)
     	    {
       	      Stash::Coordinate pop = hRunExpectedCountsMap->GetBinCoordinate(k,l);	
	      double LambdaMod = pop.GetLambda().GetDegrees();
	      double BetaMod = pop.GetBeta().GetDegrees();
	      
	      if(LambdaMod > 180)
		LambdaMod = LambdaMod-360;
	      
	      int a = hRunExpectedCountsMap->GetXaxis()->FindFixBin(-LambdaMod);
	      if(a == 0 || a > hRunExpectedCountsMap->GetNbinsX())
	      	a = hRunExpectedCountsMap->GetXaxis()->FindFixBin(-(LambdaMod+360));
	      if(a == 0 || a > hRunExpectedCountsMap->GetNbinsX())
		continue;
	      
	      int b = hRunExpectedCountsMap->GetYaxis()->FindFixBin(BetaMod);
	      double offaxis = pop.GetAngularDistance(pointingpos).GetDegrees();
	      double nth = 0;
	      
	      if(offaxis <= fPsiCut)
		nth = hRTemp->Interpolate(offaxis);
	      else
		nth = 0;

	      double valcounts = nth*LiveTime;
	      double valareatime = nth*1e8*LiveTime;

	      //double intflux = Utilities::Flux::IntFlux(1,1.,-1*fSpectralIndex,int_emin,int_emax,0);
	      //std::cout << valareatime * 1./(intflux * LiveTime) << std::endl;
	      if(pop.GetAngularDistance(pointingpos).GetDegrees() < fPsiCut){
		hRunExpectedCountsMap->SetBinContent(a,b,valcounts+hRunExpectedCountsMap->GetBinContent(a,b));
		hRunAreaTimeMap->SetBinContent(a,b,valareatime+hRunAreaTimeMap->GetBinContent(a,b));
	      }
       	    }
      	}

      //Fill final maps, only loop in the useful area
      int kmin, lmin, kmax, lmax;
      hExpectedCountsMap->FindBinPosition(hRunExpectedCountsMap->GetBinCoordinate(1,1),kmin,lmin);
      hExpectedCountsMap->FindBinPosition(hRunExpectedCountsMap->GetBinCoordinate(hRunExpectedCountsMap->GetNbinsX(),hRunExpectedCountsMap->GetNbinsY()),kmax,lmax);
      if(kmax == 0 || kmax > hExpectedCountsMap->GetNbinsX())
	kmax = hExpectedCountsMap->GetNbinsX();
      if(kmin == 0 || kmin > hExpectedCountsMap->GetNbinsX())
	kmin = 1;
      if(lmax == 0 || lmax > hExpectedCountsMap->GetNbinsY())
	lmax = hExpectedCountsMap->GetNbinsY();
      if(lmin == 0 || lmin > hExpectedCountsMap->GetNbinsY())
	lmin = 1;
      
      kmin = kmin -1;
      lmin = lmin -1;
      kmax = kmax +1;
      lmax = lmax +1;
      std::cout << kmin << " " << lmin << " " << kmax << " " << lmax << std::endl;
      for(int k = kmin; k <= kmax; ++k)
	{
       	  for(int l = lmin; l <= lmax; ++l)
     	    {
	      Stash::Coordinate pop = hExpectedCountsMap->GetBinCoordinate(k,l);	
	      
	      double val = hExpectedCountsMap->GetBinContent(k,l);
	      double valn = hRunExpectedCountsMap->GetBinContent(hRunExpectedCountsMap->FindBinPosition(pop));
	      double val_at = hAreaTimeMap->GetBinContent(k,l);
	      double valn_at = hRunAreaTimeMap->GetBinContent(hRunAreaTimeMap->FindBinPosition(pop));
	      
	      hExpectedCountsMap->SetBinContent(k,l,val+valn);
	      hAreaTimeMap->SetBinContent(k,l,val_at+valn_at);	      
	    }
	}
      delete hRunExpectedCountsMap;
      delete hRunAreaTimeMap;
      delete hRTemp;

      ofile << dsrun->GetRunNumber() << " " << pointingpos.GetLambda().GetDegrees() << " " << pointingpos.GetBeta().GetDegrees() << " " << MeanZen << " " << RelativeEfficiency << " " << LiveTime << " " <<SafeThreshold << " " << int_emin << " " << TheoricRateAtCenter*LiveTime*1e8 << std::endl;      
    }

  TCanvas *c = new TCanvas();
  c->Divide(2,1);
  c->cd(1);
  hExpectedCountsMap->Draw("colz");
  c->cd(2);
  hAreaTimeMap->Draw("colz");

  GetExpectedCountsMap() = hAreaTimeMap;
  GetExtendedExpectedCountsMap() = hAreaTimeMap;

  if(fSaveResults){
    //-----------------------------------------------------------------------------
    //Generate File Name according to analysis parameters    
    //and write in file
    std::ostringstream foutfile;
    std::ostringstream EvAndAcc;
    std::ostringstream ConfigTarget;
    std::ostringstream OverSampling;
    std::ostringstream AdaptCutAlpha;
    std::ostringstream RingMethod;
    std::ostringstream EnergyCuts;
     
    OverSampling << "_OS" << fOSRadius;
    if(fEventsAndAcceptanceFromRadial){
      EvAndAcc << "_PsiCut" << fPsiCut;
      if(fApplySafeThreshold){
	EvAndAcc << "_SafeTh" << fSafeThresholdRatioParam;
	if(fSafeThresholdFromAcceptance)
	  EvAndAcc << "_FromAcc";
	else
	  EvAndAcc << "_FromBias";
      }
    }
    else 
      EvAndAcc << "";
    
    if(!fUseConfigTarget)
      ConfigTarget << "_CustomTarget_" << fUserLambda << "_" << fUserBeta;
    else
      ConfigTarget << "";
    
    if(!fUseConfigMapParams)    
      ConfigTarget << "_MapExt_X" << fExtX << "_Y" << fExtY;
    else
      ConfigTarget << "";
            
    if(fEMin != 0 && fEMax == -1)
      EnergyCuts << "_EMin" << fEMin;
    else if(fEMin == 0 && fEMax != -1)
      EnergyCuts << "_EMax" << fEMax;
    else if(fEMin != 0 && fEMax != -1)
      EnergyCuts << "_EMin" << fEMin << "_EMax" << fEMax;
    else
      EnergyCuts << "";
    
    foutfile << OutfilePrefix << "_AreaTimeMap.root";
    
    
    TFile fileResultsOut(foutfile.str().c_str(),"RECREATE");
    //  gROOT->cd();    
    
    fileResultsOut.Append(hAreaTimeMap);    
    fileResultsOut.Write();
    fileResultsOut.Clear();
  }

}

/*
 * Create events and acceptance maps from radial lookups
 *
 */
void SurveySuite::MapMaker::CreateEventsInfosFile(const char *RadialConfig, const char *Resfile, const char *Regionsfile, std::string table_path, const char *OutfilePrefix)
{
  double fSourceSize = 0.;
  if(!GetPointLikeFluxMaps()){
    fSourceSize = 1.;
    std::cout << "Extended tables will be taken!" << std::endl;
  }

  std::cout << RadialConfig << std::endl;
  std::ostringstream f;
  f << "RadialTables/RadialTablesEnergy_" << RadialConfig << ".root";  
  TFile *RadialFile = TFile::Open(f.str().c_str());
  TH2F *hEnergy = (TH2F *)RadialFile->Get("hR2D_4Tels_0");
  double emin = fEMin;
  double emax = fEMax;
  int binemin = 0;
  int binemax = -1;
  if(fEMin != 0){
    double log10emin = TMath::Log10(fEMin);
    binemin = hEnergy->GetXaxis()->FindFixBin(log10emin);
    emin = pow(10,hEnergy->GetXaxis()->GetBinLowEdge(binemin));
  }
  if(fEMax != -1){
    double log10emax = TMath::Log10(fEMax);
    binemax = hEnergy->GetXaxis()->FindFixBin(log10emax);
    emax = pow(10,hEnergy->GetXaxis()->GetBinLowEdge(binemax)+hEnergy->GetXaxis()->GetBinWidth(binemax));
  }
  delete hEnergy;

  //
  TFile *fileResults = TFile::Open(Resfile);
  gROOT->cd();
  Sash::DataSet *results = (Sash::DataSet *)fileResults->Get("results");
  results->GetEntry(0);
  
  Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();
  ParisAnalysis::AcceptanceMap *AccTest =  hess->Handle("", (ParisAnalysis::AcceptanceMap *) 0);
  AccTest->LoadAllMembers();
  ParisAnalysis::AnalysisConfig *Config =  hess->Handle("", (ParisAnalysis::AnalysisConfig *) 0);
  Config->LoadAllMembers();
  
  Sash::DataSet *run_results = (Sash::DataSet *)fileResults->Get("run_results");
  gROOT->cd();
  
  if(fUseConfigTarget){
    std::cout << "Using target position from Results File" << std::endl;
    fUserLambda = Config->GetTargetRa();
    fUserBeta = Config->GetTargetDec();
    std::cout << "fUserLambda = " << fUserLambda << ", fUserBeta = " << fUserBeta << std::endl;
  }

  Stash::Coordinate targetpos(Stash::Lambda(fUserLambda,Stash::Angle::Degrees),
			      Stash::Beta(fUserBeta,Stash::Angle::Degrees),
			      Config->GetSystem());
  
  if(fUseConfigMapParams){
    std::cout << "Using Maps parameters from Results File" << std::endl;
    fBinSize = Config->GetMapBinSize().GetDegrees();
    fExtX  = Config->GetMapExtensionX().GetDegrees();
    fExtY  = Config->GetMapExtensionY().GetDegrees();
    std::cout << "fBinSize = " << fBinSize << ", fExtX = " << fExtX << ", fExtY = " << fExtY << std::endl;
  }

  Double_t ExtX_Run = fPsiCut;
  Double_t ExtY_Run = fPsiCut;

  //LOOP OVER RUNS TO COMPUTE ACCEPTANCE MAP FROM RADIAL LOOKUPS
  int skippedruns = 0;
  double maxExtension = sqrt(pow(fExtX,2)+pow(fExtY,2));
  double RunSelectionCut = 0;
  if(fSelectAllRunsContributingToTheMap)
    RunSelectionCut = fPsiCut+maxExtension;
  else
    RunSelectionCut = fPsiCut;

  //  TH1F *hEnergyDist = new TH1F("hE","hE",10000,0,10); 
  TH1D *hRTemp = new TH1D("hRTemp","hRTemp",100,0,9);    
  TH1D *hRTempFit = new TH1D("hRTempFit","hRTempFit",100,0,9);    
 
  std::ostringstream foutfile_all;
  foutfile_all << OutfilePrefix << "_EventsInfo_AllRuns_.txt";
  std::ofstream ofile_all(foutfile_all.str().c_str());
  
  ofile_all << "RunNumber" << " " << "time" << " " << "Lambda" << " " << "Beta" << " " << "ZenithAngle" << " " << "OffAxisAngle" << " " << "Energy" << std::endl;

  //  Display::SkyHistogram2D *hRunZenMap = 0;
  std::vector<int>::iterator it = fEntriesVector.begin();

  //  for(int i = 0; i < run_results->GetEntries(); ++i)
  for(; it != fEntriesVector.end(); ++it)
    {
      int i = *it;
      run_results->GetEntry(i);

      if(fVerbose)
	std::cout << i+1 << "/" << run_results->GetEntries() << std::endl;
      
      const ParisAnalysis::DataStorageRun* dsrun = hess->Handle("GammaFOVStorage", (ParisAnalysis::DataStorageRun *) 0);
      Stash::Coordinate pointingpos = dsrun->GetObsPos();      
      
      if(pointingpos.GetAngularDistance(targetpos).GetDegrees() > RunSelectionCut){
	//	std::cout << "skipping run" << std::endl;
	std::cout << Utilities::TextStyle::Yellow() <<  "SKIPPING RUN!!!" << Utilities::TextStyle::Reset() <<  std::endl;
	
	++skippedruns;
	continue;
      }
      
      std::ostringstream foutfile;
      foutfile << OutfilePrefix << "_EventsInfo_Run_" << dsrun->GetRunNumber() << "_PsiCut_" << fPsiCut << ".txt";
      std::ofstream ofile(foutfile.str().c_str());

      ofile << "RunNumber" << " " << "time" << " " << "Lambda" << " " << "Beta" << " " << "ZenithAngle" << " " << "OffAxisAngle" << " " << "Energy" << std::endl;
            
      for(std::vector<ParisAnalysis::EventData>::const_iterator evit = dsrun->GetEvents().begin(); evit != dsrun->GetEvents().end() ; evit++) 
	{
	  Stash::Coordinate pop(Stash::Lambda(evit->GetShowerPosX(),Stash::Angle::Degrees),
				Stash::Beta(evit->GetShowerPosY(),Stash::Angle::Degrees),
				Config->GetSystem());
	  double LambdaMod = pop.GetLambda().GetDegrees();
	  double BetaMod = pop.GetBeta().GetDegrees();
	  //std::cout << evit->GetShowerNomPosX() << " " << evit->GetShowerNomPosY() << std::endl;
	  if(LambdaMod > 180)
	    LambdaMod = LambdaMod-360;
	  
	  double time = evit->GetEventTime().GetUTC().GetModifiedJulianDate();
	  
	  ofile << dsrun->GetRunNumber() << " " << std::setprecision(12) << time << " " << LambdaMod << " " << BetaMod << " " << evit->GetZenithAngle() << " " << evit->GetOffAxisAngle() << " " << evit->GetEnergy() << std::endl;

	  ofile_all << dsrun->GetRunNumber() << " " << std::setprecision(12) << time << " " << LambdaMod << " " << BetaMod << " " << evit->GetZenithAngle() << " " << evit->GetOffAxisAngle() << " " << evit->GetEnergy() << std::endl;

	}      
    }
}

/*
 * Create events and acceptance maps from radial lookups
 *
 */
void SurveySuite::MapMaker::CreateEventsAndAcceptanceMaps_FromRadialLookups(const char *RadialConfig, const char *Resfile, const char *Regionsfile, std::string table_path, const char *OutfilePrefix)
{
  double fSourceSize = 0.;
  if(!GetPointLikeFluxMaps()){
    fSourceSize = 1.;
    std::cout << "Extended tables will be taken!" << std::endl;
  }

  std::cout << RadialConfig << std::endl;
  //
  TFile *fileResults = TFile::Open(Resfile);
  gROOT->cd();
  Sash::DataSet *results = (Sash::DataSet *)fileResults->Get("results");
  results->GetEntry(0);
  
  gStyle->SetOptStat(000);
  
  Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();
  ParisAnalysis::AcceptanceMap *AccTest =  hess->Handle("", (ParisAnalysis::AcceptanceMap *) 0);
  AccTest->LoadAllMembers();
  ParisAnalysis::AnalysisConfig *Config =  hess->Handle("", (ParisAnalysis::AnalysisConfig *) 0);
  Config->LoadAllMembers();

  ParisAnalysis::RadialAcceptance *Radial =  hess->Handle("", (ParisAnalysis::RadialAcceptance *) 0);
  Radial->LoadAllMembers();

  const TH1F * hRadialGammaAcceptance = Radial->GethRadialGammaAcceptance();
  const Display::Histogram2D * hRadialGammaAcceptanceZen = Radial->GethRadialGammaAcceptanceZen();
  
  Sash::DataSet *run_results = (Sash::DataSet *)fileResults->Get("run_results");
  gROOT->cd();
  
  if(fUseConfigTarget){
    std::cout << "Using target position from Results File" << std::endl;
    fUserLambda = Config->GetTargetRa();
    fUserBeta = Config->GetTargetDec();
    std::cout << "fUserLambda = " << fUserLambda << ", fUserBeta = " << fUserBeta << std::endl;
  }

  Stash::Coordinate targetpos(Stash::Lambda(fUserLambda,Stash::Angle::Degrees),
			      Stash::Beta(fUserBeta,Stash::Angle::Degrees),
			      Config->GetSystem());
  
  if(fUseConfigMapParams){
    std::cout << "Using Maps parameters from Results File" << std::endl;
    fBinSize = Config->GetMapBinSize().GetDegrees();
    fExtX  = Config->GetMapExtensionX().GetDegrees();
    fExtY  = Config->GetMapExtensionY().GetDegrees();
    std::cout << "fBinSize = " << fBinSize << ", fExtX = " << fExtX << ", fExtY = " << fExtY << std::endl;
  }

  Double_t ExtX_Run = fPsiCut;
  Double_t ExtY_Run = fPsiCut;

  //Global Maps ...
  Display::SkyHistogram2D *hAcceptanceMap = new Display::SkyHistogram2D("SkyAcceptance","SkyAcceptance",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
  Display::SkyHistogram2D *hEventsMap = new Display::SkyHistogram2D("SkyEvents","SkyEvents",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
  Display::SkyHistogram2D *hOffAxisMap = new Display::SkyHistogram2D("OffAxisMap","OffAxisMap",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
  Display::SkyHistogram2D *hZenithAngleMap = new Display::SkyHistogram2D("ZenithAngleMap","ZenithAngleMap",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
  Display::SkyHistogram2D *hMuonEfficiencyMap = new Display::SkyHistogram2D("MuonEfficiencyMap","MuonEfficiencyMap",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
  Display::SkyHistogram2D *hMinSafeThresholdMap = new Display::SkyHistogram2D("MinSafeThresholdMap","MinSafeThresholdMap",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
  Display::SkyHistogram2D *hAveragedSafeThresholdMap = new Display::SkyHistogram2D("AveragedSafeThresholdMap","AveragedSafeThresholdMap",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
  Display::SkyHistogram2D *hAveragedLiveTimeMap = new Display::SkyHistogram2D("AveragedLiveTimeMap","AveragedLiveTimeMap",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));

  //Create Exlusion Regions
  CreateExclusionMask(Resfile, Regionsfile);
  
  //LOOP OVER RUNS TO COMPUTE ACCEPTANCE MAP FROM RADIAL LOOKUPS
  int skippedruns = 0;
  double maxExtension = sqrt(pow(fExtX,2)+pow(fExtY,2));
  double RunSelectionCut = 0;
  if(fSelectAllRunsContributingToTheMap)
    RunSelectionCut = fPsiCut+maxExtension;
  else
    RunSelectionCut = fPsiCut;

  //  TH1F *hEnergyDist = new TH1F("hE","hE",10000,0,10); 
  TH1D *hRTemp = new TH1D("hRTemp","hRTemp",100,0,9);    
  TH1D *hRTempFit = new TH1D("hRTempFit","hRTempFit",100,0,9);    
 
  //  Display::SkyHistogram2D *hRunZenMap = 0;
  std::vector<int>::iterator it = fEntriesVector.begin();

  //  for(int i = 0; i < run_results->GetEntries(); ++i)
  for(; it != fEntriesVector.end(); ++it)
    {
      int i = *it;
      run_results->GetEntry(i);

      if(fVerbose)
	std::cout << i+1 << "/" << run_results->GetEntries() << std::endl;
      
      const ParisAnalysis::DataStorageRun* dsrun = hess->Handle("GammaFOVStorage", (ParisAnalysis::DataStorageRun *) 0);
      Stash::Coordinate pointingpos = dsrun->GetObsPos();      
      
      if(pointingpos.GetAngularDistance(targetpos).GetDegrees() > RunSelectionCut){
	//	std::cout << "skipping run" << std::endl;
	std::cout << Utilities::TextStyle::Yellow() <<  "SKIPPING RUN!!!" << Utilities::TextStyle::Reset() <<  std::endl;
		
	++skippedruns;
	continue;
      }

      Display::SkyHistogram2D *hRunExclusionMap = new Display::SkyHistogram2D("RunExclusion","RunExclusion",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
      
      Display::SkyHistogram2D *hRunAcceptanceMap = new Display::SkyHistogram2D("RunAcceptance","RunAcceptance",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
      
      Display::SkyHistogram2D *hRunAcceptanceMapCut = new Display::SkyHistogram2D("RunAcceptanceCut","RunAcceptanceCut",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
      
      Display::SkyHistogram2D *hRunEventsMap = new Display::SkyHistogram2D("RunEvents","RunEvents",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));

      Display::SkyHistogram2D *hRunEventsMapCut = new Display::SkyHistogram2D("RunEventsCut","RunEventsCut",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));

      Display::SkyHistogram2D *hRunOffAxisMap = new Display::SkyHistogram2D("RunOffAxis","RunOffAxis",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
      Display::SkyHistogram2D *hRunZenithAngleMap = new Display::SkyHistogram2D("RunZenithAngleMap","RunZenithAngleMap",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
      Display::SkyHistogram2D *hRunMuonEfficiencyMap = new Display::SkyHistogram2D("RunMuonEfficiencyMap","RunMuonEfficiencyMap",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
      
      Display::SkyHistogram2D *hRunSafeThresholdMap = new Display::SkyHistogram2D("RunSafeThresholdMap","RunSafeThresholdMap",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
      Display::SkyHistogram2D *hRunWeightedSafeThresholdMap = new Display::SkyHistogram2D("RunWeightedSafeThresholdMap","RunWeightedSafeThresholdMap",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));


      //define pointing conditions parameters
      if(fVerbose)
	std::cout << dsrun->GetMeanZenith() << " " << (dsrun->GetTelsInRun()).GetPattern() << std::endl;

      double MeanZen = dsrun->GetMeanZenith();
      double meanCosZen = cos(MeanZen*TMath::Pi()*1./180);
      double MuonEff = dsrun->GetMuonEfficiency();
      double offaxisAngle = fPsiCut;
      if(!fSafeThresholdFromPsiCut)
	offaxisAngle = fSafeThresholdFromPsiValue;
      double MeanAzimuth = dsrun->GetMeanAzimuth();
      Int_t NTels = dsrun->GetTelsInRun().size();
      Double_t SafeThreshold = 0;
      
      //if(UseZenMaps)
      TF1 *fRTemp = 0; 
      for(int i = 1; i <= 100; ++i){
	hRTemp->SetBinContent(i,hRadialGammaAcceptanceZen->GetBinContent(hRadialGammaAcceptanceZen->GetXaxis()->FindFixBin(hRTemp->GetBinCenter(i)), hRadialGammaAcceptanceZen->GetYaxis()->FindFixBin(meanCosZen)));
      }

      fRTemp = new TF1("fRTemp","pol8",0,9);
      hRTemp->Fit(fRTemp,"QN");
      
      for(int i = 1; i <= 100; ++i){
	hRTempFit->SetBinContent(i,fRTemp->Eval(hRTemp->GetBinCenter(i)));      
      }
      hRTempFit->Scale(1./hRTempFit->GetMaximum());
      
      delete fRTemp;
      
      if(fVerbose)
	std::cout << hRTemp->GetMaximum() << std::endl;
      double lastl = 0;
      double lastb = 0;
      int countdoubles = 0;
      int countdoubles2 = 0;
      
      for(std::vector<ParisAnalysis::EventData>::const_iterator evit = dsrun->GetEvents().begin(); evit != dsrun->GetEvents().end() ; evit++) 
	{
	  if(evit->GetShowerPosX() != lastl){
	    if(evit->GetShowerPosY() != lastb){
	      Stash::Coordinate pop(Stash::Lambda(evit->GetShowerPosX(),Stash::Angle::Degrees),
				    Stash::Beta(evit->GetShowerPosY(),Stash::Angle::Degrees),
				    Config->GetSystem());
	      double LambdaMod = pop.GetLambda().GetDegrees();
	      double BetaMod = pop.GetBeta().GetDegrees();
	      //std::cout << evit->GetShowerNomPosX() << " " << evit->GetShowerNomPosY() << std::endl;
	      if(LambdaMod > 180)
		LambdaMod = LambdaMod-360;
	      int a = hRunEventsMap->GetXaxis()->FindFixBin(-LambdaMod);
	      if(a == 0 || a > hRunEventsMap->GetNbinsX())
		a = hRunEventsMap->GetXaxis()->FindFixBin(-(LambdaMod+360));
	      if(a == 0 || a > hRunEventsMap->GetNbinsX())
		continue;
	      int aa = hEventsMap->GetXaxis()->FindFixBin(-LambdaMod);
	      if(aa == 0 || aa > hEventsMap->GetNbinsX())
		aa = hEventsMap->GetXaxis()->FindFixBin(-(LambdaMod+360));
	      if(aa == 0 || aa > hEventsMap->GetNbinsX())
		continue;
	      
	      int b = hRunEventsMap->GetYaxis()->FindFixBin(BetaMod);
	      int bb = hEventsMap->GetYaxis()->FindFixBin(BetaMod);
	      
	      if(pop.GetAngularDistance(pointingpos).GetDegrees() < fPsiCut){		
		hRunEventsMap->SetBinContent(a,b,1+hRunEventsMap->GetBinContent(a,b));
		hEventsMap->SetBinContent(aa,bb,1+hEventsMap->GetBinContent(aa,bb));
		
		if(GetExclusionMap()->GetBinContent(aa,bb) == 0)
		  hRunEventsMapCut->SetBinContent(a,b,1+hRunEventsMapCut->GetBinContent(a,b));
	      }
	      ++countdoubles;
	    }
	    lastl = evit->GetShowerPosX();
	    lastb = evit->GetShowerPosY();
	    ++countdoubles2;
	  }
	}      
      
      for(int k = 1; k <= hRunAcceptanceMap->GetNbinsX(); ++k)
	{
       	  for(int l = 1; l <= hRunAcceptanceMap->GetNbinsY(); ++l)
     	    {
       	      Stash::Coordinate pop = hRunAcceptanceMap->GetBinCoordinate(k,l);	
	      double LambdaMod = pop.GetLambda().GetDegrees();
	      double BetaMod = pop.GetBeta().GetDegrees();

	      if(LambdaMod > 180)
		LambdaMod = LambdaMod-360;
	      int a = GetExclusionMap()->GetXaxis()->FindFixBin(-LambdaMod);
	      if(a == 0 || a > GetExclusionMap()->GetNbinsX())
	      	a = GetExclusionMap()->GetXaxis()->FindFixBin(-(LambdaMod+360));
	      if(a == 0 || a > GetExclusionMap()->GetNbinsX())
		continue;
	      int aa = hRunAcceptanceMap->GetXaxis()->FindFixBin(-LambdaMod);
	      if(aa == 0 || aa > hRunAcceptanceMap->GetNbinsX())
	      	aa = hRunAcceptanceMap->GetXaxis()->FindFixBin(-(LambdaMod+360));
	      if(aa == 0 || aa > hRunAcceptanceMap->GetNbinsX())
		continue;

	      int b = GetExclusionMap()->GetYaxis()->FindFixBin(BetaMod);
	      int bb = hRunAcceptanceMap->GetYaxis()->FindFixBin(BetaMod);

	      double val = 0;
	      if(!fFitAcceptance)
		val = hRTemp->Interpolate(pow(pop.GetAngularDistance(pointingpos).GetDegrees(),2));
	      else
		val = hRTempFit->Interpolate(pow(pop.GetAngularDistance(pointingpos).GetDegrees(),2));

	      if(pop.GetAngularDistance(pointingpos).GetDegrees() < fPsiCut){
		hRunAcceptanceMap->SetBinContent(aa,bb,val+hRunAcceptanceMap->GetBinContent(aa,bb));
		if(GetExclusionMap()->GetBinContent(a,b) == 0)
		  hRunAcceptanceMapCut->SetBinContent(aa,bb,val+hRunAcceptanceMapCut->GetBinContent(aa,bb));
		
	      }
       	    }
      	}
      for(int k = 1; k <= hRunExclusionMap->GetNbinsX(); ++k)
	{
       	  for(int l = 1; l <= hRunExclusionMap->GetNbinsY(); ++l)
     	    {
       	      Stash::Coordinate pop = hRunExclusionMap->GetBinCoordinate(k,l);
	      double LambdaMod = pop.GetLambda().GetDegrees();
	      double BetaMod = pop.GetBeta().GetDegrees();
	      if(LambdaMod > 180)
		LambdaMod = LambdaMod-360;
	      int a = GetExclusionMap()->GetXaxis()->FindFixBin(-LambdaMod);
	      if(a == 0 || a > GetExclusionMap()->GetNbinsX())
	      	a = GetExclusionMap()->GetXaxis()->FindFixBin(-(LambdaMod+360));
	      if(a == 0 || a > GetExclusionMap()->GetNbinsX())
		continue;

	      int aa = hRunExclusionMap->GetXaxis()->FindFixBin(-LambdaMod);
	      if(aa == 0 || aa > hRunExclusionMap->GetNbinsX())
	      	aa = hRunExclusionMap->GetXaxis()->FindFixBin(-(LambdaMod+360));
	      if(aa == 0 || aa > hRunExclusionMap->GetNbinsX())
		continue;

	      int b = GetExclusionMap()->GetYaxis()->FindFixBin(BetaMod);
	      int bb = hRunExclusionMap->GetYaxis()->FindFixBin(BetaMod);
	      
	      if(pop.GetAngularDistance(pointingpos).GetDegrees() < fPsiCut){	
		if(GetExclusionMap()->GetBinContent(a,b) == 0)
		  hRunExclusionMap->SetBinContent(aa,bb,1);
		else
		  hRunExclusionMap->SetBinContent(aa,bb,0); 
	      }
	    }
	}
      if(fVerbose)
	std::cout << "Allowed Integral : events = " << hRunEventsMapCut->Integral() << ", Acceptance = " << hRunAcceptanceMapCut->Integral() << std::endl;
      if(TMath::IsNaN(hRunAcceptanceMapCut->Integral()) || hRunAcceptanceMapCut->Integral() == 0)
	{
	  std::cout << Utilities::TextStyle::Yellow() <<  "SKIPPING RUN (NAN in ACC)!!!" << Utilities::TextStyle::Reset() <<  std::endl;

	  delete hRunExclusionMap;
	  delete hRunAcceptanceMap;
	  delete hRunAcceptanceMapCut;
	  delete hRunEventsMap;
	  delete hRunEventsMapCut;
	  delete hRunOffAxisMap;
	  delete hRunZenithAngleMap;
	  delete hRunMuonEfficiencyMap;
	  delete hRunSafeThresholdMap;
	  delete hRunWeightedSafeThresholdMap;
	  hRTemp->Reset();
	  hRTempFit->Reset();
	  ++skippedruns;
	  continue;
	}
      double scalefactor = hRunEventsMapCut->Integral()*1./hRunAcceptanceMapCut->Integral();
      if(fVerbose){
	std::cout << "ScaleFactor = " << scalefactor << std::endl;
	std::cout << "Max Acceptance = " << hRunAcceptanceMap->GetMaximum() << std::endl;      
      }

      hRunAcceptanceMap->Scale(scalefactor);

      if(fVerbose)
	std::cout << "Scaled Max Acceptance = " << hRunAcceptanceMap->GetMaximum() << std::endl;      
      //Filling scale factor per run map
      fRadialIntegrationEbin[dsrun->GetRunNumber()]=std::make_pair(0,-1);      //std::make_pair(binemin,binemax);      
      fScaleFactor[dsrun->GetRunNumber()]=scalefactor;

      //FILLING POINTING CONDITIONS MAP
      for(int k = 1; k <= hRunAcceptanceMap->GetNbinsX(); ++k)
	{
       	  for(int l = 1; l <= hRunAcceptanceMap->GetNbinsY(); ++l)
     	    {
       	      Stash::Coordinate pop = hRunAcceptanceMap->GetBinCoordinate(k,l);	
	      double LambdaMod = pop.GetLambda().GetDegrees();
	      double BetaMod = pop.GetBeta().GetDegrees();

	      if(LambdaMod > 180)
		LambdaMod = LambdaMod-360;
	      int a = hRunOffAxisMap->GetXaxis()->FindFixBin(-LambdaMod);
	      if(a == 0 || a > hRunOffAxisMap->GetNbinsX())
	      	a = hRunOffAxisMap->GetXaxis()->FindFixBin(-(LambdaMod+360));
	      if(a == 0 || a > hRunOffAxisMap->GetNbinsX())
		continue;
	      int aa = hRunAcceptanceMap->GetXaxis()->FindFixBin(-LambdaMod);
	      if(aa == 0 || aa > hRunAcceptanceMap->GetNbinsX())
	      	aa = hRunAcceptanceMap->GetXaxis()->FindFixBin(-(LambdaMod+360));
	      if(aa == 0 || aa > hRunAcceptanceMap->GetNbinsX())
		continue;

	      int b = hRunOffAxisMap->GetYaxis()->FindFixBin(BetaMod);
	      int bb = hRunAcceptanceMap->GetYaxis()->FindFixBin(BetaMod);

	      if(pop.GetAngularDistance(pointingpos).GetDegrees() < fPsiCut){
		hRunOffAxisMap->SetBinContent(a,b,pop.GetAngularDistance(pointingpos).GetDegrees()*hRunAcceptanceMap->GetBinContent(aa,bb)+hRunOffAxisMap->GetBinContent(a,b));
		hRunZenithAngleMap->SetBinContent(a,b,MeanZen*hRunAcceptanceMap->GetBinContent(aa,bb)+hRunZenithAngleMap->GetBinContent(a,b));
		hRunMuonEfficiencyMap->SetBinContent(a,b,MuonEff*hRunAcceptanceMap->GetBinContent(aa,bb)+hRunMuonEfficiencyMap->GetBinContent(a,b));
		hRunSafeThresholdMap->SetBinContent(a,b,SafeThreshold+hRunSafeThresholdMap->GetBinContent(a,b));
		hRunWeightedSafeThresholdMap->SetBinContent(a,b,SafeThreshold*hRunAcceptanceMap->GetBinContent(aa,bb)+hRunWeightedSafeThresholdMap->GetBinContent(a,b));
	      }
       	    }
      	}
      

      //GO TO NOMINAL SYSTEM AND CORRECT FOR ZENITH GRADIENTS...
      if(fCorrectZenith){
	Display::SkyHistogram2D *hRunFOVAcceptanceMap = new Display::SkyHistogram2D("RunFOVAcceptance","RunFOVAcceptance",Stash::Coordinate(0,0,hess->GetNominalSystem()),Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
	Display::SkyHistogram2D *hRunFOVEventsMap = new Display::SkyHistogram2D("RunFOVEvents","RunFOVEvents",Stash::Coordinate(0,0,hess->GetNominalSystem()),Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
	Display::SkyHistogram2D *hRunFOVExclusionMap = new Display::SkyHistogram2D("RunFOVExclusion","RunFOVExclusion",Stash::Coordinate(0,0,hess->GetNominalSystem()),Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
;
 
	TH1F *hX = new TH1F("hx","hx",hRunFOVEventsMap->GetNbinsX(),hRunFOVEventsMap->GetXaxis()->GetBinLowEdge(1),hRunFOVEventsMap->GetXaxis()->GetBinLowEdge(1+hRunFOVEventsMap->GetNbinsX()));
	TH1F *hY = new TH1F("hy","hy",hRunFOVEventsMap->GetNbinsY(),hRunFOVEventsMap->GetYaxis()->GetBinLowEdge(1),hRunFOVEventsMap->GetYaxis()->GetBinLowEdge(1+hRunFOVEventsMap->GetNbinsY()));

	UInt_t nStep = 30;
	Float_t timeStep = (dsrun->GetLastEventTime() - dsrun->GetFirstEventTime()) / nStep;
	Float_t timeOffset = timeStep/2;
	
	for(UInt_t iStep = 0; iStep < nStep; ++iStep)
	  {
	    hess->SetCurrentTime(dsrun->GetFirstEventTime() + iStep * timeStep + timeOffset);    

	    Sash::NominalPointing *nomptg = hess->Handle<Sash::NominalPointing>();
	    nomptg->SetReferenceDirection(pointingpos);
	    for(Int_t fov_xbin = 1; fov_xbin <= hRunFOVAcceptanceMap->GetNbinsX(); ++fov_xbin)
	      {
		for(Int_t fov_ybin = 1; fov_ybin <= hRunFOVAcceptanceMap->GetNbinsY(); ++fov_ybin)
		  {
		    Int_t fov_bin = hRunFOVAcceptanceMap->GetBin(fov_xbin,fov_ybin);
		    // Compute bin position in nominal system
		    Stash::Coordinate nompos = hRunFOVAcceptanceMap->GetBinCoordinate(fov_xbin, fov_ybin);

		    Crash::SetRefractionType(*hess->GetHorizonSystem(),Crash::HorizonSystem::None);
		    Stash::Coordinate pos = nompos.GetCoordinate(*Config->GetSystem());
		    Int_t sky_bin = hRunAcceptanceMap->FindBinPosition(pos);
		    Int_t xb,yb,zb;
		    hRunAcceptanceMap->GetBinXYZ(sky_bin,xb,yb,zb);
		    Float_t area = hRunAcceptanceMap->GetBinSolidAngle(xb,yb);
		    hRunFOVAcceptanceMap->SetBinContent(fov_bin, hRunFOVAcceptanceMap->GetBinContent(fov_bin)+hRunAcceptanceMap->GetBinContent(sky_bin)*1./area);
		    hRunFOVEventsMap->SetBinContent(fov_bin, hRunFOVEventsMap->GetBinContent(fov_bin)+hRunEventsMap->GetBinContent(sky_bin)*1./area);
		    hRunFOVExclusionMap->SetBinContent(fov_bin, hRunFOVExclusionMap->GetBinContent(fov_bin)+hRunExclusionMap->GetBinContent(sky_bin));
		  }
	      }
	  }
      
	for(int ii = 10; ii <= hRunFOVEventsMap->GetNbinsY()-10; ++ii)
	  {
	    double sumb = 0;
	    double sumn = 0; 
	    double suma = 0; 
	    for(int jj = 10; jj <= hRunFOVEventsMap->GetNbinsX()-10; ++jj)
	      {
		if(hRunFOVExclusionMap->GetBinContent(jj,ii)>0){
		  sumb += hRunFOVExclusionMap->GetBinContent(jj,ii);
		  sumn += hRunFOVEventsMap->GetBinContent(jj,ii)*hRunFOVExclusionMap->GetBinContent(jj,ii)*1./30;
		  suma += hRunFOVAcceptanceMap->GetBinContent(jj,ii)*hRunFOVExclusionMap->GetBinContent(jj,ii)*1./30;
		}
	      }
	    if(sumb > 0)
	      hY->SetBinContent(ii,sumn*1./suma);
	  }
	
	for(int ii = 10; ii <= hRunFOVEventsMap->GetNbinsX()-10; ++ii)
	  {
	    double sumb = 0;
	    double sumn = 0; 
	    double suma = 0; 
	    for(int jj = 10; jj <= hRunFOVEventsMap->GetNbinsY()-10; ++jj)
	      {
		if(hRunFOVExclusionMap->GetBinContent(ii,jj)>0){
		  sumb += hRunFOVExclusionMap->GetBinContent(ii,jj);
		  sumn += hRunFOVEventsMap->GetBinContent(ii,jj)*hRunFOVExclusionMap->GetBinContent(ii,jj)*1./30;
		  suma += hRunFOVAcceptanceMap->GetBinContent(ii,jj)*hRunFOVExclusionMap->GetBinContent(ii,jj)*1./30;
		}
	      }
	    if(sumb > 0)
	      hX->SetBinContent(ii,sumn*1./suma);
	  }

	double sux = 0,ssx = 0;
	double xmin = -1.*(ExtX_Run-10*fBinSize);
	double xmax = ExtX_Run-10*fBinSize;
	double maxg = 0.5;
	double GradientX = 0;	
	double GradientY = 0;

	Int_t n = hX->GetNbinsX();
	for(Int_t i = 1; i <= n; ++i) {
	  //std::cout << hX->GetBinContent(i) << std::endl;
	  Double_t y1,y2,dy1,dy2;
	  Double_t x = hX->GetBinCenter(i);
	  UInt_t j = hX->GetXaxis()->FindFixBin(-x);
	  if(x < xmin || x > xmax) continue;
	  y1 = hX->GetBinContent(i);
	  y2 = hX->GetBinContent(j);
	  dy1 = hX->GetBinError(i);
	  dy2 = hX->GetBinError(j);
	  if(y1 > 0 &&  y2 > 0 && dy1 > 0 && dy2 > 0) {
	    double U = (y1-y2)/(y1 + y2);
	    double dU2 =  4 *  y1 * y2 / std::pow(y1 + y2,3);
	    sux += U * x / dU2;
	    ssx += x*x / dU2;
	  }
	}
	if(ssx > 0 && fabs(sux/ssx) < maxg) {
	  GradientX = sux/ssx;
	}
	else if(ssx > 0 && sux/ssx >  maxg) {
	  GradientX = maxg;
	}
	if(fVerbose)
	  std::cout << "Gradient X =  " << GradientX << std::endl;
	
	double suy = 0, ssy = 0;
	n = hY->GetNbinsX();
	for(Int_t i = 1; i <= n; ++i) {
	  Double_t y1,y2,dy1,dy2;
	  Double_t x = hY->GetBinCenter(i);
	  UInt_t j = hY->GetXaxis()->FindFixBin(-x);
	  if(x < xmin || x > xmax) continue;
	  y1 = hY->GetBinContent(i);
	  y2 = hY->GetBinContent(j);
	  dy1 = hY->GetBinError(i);
	  dy2 = hY->GetBinError(j);
	  if(y1 > 0 &&  y2 > 0 && dy1 > 0 && dy2 > 0) {
	    double U = (y1-y2)/(y1 + y2);
	    double dU2 =  4 *  y1 * y2 / std::pow(y1 + y2,3);
	    suy += U * x / dU2;
	    ssy += x*x / dU2;
	  }
	}
	if(ssy > 0 && fabs(suy/ssy) < maxg) {
	  GradientY = suy/ssy;
	}
	else if(ssy > 0 && suy/ssy >  maxg) {
	  GradientY = maxg;
	}
	if(fVerbose)
	  std::cout << "Gradient Y =  " << GradientY << std::endl;

	for(int ii = 1; ii <= hRunFOVAcceptanceMap->GetNbinsX(); ++ii)
	  {
	    for(int jj = 1; jj <= hRunFOVAcceptanceMap->GetNbinsY(); ++jj)
	      {
		Float_t x = hRunFOVAcceptanceMap->GetXaxis()->GetBinCenter(ii);
		Float_t y = hRunFOVAcceptanceMap->GetYaxis()->GetBinCenter(jj);
		//		std::cout << x << " " << y << std::endl;
		hRunFOVAcceptanceMap->SetBinContent(ii,jj,hRunFOVAcceptanceMap->GetBinContent(ii,jj)*(1+GradientX*x)*(1+GradientY*y));
	      }
	  }

	//REPROJECT ON THE SKY
	for(UInt_t iStep = 0; iStep < nStep; ++iStep)
	  {
	    hess->SetCurrentTime(dsrun->GetFirstEventTime() + iStep * timeStep + timeOffset);    
	    Sash::NominalPointing *nomptg = hess->Handle<Sash::NominalPointing>();
	    nomptg->SetReferenceDirection(pointingpos);
	    for(Int_t sky_xbin = 1; sky_xbin <= hRunAcceptanceMap->GetNbinsX(); ++sky_xbin)
	      {
		for(Int_t sky_ybin = 1; sky_ybin <= hRunAcceptanceMap->GetNbinsY(); ++sky_ybin)
		  {
		    Int_t sky_bin = hRunAcceptanceMap->GetBin(sky_xbin,sky_ybin);
		    Float_t area = hRunAcceptanceMap->GetBinSolidAngle(sky_xbin,sky_ybin);
		    // Compute bin position in nominal system
		    Stash::Coordinate pos = hRunAcceptanceMap->GetBinCoordinate(sky_xbin, sky_ybin);
		    Crash::SetRefractionType(*hess->GetHorizonSystem(),Crash::HorizonSystem::None);
		    Stash::Coordinate nompos = pos.GetCoordinate(*hess->GetNominalSystem());
		    Int_t fov_bin = hRunFOVAcceptanceMap->FindBinPosition(nompos);
		    hRunAcceptanceMap->SetBinContent(sky_bin,hRunAcceptanceMap->GetBinContent(sky_bin)+hRunFOVAcceptanceMap->GetBinContent(fov_bin)*area);
		  }
	      }
	  }

	Display::SkyHistogram2D *hRunAcceptanceMapCut2 = new Display::SkyHistogram2D("RunAcceptanceCut2","RunAcceptanceCut2",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));

	for(int k = 1; k <= hRunAcceptanceMap->GetNbinsX(); ++k)
	  {
	    for(int l = 1; l <= hRunAcceptanceMap->GetNbinsY(); ++l)
	      {
		Stash::Coordinate pop = hRunAcceptanceMap->GetBinCoordinate(k,l);	
		if(GetExclusionMap()->GetBinContent(GetExclusionMap()->FindBinPosition(pop)) == 0)
		  hRunAcceptanceMapCut2->SetBinContent(k,l,hRunAcceptanceMap->GetBinContent(k,l));	      
	      }
	  }
	if(fVerbose)
	  std::cout << "Allowed Integral : events = " << hRunEventsMapCut->Integral() << ", Acceptance = " << hRunAcceptanceMapCut2->Integral() << std::endl;
	hRunAcceptanceMap->Scale(hRunEventsMapCut->Integral()*1./hRunAcceptanceMapCut2->Integral());
	fScaleFactor[dsrun->GetRunNumber()]=hRunEventsMapCut->Integral()*1./hRunAcceptanceMapCut2->Integral();

	delete hRunAcceptanceMapCut2;
	delete hRunFOVAcceptanceMap;
	delete hRunFOVExclusionMap;
	delete hRunFOVEventsMap;
	delete hX;
	delete hY;
      }
      //-----------------------
      
      //only loop in the useful area
      int kmin, lmin, kmax, lmax;
      hAcceptanceMap->FindBinPosition(hRunAcceptanceMap->GetBinCoordinate(1,1),kmin,lmin);
      hAcceptanceMap->FindBinPosition(hRunAcceptanceMap->GetBinCoordinate(hRunAcceptanceMap->GetNbinsX(),hRunAcceptanceMap->GetNbinsY()),kmax,lmax);
      if(kmax == 0 || kmax > hAcceptanceMap->GetNbinsX())
	kmax = hAcceptanceMap->GetNbinsX();
      if(kmin == 0 || kmin > hAcceptanceMap->GetNbinsX())
	kmin = 1;
      if(lmax == 0 || lmax > hAcceptanceMap->GetNbinsY())
	lmax = hAcceptanceMap->GetNbinsY();
      if(lmin == 0 || lmin > hAcceptanceMap->GetNbinsY())
	lmin = 1;

      kmin = kmin -1;
      lmin = lmin -1;
      kmax = kmax +1;
      lmax = lmax +1;
      if(fVerbose)
	std::cout << kmin << " " << lmin << " " << kmax << " " << lmax << std::endl;

      double runOnTime = dsrun->GetOnLiveTime() *1./3600;

      for(int k = kmin; k <= kmax; ++k)
	{
       	  for(int l = lmin; l <= lmax; ++l)
     	    {
	      Stash::Coordinate pop = hAcceptanceMap->GetBinCoordinate(k,l);	
	      double val = hAcceptanceMap->GetBinContent(k,l);
	      double valn = hRunAcceptanceMap->GetBinContent(hRunAcceptanceMap->FindBinPosition(pop));
	      hAcceptanceMap->SetBinContent(k,l,val+valn);	      
	      double valoff = hOffAxisMap->GetBinContent(k,l);
	      double valnoff = hRunOffAxisMap->GetBinContent(hRunOffAxisMap->FindBinPosition(pop));
	      hOffAxisMap->SetBinContent(k,l,valoff+valnoff);
	      double valzen = hZenithAngleMap->GetBinContent(k,l);
	      double valnzen = hRunZenithAngleMap->GetBinContent(hRunZenithAngleMap->FindBinPosition(pop));
	      hZenithAngleMap->SetBinContent(k,l,valzen+valnzen);
	      double valeff = hMuonEfficiencyMap->GetBinContent(k,l);
	      double valneff = hRunMuonEfficiencyMap->GetBinContent(hRunMuonEfficiencyMap->FindBinPosition(pop));
	      hMuonEfficiencyMap->SetBinContent(k,l,valeff+valneff);
	     
	      double valth = hMinSafeThresholdMap->GetBinContent(k,l);
	      double valnth = hRunSafeThresholdMap->GetBinContent(hRunSafeThresholdMap->FindBinPosition(pop));
	      if(valth == 0)
		hMinSafeThresholdMap->SetBinContent(k,l,valnth);
	      else if(valnth != 0)
		hMinSafeThresholdMap->SetBinContent(k,l,TMath::Min(valth, valnth));

	      double valwth = hAveragedSafeThresholdMap->GetBinContent(k,l);
	      double valnwth = hRunWeightedSafeThresholdMap->GetBinContent(hRunWeightedSafeThresholdMap->FindBinPosition(pop));
	      hAveragedSafeThresholdMap->SetBinContent(k,l,valwth+valnwth);

	      double valtime = hAveragedLiveTimeMap->GetBinContent(k,l);
	      double valntime = runOnTime;
	      if(valn>0)
		hAveragedLiveTimeMap->SetBinContent(k,l,valtime+valntime);

	    }
	}

      delete hRunExclusionMap;
      delete hRunAcceptanceMap;
      delete hRunAcceptanceMapCut;
      delete hRunEventsMapCut;
      hRTemp->Reset();
      hRTempFit->Reset();
      delete hRunEventsMap;	
      delete hRunOffAxisMap;
      delete hRunZenithAngleMap;
      delete hRunMuonEfficiencyMap;
      delete hRunSafeThresholdMap;
      delete hRunWeightedSafeThresholdMap;

      std::cout << "Everything should be deleted" << std::endl;
    }
  
  double maxaccfinal = hAcceptanceMap->GetMaximum();

  for(int k = 1; k <= hOffAxisMap->GetNbinsX(); ++k)
    {
      for(int l = 1; l <= hOffAxisMap->GetNbinsY(); ++l)
	{
	  double ratio = hAcceptanceMap->GetBinContent(k,l)*1./maxaccfinal;
	  if(hAcceptanceMap->GetBinContent(k,l) != 0){
	    hOffAxisMap->SetBinContent(k,l,hOffAxisMap->GetBinContent(k,l)*1./hAcceptanceMap->GetBinContent(k,l));
	    hZenithAngleMap->SetBinContent(k,l,hZenithAngleMap->GetBinContent(k,l)*1./hAcceptanceMap->GetBinContent(k,l));
	    hMuonEfficiencyMap->SetBinContent(k,l,hMuonEfficiencyMap->GetBinContent(k,l)*1./hAcceptanceMap->GetBinContent(k,l));
	    hAveragedSafeThresholdMap->SetBinContent(k,l,hAveragedSafeThresholdMap->GetBinContent(k,l)*1./hAcceptanceMap->GetBinContent(k,l));
	    hAveragedLiveTimeMap->SetBinContent(k,l,hAveragedLiveTimeMap->GetBinContent(k,l)*ratio);
	  }
	}
    }
  
  TCanvas *c = new TCanvas("MeanObsCond","MeanObsCond",1000,800);
  c->Divide(3,1);
  c->cd(1);
  hOffAxisMap->Draw("colz");
  c->cd(2);
  hZenithAngleMap->Draw("colz");
  c->cd(3);
  hMuonEfficiencyMap->Draw("colz");
  
  // std::cout << fEMin << " - " << fEMax << " -> " << emin << "(bin " << binemin << ")" << " - " << emax << "(bin " << binemax << ")" << std::endl;
  std::cout << "-> " << fEntriesVector.size()-skippedruns << " runs out of " << fEntriesVector.size() << " pre-selected have been used (" << run_results->GetEntries() << " available in the result file)." << std::endl;
  GetAcceptanceMap() = hAcceptanceMap;
  GetEventsMap() = hEventsMap;
  GetOffAxisMap() = hOffAxisMap;
  GetZenithAngleMap() = hZenithAngleMap;
  GetMuonEfficiencyMap() = hMuonEfficiencyMap;
  GetMinSafeThresholdMap() = hMinSafeThresholdMap;
  GetAveragedSafeThresholdMap() = hAveragedSafeThresholdMap;
  GetAveragedLiveTimeMap() = hAveragedLiveTimeMap;

  if(fSaveResults){
    //-----------------------------------------------------------------------------
    //Generate File Name according to analysis parameters    
    //and write in file
    std::ostringstream foutfile;
    std::ostringstream canvfile;
    std::ostringstream EvAndAcc;
    std::ostringstream ConfigTarget;
    std::ostringstream OverSampling;
    std::ostringstream AdaptCutAlpha;
    std::ostringstream RingMethod;
    std::ostringstream EnergyCuts;
     
    OverSampling << "_OS" << fOSRadius;
    if(fEventsAndAcceptanceFromRadial){
      EvAndAcc << "_PsiCut" << fPsiCut;
      if(fCorrectZenith)
	EvAndAcc << "_ZenCorr";
      if(fApplySafeThreshold){
	EvAndAcc << "_SafeTh" << fSafeThresholdRatioParam;
	if(fSafeThresholdFromAcceptance)
	  EvAndAcc << "_FromAcc";
	else
	  EvAndAcc << "_FromBias";
      }
    }
    else 
      EvAndAcc << "";
    
    if(!fUseConfigTarget)
      ConfigTarget << "_CustomTarget_" << fUserLambda << "_" << fUserBeta;
    else
      ConfigTarget << "";
    
    if(!fUseConfigMapParams)    
      ConfigTarget << "_MapExt_X" << fExtX << "_Y" << fExtY;
    else
      ConfigTarget << "";
            
    if(fEMin != 0 && fEMax == -1)
      EnergyCuts << "_EMin" << fEMin;
    else if(fEMin == 0 && fEMax != -1)
      EnergyCuts << "_EMax" << fEMax;
    else if(fEMin != 0 && fEMax != -1)
      EnergyCuts << "_EMin" << fEMin << "_EMax" << fEMax;
    else
      EnergyCuts << "";
    
    foutfile << OutfilePrefix << "_RadialMaps.root";
    canvfile << OutfilePrefix << "_MeanObsCond.png";
    c->Print(canvfile.str().c_str());

    
    TFile fileResultsOut(foutfile.str().c_str(),"RECREATE");
    //  gROOT->cd();    
    fileResultsOut.Append(hAcceptanceMap);
    fileResultsOut.Append(hEventsMap);
    fileResultsOut.Append(hOffAxisMap);
    fileResultsOut.Append(hZenithAngleMap);
    fileResultsOut.Append(hMuonEfficiencyMap);
    fileResultsOut.Append(hMinSafeThresholdMap);
    fileResultsOut.Append(hAveragedSafeThresholdMap);
    fileResultsOut.Append(hAveragedLiveTimeMap);

    fileResultsOut.Write();
    fileResultsOut.Clear();
  }

}

void SurveySuite::MapMaker::CreatePSFMaps(const char *RadialConfig, const char *Resfile, const char *Regionsfile, std::string table_path, const char *OutfilePrefix)
{
  //First check that the fScaleFactor map has been filled : 
  if(fScaleFactor.size()==0){
    std::cout << "Error : the fScaleFactor map hasn't been filled -> you must run the radial acceptance generation first" << std::endl;
    return;
  }
  else if(fRadialIntegrationEbin.size()==0){
   std::cout << "Error : the fRadialIntegrationEbin map hasn't been filled -> you must run the radial acceptance generation first" << std::endl;
    return;
  } 
  else if(fRadialIntegrationEbin.size() != fScaleFactor.size()){
   std::cout << "Error : the fRadialIntegrationEbin map doesn't have the same number of entries as the fScaleFactor map" << std::endl;
   std::cout << "-> " << fRadialIntegrationEbin.size() << " " << fScaleFactor.size() << std::endl;
    return;
  } 

  std::cout << fRadialIntegrationEbin.size() << " " << fScaleFactor.size() << std::endl;

  //Load radial tables
  std::cout << RadialConfig << std::endl;
  std::ostringstream f;
  f << "RadialTables/RadialTablesEnergy_" << RadialConfig << ".root";  
  TFile *RadialFile = TFile::Open(f.str().c_str());
  
  //Opening the results file
  TFile *fileResults = TFile::Open(Resfile);
  gROOT->cd();
  Sash::DataSet *results = (Sash::DataSet *)fileResults->Get("results");
  results->GetEntry(0);
  
  Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();
  ParisAnalysis::AnalysisConfig *Config =  hess->Handle("", (ParisAnalysis::AnalysisConfig *) 0);
  Config->LoadAllMembers();
  
  Sash::DataSet *run_results = (Sash::DataSet *)fileResults->Get("run_results");
  gROOT->cd();
  
  if(fUseConfigTarget){
    std::cout << "Using target position from Results File" << std::endl;
    fUserLambda = Config->GetTargetRa();
    fUserBeta = Config->GetTargetDec();
    std::cout << "fUserLambda = " << fUserLambda << ", fUserBeta = " << fUserBeta << std::endl;
  }

  Stash::Coordinate targetpos(Stash::Lambda(fUserLambda,Stash::Angle::Degrees),
			      Stash::Beta(fUserBeta,Stash::Angle::Degrees),
			      Config->GetSystem());
  
  if(fUseConfigMapParams){
    std::cout << "Using Maps parameters from Results File" << std::endl;
    fBinSize = Config->GetMapBinSize().GetDegrees();
    fExtX  = Config->GetMapExtensionX().GetDegrees();
    fExtY  = Config->GetMapExtensionY().GetDegrees();
    std::cout << "fBinSize = " << fBinSize << ", fExtX = " << fExtX << ", fExtY = " << fExtY << std::endl;
  }

  Double_t ExtX_Run = fPsiCut;
  Double_t ExtY_Run = fPsiCut;

  //Global Maps ...
  Display::SkyHistogram2D *hPSFAmplitudeMap = new Display::SkyHistogram2D("PSFAmplitudeMap","PSFAmplitudeMap",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));

  Display::SkyHistogram2D *hPSFSigma1Map = new Display::SkyHistogram2D("PSFSigma1Map","PSFSigma1Map",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));

  Display::SkyHistogram2D *hPSFAmplitude2Map = new Display::SkyHistogram2D("PSFAmplitude2Map","PSFAmplitude2Map",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));

  Display::SkyHistogram2D *hPSFSigma2Map = new Display::SkyHistogram2D("PSFSigma2Map","PSFSigma2Map",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));

  Display::SkyHistogram2D *hPSFAmplitude3Map = new Display::SkyHistogram2D("PSFAmplitude3Map","PSFAmplitude3Map",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
  Display::SkyHistogram2D *hPSFSigma3Map = new Display::SkyHistogram2D("PSFSigma3Map","PSFSigma3Map",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));

  //Map of PSF histos
  TH3F *PSFHistoMap = new TH3F("PSFHistoMap","PSFHistoMap",hPSFAmplitudeMap->GetNbinsX(),0,hPSFAmplitudeMap->GetNbinsX(),hPSFAmplitudeMap->GetNbinsY(),0,hPSFAmplitudeMap->GetNbinsY(),50,0,0.05);

  //Retrieve morphology tables
  char *s = gSystem->ExpandPathName(table_path.c_str());
  setenv("SPECTRUM_PATH",s,1);
  delete[] s;

  std::string fAcceptanceNameBase = "Combined";
  Bool_t fAutoAzimut=false;
  Float_t fMuonEfficiency = 0.;
  Bool_t fCheckAnalysisVersion=false;
  Bool_t fVerboseTables=false;

  std::map<Int_t,std::string> fMapAzimuth;
  fMapAzimuth[0]   = "North";
  fMapAzimuth[180] = "South";
  
  std::map<Int_t, Morphology::MorphologyTableFinderMulti*> fMapTableFinder;
  std::map<Int_t, std::map<Int_t, const Morphology::AngularResolutionTableEfficiency*> > fMapTableResolution;
  
  for (std::map<Int_t,std::string>::const_iterator it = fMapAzimuth.begin(); it!=fMapAzimuth.end(); ++it) {
    std::string AcceptanceName = fAcceptanceNameBase + it->second;
    Float_t Azimuth = Float_t(it->first);
    
    Morphology::MorphologyTableFinderMulti *table = new Morphology::MorphologyTableFinderMulti(Morphology::MorphologyTableFinderMulti::Models,AcceptanceName.c_str(),Azimuth,fAutoAzimut,fCheckAnalysisVersion,fVerboseTables);
    table->Process(Sash::Folder::GetFolder("results"));
    fMapTableFinder[Azimuth] = table;
    std::cout << Utilities::TextStyle::Blue() << "MorphologyTableFinderMulti[" << Azimuth << "] = Done !" << Utilities::TextStyle::Reset() << std::endl;    
    std::cout << (void*)table << std::endl;

    for (int itel=3;itel<=4;++itel) {
      std::ostringstream oss_tabname;
      oss_tabname << AcceptanceName << "_Tels" << itel;
      fMapTableResolution[Azimuth][itel] = hess->Get<Morphology::AngularResolutionTableEfficiency>(oss_tabname.str().c_str());
      std::cout << Utilities::TextStyle::Blue() << "Retrieve AngularResolutionTableEfficiency[" << Azimuth << "][" << itel << "] Named : " << oss_tabname.str().c_str() << " (" << (void*)fMapTableResolution[Azimuth][itel] << ")" << Utilities::TextStyle::Reset() << std::endl;
      oss_tabname.str("");
    }
  }

  //LOOP OVER RUNS TO COMPUTE PSF
  int skippedruns = 0;
  double maxExtension = sqrt(pow(fExtX,2)+pow(fExtY,2));
  double RunSelectionCut = 0;
  if(fSelectAllRunsContributingToTheMap)
    RunSelectionCut = fPsiCut+maxExtension;
  else
    RunSelectionCut = fPsiCut;

  TH1D *hRTemp = new TH1D("hRTemp","hRTemp",100,0,9);    
  TH1D *hRTempFit = new TH1D("hRTempFit","hRTempFit",100,0,9);    
  int nbins_psf = int(fPsiCut*1./fBinSize);
  TH2F *hRunPSFvsOffAxis = new TH2F("PSFvsOffAxis","PSFvsOffAxis",nbins_psf,0,fPsiCut,50,0,0.05);

  Display::SkyHistogram2D *hRunAccMap;
  Display::SkyHistogram2D *hRunPSFAmplitudeMap;

  std::vector<int>::iterator it = fEntriesVector.begin();

  //  for(int i = 0; i < run_results->GetEntries(); ++i)
  for(; it != fEntriesVector.end(); ++it)
    {
      int i = *it;
      run_results->GetEntry(i);
      if(fVerbose)
	std::cout << i+1 << "/" << run_results->GetEntries() << std::endl;
      
      const ParisAnalysis::DataStorageRun* dsrun = hess->Handle("GammaFOVStorage", (ParisAnalysis::DataStorageRun *) 0);
      Stash::Coordinate pointingpos = dsrun->GetObsPos();      

      if(fScaleFactor[dsrun->GetRunNumber()] == 0)
	continue;
      
      if(pointingpos.GetAngularDistance(targetpos).GetDegrees() > RunSelectionCut){
	std::cout << Utilities::TextStyle::Yellow() <<  "SKIPPING RUN!!!" << Utilities::TextStyle::Reset() <<  std::endl;
	++skippedruns;
	continue;
      }

      ///SURVEY XCHECK
      /*if(dsrun->GetRunNumber() > 75227){
	++skippedruns;
	continue;
	}*/

      std::ostringstream namerunamp;
      std::ostringstream namerunacc;
      namerunamp << "amp" << dsrun->GetRunNumber();
      namerunacc << "acc" << dsrun->GetRunNumber();

      hRunPSFAmplitudeMap = new Display::SkyHistogram2D(namerunamp.str().c_str(),namerunamp.str().c_str(),pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));

      hRunAccMap = new Display::SkyHistogram2D(namerunacc.str().c_str(),namerunacc.str().c_str(),pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));

            
      //define pointing conditions parameters
      if(fVerbose)
	std::cout << dsrun->GetMeanZenith() << " " << (dsrun->GetTelsInRun()).GetPattern() << std::endl;

      double MeanZen = dsrun->GetMeanZenith();
      double meanCosZen = cos(MeanZen*TMath::Pi()*1./180);
      double MuonEff = dsrun->GetMuonEfficiency();
      double offaxisAngle = fPsiCut;
      if(!fSafeThresholdFromPsiCut)
	offaxisAngle = fSafeThresholdFromPsiValue;
      double MeanAzimuth = dsrun->GetMeanAzimuth();
      Int_t NTels = dsrun->GetTelsInRun().size();
      
      if (MeanAzimuth<0) {MeanAzimuth+=360.;}
      if (MeanAzimuth>360) {MeanAzimuth=MeanAzimuth-360.;}
      Int_t AzimuthCode = ( (MeanAzimuth<=90. || MeanAzimuth>270) ? 0 : 180 );
      
      std::cout << "MeanAzimuth = " << MeanAzimuth << " AzimuthCode = " << AzimuthCode << std::endl;
      //Find Table
      const Morphology::AngularResolutionTableEfficiency *ResTab = fMapTableResolution[AzimuthCode][NTels];
      if (!ResTab) {
	std::cout << Utilities::TextStyle::Yellow() << " Can't find AngularResolutionTableEfficiency for AzimuthCode : " << AzimuthCode << " NTels : " << NTels << Utilities::TextStyle::Reset() << std::endl;
	continue;
      }

      std::cout << MeanZen << " --- " << MuonEff << std::endl;
      std::ostringstream fzenrad3;
      std::ostringstream fzenrad4;
      fzenrad3 << "hR2D_3Tels_" << int(dsrun->GetMeanZenith()*1./10)*10;
      fzenrad4 << "hR2D_4Tels_" << int(dsrun->GetMeanZenith()*1./10)*10;
      TH2F *h3Tels = (TH2F *)RadialFile->Get(fzenrad3.str().c_str());
      TH2F *h4Tels = (TH2F *)RadialFile->Get(fzenrad4.str().c_str());
      
      if((dsrun->GetTelsInRun()).GetPattern() == 30)
	{
	  hRTemp = h4Tels->ProjectionY("p4",fRadialIntegrationEbin[dsrun->GetRunNumber()].first,fRadialIntegrationEbin[dsrun->GetRunNumber()].second);
	}
      else
	{
	  hRTemp = h3Tels->ProjectionY("p3",fRadialIntegrationEbin[dsrun->GetRunNumber()].first,fRadialIntegrationEbin[dsrun->GetRunNumber()].second);
	}
      delete h3Tels;
      delete h4Tels;

      TF1 *fRTemp = 0;
      if(!fFitAcceptance)      
	hRTemp->Scale(1./hRTemp->GetMaximum());
      else{
	fRTemp = new TF1("fRTemp","pol6",0,9);
	hRTemp->Fit(fRTemp,"QN");
	
	for(int u = 1; u <= 100; ++u){
	  hRTempFit->SetBinContent(u,fRTemp->Eval(hRTemp->GetBinCenter(u)));      
	}
	hRTempFit->Scale(1./hRTempFit->GetMaximum());
      }
      delete fRTemp;

      if(fVerbose)
	std::cout << hRTemp->GetMaximum() << std::endl;
      
      for(int k = 1; k <= hRunAccMap->GetNbinsX(); ++k)
	{
       	  for(int l = 1; l <= hRunAccMap->GetNbinsY(); ++l)
     	    {
       	      Stash::Coordinate pop = hRunAccMap->GetBinCoordinate(k,l);	
	      double LambdaMod = pop.GetLambda().GetDegrees();
	      double BetaMod = pop.GetBeta().GetDegrees();

	      if(LambdaMod > 180)
		LambdaMod = LambdaMod-360;

	      int aa = hRunAccMap->GetXaxis()->FindFixBin(-LambdaMod);
	      if(aa == 0 || aa > hRunAccMap->GetNbinsX())
	      	aa = hRunAccMap->GetXaxis()->FindFixBin(-(LambdaMod+360));
	      if(aa == 0 || aa > hRunAccMap->GetNbinsX())
		continue;

	      int bb = hRunAccMap->GetYaxis()->FindFixBin(BetaMod);
	      double val = 0;
	      if(!fFitAcceptance)
		val = hRTemp->Interpolate(pow(pop.GetAngularDistance(pointingpos).GetDegrees(),2));
	      else
		val = hRTempFit->Interpolate(pow(pop.GetAngularDistance(pointingpos).GetDegrees(),2));

	      if(pop.GetAngularDistance(pointingpos).GetDegrees() < fPsiCut){
		hRunAccMap->SetBinContent(aa,bb,val+hRunAccMap->GetBinContent(aa,bb));
	      }
       	    }
      	}
      hRunAccMap->Scale(fScaleFactor[dsrun->GetRunNumber()]);
      if(fVerbose)
	std::cout << "Scaled Max Acceptance = " << hRunAccMap->GetMaximum() << std::endl;      

      TAxis *xaxis = hRunPSFvsOffAxis->GetXaxis();
      TAxis *yaxis = hRunPSFvsOffAxis->GetYaxis();
      for(int ii = 1; ii <= hRunPSFvsOffAxis->GetNbinsX(); ++ii){
	for(int jj = 1; jj <= hRunPSFvsOffAxis->GetNbinsY(); ++jj){
	  double off = xaxis->GetBinCenter(ii);
	  double th2 = yaxis->GetBinCenter(jj);
	  double max = ResTab->GetPDFvalue(yaxis->GetBinCenter(1),meanCosZen,off,fSpectralIndex,MuonEff);

	  double RadialValue = 0;
	  if(!fFitAcceptance)
	    RadialValue= fScaleFactor[dsrun->GetRunNumber()]*hRTemp->Interpolate(pow(off,2));
	  else
	    RadialValue= fScaleFactor[dsrun->GetRunNumber()]*hRTempFit->Interpolate(pow(off,2));
	  // if(jj ==1)
	  // std::cout << "Rad = " << RadialValue << ", off = " << off << ", hrt = " << hRTemp->Interpolate(pow(off,2)) << std::endl;
	  hRunPSFvsOffAxis->SetBinContent(ii,jj,ResTab->GetPDFvalue(th2,meanCosZen,off,fSpectralIndex,MuonEff)*RadialValue*1./max);
	}
      }
            
      //only loop in the useful area
      int kmin, lmin, kmax, lmax;
      hPSFAmplitudeMap->FindBinPosition(hRunPSFAmplitudeMap->GetBinCoordinate(1,1),kmin,lmin);
      hPSFAmplitudeMap->FindBinPosition(hRunPSFAmplitudeMap->GetBinCoordinate(hRunPSFAmplitudeMap->GetNbinsX(),hRunPSFAmplitudeMap->GetNbinsY()),kmax,lmax);
      if(kmax == 0 || kmax > hPSFAmplitudeMap->GetNbinsX())
	kmax = hPSFAmplitudeMap->GetNbinsX();
      if(kmin == 0 || kmin > hPSFAmplitudeMap->GetNbinsX())
	kmin = 1;
      if(lmax == 0 || lmax > hPSFAmplitudeMap->GetNbinsY())
	lmax = hPSFAmplitudeMap->GetNbinsY();
      if(lmin == 0 || lmin > hPSFAmplitudeMap->GetNbinsY())
	lmin = 1;
      
      if(kmin > 1)
	kmin = kmin -1;
      if(lmin > 1)
	lmin = lmin -1;
      if(kmax < hPSFAmplitudeMap->GetNbinsX())
	kmax = kmax +1;
      if(lmax < hPSFAmplitudeMap->GetNbinsY())
	lmax = lmax +1;
      if(fVerbose)
	std::cout << kmin << " " << lmin << " " << kmax << " " << lmax << std::endl;
      
      //      TH1D *hPSFTemp = new TH1D();//"PSFTemp","PSFTemp",50,0,0.05);
      for(int k = kmin; k <= kmax; ++k)
	{
       	  for(int l = lmin; l <= lmax; ++l)
     	    {
	      Stash::Coordinate pop = hPSFAmplitudeMap->GetBinCoordinate(k,l);	
	      double PsiDist = pop.GetAngularDistance(pointingpos).GetDegrees();
	      TH1D *hPSFTemp = new TH1D();
	      int psfbin = hRunPSFvsOffAxis->GetXaxis()->FindFixBin(PsiDist);
	      hPSFTemp = hRunPSFvsOffAxis->ProjectionY("pPSF",psfbin,psfbin);
	      for(int m = 1; m<=50; ++m){
		PSFHistoMap->AddBinContent(PSFHistoMap->GetBin(k,l,m),hPSFTemp->GetBinContent(m));
	      }
	      
	      delete hPSFTemp;
	    }
	}

      hRunPSFvsOffAxis->Reset();
      hRunAccMap->Reset();
      hRunPSFAmplitudeMap->Reset();
      hRTemp->Reset();
      hRTempFit->Reset();
      std::cout << "Everything should be deleted" << std::endl;
    }
  delete hRTemp;
  delete hRTempFit;

  std::map<std::pair<int,int>,TH1D*> PSFHistoMapProjected; 
  std::cout << "Creating PSFHistoMap..." << std::endl;
  for(int i = 1; i <= PSFHistoMap->GetNbinsX(); ++i){
    for(int j = 1; j <= PSFHistoMap->GetNbinsY(); ++j){
      std::ostringstream fname_ij;
      fname_ij << "PSF_" << i << "_" << j;
      TH1D *hPSF = PSFHistoMap->ProjectionZ(fname_ij.str().c_str(),i,i,j,j);//new TH1D(fname_ij.str().c_str(),fname_ij.str().c_str(),50,0,0.05);
      PSFHistoMapProjected[std::make_pair(i,j)]=hPSF;
    }
  }
  std::cout << "...done. " << std::endl;


  //Normalize histos
  for(int i = 1; i <= hPSFAmplitudeMap->GetNbinsX(); ++i){
    for(int j = 1; j <= hPSFAmplitudeMap->GetNbinsY(); ++j){
      //std::cout << "Normalizing : " << i << " " << j << std::endl;
      if(PSFHistoMapProjected[std::make_pair(i,j)]->GetMaximum() != 0){
	PSFHistoMapProjected[std::make_pair(i,j)]->Scale(1./PSFHistoMapProjected[std::make_pair(i,j)]->GetMaximum());
      }
    }
  }

  //Fit Function
  TF1 *func = new TF1("TRIPLE","[0]*(exp(-x/(2*[1]*[1]))+[2]*exp(-x/(2*[3]*[3]))+[4]*exp(-x/(2*[5]*[5])))",0.000,0.05);
  double m = PSFHistoMapProjected[std::make_pair(1,1)]->GetMaximum();
  func->SetParameters(m,0.01,0.1,0.1,0.1,0.2);
  func->SetParLimits(0,0.001,10000);
  func->SetParLimits(1,0.0,1.0);
  func->SetParLimits(2,0.001,10000);
  func->SetParLimits(3,0.0,1.0);
  func->SetParLimits(4,0.001,10000);
  func->SetParLimits(5,0.0,1.0);
  func->SetParNames("scale","#sigma_{1}","A_{2}","#sigma_{2}","A_{3}","#sigma_{3}");
  for(int i = 1; i <= hPSFAmplitudeMap->GetNbinsX(); ++i){
    for(int j = 1; j <= hPSFAmplitudeMap->GetNbinsY(); ++j){
      std::cout << "Fitting : " << i << " " << j << std::endl;
      if(PSFHistoMapProjected[std::make_pair(i,j)]->GetMaximum() != 0){
	PSFHistoMapProjected[std::make_pair(i,j)]->Fit(func,"R0Q");
	hPSFAmplitudeMap->SetBinContent(i,j,func->GetParameter(0));
	hPSFSigma1Map->SetBinContent(i,j,func->GetParameter(1));
	hPSFAmplitude2Map->SetBinContent(i,j,func->GetParameter(2));
	hPSFSigma2Map->SetBinContent(i,j,func->GetParameter(3));
	hPSFAmplitude3Map->SetBinContent(i,j,func->GetParameter(4));
	hPSFSigma3Map->SetBinContent(i,j,func->GetParameter(5));
      }
    }
  }
  
  TCanvas *cPSF = new TCanvas("cPSF","cPSF",1000,800);
  cPSF->Divide(3,2);
  cPSF->cd(1);
  hPSFAmplitudeMap->Draw("colz");
  cPSF->cd(2);
  hPSFSigma1Map->Draw("colz");
  cPSF->cd(3);
  hPSFAmplitude2Map->Draw("colz");
  cPSF->cd(4);
  hPSFSigma2Map->Draw("colz");
  cPSF->cd(5);
  hPSFAmplitude3Map->Draw("colz");
  cPSF->cd(6);
  hPSFSigma3Map->Draw("colz");

  if(fSaveResults){
    std::ostringstream foutfile;
    foutfile << OutfilePrefix << "_PSFParams.root";

    TFile fileResultsOut(foutfile.str().c_str(),"RECREATE");
    fileResultsOut.Append(hPSFAmplitudeMap);
    fileResultsOut.Append(hPSFSigma1Map);
    fileResultsOut.Append(hPSFAmplitude2Map);
    fileResultsOut.Append(hPSFSigma2Map);
    fileResultsOut.Append(hPSFAmplitude3Map);
    fileResultsOut.Append(hPSFSigma3Map);
    fileResultsOut.Append(PSFHistoMap);    
    fileResultsOut.Write();
    fileResultsOut.Clear();
  }
}

void SurveySuite::MapMaker::CreateRingBgMaps(const char *AnalysisConfig, const char *Resfile, const char *Evtfile, const char *Accfile, const char *ExposureMapFile, const char *ExtendedExposureMapFile, const char *Regionsfile, std::string table_path, const char *OutfilePrefix)
{
  TFile *fileResults = TFile::Open(Resfile);
  gROOT->cd();
  Sash::DataSet *results = (Sash::DataSet *)fileResults->Get("results");
  results->GetEntry(0);
  Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();
  
  ParisAnalysis::AnalysisConfig *Config =  hess->Handle("", (ParisAnalysis::AnalysisConfig *) 0);
  Config->LoadAllMembers();
  
  if(fUseConfigTarget){
    std::cout << "Using target position from Results File" << std::endl;
    fUserLambda = Config->GetTargetRa();
    fUserBeta = Config->GetTargetDec();
    std::cout << "fUserLambda = " << fUserLambda << ", fUserBeta = " << fUserBeta << std::endl;
  }
  Stash::Coordinate targetpos(Stash::Lambda(fUserLambda,Stash::Angle::Degrees),
			      Stash::Beta(fUserBeta,Stash::Angle::Degrees),
			      Config->GetSystem());
  if(fUseConfigMapParams){
    std::cout << "Using Maps parameters from Results File" << std::endl;
    fBinSize = Config->GetMapBinSize().GetDegrees();
    fExtX  = Config->GetMapExtensionX().GetDegrees();
    fExtY  = Config->GetMapExtensionY().GetDegrees();
    std::cout << "fBinSize = " << fBinSize << ", fExtX = " << fExtX << ", fExtY = " << fExtY << std::endl;
  }

  //Compute excess and significance maps
  Double_t SafeExtX = fExtX + 2.0;
  Double_t SafeExtY = fExtY + 2.0;
  //  Double_t BinSize  = binsize;
  std::cout << fBinSize << std::endl;
  Double_t nSafeBins = 2.0 *1./fBinSize;

  //Events and Acceptance?
  CreateEventsAndAcceptanceMaps(AnalysisConfig, Resfile, Regionsfile, Evtfile, Accfile, table_path, OutfilePrefix);
  //Exclusion regions map
  CreateExclusionMask(Resfile, Regionsfile);
  
  Display::SkyHistogram2D * hGammaMod = (Display::SkyHistogram2D *)GetEventsMap()->Clone();    
  hGammaMod->SetFoV(SafeExtX,SafeExtY);
  
  Display::SkyHistogram2D * hGammaAccMod = (Display::SkyHistogram2D *)hGammaMod->Clone();
  hGammaAccMod->Reset();
  hGammaMod->Reset();
  
  Display::SkyHistogram2D * hExcMod = (Display::SkyHistogram2D *)hGammaMod->Clone();
  hExcMod->Reset();
  
  Display::SkyHistogram2D * hMask = (Display::SkyHistogram2D*)GetEventsMap()->Clone();
  for (int i=1;i<=hMask->GetNbinsX();++i) {
    for (int j=1;j<=hMask->GetNbinsY();++j) {
      hMask->SetBinContent(i,j,1.);
    }
  }
  
  Display::SkyHistogram2D * hMaskMod = (Display::SkyHistogram2D*)GetEventsMap()->Clone();
  hMaskMod->SetFoV(SafeExtX,SafeExtY);
  
  for (int i=1;i<=hMaskMod->GetNbinsX();++i) {
    for (int j=1;j<=hMaskMod->GetNbinsY();++j) {
      if(i > nSafeBins && j > nSafeBins && i <= hMaskMod->GetNbinsX()-nSafeBins && j <= hMaskMod->GetNbinsY()-nSafeBins){
	hMaskMod->SetBinContent(i,j,1.);
	hExcMod->SetBinContent(i,j,0.);
      }
      else{
	hMaskMod->SetBinContent(i,j,0.);
	hExcMod->SetBinContent(i,j,1.);
      }
    }
  }

  int nx0 = GetEventsMap()->GetNbinsX();
  int ny0 = GetEventsMap()->GetNbinsY();
  int nx = hGammaMod->GetNbinsX();
  int ny = hGammaMod->GetNbinsY();    
  int nxa = hGammaAccMod->GetNbinsX();
  int nya = hGammaAccMod->GetNbinsY();
  
  std::cout << nx << " " << ny << std::endl;
  std::cout << nxa << " " << nya << std::endl;
  
  Display::SkyHistogram2D *hGammaModExc = new Display::SkyHistogram2D("GammaModExc","GammaModExc",*hExcMod); 
  Display::SkyHistogram2D *hGammaAccModExc = new Display::SkyHistogram2D("GammaAccModExc","GammaAccModExc",*hGammaAccMod); 
  Display::SkyHistogram2D *hOS = new Display::SkyHistogram2D("OS","OS",*hExcMod);     
  Display::SkyHistogram2D *hRing = new Display::SkyHistogram2D("Ring","Ring",*hExcMod);     
  
  for(int i = 1; i <= nx; ++i)
    {
      for(int j = 1; j <= ny; ++j)
	{
	  hGammaMod->SetBinContent(i,j,0);
	  hGammaAccMod->SetBinContent(i,j,0);
	}
    }
  
  for(int i = 1; i <= nx0; ++i)
    {
      for(int j = 1; j <= ny0; ++j)
	{
	  Stash::Coordinate pos =  GetEventsMap()->GetBinCoordinate(i,j);
	  Int_t bin = hGammaMod->FindBinPosition(pos);
	  hGammaMod->SetBinContent(bin,GetEventsMap()->GetBinContent(i,j));
	  hExcMod->SetBinContent(bin,GetExclusionMap()->GetBinContent(i,j));
	}
    }
  for(int i = 1; i <= nxa; ++i)
    {
      for(int j = 1; j <= nya; ++j)
	{
	  Stash::Coordinate pos =  hGammaAccMod->GetBinCoordinate(i,j);
	  Int_t bin = GetAcceptanceMap()->FindBinPosition(pos);
	  hGammaAccMod->SetBinContent(i,j,GetAcceptanceMap()->GetBinContent(bin));
	}
    }
    
  for(int i = 1; i <= nx; ++i)
    {
      for(int j = 1; j <= ny; ++j)
	{
	  hGammaModExc->SetBinContent(i,j,hGammaMod->GetBinContent(i,j)*(1-hExcMod->GetBinContent(i,j)));
	  hGammaAccModExc->SetBinContent(i,j,hGammaAccMod->GetBinContent(i,j)*(1-hExcMod->GetBinContent(i,j)));
	}
    }

  hGammaMod->Multiply(hMaskMod);
  hGammaAccMod->Multiply(hMaskMod);
  hGammaModExc->Multiply(hMaskMod);
  hGammaAccModExc->Multiply(hMaskMod);
    
  //OS Kernel Generation
  for(int i = 1; i <= nx; ++i)
    {
      for(int j = 1; j <= ny; ++j)
	{
	  double xcenter = 0.5*nx+1;
	  double ycenter = 0.5*ny+1;
	  double r = TMath::Sqrt(TMath::Power(i-xcenter,2.)+TMath::Power(j-ycenter,2.));
	  if(r<=fOSRadius/fBinSize) {
	    hOS->SetBinContent(i,j,1.);
	  }
	  else {
	    hOS->SetBinContent(i,j,0.);
	  }
	}
    }

  //  Display::SkyHistogram2D *hExcess = new Display::SkyHistogram2D("Excess","Excess",*hExcMod); 
  Display::SkyHistogram2D *hExcessUncorrelated = new Display::SkyHistogram2D("ExcessUncorrelated","ExcessUncorrelated",*hExcMod); 
  Display::SkyHistogram2D *hExcessUncorrelated_Exc = new Display::SkyHistogram2D("ExcessUncorrelated_Exc","ExcessUncorrelated_Exc",*hExcMod); 
  //Display::SkyHistogram2D *hdExcess = new Display::SkyHistogram2D("dExcess","dExcess",*hExcMod);   
  Display::SkyHistogram2D *hRingON = new Display::SkyHistogram2D("RingON","RingON",*hExcMod); 
  Display::SkyHistogram2D *hRingON_Exc = new Display::SkyHistogram2D("RingON_Exc","RingON_Exc",*hExcMod); 
  Display::SkyHistogram2D *hRingONUncorrelated = new Display::SkyHistogram2D("RingONUncorrelated","RingONUncorrelated",*hExcMod); 
  Display::SkyHistogram2D *hRingONUncorrelated_Exc = new Display::SkyHistogram2D("RingONUncorrelated_Exc","RingONUncorrelated_Exc",*hExcMod); 

  Display::SkyHistogram2D *hRingOFF = new Display::SkyHistogram2D("RingOFF","RingOFF",*hExcMod); 
  Display::SkyHistogram2D *hRingOFFUncorrelated = new Display::SkyHistogram2D("RingOFFUncorrelated","RingOFFUncorrelated",*hExcMod); 
  Display::SkyHistogram2D *hRingAccON = new Display::SkyHistogram2D("RingAccON","RingAccON",*hExcMod); 
  Display::SkyHistogram2D *hRingAccON_Exc = new Display::SkyHistogram2D("RingAccON_Exc","RingAccON_Exc",*hExcMod); 
  Display::SkyHistogram2D *hRingAccOFF = new Display::SkyHistogram2D("RingAccOFF","RingAccOFF",*hExcMod);

  Display::SkyHistogram2D *hNRing = new Display::SkyHistogram2D("NR","NR",*hExcMod);     
  Display::SkyHistogram2D *hNRing2 = new Display::SkyHistogram2D("NR2","NR2",*hExcMod); 
  Display::SkyHistogram2D *hAlphaOFFUncorrelated  = new Display::SkyHistogram2D("AlphaOFFUncorrelated","AlphaOFFUncorrelated",*hExcMod);
  Display::SkyHistogram2D *hAlphaOFFUncorrelated_Exc  = new Display::SkyHistogram2D("Alphaoffuncorrelated_Exc","Alphaoffuncorrelated_Exc",*hExcMod);

  Display::SkyHistogram2D *hResOFF = new Display::SkyHistogram2D("ResOFF","ResOFF",*hGammaAccMod);
  Display::SkyHistogram2D *hResExcess = new Display::SkyHistogram2D("ResExcess","ResExcess",*hGammaAccMod);
  Display::SkyHistogram2D *hResExcess_Exc = new Display::SkyHistogram2D("ResExcess_Exc","ResExcess_Exc",*hGammaAccMod);
  //Display::SkyHistogram2D *hResAlphaOFF = new Display::SkyHistogram2D("ResAlphaOFF","ResAlphaOFF",*hGammaAccMod); 
  //Display::SkyHistogram2D *hResAlphaOFF_Exc = new Display::SkyHistogram2D("ResAlphaOFF_Exc","ResAlphaOFF_Exc",*hGammaAccMod); 
  Display::SkyHistogram2D *hResConvArea = new Display::SkyHistogram2D("ResConvArea","ResConvArea",*hGammaAccMod);
  Display::SkyHistogram2D *hResOFFExp = new Display::SkyHistogram2D("ResOFFExp","ResOFFExp",*hGammaAccMod);


  if(true)
    {
      std::cout <<  Utilities::TextStyle::Blue() << "*** Using adaptive ring method and FFT ***" << Utilities::TextStyle::Reset() << std::endl;
      //RING
      Float_t InnerRingStep = fRingStep;
      // Source and rings radii
      Int_t nSteps = 0;
      if(fAdaptFFT){
	if(fConstantArea)
	  nSteps = TMath::FloorNint((fInnerRingMax - fRingMinimalRadius)*1./fRingStep);
	else if(fConstantThickness){
	  nSteps = TMath::FloorNint(((fOuterRingMax-fRingParam_Thickness) - fRingMinimalRadius)*1./fRingStep);
	}
	else// if(fConstantInnerRadius)
	  nSteps = TMath::FloorNint((fOuterRingMax - (fRingMinimalRadius+fRingParam_Thickness))*1./fRingStep);
      }
      std::cout << "Number of Steps = " << nSteps << std::endl;
      int nrings = nSteps;
      
      std::map<int,int> binsoff;
      std::map<int,int> binsoff2;
      Display::SkyHistogram2D *hExcFrac = new Display::SkyHistogram2D("ExcFrac","ExcFrac",*hExcMod); 
      
      for(int i = 1; i <= nx; ++i)
	{
	  for(int j = 1; j <= ny; ++j)
	    {
	      if(fAdaptCut_Alpha)
		hExcFrac->SetBinContent(i,j,fRingParam_AlphaMax);
	      else
		hExcFrac->SetBinContent(i,j,1);
	      int bin = hExcFrac->GetBin(i,j);
	      binsoff[bin] = 0;
	    }
	}
	        
      for(int k = 0; k <= nrings; ++k)
	{ 	
	  //Display::SkyHistogram2D *hResRing = new Display::SkyHistogram2D("ResRing","ResRing",*hExcMod); 
	  Display::SkyHistogram2D *hResRingExp = new Display::SkyHistogram2D("ResRingExp","ResRingExp",*hGammaAccMod);
	  Display::SkyHistogram2D *hResRingAll = new Display::SkyHistogram2D("ResRingAll","ResRingAll",*hGammaAccMod); 
	  //Display::SkyHistogram2D *hRes = new Display::SkyHistogram2D("Res","Res",*hExcMod); 
	  Display::SkyHistogram2D *hResExp = new Display::SkyHistogram2D("ResExp","ResExp",*hGammaAccMod); 
	  std::cout << "InRing " << k << std::endl;      
	    
	  Float_t InnerRing = 0;
	  Float_t OuterRing = 0;
	  if(fAdaptFFT){
	    if(fConstantArea){
	      InnerRing = fRingMinimalRadius + k * fRingStep;
	      OuterRing = sqrt(fRingParam_AreaOverPi + std::pow(InnerRing,2));
	    }
	    else if(fConstantThickness){
	      InnerRing = fRingMinimalRadius + k * fRingStep;
	      OuterRing = InnerRing + fRingParam_Thickness;
	    }
	    else{// if(fConstantInnerRadius){
	      InnerRing = fRingMinimalRadius;
	      OuterRing = InnerRing + fRingParam_Thickness + k * fRingStep;
	    }
	  }
	  else{
	    InnerRing = fStandardRingRadius-fStandardRingThickness;
	    OuterRing = fStandardRingRadius+fStandardRingThickness;
	  }
	  Float_t RingRadius = 0.5*(InnerRing + OuterRing);
	  Float_t RingThickness = RingRadius - InnerRing;
	  Float_t r1 = RingRadius - RingThickness;
	  Float_t r2 = RingRadius + RingThickness;
	  std::cout << "Step :  " << k << "/" << nSteps << std::endl;
	  std::cout << "InnerRing = " << InnerRing << " - " << "OuterRing = " << OuterRing << std::endl ; 
	    
	  //Ring Kernel Generation
	  for(int i = 1; i <= nx; ++i)
	    {
	      for(int j = 1; j <= ny; ++j)
		{
		  double xcenter = 0.5*nx+1;
		  double ycenter = 0.5*ny+1;
		  double r = TMath::Sqrt(TMath::Power(i-xcenter,2.)+TMath::Power(j-ycenter,2.));
		  if(r <= r2/fBinSize && r >= r1/fBinSize){
		    hRing->SetBinContent(i,j,1);//6*TMath::Landau(dist,r1+0.05,0.1,0));
		  }
		}
	    }
	    
	  Utilities::TFftConv *convRExp = new Utilities::TFftConv(*hGammaAccModExc,*hRing);
	  Utilities::TFftConv *convRingAll = new Utilities::TFftConv(*hGammaAccMod,*hRing);
	  Utilities::TFftConv *convOExp = new Utilities::TFftConv(*hGammaAccMod,*hOS);
	    
	  convOExp->convolve(hResExp);
	  convRExp->convolve(hResRingExp);
	  convRingAll->convolve(hResRingAll);
	    
	  if(fAdaptCut_Alpha){
	    for(int i = 1; i <= nx; ++i)
	      {
		for(int j = 1; j <= ny; ++j)
		  {
		    Stash::Coordinate pos =  hGammaMod->GetBinCoordinate(i,j);
		    int bin = hExcFrac->GetBin(i,j);

		    double alpha = 0;
		    if(hResRingExp->GetBinContent(i,j) != 0 && hResRingExp->GetBinContent(i,j) > 0)
		      alpha =  hResExp->GetBinContent(i,j)*1./hResRingExp->GetBinContent(i,j);
		    
		    if((alpha < hExcFrac->GetBinContent(i,j) && hExcFrac->GetBinContent(i,j) > 0 && binsoff[bin] == 0) || k == 0)
		      {
			hExcFrac->SetBinContent(i,j,alpha);
			
			if(i > nSafeBins && j > nSafeBins && i <= nx-nSafeBins+1 && j <= ny-nSafeBins+1)
			  {
			    hNRing->SetBinContent(i,j,k);
			    hExcFrac->SetBinContent(i,j,alpha);			      
			    if(alpha < fRingParam_AlphaMax && alpha > 0){
			      binsoff[bin] = 1;
			    }
			  }
		      }
		    if(hExcFrac->GetBinContent(i,j) <= 0)
		      hExcFrac->SetBinContent(i,j,1);
		  }		    		    
	      }
	  }
	  else{
	    for(int i = 1; i <= nx; ++i)
	      {
		for(int j = 1; j <= ny; ++j)
		  {

		    Stash::Coordinate pos =  hGammaMod->GetBinCoordinate(i,j);
		    int bin = hExcFrac->GetBin(i,j);
		    double vallRing = hResRingAll->GetBinContent(i,j);
		    double vaccRing = hResRingExp->GetBinContent(i,j);
		    double excfrac = (vallRing-vaccRing)*1./vallRing;
		    
		    if(excfrac < hExcFrac->GetBinContent(i,j) && binsoff[bin] == 0)
		      {
			hExcFrac->SetBinContent(i,j,excfrac);
			
			if(i > nSafeBins && j > nSafeBins && i <= nx-nSafeBins+1 && j <= ny-nSafeBins+1)
			  {
			    hNRing->SetBinContent(i,j,k);
			    hExcFrac->SetBinContent(i,j,excfrac);
			    if(excfrac < fRingParam_ExcFracMax){
				binsoff[bin] = 1;
			    }
			  }
		      }
		  }
	      }
	  }
	  
	  convOExp->Delete();
	  convRExp->Delete();
	  convRingAll->Delete();      
	  
	  hRing->Reset();
	  hResRingExp->Delete();
	  hResRingAll->Delete();
	  hResExp->Delete();
	}
	
      std::cout << hNRing->GetMaximum() << std::endl;
      Int_t NRingsMax = (int)hNRing->GetMaximum();
      
      for(int i = 1; i <= nx; ++i)
	{
	  for(int j = 1; j <= ny; ++j)
	    {
	      int val = (int)hNRing->GetBinContent(i,j);
	      int bin = hNRing->GetBin(i,j);
	      hNRing2->SetBinContent(i,j,val);
	      binsoff2[bin] = val;
	    }
	}
      
      if(fSmoothRings){
	while(NRingsMax > 0)
	  {
	    std::cout << NRingsMax << std::endl;
	    for(int i = 1; i <= nx; ++i)
	      {
		for(int j = 1; j <= ny; ++j)
		  {
		    double val = hNRing2->GetBinContent(i,j);
		    int bin = hNRing2->GetBin(i,j);
		    double nextval = 0;
			
		    if(binsoff2[bin] == NRingsMax && val == NRingsMax)
		      {	  
			for(int ix = i-1; ix <= i+1; ++ix)
			  {
			    for(int iy = j-1; iy <= j+1; ++iy)
			      {
				int numx = ix-i;
				int numy = iy-j;
				int nextBin = hNRing2->GetBin(ix,iy);
				if(binsoff2[nextBin] < val)
				  {
				    if(numx != 0 || numy != 0) 
				      { 
					nextval = hNRing2->GetBinContent(ix,iy);
					hNRing2->SetBinContent(ix,iy,val-1);
					binsoff2[nextBin] = NRingsMax-1;
				      }
				  }
			      }
			  }
		      } 
		  }
	      }
	    NRingsMax--;
	  }
      }
      //Second adapt ring knowing ring numbers for each pixel
      Int_t MaxRingNumber = (int)hNRing->GetMaximum();
      std::cout << "MaxRingNumber = " << MaxRingNumber << std::endl;
      for(int k = 0; k <= MaxRingNumber; ++k)
	{ 
	  Display::SkyHistogram2D *hResRing = new Display::SkyHistogram2D("ResRing","ResRing",*hExcMod); 
	  Display::SkyHistogram2D *hResRingExp = new Display::SkyHistogram2D("ResRingExp","ResRingExp",*hExcMod);
	  Display::SkyHistogram2D *hRes = new Display::SkyHistogram2D("Res","Res",*hExcMod); 
	  Display::SkyHistogram2D *hResExp = new Display::SkyHistogram2D("ResExp","ResExp",*hExcMod); 
	  Display::SkyHistogram2D *hResExc = new Display::SkyHistogram2D("ResExc","ResExc",*hExcMod); 
	  Display::SkyHistogram2D *hResExpExc = new Display::SkyHistogram2D("ResExpExc","ResExpExc",*hExcMod); 

	  std::cout << "InRing " << k << std::endl;      
	  
	  Float_t InnerRing = 0;
	  Float_t OuterRing = 0;
	  if(fAdaptFFT){
	    if(fConstantArea){
	      InnerRing = fRingMinimalRadius + k * fRingStep;
	      OuterRing = sqrt(fRingParam_AreaOverPi + std::pow(InnerRing,2));
	    }
	    else if(fConstantThickness){
	      InnerRing = fRingMinimalRadius + k * fRingStep;
	      OuterRing = InnerRing + fRingParam_Thickness;
	    }
	    else{// if(fConstantInnerRadius){
	      InnerRing = fRingMinimalRadius;
	      OuterRing = InnerRing + fRingParam_Thickness + k * fRingStep;
	    }
	  }
	  else{
	    InnerRing = fStandardRingRadius-fStandardRingThickness;
	    OuterRing = fStandardRingRadius+fStandardRingThickness;
	  }
	  Float_t RingRadius = 0.5*(InnerRing + OuterRing);
	  Float_t RingThickness = RingRadius - InnerRing;
	  Float_t r1 = RingRadius - RingThickness;
	  Float_t r2 = RingRadius + RingThickness;
	  std::cout << "Pass 2 / Step : " << k << "/" << nSteps << std::endl;
	  std::cout << "InnerRing = " << InnerRing << " / " << "OuterRing = " << OuterRing << std::endl ; 
	  
	  //Ring Kernel Generation
	  for(int i = 1; i <= nx; ++i)
	    {
	      for(int j = 1; j <= ny; ++j)
		{
		  double xcenter = 0.5*nx+1;
		  double ycenter = 0.5*ny+1;
		  double r = TMath::Sqrt(TMath::Power(i-xcenter,2.)+TMath::Power(j-ycenter,2.));
		  if(r <= r2/fBinSize && r >= r1/fBinSize){
		    hRing->SetBinContent(i,j,1);
		  }
		}
	    }
	  
	  Utilities::TFftConv *convR = new Utilities::TFftConv(*hGammaModExc,*hRing);
	  Utilities::TFftConv *convRExp = new Utilities::TFftConv(*hGammaAccModExc,*hRing);
	  Utilities::TFftConv *convO = new Utilities::TFftConv(*hGammaMod,*hOS);
	  Utilities::TFftConv *convOExp = new Utilities::TFftConv(*hGammaAccMod,*hOS);
	  Utilities::TFftConv *convOExc = new Utilities::TFftConv(*hGammaModExc,*hOS);
	  Utilities::TFftConv *convOExpExc = new Utilities::TFftConv(*hGammaAccModExc,*hOS);
	  
	  convO->convolve(hRes);
	  convOExp->convolve(hResExp);
	  convOExc->convolve(hResExc);
	  convOExpExc->convolve(hResExpExc);
	  convR->convolve(hResRing);  
	  convRExp->convolve(hResRingExp);
	  
	  for(int i = 1; i <= nx; ++i)
	    {
	      for(int j = 1; j <= ny; ++j)
		{
		  int bin = hExcFrac->GetBin(i,j);
		  
		  if(k == binsoff2[bin])
		    {
		      double vON = hRes->GetBinContent(i,j);
		      double vONExc = hResExc->GetBinContent(i,j);		      
		      double vONUncorrelated = hGammaMod->GetBinContent(i,j);
		      double vONUncorrelated_Exc = hGammaModExc->GetBinContent(i,j);
		      
		      double vOFF = hResRing->GetBinContent(i,j);
		      double vOFFUncorrelated = hResRing->GetBinContent(i,j);

		      double alphaOFFUncorrelated =  0;
		      double alphaOFFUncorrelated_Exc =  0;

		      if(hResRingExp->GetBinContent(i,j) != 0){
			alphaOFFUncorrelated =  vOFFUncorrelated* (hGammaAccMod->GetBinContent(i,j)*1./hResRingExp->GetBinContent(i,j));
			alphaOFFUncorrelated_Exc =  vOFFUncorrelated* (hGammaAccModExc->GetBinContent(i,j)*1./hResRingExp->GetBinContent(i,j));
		      }
		      double AON = hResExp->GetBinContent(i,j);
		      double AONExc = hResExpExc->GetBinContent(i,j);
		      
		      if(i > nSafeBins && j > nSafeBins && i <= nx-nSafeBins+1 && j <= ny-nSafeBins+1)
			{
			  hExcessUncorrelated->SetBinContent(i,j,vONUncorrelated-alphaOFFUncorrelated);
			  hExcessUncorrelated_Exc->SetBinContent(i,j,vONUncorrelated_Exc-alphaOFFUncorrelated_Exc);
			  hRingON->SetBinContent(i,j,vON);
			  hRingON_Exc->SetBinContent(i,j,vONExc);
			  hRingONUncorrelated->SetBinContent(i,j,vONUncorrelated);
			  hRingONUncorrelated_Exc->SetBinContent(i,j,vONUncorrelated_Exc);
			  hRingOFF->SetBinContent(i,j,vOFF);
			  hRingOFFUncorrelated->SetBinContent(i,j,vOFFUncorrelated);
			  hRingAccON->SetBinContent(i,j,AON);
			  hRingAccON_Exc->SetBinContent(i,j,AONExc);
			  hRingAccOFF->SetBinContent(i,j,hResRingExp->GetBinContent(i,j));
			  hAlphaOFFUncorrelated->SetBinContent(i,j,alphaOFFUncorrelated);
			  hAlphaOFFUncorrelated_Exc->SetBinContent(i,j,alphaOFFUncorrelated_Exc);
			}
		    }
		}
	    }
	  
	  convO->Delete();
	  convOExp->Delete();
	  convOExc->Delete();
	  convOExpExc->Delete();
	  convR->Delete();
	  convRExp->Delete();
	  
	  hRing->Reset();
	  hResRing->Delete();
	  hRes->Delete();
	  hResExc->Delete();
	  hResRingExp->Delete();
	  hResExp->Delete();	  
	  hResExpExc->Delete();	  
	}
    
      hExcFrac->Delete();
      hNRing->Delete();
      hNRing2->Delete();
      
      std::cout << "LASTMAP" << std::endl;
      Utilities::TFftConv *convOFF = new Utilities::TFftConv(*hRingOFFUncorrelated,*hOS);
      Utilities::TFftConv *convArea = new Utilities::TFftConv(*hMaskMod,*hOS);
      //
      Utilities::TFftConv *convOFFExp = new Utilities::TFftConv(*hRingAccOFF,*hOS);
      Utilities::TFftConv *convExcess = new Utilities::TFftConv(*hExcessUncorrelated,*hOS);
      Utilities::TFftConv *convExcess_Exc = new Utilities::TFftConv(*hExcessUncorrelated_Exc,*hOS);
      //Utilities::TFftConv *convAlphaOFF = new Utilities::TFftConv(*hAlphaOFFUncorrelated,*hOS);
      //Utilities::TFftConv *convAlphaOFF_Exc = new Utilities::TFftConv(*hAlphaOFFUncorrelated_Exc,*hOS);
      
      convOFFExp->convolve(hResOFFExp);
      convExcess->convolve(hResExcess);
      convExcess_Exc->convolve(hResExcess_Exc);
      convOFF->convolve(hResOFF);
      //convAlphaOFF->convolve(hResAlphaOFF);
      //convAlphaOFF_Exc->convolve(hResAlphaOFF_Exc);
      convArea->convolve(hResConvArea);

      hResOFFExp->Multiply(hMaskMod);
      hResExcess->Multiply(hMaskMod);
      hResExcess_Exc->Multiply(hMaskMod);
      hResOFF->Multiply(hMaskMod);
      //hResAlphaOFF->Multiply(hMaskMod);
      //hResAlphaOFF_Exc->Multiply(hMaskMod);
      hResConvArea->Multiply(hMaskMod);
      
    }

  //UnCorrelated
  Display::SkyHistogram2D *ON =  new Display::SkyHistogram2D("ON","ON",*GetEventsMap()); 
  Display::SkyHistogram2D *ONExposure =  new Display::SkyHistogram2D("ONExposure","ONExposure",*GetEventsMap()); 
  Display::SkyHistogram2D *OFF =  new Display::SkyHistogram2D("OFF","OFF",*GetEventsMap()); 
  Display::SkyHistogram2D *OFFExposure =  new Display::SkyHistogram2D("OFFExposure","OFFExposure",*GetEventsMap()); 
  Display::SkyHistogram2D *Alpha =  new Display::SkyHistogram2D("Alpha","Alpha",*GetEventsMap()); 
  Display::SkyHistogram2D *AlphaOff =  new Display::SkyHistogram2D("AlphaOff","AlphaOff",*GetEventsMap());
  Display::SkyHistogram2D *Excess =  new Display::SkyHistogram2D("Excess","Excess",*GetEventsMap()); 
  //Correlated
  Display::SkyHistogram2D *ONCorrelated =  new Display::SkyHistogram2D("ONCorrelated","ONCorrelated",*GetEventsMap()); 
  Display::SkyHistogram2D *ONExposureCorrelated =  new Display::SkyHistogram2D("ONExposureCorrelated","ONExposureCorrelated",*GetEventsMap()); 
  Display::SkyHistogram2D *OFFCorrelated =  new Display::SkyHistogram2D("OFFCorrelated","OFFCorrelated",*GetEventsMap()); 
  Display::SkyHistogram2D *OFFExposureCorrelated =  new Display::SkyHistogram2D("OFFExposureCorrelated","OFFExposureCorrelated",*GetEventsMap()); 
  Display::SkyHistogram2D *AlphaCorrelated =  new Display::SkyHistogram2D("AlphaCorrelated","AlphaCorrelated",*GetEventsMap()); 
  Display::SkyHistogram2D *AlphaOffCorrelated =  new Display::SkyHistogram2D("AlphaOffCorrelated","AlphaOffCorrelated",*GetEventsMap());
  Display::SkyHistogram2D *ExcessCorrelated =  new Display::SkyHistogram2D("ExcessCorrelated","ExcessCorrelated",*GetEventsMap()); 
  //Significances
  Display::SkyHistogram2D *SignificanceMap =  new Display::SkyHistogram2D("SignificanceMap","SignificanceMap",*GetEventsMap());
  Display::SkyHistogram2D *AllowedSignificanceMap =  new Display::SkyHistogram2D("AllowedSignificanceMap","AllowedSignificanceMap",*GetEventsMap());
  TH1F *SkySigmaDist = new TH1F("SigmaDist","SigmaDist",100,-10,10);    
  TH1F *AllowedRegionsSkySigmaDist = new TH1F("AllowedSigmaDist","AllowedSigmaDist",100,-10,10);

  //Loop to fill global map
  std::cout << "Filling the global maps " << std::endl;

  for(int i = 1; i <= nx; ++i)
    {
      for(int j = 1; j <= ny; ++j)
	  {
	    Stash::Coordinate pos =  hExcMod->GetBinCoordinate(i,j);
	    Int_t bin = GetEventsMap()->FindBinPosition(pos);
	    Int_t binacc = GetAcceptanceMap()->FindBinPosition(pos);
	    
	    if(GetAcceptanceMap()->GetBinContent(binacc) > 0)
	      {
		//Uncorrelated
		double on = hGammaMod->GetBinContent(i,j);
		double off = hRingOFFUncorrelated->GetBinContent(i,j);
		double excess = hExcessUncorrelated->GetBinContent(i,j);
		double alpha = 0;
		if(off > 0)
		  alpha = (on-excess)*1./off;
		double aon = hGammaAccMod->GetBinContent(i,j);
		double aoff = aon*1./alpha;
		
		ON->SetBinContent(bin,on);
		ON->SetBinError(bin,sqrt(on));
		ONExposure->SetBinContent(bin,aon);
		OFF->SetBinContent(bin,off);
		OFF->SetBinError(bin,sqrt(off));
		OFFExposure->SetBinContent(bin,aoff);
		Excess->SetBinContent(bin,on-alpha*off);
		Excess->SetBinError(bin,Utilities::Statistics::LiMa_dExcess_Up(on,off,alpha));
		Alpha->SetBinContent(bin,alpha);
		Alpha->SetBinError(bin,1);
		AlphaOff->SetBinContent(bin,alpha*off);
		AlphaOff->SetBinError(bin,1);
		
		//Correlated
		double on_corr = hRingON->GetBinContent(i,j);
		double off_corr = 0;
		if(fAverageOff)
		  off_corr = hResOFF->GetBinContent(i,j)*1./hResConvArea->GetBinContent(i,j);
		else
		  off_corr = hRingOFF->GetBinContent(i,j);
		double excess_corr = hResExcess->GetBinContent(i,j);
		
		if(off_corr > 5)
		  {
		    double alpha_corr = 0;
		    double onexp_corr = 0;
		    double offexp_corr = 0;
		    if(fAverageOff){
		      //alpha_corr = hResAlphaOFF->GetBinContent(i,j)*1./off_corr;
		      onexp_corr = hRingAccON->GetBinContent(i,j);
		      offexp_corr = hResOFFExp->GetBinContent(i,j)*1./hResConvArea->GetBinContent(i,j);
		      alpha_corr = onexp_corr*1./offexp_corr;
		    }
		    else
		      alpha_corr = (on_corr-excess_corr)*1./off_corr;
		    
		    if(alpha_corr*off_corr > 1)
		      {
			ONCorrelated->SetBinContent(bin,hRingON->GetBinContent(i,j));
			ONCorrelated->SetBinError(bin,sqrt(hRingON->GetBinContent(i,j)));
			ONExposureCorrelated->SetBinContent(bin,onexp_corr);//hRingAccON->GetBinContent(i,j));
			OFFCorrelated->SetBinContent(bin,off_corr);
			OFFCorrelated->SetBinError(bin,sqrt(off_corr));
			OFFExposureCorrelated->SetBinContent(bin,offexp_corr);//1./alpha_corr*hRingAccON->GetBinContent(i,j));
			double significance = Utilities::Statistics::LiMa(on_corr,off_corr,alpha_corr);
			SignificanceMap->SetBinContent(bin,significance);
			SignificanceMap->SetBinError(bin,1);	  
			
			ExcessCorrelated->SetBinContent(bin,on_corr-alpha_corr*off_corr);
			ExcessCorrelated->SetBinError(bin,Utilities::Statistics::LiMa_dExcess_Up(on_corr,off_corr,alpha_corr));
			AlphaCorrelated->SetBinContent(bin,alpha_corr);
			AlphaCorrelated->SetBinError(bin,1);
			AlphaOffCorrelated->SetBinContent(bin,alpha_corr*off_corr);
			AlphaOffCorrelated->SetBinError(bin,1);
	
			if(significance != 0)
			  SkySigmaDist->Fill(significance);
			
			if(GetExclusionMap()->GetBinContent(GetExclusionMap()->FindBinPosition(pos)) == 0 && significance != 0)
			  {
			    double onexc = hRingON_Exc->GetBinContent(i,j);
			    double offexc = 0;
			    if(fAverageOff)
			      offexc = hResOFF->GetBinContent(i,j)*1./hResConvArea->GetBinContent(i,j);
			    else
			      offexc = hRingOFF->GetBinContent(i,j);
			    
			    double alphaexc = 0;
			    double onexpexc = 0;
			    double offexpexc = 0;
			    if(fAverageOff){
			      onexpexc = hRingAccON_Exc->GetBinContent(i,j);
			      offexpexc = hResOFFExp->GetBinContent(i,j)*1./hResConvArea->GetBinContent(i,j);;
			      alphaexc = onexpexc*1./offexpexc;//hResAlphaOFF_Exc->GetBinContent(i,j)*1./offexc;
			    }
			    else
			      alphaexc = (onexc-hResExcess_Exc->GetBinContent(i,j))*1./hRingOFF->GetBinContent(i,j);
			    
			    double significance_exc = Utilities::Statistics::LiMa(onexc,offexc,alphaexc);
			    AllowedRegionsSkySigmaDist->Fill(significance_exc);
			    AllowedSignificanceMap->SetBinContent(bin,significance_exc);
			    AllowedSignificanceMap->SetBinError(bin,1);	  
			  }
		      }
		  }
	      }
	  }
      }
    
  /*TCanvas *c0 = new TCanvas("FFTMaps","FFTMaps",1000,500);
  c0->Divide(3,2);
  c0->cd(1);
  hGammaMod->Draw("colz");
  c0->cd(2);
  hGammaAccModExc->Draw("colz");
  c0->cd(3);
  hGammaAccMod->Draw("colz");
  c0->cd(4);
  hGammaModExc->Draw("colz");
  c0->cd(5);
  hExcMod->Draw("colz");
  c0->cd(6);
  hOS->Draw("colz");*/      
  
  /*  TCanvas *c0 = new TCanvas("NRMaps","NRMaps",1000,500);
  c0->Divide(2,1);
  c0->cd(1);
  hNRing->Draw("colz");
  c0->cd(2);
  hNRing2->Draw("colz");*/

  TCanvas *cCorr = new TCanvas("BgMapsCorrelated","BgMapsCorrelated",1000,700);
  cCorr->Divide(3,2);
  cCorr->cd(1);
  ExcessCorrelated->Draw("colz");
  cCorr->cd(2);
  SignificanceMap->Draw("colz");
  cCorr->cd(3);
  AllowedSignificanceMap->Draw("colz");
  cCorr->cd(4);
  AlphaCorrelated->Draw("colz");
  cCorr->cd(5);
  AlphaOffCorrelated->Draw("colz");
  cCorr->cd(6)->SetLogy();
  AllowedRegionsSkySigmaDist->SetLineColor(2);
  AllowedRegionsSkySigmaDist->Draw();
  SkySigmaDist->Draw("same");
  
  GetExcessMap() = ExcessCorrelated;
  GetAlphaMap() = AlphaCorrelated;
  GetOffMap() = OFFCorrelated;

  if(fSaveResults){
    //-----------------------------------------------------------------------------
    //Generate File Name according to analysis parameters    
    //and write in file
    std::ostringstream foutfile;
    std::ostringstream canvfile;

    std::ostringstream EvAndAcc;
    std::ostringstream ConfigTarget;
    std::ostringstream OverSampling;
    std::ostringstream AdaptCutAlpha;
    std::ostringstream RingMethod;
    std::ostringstream EnergyCuts;
     
    OverSampling << "_OS" << fOSRadius;
    if(fEventsAndAcceptanceFromRadial){
      EvAndAcc << "_RadialLookups_PsiCut" << fPsiCut;
      if(fCorrectZenith)
	EvAndAcc << "_ZenCorr";
      if(fApplySafeThreshold){
	EvAndAcc << "_SafeTh" << fSafeThresholdRatioParam;
	if(fSafeThresholdFromAcceptance)
	  EvAndAcc << "_FromAcc";
	else
	  EvAndAcc << "_FromBias";
      }
    }
    else 
      EvAndAcc << "";
    
    if(!fUseConfigTarget)
      ConfigTarget << "_CustomTarget_" << fUserLambda << "_" << fUserBeta;
    else
      ConfigTarget << "";
    
    if(!fUseConfigMapParams)    
      ConfigTarget << "_MapExt_X" << fExtX << "_Y" << fExtY;
    else
      ConfigTarget << "";
    
    if(fAdaptFFT){
      RingMethod << "_AdaptRing";
      if(fAdaptCut_Alpha)
	RingMethod << "_AlphaCut" << fRingParam_AlphaMax;
      else
	RingMethod << "_AreaCut" << fRingParam_ExcFracMax;
      
      if(fConstantArea)
	RingMethod << "_CstArea" << fRingParam_AreaOverPi << "_RInMin" << fRingMinimalRadius << "_RInMax" << fInnerRingMax << "_RStep" << fRingStep;
      if(fConstantThickness)
	RingMethod << "_CstThick" << fRingParam_Thickness << "_RInMin" << fRingMinimalRadius << "_RInMax" << fInnerRingMax << "_RStep" << fRingStep;
      else//if(fConstantInnerRadius) 
	RingMethod << "_CstInner" << fRingMinimalRadius << "_ROutMin" << fRingMinimalRadius+fRingParam_Thickness << "_ROutMax" << fOuterRingMax << "_RStep" << fRingStep;
    }
    else
      RingMethod << "_StdRing_R" << fStandardRingRadius << "_Thick" << fStandardRingThickness;
    
    if(fEMin != 0 && fEMax == -1)
      EnergyCuts << "_EMin" << fEMin;
    else if(fEMin == 0 && fEMax != -1)
      EnergyCuts << "_EMax" << fEMax;
    else if(fEMin != 0 && fEMax != -1)
      EnergyCuts << "_EMin" << fEMin << "_EMax" << fEMax;
    else
      EnergyCuts << "";
    
        
    foutfile << OutfilePrefix << "_RingBgMaps.root";
    canvfile << OutfilePrefix << "_BgMapsCorrelated.png";
    cCorr->Print(canvfile.str().c_str());
    
    TFile fileResultsOut(foutfile.str().c_str(),"RECREATE");
    //  gROOT->cd();    
    
    fileResultsOut.Append(ON);
    fileResultsOut.Append(OFF);
    fileResultsOut.Append(ONExposure);
    fileResultsOut.Append(OFFExposure);
    fileResultsOut.Append(Alpha);
    fileResultsOut.Append(AlphaOff);
    fileResultsOut.Append(Excess);
    fileResultsOut.Append(ONCorrelated);
    fileResultsOut.Append(OFFCorrelated);
    fileResultsOut.Append(ONExposureCorrelated);
    fileResultsOut.Append(OFFExposureCorrelated);
    fileResultsOut.Append(AlphaCorrelated);
    fileResultsOut.Append(AlphaOffCorrelated);
    fileResultsOut.Append(ExcessCorrelated);
    
    fileResultsOut.Append(SignificanceMap);
    fileResultsOut.Append(AllowedSignificanceMap);
    fileResultsOut.Append(SkySigmaDist);
    fileResultsOut.Append(AllowedRegionsSkySigmaDist);
    
    fileResultsOut.Append(GetAcceptanceMap());
    fileResultsOut.Append(GetEventsMap());
    fileResultsOut.Append(GetOffAxisMap());
    fileResultsOut.Append(GetZenithAngleMap());
    fileResultsOut.Append(GetMuonEfficiencyMap());
    if(fApplySafeThreshold){
      fileResultsOut.Append(GetMinSafeThresholdMap());
      fileResultsOut.Append(GetAveragedSafeThresholdMap());
    }
    fileResultsOut.Append(GetAveragedLiveTimeMap());
    fileResultsOut.Write();
    fileResultsOut.Clear();

    if(fProduceFluxProducts){
      CreateFluxMaps_FromResultFile(AnalysisConfig, foutfile.str().c_str(), ExposureMapFile, ExtendedExposureMapFile, OutfilePrefix);
    }
    
  }

}

void SurveySuite::MapMaker::CreateFluxMaps_FromResultFile(const char *AnalysisConfig, const char *ResultFile, const char *ExposureMapFile, const char *ExtendedExposureMapFile, const char *OutfilePrefix)
{
  TFile *fileResults = TFile::Open(ResultFile);
  gROOT->cd();
  Display::SkyHistogram2D *OFFMap = (Display::SkyHistogram2D *)fileResults->Get("OFFCorrelated");
  Display::SkyHistogram2D *AlphaMap = (Display::SkyHistogram2D *)fileResults->Get("AlphaCorrelated");
  Display::SkyHistogram2D *ExcessMap = (Display::SkyHistogram2D *)fileResults->Get("ExcessCorrelated");
  Display::SkyHistogram2D *UncorrelatedExcessMap = (Display::SkyHistogram2D *)fileResults->Get("Excess");

  std::cout << "WARNING : Make sure all your maps match!!! " << std::endl;
  std::cout << "Create exposure" << std::endl;
  if(!GetExpectedCountsMap()){
    if(fExposureMapsFromFits){
      CreateExposureMaps_FromFITS(ExposureMapFile);
      CreateExtendedExposureMaps_FromFITS(ExtendedExposureMapFile);
    }
    else{
      CreateExposureMaps_FromFile(ExposureMapFile);
      CreateExtendedExposureMaps_FromFile(ExtendedExposureMapFile);
      }
  }
  else{
    
  }

  std::cout << GetExpectedCountsMap()->GetBinContent(3000,200) << " " << GetExtendedExpectedCountsMap()->GetBinContent(3000,200) << std::endl;  

  GetUncorrelatedExcessMap() = UncorrelatedExcessMap;
  GetExcessMap() = ExcessMap;
  GetAlphaMap() = AlphaMap;
  GetOffMap() = OFFMap;
  
  Display::SkyHistogram2D *hFluxMap =  new Display::SkyHistogram2D("FluxMap","Differential Flux Map (1 TeV)",*GetExcessMap());   
  Display::SkyHistogram2D *hFluxMap_Uncorrelated =  new Display::SkyHistogram2D("FluxMap_Uncorrelated","Uncorrelated Differential Flux Map (1 TeV)",*GetExcessMap());   
  Display::SkyHistogram2D *hFluxErrorMap =  new Display::SkyHistogram2D("FluxErrorMap","Differential Flux Error Map (1 TeV)",*GetExpectedCountsMap());   
  Display::SkyHistogram2D *hFluxErrorMap_Uncorrelated =  new Display::SkyHistogram2D("FluxErrorMap_Uncorrelated","Uncorrelated Differential Flux Error Map (1 TeV)",*GetExpectedCountsMap());   
  Display::SkyHistogram2D *hIntFluxMap =  new Display::SkyHistogram2D("IntFluxMap","Integral Flux Map (>1 TeV)",*GetExpectedCountsMap());
  Display::SkyHistogram2D *hIntFluxMap_Uncorrelated =  new Display::SkyHistogram2D("IntFluxMap_Uncorrelated","Uncorrelated Integral Flux Map (>1 TeV)",*GetExpectedCountsMap());
  Display::SkyHistogram2D *hIntFluxErrorMap =  new Display::SkyHistogram2D("IntFluxErrorMap","Integral Flux Error Map (>1 TeV)",*GetExpectedCountsMap());
  Display::SkyHistogram2D *hIntFluxErrorMap_Uncorrelated =  new Display::SkyHistogram2D("IntFluxErrorMap_Uncorrelated","Uncorrelated Integral Flux Error Map (>1 TeV)",*GetExpectedCountsMap());
  Display::SkyHistogram2D *hFluxULMap =  new Display::SkyHistogram2D("FluxULMap","Upper limit Flux Map (CL 0.95, 1 TeV)",*GetExpectedCountsMap());
  Display::SkyHistogram2D *hIntFluxULMap =  new Display::SkyHistogram2D("IntFluxULMap","Upper limit Int Flux Map (CL 0.95, >1 TeV)",*GetExpectedCountsMap());
  Display::SkyHistogram2D *hSensitivityMap =  new Display::SkyHistogram2D("FluxSensitivityMap","Sensitivity Map (1 TeV)",*GetExpectedCountsMap());
  Display::SkyHistogram2D *hIntSensitivityMap =  new Display::SkyHistogram2D("IntFluxSensitivityMap","Integral Sensitivity Map (>1 TeV)",*GetExpectedCountsMap());
  Display::SkyHistogram2D *hSensitivityMapCrab =  new Display::SkyHistogram2D("FluxSensitivityMapCrab","Sensitivity Map (percent Crab, 1 TeV)",*GetExpectedCountsMap());
  Display::SkyHistogram2D *hIntSensitivityMapCrab =  new Display::SkyHistogram2D("IntFluxSensitivityMapCrab","Sensitivity Map (percent Crab, > 1 TeV)",*GetExpectedCountsMap());
  Display::SkyHistogram2D *hSurfaceBrightnessMap =  new Display::SkyHistogram2D("SurfaceBrightnessMap","Surface Brightness Map (1 TeV)",*GetExpectedCountsMap());
  Display::SkyHistogram2D *hSurfaceBrightnessErrorMap =  new Display::SkyHistogram2D("SurfaceBrightnessErrorMap","Surface Brightness Error Map (1 TeV)",*GetExpectedCountsMap());
  Display::SkyHistogram2D *hIntSurfaceBrightnessMap =  new Display::SkyHistogram2D("IntSurfaceBrightnessMap","Surface Brightness Map (>1 TeV)",*GetExpectedCountsMap());
  Display::SkyHistogram2D *hIntSurfaceBrightnessErrorMap =  new Display::SkyHistogram2D("IntSurfaceBrightnessErrorMap","Surface Brightness Error Map (>1 TeV)",*GetExpectedCountsMap());

  Display::SkyHistogram2D *hFluxMapCleaned =  new Display::SkyHistogram2D("FluxMapCleaned","Differential Flux Map (1 TeV) - FluxSens < 3pc Crab",*GetExcessMap());   
  Display::SkyHistogram2D *hFluxErrorMapCleaned =  new Display::SkyHistogram2D("FluxErrorMapCleaned","Differential Flux Error Map (1 TeV) - FluxSens < 3pc Crab",*GetExpectedCountsMap());   
  Display::SkyHistogram2D *hIntFluxMapCleaned =  new Display::SkyHistogram2D("IntFluxMapCleaned","Integral Flux Map (>1 TeV)",*GetExpectedCountsMap());
  Display::SkyHistogram2D *hIntFluxErrorMapCleaned =  new Display::SkyHistogram2D("IntFluxErrorMapCleaned","Integral Flux Error Map (>1 TeV) - FluxSens < 3pc Crab",*GetExpectedCountsMap());
  Display::SkyHistogram2D *hFluxULMapCleaned =  new Display::SkyHistogram2D("FluxULMapCleaned","Upper limit Flux Map (CL 0.95, 1 TeV) - FluxSens < 3pc Crab",*GetExpectedCountsMap());
  Display::SkyHistogram2D *hIntFluxULMapCleaned =  new Display::SkyHistogram2D("IntFluxULMapCleaned","Upper limit Int Flux Map (CL 0.95, >1 TeV) - FluxSens < 3pc Crab",*GetExpectedCountsMap());

  
  TRolke* rolke = 0;
  rolke = new TRolke(0.95);
  rolke->SetBounding(true); //use method 2 of Rolke et al (2005)

  TRolke* rolkeErr = 0;
  rolkeErr = new TRolke(0.68);
  rolkeErr->SetBounding(true); //use method 2 of Rolke et al (2005)

  for(int i = 1; i <= hFluxMap->GetNbinsX(); ++i)
    {
      for(int j = 1; j <= hFluxMap->GetNbinsY(); ++j)
	{
	  Stash::Coordinate pos =  hFluxMap->GetBinCoordinate(i,j);
	  Int_t bin = GetExpectedCountsMap()->FindBinPosition(pos);

	  double nexp_pl = GetExpectedCountsMap()->GetCircularAverage(pos,Stash::Lambda(fOSRadius,Stash::Angle::Degrees));//GetBinContent(bin);
	  double nexp_ext = GetExtendedExpectedCountsMap()->GetCircularAverage(pos,Stash::Lambda(fOSRadius,Stash::Angle::Degrees));//GetBinContent(bin);

	  double nexp_pl_uncorr = GetExpectedCountsMap()->GetBinContent(bin);
	  double nexp_ext_uncorr = GetExtendedExpectedCountsMap()->GetBinContent(bin);

	  double nexp = 0;
	  double nexp_uncorr = 0;
	  if(fPointLikeFluxMaps){
	    nexp = nexp_pl;
	    nexp_uncorr = nexp_pl_uncorr;
	  }
	  else{
	    nexp = nexp_ext;
	    nexp_uncorr = nexp_ext_uncorr;
	  }

	  double nobs_uncorr = GetUncorrelatedExcessMap()->GetBinContent(bin);	  
	  double nobs_err_uncorr = GetUncorrelatedExcessMap()->GetBinError(bin);	  
	  double nobs = GetExcessMap()->GetBinContent(bin);	  
	  double noff = GetOffMap()->GetBinContent(bin);
	  double alpha = GetAlphaMap()->GetBinContent(bin);
	  double nobs_err = GetExcessMap()->GetBinError(bin);
	  
	  if(nexp > 0){

	    //SENSITIVITY
	    int n5sigma = 0;//Utilities::Statistics::LiMa_requiredExcess(5,noff,alpha);
	    //  try{
	      n5sigma = Utilities::Statistics::LiMa_requiredExcess( 5, noff, alpha);
	      /*}
	    catch(Utilities::Statistics::InvalidResult &res) {
	      n5sigma = 0.;
	    }
	    catch(...){
	      n5sigma = 0.;
	      }*/
	    
	    double flux1TeVSensitivity = n5sigma*1./nexp;
	    double intFlux1TeVSensitivity = Utilities::Flux::IntFlux(flux1TeVSensitivity,1.,-1*fSpectralIndex,1.,-1,0);
	    double flux1TeVSensitivity_PercentCrabUnits = flux1TeVSensitivity * 100./(3.76e-7); //Crab units
	    double intFlux1TeVSensitivity_PercentCrabUnits = intFlux1TeVSensitivity * 100./(2.27e-7); //Crab units
	    hSensitivityMap->SetBinContent(i,j,flux1TeVSensitivity);
	    hIntSensitivityMap->SetBinContent(i,j,intFlux1TeVSensitivity);
	    hSensitivityMapCrab->SetBinContent(i,j,flux1TeVSensitivity_PercentCrabUnits);
	    hIntSensitivityMapCrab->SetBinContent(i,j,intFlux1TeVSensitivity_PercentCrabUnits);

	    //DIFF FLUX
	    double flux1TeV = nobs*1./nexp;
	    hFluxMap->SetBinContent(i,j,flux1TeV);
	    //...uncorrelated
	    double flux1TeV_uncorr = nobs_uncorr*1./nexp_uncorr;
	    hFluxMap_Uncorrelated->SetBinContent(i,j,flux1TeV_uncorr);
	    //INT FLUX
	    double intFlux1TeV = Utilities::Flux::IntFlux(flux1TeV,1.,-1*fSpectralIndex,1.,-1,0);
	    hIntFluxMap->SetBinContent(i,j,intFlux1TeV);
	    //...uncorrelated
	    double intFlux1TeV_uncorr = Utilities::Flux::IntFlux(flux1TeV_uncorr,1.,-1*fSpectralIndex,1.,-1,0);
	    hIntFluxMap_Uncorrelated->SetBinContent(i,j,intFlux1TeV_uncorr);
	    //ERRORS
	    /*	    double nerr_up, nerr_dn, nerr_mn;
	    if(alpha > 0.0){
	      double non = nobs + alpha*noff;
	      rolkeErr->SetPoissonBkgKnownEff(non, noff, 1.0 *1./alpha, 1);
	      rolkeErr->GetLimits(nerr_dn, nerr_up);
	      nerr_dn = nobs-nerr_dn;
	      nerr_up = nerr_up-nobs;
	      }*/
	    double fluxErr1TeV = nobs_err*1./nexp;
	    double intFluxErr1TeV = Utilities::Flux::IntFlux(fluxErr1TeV,1.,-1*fSpectralIndex,1.,-1,0);
	    hFluxErrorMap->SetBinContent(i,j,fluxErr1TeV);
	    hIntFluxErrorMap->SetBinContent(i,j,intFluxErr1TeV);
	    //...uncorrelated
	    double fluxErr1TeV_uncorr = nobs_err_uncorr*1./nexp_uncorr;
	    double intFluxErr1TeV_uncorr = Utilities::Flux::IntFlux(fluxErr1TeV_uncorr,1.,-1*fSpectralIndex,1.,-1,0);
	    hFluxErrorMap_Uncorrelated->SetBinContent(i,j,fluxErr1TeV_uncorr);
	    hIntFluxErrorMap_Uncorrelated->SetBinContent(i,j,intFluxErr1TeV_uncorr);	    

	    //UL
	    double nSignalUL = 0;	    
	    if(alpha > 0.0){
	      double excess_down, excess_up;
	      double non = nobs + alpha*noff;
	      rolke->SetPoissonBkgKnownEff(non, noff, 1.0 *1./alpha, 1);
	      rolke->GetLimits(excess_down, excess_up);
	      nSignalUL = excess_up;
	    }
	    
	    double flux1TeVUL = nSignalUL*1./nexp;
	    hFluxULMap->SetBinContent(i,j,flux1TeVUL);
	    double intFlux1TeVUL = Utilities::Flux::IntFlux(flux1TeVUL,1.,-1*fSpectralIndex,1.,-1,0);
	    hIntFluxULMap->SetBinContent(i,j,intFlux1TeVUL);
	    

	    if(flux1TeVSensitivity_PercentCrabUnits < 3){
	      hFluxMapCleaned->SetBinContent(i,j,flux1TeV);
	      hFluxErrorMapCleaned->SetBinContent(i,j,fluxErr1TeV);	      
	      hIntFluxMapCleaned->SetBinContent(i,j,intFlux1TeV);
	      hIntFluxErrorMapCleaned->SetBinContent(i,j,intFluxErr1TeV);
	      hFluxULMapCleaned->SetBinContent(i,j,flux1TeVUL);
	      hIntFluxULMapCleaned->SetBinContent(i,j,intFlux1TeVUL);
	    }

	    //SURFACE BRIGHTNESS - Exposure map must be from extended config!!!
	    double intsb = 1./(TMath::Pi()*pow(fOSRadius,2))*Utilities::Flux::IntFlux(nobs*1./nexp_ext,1.,-1*fSpectralIndex,1.,-1,0);
	    double intsberr = 1./(TMath::Pi()*pow(fOSRadius,2))*Utilities::Flux::IntFlux(nobs_err*1./nexp_ext,1.,-1*fSpectralIndex,1.,-1,0);
	    hIntSurfaceBrightnessMap->SetBinContent(i,j,intsb);
	    hIntSurfaceBrightnessErrorMap->SetBinContent(i,j,intsberr);
	    
	    double sb = 1./(TMath::Pi()*pow(fOSRadius,2))*nobs*1./nexp_ext;
	    double sberr = 1./(TMath::Pi()*pow(fOSRadius,2))*nobs_err*1./nexp_ext;
	    hSurfaceBrightnessMap->SetBinContent(i,j,sb);
	    hSurfaceBrightnessErrorMap->SetBinContent(i,j,sberr);
	  }
	}
    }
  
  TCanvas *cFlux = new TCanvas("FluxMaps","FluxMaps",600,500);
  hIntFluxMap->Draw("colz");
  
  GetExpectedCountsMap()->SetNameTitle("AreaTimeMap","AreaTimeMap");
  GetExtendedExpectedCountsMap()->SetNameTitle("AreaTimeMap_Extended","AreaTimeMap_Extended");
  std::ostringstream foutfile;
  std::ostringstream EnergyCuts;

  if(fEMin != 0 && fEMax == -1)
    EnergyCuts << "_EMin" << fEMin;
  else if(fEMin == 0 && fEMax != -1)
    EnergyCuts << "_EMax" << fEMax;
  else if(fEMin != 0 && fEMax != -1)
    EnergyCuts << "_EMin" << fEMin << "_EMax" << fEMax;
  else
    EnergyCuts << "";
  
  foutfile << OutfilePrefix << "_FluxMaps.root";
  
  TFile fileResultsOut(foutfile.str().c_str(),"RECREATE");

  fileResultsOut.Append(GetExcessMap());
  fileResultsOut.Append(GetAlphaMap());
  fileResultsOut.Append(GetOffMap());
  fileResultsOut.Append(GetExpectedCountsMap());
  fileResultsOut.Append(GetExtendedExpectedCountsMap());
  fileResultsOut.Append(hFluxMap);
  fileResultsOut.Append(hFluxErrorMap);
  fileResultsOut.Append(hIntFluxMap);
  fileResultsOut.Append(hIntFluxErrorMap);
  fileResultsOut.Append(hFluxMap_Uncorrelated);
  fileResultsOut.Append(hFluxErrorMap_Uncorrelated);
  fileResultsOut.Append(hIntFluxMap_Uncorrelated);
  fileResultsOut.Append(hIntFluxErrorMap_Uncorrelated);
  fileResultsOut.Append(hFluxULMap);
  fileResultsOut.Append(hIntFluxULMap);
  fileResultsOut.Append(hSensitivityMap);
  fileResultsOut.Append(hIntSensitivityMap);
  fileResultsOut.Append(hSensitivityMapCrab);
  fileResultsOut.Append(hIntSensitivityMapCrab);    
  fileResultsOut.Append(hSurfaceBrightnessMap);
  fileResultsOut.Append(hSurfaceBrightnessErrorMap);
  fileResultsOut.Append(hIntSurfaceBrightnessMap);
  fileResultsOut.Append(hIntSurfaceBrightnessErrorMap);

  fileResultsOut.Append(hFluxMapCleaned);
  fileResultsOut.Append(hFluxErrorMapCleaned);
  fileResultsOut.Append(hIntFluxMapCleaned);
  fileResultsOut.Append(hIntFluxErrorMapCleaned);
  fileResultsOut.Append(hFluxULMapCleaned);
  fileResultsOut.Append(hIntFluxULMapCleaned);

  fileResultsOut.Write();
  fileResultsOut.Clear();
}

void SurveySuite::MapMaker::CreateEventsAndAcceptanceMaps_FromRadialLookups_PerRun(const char *RadialConfig, const char *Resfile, const char *Regionsfile, std::string table_path, const char *OutfilePrefix)
{
  std::vector<int>::iterator it = fEntriesVector.begin();
  for(; it != fEntriesVector.end(); ++it)
    {
      int i = *it;
      CreateEventsAndAcceptanceMaps_FromRadialLookups_OneEntry(RadialConfig, Resfile, Regionsfile, table_path, OutfilePrefix, i);
    }
}

void SurveySuite::MapMaker::CreateBgMaps_PerRun(const char *AnalysisConfig, const char *Resfile, const char *Evtfile, const char *Accfile, const char *ExposureMapFile, const char *ExtendedExposureMapFile, const char *Regionsfile, std::string table_path, const char *OutfilePrefix)
{
  /*  std::map<int,Display::SkyHistogram2D *> EventsPerRuns;
  std::map<int,Display::SkyHistogram2D *> AccPerRuns;
  std::map<int,Display::SkyHistogram2D *> ExcessPerRuns;
  std::map<int,Display::SkyHistogram2D *> SignificancePerRuns;*/
  int last = 0;

  std::ostringstream foutfile;
  std::ostringstream EvAndAcc;
  std::ostringstream ConfigTarget;
  std::ostringstream OverSampling;
  std::ostringstream AdaptCutAlpha;
  std::ostringstream RingMethod;
  std::ostringstream EnergyCuts;
  
  if(fSaveResults){   
    OverSampling << "_OS" << fOSRadius;
    if(fEventsAndAcceptanceFromRadial){
      EvAndAcc << "_RadialLookups_PsiCut" << fPsiCut;
      if(fCorrectZenith)
  	EvAndAcc << "_ZenCorr";
      if(fApplySafeThreshold){
  	EvAndAcc << "_SafeTh" << fSafeThresholdRatioParam;
  	if(fSafeThresholdFromAcceptance)
  	  EvAndAcc << "_FromAcc";
  	else
  	  EvAndAcc << "_FromBias";
      }
    }
    else 
      EvAndAcc << "";
    
    if(!fUseConfigTarget)
      ConfigTarget << "_CustomTarget_" << fUserLambda << "_" << fUserBeta;
    else
      ConfigTarget << "";
    
    if(!fUseConfigMapParams)    
      ConfigTarget << "_MapExt_X" << fExtX << "_Y" << fExtY;
    else
      ConfigTarget << "";
    
    if(fAdaptFFT){
      RingMethod << "_AdaptRing";
      if(fAdaptCut_Alpha)
  	RingMethod << "_AlphaCut" << fRingParam_AlphaMax;
      else
  	RingMethod << "_AreaCut" << fRingParam_ExcFracMax;
      
      if(fConstantArea)
  	RingMethod << "_CstArea" << fRingParam_AreaOverPi << "_RInMin" << fRingMinimalRadius << "_RInMax" << fInnerRingMax << "_RStep" << fRingStep;
      if(fConstantThickness)
  	RingMethod << "_CstThick" << fRingParam_Thickness << "_RInMin" << fRingMinimalRadius << "_RInMax" << fInnerRingMax << "_RStep" << fRingStep;
      else//if(fConstantInnerRadius) 
  	RingMethod << "_CstInner" << fRingMinimalRadius << "_ROutMin" << fRingMinimalRadius+fRingParam_Thickness << "_ROutMax" << fOuterRingMax << "_RStep" << fRingStep;
    }
    else
      RingMethod << "_StdRing_R" << fStandardRingRadius << "_Thick" << fStandardRingThickness;
    
    if(fEMin != 0 && fEMax == -1)
      EnergyCuts << "_EMin" << fEMin;
    else if(fEMin == 0 && fEMax != -1)
      EnergyCuts << "_EMax" << fEMax;
    else if(fEMin != 0 && fEMax != -1)
      EnergyCuts << "_EMin" << fEMin << "_EMax" << fEMax;
    else
      EnergyCuts << "";
    
    foutfile << OutfilePrefix << "_PerRunMaps.root";
  }

  TFile fileResultsOut(foutfile.str().c_str(),"RECREATE");
  //  gROOT->cd();    

  std::vector<int>::iterator it = fEntriesVector.begin();
  for(; it != fEntriesVector.end(); ++it)
    {
      int i = *it;
      CreateRingBgMaps_OneEntry(AnalysisConfig, Resfile, Evtfile, Accfile, ExposureMapFile, ExtendedExposureMapFile, Regionsfile, table_path, OutfilePrefix, i);

      Display::SkyHistogram2D * hEvents = (Display::SkyHistogram2D *)GetEventsMap()->Clone();    
      Display::SkyHistogram2D * hAcceptance = (Display::SkyHistogram2D *)GetAcceptanceMap()->Clone();    
      Display::SkyHistogram2D * hExcess = (Display::SkyHistogram2D *)GetExcessMap()->Clone();    
      Display::SkyHistogram2D * hSignificance = (Display::SkyHistogram2D *)GetSignificanceMap()->Clone();    

      std::ostringstream fname_Events;
      fname_Events << "Events_" << i;
      std::ostringstream fname_Acceptance;
      fname_Acceptance << "Acceptance_" << i;
      std::ostringstream fname_Excess;
      fname_Excess << "Excess_" << i;
      std::ostringstream fname_Significance;
      fname_Significance << "Significance_" << i;

      hEvents->SetNameTitle(fname_Events.str().c_str(),fname_Events.str().c_str());
      hAcceptance->SetNameTitle(fname_Acceptance.str().c_str(),fname_Acceptance.str().c_str());
      hExcess->SetNameTitle(fname_Excess.str().c_str(),fname_Excess.str().c_str());
      hSignificance->SetNameTitle(fname_Significance.str().c_str(),fname_Significance.str().c_str());
      
      /*EventsPerRuns[i] = GetEventsMap();
      AccPerRuns[i] = GetAcceptanceMap();
      ExcessPerRuns[i] = GetExcessMap();
      SignificancePerRuns[i] = GetSignificanceMap();
      */
      fileResultsOut.Append(hEvents);
      fileResultsOut.Append(hAcceptance);
      fileResultsOut.Append(hExcess);
      fileResultsOut.Append(hSignificance);

      last = i;      
    }

  /*  TCanvas *cMaps = new TCanvas("cc","cc",1200,600);
  cMaps->Divide(2,2);
  cMaps->cd(1);
  EventsPerRuns[last]->Draw("colz");
  cMaps->cd(2);
  AccPerRuns[last]->Draw("colz");
  cMaps->cd(3);
  ExcessPerRuns[last]->Draw("colz");
  cMaps->cd(4);
  SignificancePerRuns[last]->Draw("colz");*/
  
  fileResultsOut.Write();
  fileResultsOut.Clear();
 
}

void SurveySuite::MapMaker::CreateONOFFTestMaps(const char *AnalysisConfig, const char *Resfile, const char *Evtfile, const char *Accfile, const char *ExposureMapFile, const char *ExtendedExposureMapFile, const char *Regionsfile, std::string table_path, const char *OutfilePrefix)
{
  //Opening the results file
  TFile *fileResults = TFile::Open(Resfile);
  gROOT->cd();
  Sash::DataSet *results = (Sash::DataSet *)fileResults->Get("results");
  results->GetEntry(0);
  
  Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();
  ParisAnalysis::AnalysisConfig *Config =  hess->Handle("", (ParisAnalysis::AnalysisConfig *) 0);
  Config->LoadAllMembers();
  
  Sash::DataSet *run_results = (Sash::DataSet *)fileResults->Get("run_results");
  gROOT->cd();
  
  gStyle->SetOptStat(000);

  Stash::Coordinate targetpos(Stash::Lambda(Config->GetTargetRa(),Stash::Angle::Degrees),
			      Stash::Beta(Config->GetTargetDec(),Stash::Angle::Degrees),
			      Config->GetSystem());

  int count = 0;
  std::map<std::pair<int,int>,std::map <int,double> > vAccs, vEvents, vAllowed, vDate;

  // std::map<int,Display::SkyHistogram2D *> EventsPerRuns;
  // std::map<int,Display::SkyHistogram2D *> AccPerRuns;
  // std::map<int,Display::SkyHistogram2D *> ExcessPerRuns;
  // std::map<int,Display::SkyHistogram2D *> SignificancePerRuns;
  // int last = 0;

  TH1F *hSignificanceDistribution = new TH1F("ONOFFSignifDist","ONOFFSignifDist",100,-10,10);
  TH1F *hSignificanceDistribution_Allowed = new TH1F("ONOFFSignifDist_Allowed","ONOFFSignifDist_Allowed",100,-10,10);
  
  //FILE SAVING
  std::ostringstream foutfile;
  std::ostringstream EvAndAcc;
  std::ostringstream ConfigTarget;
  std::ostringstream OverSampling;
  std::ostringstream AdaptCutAlpha;
  std::ostringstream RingMethod;
  std::ostringstream EnergyCuts;

  if(fSaveResults){   
    OverSampling << "_OS" << fOSRadius;
    if(fEventsAndAcceptanceFromRadial){
      EvAndAcc << "_RadialLookups_PsiCut" << fPsiCut;
      if(fCorrectZenith)
  	EvAndAcc << "_ZenCorr";
      if(fApplySafeThreshold){
  	EvAndAcc << "_SafeTh" << fSafeThresholdRatioParam;
  	if(fSafeThresholdFromAcceptance)
  	  EvAndAcc << "_FromAcc";
  	else
  	  EvAndAcc << "_FromBias";
      }
    }
    else 
      EvAndAcc << "";
    
    if(!fUseConfigTarget)
      ConfigTarget << "_CustomTarget_" << fUserLambda << "_" << fUserBeta;
    else
      ConfigTarget << "";
    
    if(!fUseConfigMapParams)    
      ConfigTarget << "_MapExt_X" << fExtX << "_Y" << fExtY;
    else
      ConfigTarget << "";
    
    if(fAdaptFFT){
      RingMethod << "_AdaptRing";
      if(fAdaptCut_Alpha)
  	RingMethod << "_AlphaCut" << fRingParam_AlphaMax;
      else
  	RingMethod << "_AreaCut" << fRingParam_ExcFracMax;
      
      if(fConstantArea)
  	RingMethod << "_CstArea" << fRingParam_AreaOverPi << "_RInMin" << fRingMinimalRadius << "_RInMax" << fInnerRingMax << "_RStep" << fRingStep;
      if(fConstantThickness)
  	RingMethod << "_CstThick" << fRingParam_Thickness << "_RInMin" << fRingMinimalRadius << "_RInMax" << fInnerRingMax << "_RStep" << fRingStep;
      else//if(fConstantInnerRadius) 
  	RingMethod << "_CstInner" << fRingMinimalRadius << "_ROutMin" << fRingMinimalRadius+fRingParam_Thickness << "_ROutMax" << fOuterRingMax << "_RStep" << fRingStep;
    }
    else
      RingMethod << "_StdRing_R" << fStandardRingRadius << "_Thick" << fStandardRingThickness;
    
    if(fEMin != 0 && fEMax == -1)
      EnergyCuts << "_EMin" << fEMin;
    else if(fEMin == 0 && fEMax != -1)
      EnergyCuts << "_EMax" << fEMax;
    else if(fEMin != 0 && fEMax != -1)
      EnergyCuts << "_EMin" << fEMin << "_EMax" << fEMax;
    else
      EnergyCuts << "";
    
    foutfile << OutfilePrefix << "_ONOFFTestMaps.root";
  }

  TFile fileResultsOut(foutfile.str().c_str(),"RECREATE");
  //

  std::vector<int>::iterator it = fEntriesVector.begin();
  for(; it != fEntriesVector.end(); ++it)
    {
      int i = *it;
      run_results->GetEntry(i);
      const ParisAnalysis::DataStorageRun* dsrun = hess->Handle("GammaFOVStorage", (ParisAnalysis::DataStorageRun *) 0);
     
      CreateRingBgMaps_OneEntry(AnalysisConfig, Resfile, Evtfile, Accfile, ExposureMapFile, ExtendedExposureMapFile, Regionsfile, table_path, OutfilePrefix, i);
      // EventsPerRuns[i] = GetEventsMap();
      // AccPerRuns[i] = GetAcceptanceMap();
      // ExcessPerRuns[i] = GetExcessMap();
      // SignificancePerRuns[i] = GetSignificanceMap();
      // last = i;      

      Display::SkyHistogram2D * hEvents = (Display::SkyHistogram2D *)GetEventsMap()->Clone();    
      Display::SkyHistogram2D * hAcceptance = (Display::SkyHistogram2D *)GetAcceptanceMap()->Clone();    
      Display::SkyHistogram2D * hExcess = (Display::SkyHistogram2D *)GetExcessMap()->Clone();    
      Display::SkyHistogram2D * hSignificance = (Display::SkyHistogram2D *)GetSignificanceMap()->Clone();    

      std::ostringstream fname_Events;
      fname_Events << "Events_" << i;
      std::ostringstream fname_Acceptance;
      fname_Acceptance << "Acceptance_" << i;
      std::ostringstream fname_Excess;
      fname_Excess << "Excess_" << i;
      std::ostringstream fname_Significance;
      fname_Significance << "Significance_" << i;

      hEvents->SetNameTitle(fname_Events.str().c_str(),fname_Events.str().c_str());
      hAcceptance->SetNameTitle(fname_Acceptance.str().c_str(),fname_Acceptance.str().c_str());
      hExcess->SetNameTitle(fname_Excess.str().c_str(),fname_Excess.str().c_str());
      hSignificance->SetNameTitle(fname_Significance.str().c_str(),fname_Significance.str().c_str());
      
      fileResultsOut.Append(hEvents);
      fileResultsOut.Append(hAcceptance);
      fileResultsOut.Append(hExcess);
      fileResultsOut.Append(hSignificance);
      
      if(GetEventsMap()->Integral() == 0){
	continue;
      }

      for(int ii = 1; ii <= GetEventsMap()->GetNbinsX(); ++ii)
	{
	  for(int jj = 1; jj <= GetEventsMap()->GetNbinsY(); ++jj)
	    {
	      Stash::Coordinate binpos =  GetEventsMap()->GetBinCoordinate(ii,jj);
	      double RunNumber = dsrun->GetRunNumber();
	      double Events = GetEventsMap()->GetCircularIntegral(binpos,Stash::Lambda(fOSRadius,Stash::Angle::Degrees));
	      double Exposure = GetAcceptanceMap()->GetCircularIntegral(binpos,Stash::Lambda(fOSRadius,Stash::Angle::Degrees));//GetBinContent(ii,jj);
	      double MJD_Start = dsrun->GetFirstEventTime().GetUTC().GetModifiedJulianDate();
	      double MJD_End = dsrun->GetLastEventTime().GetUTC().GetModifiedJulianDate();
	      
	      //std::cout << "Run " << RunNumber << std::setprecision(8) << ":  (" << MJD_Start << ", " << MJD_End << ") -> " << Exposure << " " << " +/- " << Events << std::endl; 
	      
	      if(Exposure > 0){
		vDate[std::make_pair(ii,jj)][count] = (MJD_Start+MJD_End)*0.5;
		vAccs[std::make_pair(ii,jj)][count] = Exposure;
		vEvents[std::make_pair(ii,jj)][count] = Events;
		vAllowed[std::make_pair(ii,jj)][count] = 1;
		++count;
	      }
	    }
	}
      // delete hEvents;
      // delete hAcceptance;
      // delete hExcess;
      // delete hSignificance;

    }

  Display::SkyHistogram2D *hEvents = GetEventsMap();//EventsPerRuns[last];
  Display::SkyHistogram2D *hONOFFTestMap = new Display::SkyHistogram2D("ONOFFTest","ONOFFTest",*hEvents);
  Display::SkyHistogram2D *hONOFFTestMap_Post = new Display::SkyHistogram2D("ONOFFTest_Post","ONOFFTest_Post",*hEvents);
  int nx = hEvents->GetNbinsX();
  int last_step = 0;  
for(int ii = 1; ii <= hEvents->GetNbinsX(); ++ii)
    {
      if(((ii-1)*20 / (nx-1)) > last_step) {
	last_step = ((ii-1)*20 / (nx-1));
	std::cout << "i>" << ii - 1 << " (-> " << nx - 1 << ")" << std::endl;
	std::cout << "." << std::flush;
      }
      
      for(int jj = 1; jj <= hEvents->GetNbinsY(); ++jj)
	{

	  std::map<int,double> vDate_In = vDate[std::make_pair(ii,jj)];
	  std::map<int,double> vAccs_In = vAccs[std::make_pair(ii,jj)];
	  std::map<int,double> vEvents_In = vEvents[std::make_pair(ii,jj)];
	  std::map<int,double> vAllowed_In = vAllowed[std::make_pair(ii,jj)];
  

	  //DO THE TEST
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
    
	    std::map<int,double>::iterator itD = vDate_In.begin();
	    std::map<int,double>::iterator itA = vAccs_In.begin();
	    std::map<int,double>::iterator itE = vEvents_In.begin();
	    std::map<int,double>::iterator itX = vAllowed_In.begin();
    
	    for(; itD != vDate_In.end(); ++itD, ++itA, ++itE, ++itX){
    
	      //      std::cout << itD->first << " " << itD->second << std::endl;
    
	      double this_entry = itD->first;
	      vON[itD->first]= itE->second;
	      vAlphaON[itD->first]= itA->second;
	      std::map<int,double>::iterator itDIn = vDate_In.begin();
	      std::map<int,double>::iterator itAIn = vAccs_In.begin();
	      std::map<int,double>::iterator itEIn = vEvents_In.begin();
	      std::map<int,double>::iterator itXIn = vAllowed_In.begin();
	      double sumAcc = 0;
	      double sumEv = 0;
    
	      for(; itDIn != vDate_In.end(); ++itDIn, ++itAIn, ++itEIn, ++itXIn){
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
	    std::map<int,double>::iterator itX2 = vAllowed_In.begin();

	    for(; itON != vON.end(); ++itON, ++itOFF, ++itAON, ++itAOFF, ++itX2){
	      int this_entry = itON->first;
	      double non = itON->second;
	      double noff = itOFF->second;
	      double alpha = itAON->second*1./itAOFF->second;
	      double sig = Utilities::Statistics::Significance(non,noff,alpha);

	      if(sig > 4. && itX2->second == 1){
		vAllowed_In[itX2->first-2] = 0;
		vAllowed_In[itX2->first-1] = 0;
		vAllowed_In[itX2->first] = 0;
		vAllowed_In[itX2->first+1] = 0;
		vAllowed_In[itX2->first+2] = 0;
		
		++countchanges;
	      }
	      //    std::cout << non-alpha*noff << " " << sig << " " << countchanges << std::endl;
	    }
	    //std::cout << iterations << " " << countchanges << std::endl;
	    changes = countchanges;
	    ++iterations;
	  }

	  //std::cout << "Size " << vON.size() << std::endl;
	  double vecsize = vON.size();
	  std::map<int,double>::iterator itON = vON.begin();
	  std::map<int,double>::iterator itOFF = vOFF.begin();
	  std::map<int,double>::iterator itAON = vAlphaON.begin();
	  std::map<int,double>::iterator itAOFF = vAlphaOFF.begin();
	  std::map<int,double>::iterator itX2 = vAllowed_In.begin();  
	  std::map<int,double>::iterator itONBegin = vON.begin();
	  std::map<int,double>::iterator itONEnd = vON.end();
  
	  double maxsig = 0;
	  for(; itON != vON.end(); ++itON, ++itOFF, ++itAON, ++itAOFF, ++itX2){
	    int this_entry = itON->first;
	    double non = itON->second;
	    double noff = itOFF->second;
	    double alpha = itAON->second*1./itAOFF->second;
	    double sig = Utilities::Statistics::Significance(non,noff,alpha);
	    double allowed = itX2->second;
	    //    std::cout << non-alpha*noff << " " << sig  << std::endl;
	    hSignificanceDistribution->Fill(sig);
	    if(allowed)
	      hSignificanceDistribution_Allowed->Fill(sig);

	    if(sig > maxsig)
	      maxsig = sig;

	  }

	  double maxsig_posttrials = 0;
	  if(maxsig != 0){
	    if(maxsig > 8.3)
	      maxsig_posttrials = TMath::ErfcInverse(TMath::Erfc(maxsig*1./sqrt(2))*vecsize)*sqrt(2);
	    else
	      maxsig_posttrials = TMath::ErfcInverse(1-pow(1-TMath::Erfc(maxsig*1./sqrt(2)),vecsize))*sqrt(2);
	  }
	  
	  hONOFFTestMap->SetBinContent(ii,jj,maxsig);
	  hONOFFTestMap_Post->SetBinContent(ii,jj,maxsig_posttrials);

	}
    }

//if(!SaveResults){
    TCanvas *c = new TCanvas("ONOFF","ONOFF",1000,500);
    c->Divide(2,1);
    c->cd(1);
    hONOFFTestMap->Draw("colz");
    c->cd(2);
    hONOFFTestMap_Post->Draw("colz");
    //}
  std::ostringstream canvfile;
  canvfile << OutfilePrefix << "_ONOFF.png";
  c->Print(canvfile.str().c_str());

  fileResultsOut.Append(hONOFFTestMap);
  fileResultsOut.Append(hONOFFTestMap_Post);
  fileResultsOut.Append(hSignificanceDistribution);
  fileResultsOut.Append(hSignificanceDistribution_Allowed);
      
  fileResultsOut.Write();
  fileResultsOut.Clear();
}

/*
 * Create events and acceptance maps from radial lookups
 *
 */
bool SurveySuite::MapMaker::CreateEventsAndAcceptanceMaps_FromRadialLookups_OneEntry(const char *RadialConfig, const char *Resfile, const char *Regionsfile, std::string table_path, const char *OutfilePrefix, int EntryNumber)
{
  double fSourceSize = 0.;
  if(!GetPointLikeFluxMaps()){
    fSourceSize = 1.;
    std::cout << "Extended tables will be taken!" << std::endl;
  }

  std::cout << RadialConfig << std::endl;
  //
  TFile *fileResults = TFile::Open(Resfile);
  gROOT->cd();
  Sash::DataSet *results = (Sash::DataSet *)fileResults->Get("results");
  results->GetEntry(0);
  
  Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();
  ParisAnalysis::AcceptanceMap *AccTest =  hess->Handle("", (ParisAnalysis::AcceptanceMap *) 0);
  AccTest->LoadAllMembers();
  ParisAnalysis::AnalysisConfig *Config =  hess->Handle("", (ParisAnalysis::AnalysisConfig *) 0);
  Config->LoadAllMembers();
  
  ParisAnalysis::RadialAcceptance *Radial =  hess->Handle("", (ParisAnalysis::RadialAcceptance *) 0);
  Radial->LoadAllMembers();

  const TH1F * hRadialGammaAcceptance = Radial->GethRadialGammaAcceptance();
  const Display::Histogram2D * hRadialGammaAcceptanceZen = Radial->GethRadialGammaAcceptanceZen();

  Sash::DataSet *run_results = (Sash::DataSet *)fileResults->Get("run_results");
  gROOT->cd();
  
  if(fUseConfigTarget){
    std::cout << "Using target position from Results File" << std::endl;
    fUserLambda = Config->GetTargetRa();
    fUserBeta = Config->GetTargetDec();
    std::cout << "fUserLambda = " << fUserLambda << ", fUserBeta = " << fUserBeta << std::endl;
  }

  Stash::Coordinate targetpos(Stash::Lambda(fUserLambda,Stash::Angle::Degrees),
			      Stash::Beta(fUserBeta,Stash::Angle::Degrees),
			      Config->GetSystem());
  
  if(fUseConfigMapParams){
    std::cout << "Using Maps parameters from Results File" << std::endl;
    fBinSize = Config->GetMapBinSize().GetDegrees();
    fExtX  = Config->GetMapExtensionX().GetDegrees();
    fExtY  = Config->GetMapExtensionY().GetDegrees();
    std::cout << "fBinSize = " << fBinSize << ", fExtX = " << fExtX << ", fExtY = " << fExtY << std::endl;
  }

  Double_t ExtX_Run = fPsiCut;
  Double_t ExtY_Run = fPsiCut;

  //Global Maps ...
  Display::SkyHistogram2D *hAcceptanceMap = new Display::SkyHistogram2D("SkyAcceptance","SkyAcceptance",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
  Display::SkyHistogram2D *hEventsMap = new Display::SkyHistogram2D("SkyEvents","SkyEvents",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
  Display::SkyHistogram2D *hOffAxisMap = new Display::SkyHistogram2D("OffAxisMap","OffAxisMap",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
  Display::SkyHistogram2D *hZenithAngleMap = new Display::SkyHistogram2D("ZenithAngleMap","ZenithAngleMap",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
  Display::SkyHistogram2D *hMuonEfficiencyMap = new Display::SkyHistogram2D("MuonEfficiencyMap","MuonEfficiencyMap",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
  Display::SkyHistogram2D *hMinSafeThresholdMap = new Display::SkyHistogram2D("MinSafeThresholdMap","MinSafeThresholdMap",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
  Display::SkyHistogram2D *hAveragedSafeThresholdMap = new Display::SkyHistogram2D("AveragedSafeThresholdMap","AveragedSafeThresholdMap",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
  Display::SkyHistogram2D *hAveragedLiveTimeMap = new Display::SkyHistogram2D("AveragedLiveTimeMap","AveragedLiveTimeMap",targetpos,Stash::Lambda(fExtX,Stash::Angle::Degrees),Stash::Lambda(fExtY,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));

  //Create Exlusion Regions
  CreateExclusionMask(Resfile, Regionsfile);
  
  //LOOP OVER RUNS TO COMPUTE ACCEPTANCE MAP FROM RADIAL LOOKUPS
  int skippedruns = 0;
  double maxExtension = sqrt(pow(fExtX,2)+pow(fExtY,2));
  double RunSelectionCut = 0;
  if(fSelectAllRunsContributingToTheMap)
    RunSelectionCut = fPsiCut+maxExtension;
  else
    RunSelectionCut = fPsiCut;

  //  TH1F *hEnergyDist = new TH1F("hE","hE",10000,0,10); 
  TH1D *hRTemp = new TH1D("hRTemp","hRTemp",100,0,9);    
  TH1D *hRTempFit = new TH1D("hRTempFit","hRTempFit",100,0,9);    
 
  //  Display::SkyHistogram2D *hRunZenMap = 0;
  std::vector<int>::iterator it = fEntriesVector.begin();

  //  for(int i = 0; i < run_results->GetEntries(); ++i)
  if(true)//for(; it != fEntriesVector.end(); ++it)
    {
      int i = EntryNumber;//*it;
      run_results->GetEntry(i);

      if(fVerbose)
	std::cout << i+1 << "/" << run_results->GetEntries() << std::endl;
      
      const ParisAnalysis::DataStorageRun* dsrun = hess->Handle("GammaFOVStorage", (ParisAnalysis::DataStorageRun *) 0);
      Stash::Coordinate pointingpos = dsrun->GetObsPos();      
      
      if(pointingpos.GetAngularDistance(targetpos).GetDegrees() > RunSelectionCut){
	//	std::cout << "skipping run" << std::endl;
	std::cout << Utilities::TextStyle::Yellow() <<  "SKIPPING RUN!!!" << Utilities::TextStyle::Reset() <<  std::endl;
		
	++skippedruns;
	return 0;
      }

      ///SURVEY XCHECK
      /*if(dsrun->GetRunNumber() > 75227){
	++skippedruns;
	continue;
	}*/

      Display::SkyHistogram2D *hRunExclusionMap = new Display::SkyHistogram2D("RunExclusion","RunExclusion",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
      
      Display::SkyHistogram2D *hRunAcceptanceMap = new Display::SkyHistogram2D("RunAcceptance","RunAcceptance",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
      
      Display::SkyHistogram2D *hRunAcceptanceMapCut = new Display::SkyHistogram2D("RunAcceptanceCut","RunAcceptanceCut",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
      
      Display::SkyHistogram2D *hRunEventsMap = new Display::SkyHistogram2D("RunEvents","RunEvents",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));

      Display::SkyHistogram2D *hRunEventsMapCut = new Display::SkyHistogram2D("RunEventsCut","RunEventsCut",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));

      Display::SkyHistogram2D *hRunOffAxisMap = new Display::SkyHistogram2D("RunOffAxis","RunOffAxis",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
      Display::SkyHistogram2D *hRunZenithAngleMap = new Display::SkyHistogram2D("RunZenithAngleMap","RunZenithAngleMap",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
      Display::SkyHistogram2D *hRunMuonEfficiencyMap = new Display::SkyHistogram2D("RunMuonEfficiencyMap","RunMuonEfficiencyMap",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
      
      Display::SkyHistogram2D *hRunSafeThresholdMap = new Display::SkyHistogram2D("RunSafeThresholdMap","RunSafeThresholdMap",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
      Display::SkyHistogram2D *hRunWeightedSafeThresholdMap = new Display::SkyHistogram2D("RunWeightedSafeThresholdMap","RunWeightedSafeThresholdMap",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));


      //define pointing conditions parameters
      if(fVerbose)
	std::cout << dsrun->GetMeanZenith() << " " << (dsrun->GetTelsInRun()).GetPattern() << std::endl;

      double MeanZen = dsrun->GetMeanZenith();
      double meanCosZen = cos(MeanZen*TMath::Pi()*1./180);
      double MuonEff = dsrun->GetMuonEfficiency();
      double offaxisAngle = fPsiCut;
      if(!fSafeThresholdFromPsiCut)
	offaxisAngle = fSafeThresholdFromPsiValue;
      double MeanAzimuth = dsrun->GetMeanAzimuth();
      Int_t NTels = dsrun->GetTelsInRun().size();
      Double_t SafeThreshold = 0;

      TF1 *fRTemp = 0; 
      for(int i = 1; i <= 100; ++i){
	hRTemp->SetBinContent(i,hRadialGammaAcceptanceZen->GetBinContent(hRadialGammaAcceptanceZen->GetXaxis()->FindFixBin(hRTemp->GetBinCenter(i)), hRadialGammaAcceptanceZen->GetYaxis()->FindFixBin(meanCosZen)));
      }

      fRTemp = new TF1("fRTemp","pol8",0,9);
      hRTemp->Fit(fRTemp,"QN");
      
      for(int i = 1; i <= 100; ++i){
	hRTempFit->SetBinContent(i,fRTemp->Eval(hRTemp->GetBinCenter(i)));      
      }
      hRTempFit->Scale(1./hRTempFit->GetMaximum());
      
      delete fRTemp;

      if(fVerbose)
	std::cout << hRTemp->GetMaximum() << std::endl;
      double lastl = 0;
      double lastb = 0;
      int countdoubles = 0;
      int countdoubles2 = 0;
      
      for(std::vector<ParisAnalysis::EventData>::const_iterator evit = dsrun->GetEvents().begin(); evit != dsrun->GetEvents().end() ; evit++) 
	{
	  if(evit->GetShowerPosX() != lastl){
	    if(evit->GetShowerPosY() != lastb){
	      Stash::Coordinate pop(Stash::Lambda(evit->GetShowerPosX(),Stash::Angle::Degrees),
				    Stash::Beta(evit->GetShowerPosY(),Stash::Angle::Degrees),
				    Config->GetSystem());
	      double LambdaMod = pop.GetLambda().GetDegrees();
	      double BetaMod = pop.GetBeta().GetDegrees();
	      //std::cout << evit->GetShowerNomPosX() << " " << evit->GetShowerNomPosY() << std::endl;
	      if(LambdaMod > 180)
		LambdaMod = LambdaMod-360;
	      int a = hRunEventsMap->GetXaxis()->FindFixBin(-LambdaMod);
	      if(a == 0 || a > hRunEventsMap->GetNbinsX())
		a = hRunEventsMap->GetXaxis()->FindFixBin(-(LambdaMod+360));
	      if(a == 0 || a > hRunEventsMap->GetNbinsX())
		continue;
	      int aa = hEventsMap->GetXaxis()->FindFixBin(-LambdaMod);
	      if(aa == 0 || aa > hEventsMap->GetNbinsX())
		aa = hEventsMap->GetXaxis()->FindFixBin(-(LambdaMod+360));
	      if(aa == 0 || aa > hEventsMap->GetNbinsX())
		continue;
	      
	      int b = hRunEventsMap->GetYaxis()->FindFixBin(BetaMod);
	      int bb = hEventsMap->GetYaxis()->FindFixBin(BetaMod);
	      
	      if(pop.GetAngularDistance(pointingpos).GetDegrees() < fPsiCut){		
		hRunEventsMap->SetBinContent(a,b,1+hRunEventsMap->GetBinContent(a,b));
		hEventsMap->SetBinContent(aa,bb,1+hEventsMap->GetBinContent(aa,bb));
		
		if(GetExclusionMap()->GetBinContent(aa,bb) == 0)
		  hRunEventsMapCut->SetBinContent(a,b,1+hRunEventsMapCut->GetBinContent(a,b));
	      }
	      ++countdoubles;
	    }
	    lastl = evit->GetShowerPosX();
	    lastb = evit->GetShowerPosY();
	    ++countdoubles2;
	  }
	}      
      
      for(int k = 1; k <= hRunAcceptanceMap->GetNbinsX(); ++k)
	{
       	  for(int l = 1; l <= hRunAcceptanceMap->GetNbinsY(); ++l)
     	    {
       	      Stash::Coordinate pop = hRunAcceptanceMap->GetBinCoordinate(k,l);	
	      double LambdaMod = pop.GetLambda().GetDegrees();
	      double BetaMod = pop.GetBeta().GetDegrees();

	      if(LambdaMod > 180)
		LambdaMod = LambdaMod-360;
	      int a = GetExclusionMap()->GetXaxis()->FindFixBin(-LambdaMod);
	      if(a == 0 || a > GetExclusionMap()->GetNbinsX())
	      	a = GetExclusionMap()->GetXaxis()->FindFixBin(-(LambdaMod+360));
	      if(a == 0 || a > GetExclusionMap()->GetNbinsX())
		continue;
	      int aa = hRunAcceptanceMap->GetXaxis()->FindFixBin(-LambdaMod);
	      if(aa == 0 || aa > hRunAcceptanceMap->GetNbinsX())
	      	aa = hRunAcceptanceMap->GetXaxis()->FindFixBin(-(LambdaMod+360));
	      if(aa == 0 || aa > hRunAcceptanceMap->GetNbinsX())
		continue;

	      int b = GetExclusionMap()->GetYaxis()->FindFixBin(BetaMod);
	      int bb = hRunAcceptanceMap->GetYaxis()->FindFixBin(BetaMod);

	      double val = 0;
	      if(!fFitAcceptance)
		val = hRTemp->Interpolate(pow(pop.GetAngularDistance(pointingpos).GetDegrees(),2));
	      else
		val = hRTempFit->Interpolate(pow(pop.GetAngularDistance(pointingpos).GetDegrees(),2));

	      if(pop.GetAngularDistance(pointingpos).GetDegrees() < fPsiCut){
		hRunAcceptanceMap->SetBinContent(aa,bb,val+hRunAcceptanceMap->GetBinContent(aa,bb));
		if(GetExclusionMap()->GetBinContent(a,b) == 0)
		  hRunAcceptanceMapCut->SetBinContent(aa,bb,val+hRunAcceptanceMapCut->GetBinContent(aa,bb));
		
	      }
       	    }
      	}
      for(int k = 1; k <= hRunExclusionMap->GetNbinsX(); ++k)
	{
       	  for(int l = 1; l <= hRunExclusionMap->GetNbinsY(); ++l)
     	    {
       	      Stash::Coordinate pop = hRunExclusionMap->GetBinCoordinate(k,l);
	      double LambdaMod = pop.GetLambda().GetDegrees();
	      double BetaMod = pop.GetBeta().GetDegrees();
	      if(LambdaMod > 180)
		LambdaMod = LambdaMod-360;
	      int a = GetExclusionMap()->GetXaxis()->FindFixBin(-LambdaMod);
	      if(a == 0 || a > GetExclusionMap()->GetNbinsX())
	      	a = GetExclusionMap()->GetXaxis()->FindFixBin(-(LambdaMod+360));
	      if(a == 0 || a > GetExclusionMap()->GetNbinsX())
		continue;

	      int aa = hRunExclusionMap->GetXaxis()->FindFixBin(-LambdaMod);
	      if(aa == 0 || aa > hRunExclusionMap->GetNbinsX())
	      	aa = hRunExclusionMap->GetXaxis()->FindFixBin(-(LambdaMod+360));
	      if(aa == 0 || aa > hRunExclusionMap->GetNbinsX())
		continue;

	      int b = GetExclusionMap()->GetYaxis()->FindFixBin(BetaMod);
	      int bb = hRunExclusionMap->GetYaxis()->FindFixBin(BetaMod);
	      
	      if(pop.GetAngularDistance(pointingpos).GetDegrees() < fPsiCut){	
		if(GetExclusionMap()->GetBinContent(a,b) == 0)
		  hRunExclusionMap->SetBinContent(aa,bb,1);
		else
		  hRunExclusionMap->SetBinContent(aa,bb,0); 
	      }
	    }
	}
      if(fVerbose)
	std::cout << "Allowed Integral : events = " << hRunEventsMapCut->Integral() << ", Acceptance = " << hRunAcceptanceMapCut->Integral() << std::endl;
      if(TMath::IsNaN(hRunAcceptanceMapCut->Integral()) || hRunAcceptanceMapCut->Integral() == 0)
	{
	  std::cout << Utilities::TextStyle::Yellow() <<  "SKIPPING RUN (NAN in ACC)!!!" << Utilities::TextStyle::Reset() <<  std::endl;

	  delete hRunExclusionMap;
	  delete hRunAcceptanceMap;
	  delete hRunAcceptanceMapCut;
	  delete hRunEventsMap;
	  delete hRunEventsMapCut;
	  delete hRunOffAxisMap;
	  delete hRunZenithAngleMap;
	  delete hRunMuonEfficiencyMap;
	  delete hRunSafeThresholdMap;
	  delete hRunWeightedSafeThresholdMap;
	  hRTemp->Reset();
	  hRTempFit->Reset();
	  ++skippedruns;
	  return 0;
	}
      double scalefactor = hRunEventsMapCut->Integral()*1./hRunAcceptanceMapCut->Integral();
      if(fVerbose){
	std::cout << "ScaleFactor = " << scalefactor << std::endl;
	std::cout << "Max Acceptance = " << hRunAcceptanceMap->GetMaximum() << std::endl;      
      }

      hRunAcceptanceMap->Scale(scalefactor);

      if(fVerbose)
	std::cout << "Scaled Max Acceptance = " << hRunAcceptanceMap->GetMaximum() << std::endl;      
      //Filling scale factor per run map
      fRadialIntegrationEbin[dsrun->GetRunNumber()]=std::make_pair(0,-1);//binemin,binemax);      
      fScaleFactor[dsrun->GetRunNumber()]=scalefactor;

      //FILLING POINTING CONDITIONS MAP
      for(int k = 1; k <= hRunAcceptanceMap->GetNbinsX(); ++k)
	{
       	  for(int l = 1; l <= hRunAcceptanceMap->GetNbinsY(); ++l)
     	    {
       	      Stash::Coordinate pop = hRunAcceptanceMap->GetBinCoordinate(k,l);	
	      double LambdaMod = pop.GetLambda().GetDegrees();
	      double BetaMod = pop.GetBeta().GetDegrees();

	      if(LambdaMod > 180)
		LambdaMod = LambdaMod-360;
	      int a = hRunOffAxisMap->GetXaxis()->FindFixBin(-LambdaMod);
	      if(a == 0 || a > hRunOffAxisMap->GetNbinsX())
	      	a = hRunOffAxisMap->GetXaxis()->FindFixBin(-(LambdaMod+360));
	      if(a == 0 || a > hRunOffAxisMap->GetNbinsX())
		continue;
	      int aa = hRunAcceptanceMap->GetXaxis()->FindFixBin(-LambdaMod);
	      if(aa == 0 || aa > hRunAcceptanceMap->GetNbinsX())
	      	aa = hRunAcceptanceMap->GetXaxis()->FindFixBin(-(LambdaMod+360));
	      if(aa == 0 || aa > hRunAcceptanceMap->GetNbinsX())
		continue;

	      int b = hRunOffAxisMap->GetYaxis()->FindFixBin(BetaMod);
	      int bb = hRunAcceptanceMap->GetYaxis()->FindFixBin(BetaMod);

	      if(pop.GetAngularDistance(pointingpos).GetDegrees() < fPsiCut){
		hRunOffAxisMap->SetBinContent(a,b,pop.GetAngularDistance(pointingpos).GetDegrees()*hRunAcceptanceMap->GetBinContent(aa,bb)+hRunOffAxisMap->GetBinContent(a,b));
		hRunZenithAngleMap->SetBinContent(a,b,MeanZen*hRunAcceptanceMap->GetBinContent(aa,bb)+hRunZenithAngleMap->GetBinContent(a,b));
		hRunMuonEfficiencyMap->SetBinContent(a,b,MuonEff*hRunAcceptanceMap->GetBinContent(aa,bb)+hRunMuonEfficiencyMap->GetBinContent(a,b));
		hRunSafeThresholdMap->SetBinContent(a,b,SafeThreshold+hRunSafeThresholdMap->GetBinContent(a,b));
		hRunWeightedSafeThresholdMap->SetBinContent(a,b,SafeThreshold*hRunAcceptanceMap->GetBinContent(aa,bb)+hRunWeightedSafeThresholdMap->GetBinContent(a,b));
	      }
       	    }
      	}
      

      //GO TO NOMINAL SYSTEM AND CORRECT FOR ZENITH GRADIENTS...
      if(fCorrectZenith){
	Display::SkyHistogram2D *hRunFOVAcceptanceMap = new Display::SkyHistogram2D("RunFOVAcceptance","RunFOVAcceptance",Stash::Coordinate(0,0,hess->GetNominalSystem()),Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
	Display::SkyHistogram2D *hRunFOVEventsMap = new Display::SkyHistogram2D("RunFOVEvents","RunFOVEvents",Stash::Coordinate(0,0,hess->GetNominalSystem()),Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
	Display::SkyHistogram2D *hRunFOVExclusionMap = new Display::SkyHistogram2D("RunFOVExclusion","RunFOVExclusion",Stash::Coordinate(0,0,hess->GetNominalSystem()),Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));
;
 
	TH1F *hX = new TH1F("hx","hx",hRunFOVEventsMap->GetNbinsX(),hRunFOVEventsMap->GetXaxis()->GetBinLowEdge(1),hRunFOVEventsMap->GetXaxis()->GetBinLowEdge(1+hRunFOVEventsMap->GetNbinsX()));
	TH1F *hY = new TH1F("hy","hy",hRunFOVEventsMap->GetNbinsY(),hRunFOVEventsMap->GetYaxis()->GetBinLowEdge(1),hRunFOVEventsMap->GetYaxis()->GetBinLowEdge(1+hRunFOVEventsMap->GetNbinsY()));

	UInt_t nStep = 30;
	Float_t timeStep = (dsrun->GetLastEventTime() - dsrun->GetFirstEventTime()) / nStep;
	Float_t timeOffset = timeStep/2;
	
	for(UInt_t iStep = 0; iStep < nStep; ++iStep)
	  {
	    hess->SetCurrentTime(dsrun->GetFirstEventTime() + iStep * timeStep + timeOffset);    

	    Sash::NominalPointing *nomptg = hess->Handle<Sash::NominalPointing>();
	    nomptg->SetReferenceDirection(pointingpos);
	    for(Int_t fov_xbin = 1; fov_xbin <= hRunFOVAcceptanceMap->GetNbinsX(); ++fov_xbin)
	      {
		for(Int_t fov_ybin = 1; fov_ybin <= hRunFOVAcceptanceMap->GetNbinsY(); ++fov_ybin)
		  {
		    Int_t fov_bin = hRunFOVAcceptanceMap->GetBin(fov_xbin,fov_ybin);
		    // Compute bin position in nominal system
		    Stash::Coordinate nompos = hRunFOVAcceptanceMap->GetBinCoordinate(fov_xbin, fov_ybin);

		    Crash::SetRefractionType(*hess->GetHorizonSystem(),Crash::HorizonSystem::None);
		    Stash::Coordinate pos = nompos.GetCoordinate(*Config->GetSystem());
		    Int_t sky_bin = hRunAcceptanceMap->FindBinPosition(pos);
		    Int_t xb,yb,zb;
		    hRunAcceptanceMap->GetBinXYZ(sky_bin,xb,yb,zb);
		    Float_t area = hRunAcceptanceMap->GetBinSolidAngle(xb,yb);
		    hRunFOVAcceptanceMap->SetBinContent(fov_bin, hRunFOVAcceptanceMap->GetBinContent(fov_bin)+hRunAcceptanceMap->GetBinContent(sky_bin)*1./area);
		    hRunFOVEventsMap->SetBinContent(fov_bin, hRunFOVEventsMap->GetBinContent(fov_bin)+hRunEventsMap->GetBinContent(sky_bin)*1./area);
		    hRunFOVExclusionMap->SetBinContent(fov_bin, hRunFOVExclusionMap->GetBinContent(fov_bin)+hRunExclusionMap->GetBinContent(sky_bin));
		  }
	      }
	  }
      
	for(int ii = 10; ii <= hRunFOVEventsMap->GetNbinsY()-10; ++ii)
	  {
	    double sumb = 0;
	    double sumn = 0; 
	    double suma = 0; 
	    for(int jj = 10; jj <= hRunFOVEventsMap->GetNbinsX()-10; ++jj)
	      {
		if(hRunFOVExclusionMap->GetBinContent(jj,ii)>0){
		  sumb += hRunFOVExclusionMap->GetBinContent(jj,ii);
		  sumn += hRunFOVEventsMap->GetBinContent(jj,ii)*hRunFOVExclusionMap->GetBinContent(jj,ii)*1./30;
		  suma += hRunFOVAcceptanceMap->GetBinContent(jj,ii)*hRunFOVExclusionMap->GetBinContent(jj,ii)*1./30;
		}
	      }
	    if(sumb > 0)
	      hY->SetBinContent(ii,sumn*1./suma);
	  }
	
	for(int ii = 10; ii <= hRunFOVEventsMap->GetNbinsX()-10; ++ii)
	  {
	    double sumb = 0;
	    double sumn = 0; 
	    double suma = 0; 
	    for(int jj = 10; jj <= hRunFOVEventsMap->GetNbinsY()-10; ++jj)
	      {
		if(hRunFOVExclusionMap->GetBinContent(ii,jj)>0){
		  sumb += hRunFOVExclusionMap->GetBinContent(ii,jj);
		  sumn += hRunFOVEventsMap->GetBinContent(ii,jj)*hRunFOVExclusionMap->GetBinContent(ii,jj)*1./30;
		  suma += hRunFOVAcceptanceMap->GetBinContent(ii,jj)*hRunFOVExclusionMap->GetBinContent(ii,jj)*1./30;
		}
	      }
	    if(sumb > 0)
	      hX->SetBinContent(ii,sumn*1./suma);
	  }

	double sux = 0,ssx = 0;
	double xmin = -1.*(ExtX_Run-10*fBinSize);
	double xmax = ExtX_Run-10*fBinSize;
	double maxg = 0.5;
	double GradientX = 0;	
	double GradientY = 0;

	Int_t n = hX->GetNbinsX();
	for(Int_t i = 1; i <= n; ++i) {
	  //std::cout << hX->GetBinContent(i) << std::endl;
	  Double_t y1,y2,dy1,dy2;
	  Double_t x = hX->GetBinCenter(i);
	  UInt_t j = hX->GetXaxis()->FindFixBin(-x);
	  if(x < xmin || x > xmax) continue;
	  y1 = hX->GetBinContent(i);
	  y2 = hX->GetBinContent(j);
	  dy1 = hX->GetBinError(i);
	  dy2 = hX->GetBinError(j);
	  if(y1 > 0 &&  y2 > 0 && dy1 > 0 && dy2 > 0) {
	    double U = (y1-y2)/(y1 + y2);
	    double dU2 =  4 *  y1 * y2 / std::pow(y1 + y2,3);
	    sux += U * x / dU2;
	    ssx += x*x / dU2;
	  }
	}
	if(ssx > 0 && fabs(sux/ssx) < maxg) {
	  GradientX = sux/ssx;
	}
	else if(ssx > 0 && sux/ssx >  maxg) {
	  GradientX = maxg;
	}
	if(fVerbose)
	  std::cout << "Gradient X =  " << GradientX << std::endl;
	
	double suy = 0, ssy = 0;
	n = hY->GetNbinsX();
	for(Int_t i = 1; i <= n; ++i) {
	  Double_t y1,y2,dy1,dy2;
	  Double_t x = hY->GetBinCenter(i);
	  UInt_t j = hY->GetXaxis()->FindFixBin(-x);
	  if(x < xmin || x > xmax) continue;
	  y1 = hY->GetBinContent(i);
	  y2 = hY->GetBinContent(j);
	  dy1 = hY->GetBinError(i);
	  dy2 = hY->GetBinError(j);
	  if(y1 > 0 &&  y2 > 0 && dy1 > 0 && dy2 > 0) {
	    double U = (y1-y2)/(y1 + y2);
	    double dU2 =  4 *  y1 * y2 / std::pow(y1 + y2,3);
	    suy += U * x / dU2;
	    ssy += x*x / dU2;
	  }
	}
	if(ssy > 0 && fabs(suy/ssy) < maxg) {
	  GradientY = suy/ssy;
	}
	else if(ssy > 0 && suy/ssy >  maxg) {
	  GradientY = maxg;
	}
	if(fVerbose)
	  std::cout << "Gradient Y =  " << GradientY << std::endl;

	for(int ii = 1; ii <= hRunFOVAcceptanceMap->GetNbinsX(); ++ii)
	  {
	    for(int jj = 1; jj <= hRunFOVAcceptanceMap->GetNbinsY(); ++jj)
	      {
		Float_t x = hRunFOVAcceptanceMap->GetXaxis()->GetBinCenter(ii);
		Float_t y = hRunFOVAcceptanceMap->GetYaxis()->GetBinCenter(jj);
		//		std::cout << x << " " << y << std::endl;
		hRunFOVAcceptanceMap->SetBinContent(ii,jj,hRunFOVAcceptanceMap->GetBinContent(ii,jj)*(1+GradientX*x)*(1+GradientY*y));
	      }
	  }

	//REPROJECT ON THE SKY
	for(UInt_t iStep = 0; iStep < nStep; ++iStep)
	  {
	    hess->SetCurrentTime(dsrun->GetFirstEventTime() + iStep * timeStep + timeOffset);    
	    Sash::NominalPointing *nomptg = hess->Handle<Sash::NominalPointing>();
	    nomptg->SetReferenceDirection(pointingpos);
	    for(Int_t sky_xbin = 1; sky_xbin <= hRunAcceptanceMap->GetNbinsX(); ++sky_xbin)
	      {
		for(Int_t sky_ybin = 1; sky_ybin <= hRunAcceptanceMap->GetNbinsY(); ++sky_ybin)
		  {
		    Int_t sky_bin = hRunAcceptanceMap->GetBin(sky_xbin,sky_ybin);
		    Float_t area = hRunAcceptanceMap->GetBinSolidAngle(sky_xbin,sky_ybin);
		    // Compute bin position in nominal system
		    Stash::Coordinate pos = hRunAcceptanceMap->GetBinCoordinate(sky_xbin, sky_ybin);
		    Crash::SetRefractionType(*hess->GetHorizonSystem(),Crash::HorizonSystem::None);
		    Stash::Coordinate nompos = pos.GetCoordinate(*hess->GetNominalSystem());
		    Int_t fov_bin = hRunFOVAcceptanceMap->FindBinPosition(nompos);
		    hRunAcceptanceMap->SetBinContent(sky_bin,hRunAcceptanceMap->GetBinContent(sky_bin)+hRunFOVAcceptanceMap->GetBinContent(fov_bin)*area);
		  }
	      }
	  }

	Display::SkyHistogram2D *hRunAcceptanceMapCut2 = new Display::SkyHistogram2D("RunAcceptanceCut2","RunAcceptanceCut2",pointingpos,Stash::Lambda(ExtX_Run,Stash::Angle::Degrees),Stash::Lambda(ExtY_Run,Stash::Angle::Degrees),Stash::SmallAngle(fBinSize,Stash::Angle::Degrees));

	for(int k = 1; k <= hRunAcceptanceMap->GetNbinsX(); ++k)
	  {
	    for(int l = 1; l <= hRunAcceptanceMap->GetNbinsY(); ++l)
	      {
		Stash::Coordinate pop = hRunAcceptanceMap->GetBinCoordinate(k,l);	
		if(GetExclusionMap()->GetBinContent(GetExclusionMap()->FindBinPosition(pop)) == 0)
		  hRunAcceptanceMapCut2->SetBinContent(k,l,hRunAcceptanceMap->GetBinContent(k,l));	      
	      }
	  }
	if(fVerbose)
	  std::cout << "Allowed Integral : events = " << hRunEventsMapCut->Integral() << ", Acceptance = " << hRunAcceptanceMapCut2->Integral() << std::endl;
	hRunAcceptanceMap->Scale(hRunEventsMapCut->Integral()*1./hRunAcceptanceMapCut2->Integral());
	fScaleFactor[dsrun->GetRunNumber()]=hRunEventsMapCut->Integral()*1./hRunAcceptanceMapCut2->Integral();

	delete hRunAcceptanceMapCut2;
	delete hRunFOVAcceptanceMap;
	delete hRunFOVExclusionMap;
	delete hRunFOVEventsMap;
	delete hX;
	delete hY;
      }
      //-----------------------
      
      //only loop in the useful area
      int kmin, lmin, kmax, lmax;
      hAcceptanceMap->FindBinPosition(hRunAcceptanceMap->GetBinCoordinate(1,1),kmin,lmin);
      hAcceptanceMap->FindBinPosition(hRunAcceptanceMap->GetBinCoordinate(hRunAcceptanceMap->GetNbinsX(),hRunAcceptanceMap->GetNbinsY()),kmax,lmax);
      if(kmax == 0 || kmax > hAcceptanceMap->GetNbinsX())
	kmax = hAcceptanceMap->GetNbinsX();
      if(kmin == 0 || kmin > hAcceptanceMap->GetNbinsX())
	kmin = 1;
      if(lmax == 0 || lmax > hAcceptanceMap->GetNbinsY())
	lmax = hAcceptanceMap->GetNbinsY();
      if(lmin == 0 || lmin > hAcceptanceMap->GetNbinsY())
	lmin = 1;

      kmin = kmin -1;
      lmin = lmin -1;
      kmax = kmax +1;
      lmax = lmax +1;
      if(fVerbose)
	std::cout << kmin << " " << lmin << " " << kmax << " " << lmax << std::endl;

      double runOnTime = dsrun->GetOnLiveTime() *1./3600;

      for(int k = kmin; k <= kmax; ++k)
	{
       	  for(int l = lmin; l <= lmax; ++l)
     	    {
	      Stash::Coordinate pop = hAcceptanceMap->GetBinCoordinate(k,l);	
	      double val = hAcceptanceMap->GetBinContent(k,l);
	      double valn = hRunAcceptanceMap->GetBinContent(hRunAcceptanceMap->FindBinPosition(pop));
	      hAcceptanceMap->SetBinContent(k,l,val+valn);	      
	      double valoff = hOffAxisMap->GetBinContent(k,l);
	      double valnoff = hRunOffAxisMap->GetBinContent(hRunOffAxisMap->FindBinPosition(pop));
	      hOffAxisMap->SetBinContent(k,l,valoff+valnoff);
	      double valzen = hZenithAngleMap->GetBinContent(k,l);
	      double valnzen = hRunZenithAngleMap->GetBinContent(hRunZenithAngleMap->FindBinPosition(pop));
	      hZenithAngleMap->SetBinContent(k,l,valzen+valnzen);
	      double valeff = hMuonEfficiencyMap->GetBinContent(k,l);
	      double valneff = hRunMuonEfficiencyMap->GetBinContent(hRunMuonEfficiencyMap->FindBinPosition(pop));
	      hMuonEfficiencyMap->SetBinContent(k,l,valeff+valneff);
	     
	      double valth = hMinSafeThresholdMap->GetBinContent(k,l);
	      double valnth = hRunSafeThresholdMap->GetBinContent(hRunSafeThresholdMap->FindBinPosition(pop));
	      if(valth == 0)
		hMinSafeThresholdMap->SetBinContent(k,l,valnth);
	      else if(valnth != 0)
		hMinSafeThresholdMap->SetBinContent(k,l,TMath::Min(valth, valnth));

	      double valwth = hAveragedSafeThresholdMap->GetBinContent(k,l);
	      double valnwth = hRunWeightedSafeThresholdMap->GetBinContent(hRunWeightedSafeThresholdMap->FindBinPosition(pop));
	      hAveragedSafeThresholdMap->SetBinContent(k,l,valwth+valnwth);

	      double valtime = hAveragedLiveTimeMap->GetBinContent(k,l);
	      double valntime = runOnTime;
	      if(valn>0)
		hAveragedLiveTimeMap->SetBinContent(k,l,valtime+valntime);

	    }
	}

      delete hRunExclusionMap;
      delete hRunAcceptanceMap;
      delete hRunAcceptanceMapCut;
      delete hRunEventsMapCut;
      hRTemp->Reset();
      hRTempFit->Reset();
      delete hRunEventsMap;	
      delete hRunOffAxisMap;
      delete hRunZenithAngleMap;
      delete hRunMuonEfficiencyMap;
      delete hRunSafeThresholdMap;
      delete hRunWeightedSafeThresholdMap;

      std::cout << "Everything should be deleted" << std::endl;
    }
  
  double maxaccfinal = hAcceptanceMap->GetMaximum();

  for(int k = 1; k <= hOffAxisMap->GetNbinsX(); ++k)
    {
      for(int l = 1; l <= hOffAxisMap->GetNbinsY(); ++l)
	{
	  double ratio = hAcceptanceMap->GetBinContent(k,l)*1./maxaccfinal;
	  if(hAcceptanceMap->GetBinContent(k,l) != 0){
	    hOffAxisMap->SetBinContent(k,l,hOffAxisMap->GetBinContent(k,l)*1./hAcceptanceMap->GetBinContent(k,l));
	    hZenithAngleMap->SetBinContent(k,l,hZenithAngleMap->GetBinContent(k,l)*1./hAcceptanceMap->GetBinContent(k,l));
	    hMuonEfficiencyMap->SetBinContent(k,l,hMuonEfficiencyMap->GetBinContent(k,l)*1./hAcceptanceMap->GetBinContent(k,l));
	    hAveragedSafeThresholdMap->SetBinContent(k,l,hAveragedSafeThresholdMap->GetBinContent(k,l)*1./hAcceptanceMap->GetBinContent(k,l));
	    hAveragedLiveTimeMap->SetBinContent(k,l,hAveragedLiveTimeMap->GetBinContent(k,l)*ratio);
	  }
	}
    }
  
  //std::cout << fEMin << " - " << fEMax << " -> " << emin << "(bin " << binemin << ")" << " - " << emax << "(bin " << binemax << ")" << std::endl;
  std::cout << "-> " << fEntriesVector.size()-skippedruns << " runs out of " << fEntriesVector.size() << " pre-selected have been used (" << run_results->GetEntries() << " available in the result file)." << std::endl;
  GetAcceptanceMap() = hAcceptanceMap;
  GetEventsMap() = hEventsMap;
  GetOffAxisMap() = hOffAxisMap;
  GetZenithAngleMap() = hZenithAngleMap;
  GetMuonEfficiencyMap() = hMuonEfficiencyMap;
  GetMinSafeThresholdMap() = hMinSafeThresholdMap;
  GetAveragedSafeThresholdMap() = hAveragedSafeThresholdMap;
  GetAveragedLiveTimeMap() = hAveragedLiveTimeMap;

  //delete hAcceptanceMap
  //delete hEventsMap
  delete hOffAxisMap;
  delete hZenithAngleMap;
  delete hMuonEfficiencyMap;
  delete hMinSafeThresholdMap;
  delete hAveragedSafeThresholdMap;
  delete hAveragedLiveTimeMap;

  fileResults->Close();
  //  RadialFile->Close();
  
  return true;
}



bool SurveySuite::MapMaker::CreateRingBgMaps_OneEntry(const char *AnalysisConfig, const char *Resfile, const char *Evtfile, const char *Accfile, const char *ExposureMapFile, const char *ExtendedExposureMapFile, const char *Regionsfile, std::string table_path, const char *OutfilePrefix, int EntryNumber)
{
  TFile *fileResults = TFile::Open(Resfile);
  gROOT->cd();
  Sash::DataSet *results = (Sash::DataSet *)fileResults->Get("results");
  results->GetEntry(0);
  Sash::HESSArray *hess = Sash::Folder::GetFolder(0)->GetHESSArray();
  
  ParisAnalysis::AnalysisConfig *Config =  hess->Handle("", (ParisAnalysis::AnalysisConfig *) 0);
  Config->LoadAllMembers();
  
  if(fUseConfigTarget){
    std::cout << "Using target position from Results File" << std::endl;
    fUserLambda = Config->GetTargetRa();
    fUserBeta = Config->GetTargetDec();
    std::cout << "fUserLambda = " << fUserLambda << ", fUserBeta = " << fUserBeta << std::endl;
  }
  Stash::Coordinate targetpos(Stash::Lambda(fUserLambda,Stash::Angle::Degrees),
			      Stash::Beta(fUserBeta,Stash::Angle::Degrees),
			      Config->GetSystem());
  if(fUseConfigMapParams){
    std::cout << "Using Maps parameters from Results File" << std::endl;
    fBinSize = Config->GetMapBinSize().GetDegrees();
    fExtX  = Config->GetMapExtensionX().GetDegrees();
    fExtY  = Config->GetMapExtensionY().GetDegrees();
    std::cout << "fBinSize = " << fBinSize << ", fExtX = " << fExtX << ", fExtY = " << fExtY << std::endl;
  }

  //Compute excess and significance maps
  Double_t SafeExtX = fExtX + 2.0;
  Double_t SafeExtY = fExtY + 2.0;
  //  Double_t BinSize  = binsize;
  std::cout << fBinSize << std::endl;
  Double_t nSafeBins = 2.0 *1./fBinSize;

  //Events and Acceptance?
  //  CreateEventsAndAcceptanceMaps(AnalysisConfig, Resfile, Regionsfile, Evtfile, Accfile, table_path, OutfilePrefix);
  CreateEventsAndAcceptanceMaps_FromRadialLookups_OneEntry(AnalysisConfig, Resfile, Regionsfile, table_path, OutfilePrefix, EntryNumber);
  //Exclusion regions map
  CreateExclusionMask(Resfile, Regionsfile);
  
  Display::SkyHistogram2D * hGammaMod = (Display::SkyHistogram2D *)GetEventsMap()->Clone();    
  hGammaMod->SetFoV(SafeExtX,SafeExtY);
  
  Display::SkyHistogram2D * hGammaAccMod = (Display::SkyHistogram2D *)hGammaMod->Clone();
  hGammaAccMod->Reset();
  hGammaMod->Reset();
  
  Display::SkyHistogram2D * hExcMod = (Display::SkyHistogram2D *)hGammaMod->Clone();
  hExcMod->Reset();
  
  Display::SkyHistogram2D * hMask = (Display::SkyHistogram2D*)GetEventsMap()->Clone();
  for (int i=1;i<=hMask->GetNbinsX();++i) {
    for (int j=1;j<=hMask->GetNbinsY();++j) {
      hMask->SetBinContent(i,j,1.);
    }
  }
  
  Display::SkyHistogram2D * hMaskMod = (Display::SkyHistogram2D*)GetEventsMap()->Clone();
  hMaskMod->SetFoV(SafeExtX,SafeExtY);
  
  for (int i=1;i<=hMaskMod->GetNbinsX();++i) {
    for (int j=1;j<=hMaskMod->GetNbinsY();++j) {
      if(i > nSafeBins && j > nSafeBins && i <= hMaskMod->GetNbinsX()-nSafeBins && j <= hMaskMod->GetNbinsY()-nSafeBins){
	hMaskMod->SetBinContent(i,j,1.);
	hExcMod->SetBinContent(i,j,0.);
      }
      else{
	hMaskMod->SetBinContent(i,j,0.);
	hExcMod->SetBinContent(i,j,1.);
      }
    }
  }

  int nx0 = GetEventsMap()->GetNbinsX();
  int ny0 = GetEventsMap()->GetNbinsY();
  int nx = hGammaMod->GetNbinsX();
  int ny = hGammaMod->GetNbinsY();    
  int nxa = hGammaAccMod->GetNbinsX();
  int nya = hGammaAccMod->GetNbinsY();
  
  std::cout << nx << " " << ny << std::endl;
  std::cout << nxa << " " << nya << std::endl;
  
  Display::SkyHistogram2D *hGammaModExc = new Display::SkyHistogram2D("GammaModExc","GammaModExc",*hExcMod); 
  Display::SkyHistogram2D *hGammaAccModExc = new Display::SkyHistogram2D("GammaAccModExc","GammaAccModExc",*hGammaAccMod); 
  Display::SkyHistogram2D *hOS = new Display::SkyHistogram2D("OS","OS",*hExcMod);     
  Display::SkyHistogram2D *hRing = new Display::SkyHistogram2D("Ring","Ring",*hExcMod);     
  
  for(int i = 1; i <= nx; ++i)
    {
      for(int j = 1; j <= ny; ++j)
	{
	  hGammaMod->SetBinContent(i,j,0);
	  hGammaAccMod->SetBinContent(i,j,0);
	}
    }
  
  for(int i = 1; i <= nx0; ++i)
    {
      for(int j = 1; j <= ny0; ++j)
	{
	  Stash::Coordinate pos =  GetEventsMap()->GetBinCoordinate(i,j);
	  Int_t bin = hGammaMod->FindBinPosition(pos);
	  hGammaMod->SetBinContent(bin,GetEventsMap()->GetBinContent(i,j));
	  hExcMod->SetBinContent(bin,GetExclusionMap()->GetBinContent(i,j));
	}
    }
  for(int i = 1; i <= nxa; ++i)
    {
      for(int j = 1; j <= nya; ++j)
	{
	  Stash::Coordinate pos =  hGammaAccMod->GetBinCoordinate(i,j);
	  Int_t bin = GetAcceptanceMap()->FindBinPosition(pos);
	  hGammaAccMod->SetBinContent(i,j,GetAcceptanceMap()->GetBinContent(bin));
	}
    }
    
  for(int i = 1; i <= nx; ++i)
    {
      for(int j = 1; j <= ny; ++j)
	{
	  hGammaModExc->SetBinContent(i,j,hGammaMod->GetBinContent(i,j)*(1-hExcMod->GetBinContent(i,j)));
	  hGammaAccModExc->SetBinContent(i,j,hGammaAccMod->GetBinContent(i,j)*(1-hExcMod->GetBinContent(i,j)));
	}
    }

  hGammaMod->Multiply(hMaskMod);
  hGammaAccMod->Multiply(hMaskMod);
  hGammaModExc->Multiply(hMaskMod);
  hGammaAccModExc->Multiply(hMaskMod);
    
  //OS Kernel Generation
  for(int i = 1; i <= nx; ++i)
    {
      for(int j = 1; j <= ny; ++j)
	{
	  double xcenter = 0.5*nx+1;
	  double ycenter = 0.5*ny+1;
	  double r = TMath::Sqrt(TMath::Power(i-xcenter,2.)+TMath::Power(j-ycenter,2.));
	  if(r<=fOSRadius/fBinSize) {
	    hOS->SetBinContent(i,j,1.);
	  }
	  else {
	    hOS->SetBinContent(i,j,0.);
	  }
	}
    }

  //  Display::SkyHistogram2D *hExcess = new Display::SkyHistogram2D("Excess","Excess",*hExcMod); 
  Display::SkyHistogram2D *hExcessUncorrelated = new Display::SkyHistogram2D("ExcessUncorrelated","ExcessUncorrelated",*hExcMod); 
  Display::SkyHistogram2D *hExcessUncorrelated_Exc = new Display::SkyHistogram2D("ExcessUncorrelated_Exc","ExcessUncorrelated_Exc",*hExcMod); 
  //Display::SkyHistogram2D *hdExcess = new Display::SkyHistogram2D("dExcess","dExcess",*hExcMod);   
  Display::SkyHistogram2D *hRingON = new Display::SkyHistogram2D("RingON","RingON",*hExcMod); 
  Display::SkyHistogram2D *hRingON_Exc = new Display::SkyHistogram2D("RingON_Exc","RingON_Exc",*hExcMod); 
  Display::SkyHistogram2D *hRingONUncorrelated = new Display::SkyHistogram2D("RingONUncorrelated","RingONUncorrelated",*hExcMod); 
  Display::SkyHistogram2D *hRingONUncorrelated_Exc = new Display::SkyHistogram2D("RingONUncorrelated_Exc","RingONUncorrelated_Exc",*hExcMod); 

  Display::SkyHistogram2D *hRingOFF = new Display::SkyHistogram2D("RingOFF","RingOFF",*hExcMod); 
  Display::SkyHistogram2D *hRingOFFUncorrelated = new Display::SkyHistogram2D("RingOFFUncorrelated","RingOFFUncorrelated",*hExcMod); 
  Display::SkyHistogram2D *hRingAccON = new Display::SkyHistogram2D("RingAccON","RingAccON",*hExcMod); 
  Display::SkyHistogram2D *hRingAccON_Exc = new Display::SkyHistogram2D("RingAccON_Exc","RingAccON_Exc",*hExcMod); 
  Display::SkyHistogram2D *hRingAccOFF = new Display::SkyHistogram2D("RingAccOFF","RingAccOFF",*hExcMod);

  Display::SkyHistogram2D *hNRing = new Display::SkyHistogram2D("NR","NR",*hExcMod);     
  Display::SkyHistogram2D *hNRing2 = new Display::SkyHistogram2D("NR2","NR2",*hExcMod); 
  Display::SkyHistogram2D *hAlphaOFFUncorrelated  = new Display::SkyHistogram2D("AlphaOFFUncorrelated","AlphaOFFUncorrelated",*hExcMod);
  Display::SkyHistogram2D *hAlphaOFFUncorrelated_Exc  = new Display::SkyHistogram2D("Alphaoffuncorrelated_Exc","Alphaoffuncorrelated_Exc",*hExcMod);

  Display::SkyHistogram2D *hResOFF = new Display::SkyHistogram2D("ResOFF","ResOFF",*hGammaAccMod);
  Display::SkyHistogram2D *hResExcess = new Display::SkyHistogram2D("ResExcess","ResExcess",*hGammaAccMod);
  Display::SkyHistogram2D *hResExcess_Exc = new Display::SkyHistogram2D("ResExcess_Exc","ResExcess_Exc",*hGammaAccMod);
  //Display::SkyHistogram2D *hResAlphaOFF = new Display::SkyHistogram2D("ResAlphaOFF","ResAlphaOFF",*hGammaAccMod); 
  //Display::SkyHistogram2D *hResAlphaOFF_Exc = new Display::SkyHistogram2D("ResAlphaOFF_Exc","ResAlphaOFF_Exc",*hGammaAccMod); 
  Display::SkyHistogram2D *hResConvArea = new Display::SkyHistogram2D("ResConvArea","ResConvArea",*hGammaAccMod);
  Display::SkyHistogram2D *hResOFFExp = new Display::SkyHistogram2D("ResOFFExp","ResOFFExp",*hGammaAccMod);


  if(true)
    {
      std::cout <<  Utilities::TextStyle::Blue() << "*** Using adaptive ring method and FFT ***" << Utilities::TextStyle::Reset() << std::endl;
      //RING
      Float_t InnerRingStep = fRingStep;
      // Source and rings radii
      Int_t nSteps = 0;
      if(fAdaptFFT){
	if(fConstantArea)
	  nSteps = TMath::FloorNint((fInnerRingMax - fRingMinimalRadius)*1./fRingStep);
	else if(fConstantThickness){
	  nSteps = TMath::FloorNint(((fOuterRingMax-fRingParam_Thickness) - fRingMinimalRadius)*1./fRingStep);
	}
	else// if(fConstantInnerRadius)
	  nSteps = TMath::FloorNint((fOuterRingMax - (fRingMinimalRadius+fRingParam_Thickness))*1./fRingStep);
      }
      std::cout << "Number of Steps = " << nSteps << std::endl;
      int nrings = nSteps;
      
      std::map<int,int> binsoff;
      std::map<int,int> binsoff2;
      Display::SkyHistogram2D *hExcFrac = new Display::SkyHistogram2D("ExcFrac","ExcFrac",*hExcMod); 
      
      for(int i = 1; i <= nx; ++i)
	{
	  for(int j = 1; j <= ny; ++j)
	    {
	      if(fAdaptCut_Alpha)
		hExcFrac->SetBinContent(i,j,fRingParam_AlphaMax);
	      else
		hExcFrac->SetBinContent(i,j,1);
	      int bin = hExcFrac->GetBin(i,j);
	      binsoff[bin] = 0;
	    }
	}
	        
      for(int k = 0; k <= nrings; ++k)
	{ 	
	  //Display::SkyHistogram2D *hResRing = new Display::SkyHistogram2D("ResRing","ResRing",*hExcMod); 
	  Display::SkyHistogram2D *hResRingExp = new Display::SkyHistogram2D("ResRingExp","ResRingExp",*hGammaAccMod);
	  Display::SkyHistogram2D *hResRingAll = new Display::SkyHistogram2D("ResRingAll","ResRingAll",*hGammaAccMod); 
	  //Display::SkyHistogram2D *hRes = new Display::SkyHistogram2D("Res","Res",*hExcMod); 
	  Display::SkyHistogram2D *hResExp = new Display::SkyHistogram2D("ResExp","ResExp",*hGammaAccMod); 
	  std::cout << "InRing " << k << std::endl;      
	    
	  Float_t InnerRing = 0;
	  Float_t OuterRing = 0;
	  if(fAdaptFFT){
	    if(fConstantArea){
	      InnerRing = fRingMinimalRadius + k * fRingStep;
	      OuterRing = sqrt(fRingParam_AreaOverPi + std::pow(InnerRing,2));
	    }
	    else if(fConstantThickness){
	      InnerRing = fRingMinimalRadius + k * fRingStep;
	      OuterRing = InnerRing + fRingParam_Thickness;
	    }
	    else{// if(fConstantInnerRadius){
	      InnerRing = fRingMinimalRadius;
	      OuterRing = InnerRing + fRingParam_Thickness + k * fRingStep;
	    }
	  }
	  else{
	    InnerRing = fStandardRingRadius-fStandardRingThickness;
	    OuterRing = fStandardRingRadius+fStandardRingThickness;
	  }
	  Float_t RingRadius = 0.5*(InnerRing + OuterRing);
	  Float_t RingThickness = RingRadius - InnerRing;
	  Float_t r1 = RingRadius - RingThickness;
	  Float_t r2 = RingRadius + RingThickness;
	  std::cout << "Step :  " << k << "/" << nSteps << std::endl;
	  std::cout << "InnerRing = " << InnerRing << " - " << "OuterRing = " << OuterRing << std::endl ; 
	    
	  //Ring Kernel Generation
	  for(int i = 1; i <= nx; ++i)
	    {
	      for(int j = 1; j <= ny; ++j)
		{
		  double xcenter = 0.5*nx+1;
		  double ycenter = 0.5*ny+1;
		  double r = TMath::Sqrt(TMath::Power(i-xcenter,2.)+TMath::Power(j-ycenter,2.));
		  if(r <= r2/fBinSize && r >= r1/fBinSize){
		    hRing->SetBinContent(i,j,1);//6*TMath::Landau(dist,r1+0.05,0.1,0));
		  }
		}
	    }
	    
	  Utilities::TFftConv *convRExp = new Utilities::TFftConv(*hGammaAccModExc,*hRing);
	  Utilities::TFftConv *convRingAll = new Utilities::TFftConv(*hGammaAccMod,*hRing);
	  Utilities::TFftConv *convOExp = new Utilities::TFftConv(*hGammaAccMod,*hOS);
	    
	  convOExp->convolve(hResExp);
	  convRExp->convolve(hResRingExp);
	  convRingAll->convolve(hResRingAll);
	    
	  if(fAdaptCut_Alpha){
	    for(int i = 1; i <= nx; ++i)
	      {
		for(int j = 1; j <= ny; ++j)
		  {
		    Stash::Coordinate pos =  hGammaMod->GetBinCoordinate(i,j);
		    int bin = hExcFrac->GetBin(i,j);

		    double alpha = 0;
		    if(hResRingExp->GetBinContent(i,j) != 0 && hResRingExp->GetBinContent(i,j) > 0)
		      alpha =  hResExp->GetBinContent(i,j)*1./hResRingExp->GetBinContent(i,j);
		    
		    if((alpha < hExcFrac->GetBinContent(i,j) && hExcFrac->GetBinContent(i,j) > 0 && binsoff[bin] == 0) || k == 0)
		      {
			hExcFrac->SetBinContent(i,j,alpha);
			
			if(i > nSafeBins && j > nSafeBins && i <= nx-nSafeBins+1 && j <= ny-nSafeBins+1)
			  {
			    hNRing->SetBinContent(i,j,k);
			    hExcFrac->SetBinContent(i,j,alpha);			      
			    if(alpha < fRingParam_AlphaMax && alpha > 0){
			      binsoff[bin] = 1;
			    }
			  }
		      }
		    if(hExcFrac->GetBinContent(i,j) <= 0)
		      hExcFrac->SetBinContent(i,j,1);
		  }		    		    
	      }
	  }
	  else{
	    for(int i = 1; i <= nx; ++i)
	      {
		for(int j = 1; j <= ny; ++j)
		  {

		    Stash::Coordinate pos =  hGammaMod->GetBinCoordinate(i,j);
		    int bin = hExcFrac->GetBin(i,j);
		    double vallRing = hResRingAll->GetBinContent(i,j);
		    double vaccRing = hResRingExp->GetBinContent(i,j);
		    double excfrac = (vallRing-vaccRing)*1./vallRing;
		    
		    if(excfrac < hExcFrac->GetBinContent(i,j) && binsoff[bin] == 0)
		      {
			hExcFrac->SetBinContent(i,j,excfrac);
			
			if(i > nSafeBins && j > nSafeBins && i <= nx-nSafeBins+1 && j <= ny-nSafeBins+1)
			  {
			    hNRing->SetBinContent(i,j,k);
			    hExcFrac->SetBinContent(i,j,excfrac);
			    if(excfrac < fRingParam_ExcFracMax){
				binsoff[bin] = 1;
			    }
			  }
		      }
		  }
	      }
	  }
	  
	  convOExp->Delete();
	  convRExp->Delete();
	  convRingAll->Delete();      
	  
	  hRing->Reset();
	  hResRingExp->Delete();
	  hResRingAll->Delete();
	  hResExp->Delete();
	}
	
      std::cout << hNRing->GetMaximum() << std::endl;
      Int_t NRingsMax = (int)hNRing->GetMaximum();
      
      for(int i = 1; i <= nx; ++i)
	{
	  for(int j = 1; j <= ny; ++j)
	    {
	      int val = (int)hNRing->GetBinContent(i,j);
	      int bin = hNRing->GetBin(i,j);
	      hNRing2->SetBinContent(i,j,val);
	      binsoff2[bin] = val;
	    }
	}
      
      if(fSmoothRings){
	while(NRingsMax > 0)
	  {
	    std::cout << NRingsMax << std::endl;
	    for(int i = 1; i <= nx; ++i)
	      {
		for(int j = 1; j <= ny; ++j)
		  {
		    double val = hNRing2->GetBinContent(i,j);
		    int bin = hNRing2->GetBin(i,j);
		    double nextval = 0;
			
		    if(binsoff2[bin] == NRingsMax && val == NRingsMax)
		      {	  
			for(int ix = i-1; ix <= i+1; ++ix)
			  {
			    for(int iy = j-1; iy <= j+1; ++iy)
			      {
				int numx = ix-i;
				int numy = iy-j;
				int nextBin = hNRing2->GetBin(ix,iy);
				if(binsoff2[nextBin] < val)
				  {
				    if(numx != 0 || numy != 0) 
				      { 
					nextval = hNRing2->GetBinContent(ix,iy);
					hNRing2->SetBinContent(ix,iy,val-1);
					binsoff2[nextBin] = NRingsMax-1;
				      }
				  }
			      }
			  }
		      } 
		  }
	      }
	    NRingsMax--;
	  }
      }
      //Second adapt ring knowing ring numbers for each pixel
      Int_t MaxRingNumber = (int)hNRing->GetMaximum();
      std::cout << "MaxRingNumber = " << MaxRingNumber << std::endl;
      for(int k = 0; k <= MaxRingNumber; ++k)
	{ 
	  Display::SkyHistogram2D *hResRing = new Display::SkyHistogram2D("ResRing","ResRing",*hExcMod); 
	  Display::SkyHistogram2D *hResRingExp = new Display::SkyHistogram2D("ResRingExp","ResRingExp",*hExcMod);
	  Display::SkyHistogram2D *hRes = new Display::SkyHistogram2D("Res","Res",*hExcMod); 
	  Display::SkyHistogram2D *hResExp = new Display::SkyHistogram2D("ResExp","ResExp",*hExcMod); 
	  Display::SkyHistogram2D *hResExc = new Display::SkyHistogram2D("ResExc","ResExc",*hExcMod); 
	  Display::SkyHistogram2D *hResExpExc = new Display::SkyHistogram2D("ResExpExc","ResExpExc",*hExcMod); 

	  std::cout << "InRing " << k << std::endl;      
	  
	  Float_t InnerRing = 0;
	  Float_t OuterRing = 0;
	  if(fAdaptFFT){
	    if(fConstantArea){
	      InnerRing = fRingMinimalRadius + k * fRingStep;
	      OuterRing = sqrt(fRingParam_AreaOverPi + std::pow(InnerRing,2));
	    }
	    else if(fConstantThickness){
	      InnerRing = fRingMinimalRadius + k * fRingStep;
	      OuterRing = InnerRing + fRingParam_Thickness;
	    }
	    else{// if(fConstantInnerRadius){
	      InnerRing = fRingMinimalRadius;
	      OuterRing = InnerRing + fRingParam_Thickness + k * fRingStep;
	    }
	  }
	  else{
	    InnerRing = fStandardRingRadius-fStandardRingThickness;
	    OuterRing = fStandardRingRadius+fStandardRingThickness;
	  }
	  Float_t RingRadius = 0.5*(InnerRing + OuterRing);
	  Float_t RingThickness = RingRadius - InnerRing;
	  Float_t r1 = RingRadius - RingThickness;
	  Float_t r2 = RingRadius + RingThickness;
	  std::cout << "Pass 2 / Step : " << k << "/" << nSteps << std::endl;
	  std::cout << "InnerRing = " << InnerRing << " / " << "OuterRing = " << OuterRing << std::endl ; 
	  
	  //Ring Kernel Generation
	  for(int i = 1; i <= nx; ++i)
	    {
	      for(int j = 1; j <= ny; ++j)
		{
		  double xcenter = 0.5*nx+1;
		  double ycenter = 0.5*ny+1;
		  double r = TMath::Sqrt(TMath::Power(i-xcenter,2.)+TMath::Power(j-ycenter,2.));
		  if(r <= r2/fBinSize && r >= r1/fBinSize){
		    hRing->SetBinContent(i,j,1);
		  }
		}
	    }
	  
	  Utilities::TFftConv *convR = new Utilities::TFftConv(*hGammaModExc,*hRing);
	  Utilities::TFftConv *convRExp = new Utilities::TFftConv(*hGammaAccModExc,*hRing);
	  Utilities::TFftConv *convO = new Utilities::TFftConv(*hGammaMod,*hOS);
	  Utilities::TFftConv *convOExp = new Utilities::TFftConv(*hGammaAccMod,*hOS);
	  Utilities::TFftConv *convOExc = new Utilities::TFftConv(*hGammaModExc,*hOS);
	  Utilities::TFftConv *convOExpExc = new Utilities::TFftConv(*hGammaAccModExc,*hOS);
	  
	  convO->convolve(hRes);
	  convOExp->convolve(hResExp);
	  convOExc->convolve(hResExc);
	  convOExpExc->convolve(hResExpExc);
	  convR->convolve(hResRing);  
	  convRExp->convolve(hResRingExp);
	  
	  for(int i = 1; i <= nx; ++i)
	    {
	      for(int j = 1; j <= ny; ++j)
		{
		  int bin = hExcFrac->GetBin(i,j);
		  
		  if(k == binsoff2[bin])
		    {
		      double vON = hRes->GetBinContent(i,j);
		      double vONExc = hResExc->GetBinContent(i,j);		      
		      double vONUncorrelated = hGammaMod->GetBinContent(i,j);
		      double vONUncorrelated_Exc = hGammaModExc->GetBinContent(i,j);
		      
		      double vOFF = hResRing->GetBinContent(i,j);
		      double vOFFUncorrelated = hResRing->GetBinContent(i,j);

		      double alphaOFFUncorrelated =  0;
		      double alphaOFFUncorrelated_Exc =  0;

		      if(hResRingExp->GetBinContent(i,j) != 0){
			alphaOFFUncorrelated =  vOFFUncorrelated* (hGammaAccMod->GetBinContent(i,j)*1./hResRingExp->GetBinContent(i,j));
			alphaOFFUncorrelated_Exc =  vOFFUncorrelated* (hGammaAccModExc->GetBinContent(i,j)*1./hResRingExp->GetBinContent(i,j));
		      }
		      double AON = hResExp->GetBinContent(i,j);
		      double AONExc = hResExpExc->GetBinContent(i,j);
		      
		      if(i > nSafeBins && j > nSafeBins && i <= nx-nSafeBins+1 && j <= ny-nSafeBins+1)
			{
			  hExcessUncorrelated->SetBinContent(i,j,vONUncorrelated-alphaOFFUncorrelated);
			  hExcessUncorrelated_Exc->SetBinContent(i,j,vONUncorrelated_Exc-alphaOFFUncorrelated_Exc);
			  hRingON->SetBinContent(i,j,vON);
			  hRingON_Exc->SetBinContent(i,j,vONExc);
			  hRingONUncorrelated->SetBinContent(i,j,vONUncorrelated);
			  hRingONUncorrelated_Exc->SetBinContent(i,j,vONUncorrelated_Exc);
			  hRingOFF->SetBinContent(i,j,vOFF);
			  hRingOFFUncorrelated->SetBinContent(i,j,vOFFUncorrelated);
			  hRingAccON->SetBinContent(i,j,AON);
			  hRingAccON_Exc->SetBinContent(i,j,AONExc);
			  hRingAccOFF->SetBinContent(i,j,hResRingExp->GetBinContent(i,j));
			  hAlphaOFFUncorrelated->SetBinContent(i,j,alphaOFFUncorrelated);
			  hAlphaOFFUncorrelated_Exc->SetBinContent(i,j,alphaOFFUncorrelated_Exc);
			}
		    }
		}
	    }
	  
	  convO->Delete();
	  convOExp->Delete();
	  convOExc->Delete();
	  convOExpExc->Delete();
	  convR->Delete();
	  convRExp->Delete();
	  
	  hRing->Reset();
	  hResRing->Delete();
	  hRes->Delete();
	  hResExc->Delete();
	  hResRingExp->Delete();
	  hResExp->Delete();	  
	  hResExpExc->Delete();	  
	}
    
      hExcFrac->Delete();
      hNRing->Delete();
      hNRing2->Delete();
      
      std::cout << "LASTMAP" << std::endl;
      Utilities::TFftConv *convOFF = new Utilities::TFftConv(*hRingOFFUncorrelated,*hOS);
      Utilities::TFftConv *convArea = new Utilities::TFftConv(*hMaskMod,*hOS);
      //
      Utilities::TFftConv *convOFFExp = new Utilities::TFftConv(*hRingAccOFF,*hOS);
      Utilities::TFftConv *convExcess = new Utilities::TFftConv(*hExcessUncorrelated,*hOS);
      Utilities::TFftConv *convExcess_Exc = new Utilities::TFftConv(*hExcessUncorrelated_Exc,*hOS);
      //Utilities::TFftConv *convAlphaOFF = new Utilities::TFftConv(*hAlphaOFFUncorrelated,*hOS);
      //Utilities::TFftConv *convAlphaOFF_Exc = new Utilities::TFftConv(*hAlphaOFFUncorrelated_Exc,*hOS);
      
      convOFFExp->convolve(hResOFFExp);
      convExcess->convolve(hResExcess);
      convExcess_Exc->convolve(hResExcess_Exc);
      convOFF->convolve(hResOFF);
      //convAlphaOFF->convolve(hResAlphaOFF);
      //convAlphaOFF_Exc->convolve(hResAlphaOFF_Exc);
      convArea->convolve(hResConvArea);

      hResOFFExp->Multiply(hMaskMod);
      hResExcess->Multiply(hMaskMod);
      hResExcess_Exc->Multiply(hMaskMod);
      hResOFF->Multiply(hMaskMod);
      //hResAlphaOFF->Multiply(hMaskMod);
      //hResAlphaOFF_Exc->Multiply(hMaskMod);
      hResConvArea->Multiply(hMaskMod);
      
    }

  //UnCorrelated
  Display::SkyHistogram2D *ON =  new Display::SkyHistogram2D("ON","ON",*GetEventsMap()); 
  Display::SkyHistogram2D *ONExposure =  new Display::SkyHistogram2D("ONExposure","ONExposure",*GetEventsMap()); 
  Display::SkyHistogram2D *OFF =  new Display::SkyHistogram2D("OFF","OFF",*GetEventsMap()); 
  Display::SkyHistogram2D *OFFExposure =  new Display::SkyHistogram2D("OFFExposure","OFFExposure",*GetEventsMap()); 
  Display::SkyHistogram2D *Alpha =  new Display::SkyHistogram2D("Alpha","Alpha",*GetEventsMap()); 
  Display::SkyHistogram2D *AlphaOff =  new Display::SkyHistogram2D("AlphaOff","AlphaOff",*GetEventsMap());
  Display::SkyHistogram2D *Excess =  new Display::SkyHistogram2D("Excess","Excess",*GetEventsMap()); 
  //Correlated
  Display::SkyHistogram2D *ONCorrelated =  new Display::SkyHistogram2D("ONCorrelated","ONCorrelated",*GetEventsMap()); 
  Display::SkyHistogram2D *ONExposureCorrelated =  new Display::SkyHistogram2D("ONExposureCorrelated","ONExposureCorrelated",*GetEventsMap()); 
  Display::SkyHistogram2D *OFFCorrelated =  new Display::SkyHistogram2D("OFFCorrelated","OFFCorrelated",*GetEventsMap()); 
  Display::SkyHistogram2D *OFFExposureCorrelated =  new Display::SkyHistogram2D("OFFExposureCorrelated","OFFExposureCorrelated",*GetEventsMap()); 
  Display::SkyHistogram2D *AlphaCorrelated =  new Display::SkyHistogram2D("AlphaCorrelated","AlphaCorrelated",*GetEventsMap()); 
  Display::SkyHistogram2D *AlphaOffCorrelated =  new Display::SkyHistogram2D("AlphaOffCorrelated","AlphaOffCorrelated",*GetEventsMap());
  Display::SkyHistogram2D *ExcessCorrelated =  new Display::SkyHistogram2D("ExcessCorrelated","ExcessCorrelated",*GetEventsMap()); 
  //Significances
  Display::SkyHistogram2D *SignificanceMap =  new Display::SkyHistogram2D("SignificanceMap","SignificanceMap",*GetEventsMap());
  Display::SkyHistogram2D *AllowedSignificanceMap =  new Display::SkyHistogram2D("AllowedSignificanceMap","AllowedSignificanceMap",*GetEventsMap());
  TH1F *SkySigmaDist = new TH1F("SigmaDist","SigmaDist",100,-10,10);    
  TH1F *AllowedRegionsSkySigmaDist = new TH1F("AllowedSigmaDist","AllowedSigmaDist",100,-10,10);

  //Loop to fill global map
  std::cout << "Filling the global maps " << std::endl;

  for(int i = 1; i <= nx; ++i)
    {
      for(int j = 1; j <= ny; ++j)
	  {
	    Stash::Coordinate pos =  hExcMod->GetBinCoordinate(i,j);
	    Int_t bin = GetEventsMap()->FindBinPosition(pos);
	    Int_t binacc = GetAcceptanceMap()->FindBinPosition(pos);
	    
	    if(GetAcceptanceMap()->GetBinContent(binacc) > 0)
	      {
		//Uncorrelated
		double on = hGammaMod->GetBinContent(i,j);
		double off = hRingOFFUncorrelated->GetBinContent(i,j);
		double excess = hExcessUncorrelated->GetBinContent(i,j);
		double alpha = 0;
		if(off > 0)
		  alpha = (on-excess)*1./off;
		double aon = hGammaAccMod->GetBinContent(i,j);
		double aoff = aon*1./alpha;
		
		ON->SetBinContent(bin,on);
		ON->SetBinError(bin,sqrt(on));
		ONExposure->SetBinContent(bin,aon);
		OFF->SetBinContent(bin,off);
		OFF->SetBinError(bin,sqrt(off));
		OFFExposure->SetBinContent(bin,aoff);
		Excess->SetBinContent(bin,on-alpha*off);
		Excess->SetBinError(bin,Utilities::Statistics::LiMa_dExcess_Up(on,off,alpha));
		Alpha->SetBinContent(bin,alpha);
		Alpha->SetBinError(bin,1);
		AlphaOff->SetBinContent(bin,alpha*off);
		AlphaOff->SetBinError(bin,1);
		
		//Correlated
		double on_corr = hRingON->GetBinContent(i,j);
		double off_corr = 0;
		if(fAverageOff)
		  off_corr = hResOFF->GetBinContent(i,j)*1./hResConvArea->GetBinContent(i,j);
		else
		  off_corr = hRingOFF->GetBinContent(i,j);
		double excess_corr = hResExcess->GetBinContent(i,j);
		
		if(off_corr > 5)
		  {
		    double alpha_corr = 0;
		    double onexp_corr = 0;
		    double offexp_corr = 0;
		    if(fAverageOff){
		      //alpha_corr = hResAlphaOFF->GetBinContent(i,j)*1./off_corr;
		      onexp_corr = hRingAccON->GetBinContent(i,j);
		      offexp_corr = hResOFFExp->GetBinContent(i,j)*1./hResConvArea->GetBinContent(i,j);
		      alpha_corr = onexp_corr*1./offexp_corr;
		    }
		    else
		      alpha_corr = (on_corr-excess_corr)*1./off_corr;
		    
		    if(alpha_corr*off_corr > 1)
		      {
			ONCorrelated->SetBinContent(bin,hRingON->GetBinContent(i,j));
			ONCorrelated->SetBinError(bin,sqrt(hRingON->GetBinContent(i,j)));
			ONExposureCorrelated->SetBinContent(bin,onexp_corr);//hRingAccON->GetBinContent(i,j));
			OFFCorrelated->SetBinContent(bin,off_corr);
			OFFCorrelated->SetBinError(bin,sqrt(off_corr));
			OFFExposureCorrelated->SetBinContent(bin,offexp_corr);//1./alpha_corr*hRingAccON->GetBinContent(i,j));
			double significance = Utilities::Statistics::LiMa(on_corr,off_corr,alpha_corr);
			SignificanceMap->SetBinContent(bin,significance);
			SignificanceMap->SetBinError(bin,1);	  
			
			ExcessCorrelated->SetBinContent(bin,on_corr-alpha_corr*off_corr);
			ExcessCorrelated->SetBinError(bin,Utilities::Statistics::LiMa_dExcess_Up(on_corr,off_corr,alpha_corr));
			AlphaCorrelated->SetBinContent(bin,alpha_corr);
			AlphaCorrelated->SetBinError(bin,1);
			AlphaOffCorrelated->SetBinContent(bin,alpha_corr*off_corr);
			AlphaOffCorrelated->SetBinError(bin,1);
	
			if(significance != 0)
			  SkySigmaDist->Fill(significance);
			
			if(GetExclusionMap()->GetBinContent(GetExclusionMap()->FindBinPosition(pos)) == 0 && significance != 0)
			  {
			    double onexc = hRingON_Exc->GetBinContent(i,j);
			    double offexc = 0;
			    if(fAverageOff)
			      offexc = hResOFF->GetBinContent(i,j)*1./hResConvArea->GetBinContent(i,j);
			    else
			      offexc = hRingOFF->GetBinContent(i,j);
			    
			    double alphaexc = 0;
			    double onexpexc = 0;
			    double offexpexc = 0;
			    if(fAverageOff){
			      onexpexc = hRingAccON_Exc->GetBinContent(i,j);
			      offexpexc = hResOFFExp->GetBinContent(i,j)*1./hResConvArea->GetBinContent(i,j);;
			      alphaexc = onexpexc*1./offexpexc;//hResAlphaOFF_Exc->GetBinContent(i,j)*1./offexc;
			    }
			    else
			      alphaexc = (onexc-hResExcess_Exc->GetBinContent(i,j))*1./hRingOFF->GetBinContent(i,j);
			    
			    double significance_exc = Utilities::Statistics::LiMa(onexc,offexc,alphaexc);
			    AllowedRegionsSkySigmaDist->Fill(significance_exc);
			    AllowedSignificanceMap->SetBinContent(bin,significance_exc);
			    AllowedSignificanceMap->SetBinError(bin,1);	  
			  }
		      }
		  }
	      }
	  }
      }
      
  GetSignificanceMap() = SignificanceMap;
  GetExcessMap() = ExcessCorrelated;
  GetAlphaMap() = AlphaCorrelated;
  GetOffMap() = OFFCorrelated;

  delete hGammaModExc;
  delete hGammaAccModExc;
  delete hOS;
  delete hRing;
  delete hExcessUncorrelated;
  delete hExcessUncorrelated_Exc;
  delete hRingON;
  delete hRingON_Exc;
  delete hRingONUncorrelated;
  delete hRingONUncorrelated_Exc;
  delete hRingOFF;
  delete hRingOFFUncorrelated;
  delete hRingAccON;
  delete hRingAccON_Exc;
  delete hRingAccOFF;
  delete hAlphaOFFUncorrelated;
  delete hAlphaOFFUncorrelated_Exc;
  delete hResOFF;
  delete hResExcess;
  delete hResExcess_Exc;
  delete hResConvArea;
  delete hResOFFExp;


  delete ON;
  delete ONExposure;
  delete OFF;
  delete OFFExposure;
  delete Alpha;
  delete AlphaOff;
  delete Excess;
  delete ONCorrelated;
  delete ONExposureCorrelated;
  delete OFFExposureCorrelated;
  delete AlphaOffCorrelated;
  delete AllowedSignificanceMap;
  delete SkySigmaDist;
  delete AllowedRegionsSkySigmaDist;

  fileResults->Close();
  
  return true;
}

