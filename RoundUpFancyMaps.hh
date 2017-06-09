
#include <TNamed.h>

#include <TGraph.h>
#include <TString.h>
#include <TNtuple.h>
#include <TNtupleD.h>
#include <TH2F.h>

#include <vector>
#include <iostream>

#include <sash/MonitorBase.hh>
#include <display/CanvasHolder.hh>
#include <display/SkyHistogram2D.hh>

namespace SurveySuite {
  
  class MapMaker : public Sash::MonitorBase,
		   public Display::CanvasHolder
  {
  public:
    typedef TH1F* pTH1F;
    typedef Display::SkyHistogram2D * pSkyHisto;

    MapMaker();//Sash::HESSArray & hess
	     //#if defined(__CINT__) || defined(CINTOBJECT)
	//	   = hessdummy
	     //#endif
    //	   ,const char *name = "");
    
    virtual ~MapMaker();

    GETMEMBER(EventsMap,pSkyHisto);
    GETMEMBER(AcceptanceMap,pSkyHisto);
    GETMEMBER(ExclusionMap,pSkyHisto);
    GETMEMBER(ExpectedCountsMap,pSkyHisto);
    GETMEMBER(ExtendedExpectedCountsMap,pSkyHisto);
    GETMEMBER(ExcessMap,pSkyHisto);
    GETMEMBER(UncorrelatedExcessMap,pSkyHisto);
    GETMEMBER(SignificanceMap,pSkyHisto);
    GETMEMBER(AlphaMap,pSkyHisto);
    GETMEMBER(OffMap,pSkyHisto);
    GETMEMBER(OffAxisMap,pSkyHisto);
    GETMEMBER(ZenithAngleMap,pSkyHisto);
    GETMEMBER(MuonEfficiencyMap,pSkyHisto);
    GETMEMBER(MinSafeThresholdMap,pSkyHisto);
    GETMEMBER(AveragedSafeThresholdMap,pSkyHisto);
    GETMEMBER(AveragedLiveTimeMap,pSkyHisto);

    GETMEMBER(EventsAndAcceptanceFromFile,bool);
    GETMEMBER(EventsAndAcceptanceFromRadial,bool);
    GETMEMBER(EventsAndAcceptanceFromRadialFile,bool);
    GETMEMBER(ExclusionFromFits, bool); 
    GETMEMBER(ExclusionFromRegionFile, bool);
    GETMEMBER(UseConfigTarget, bool); 
    GETMEMBER(UseConfigMapParams, bool);
    GETMEMBER(UserLambda,double);
    GETMEMBER(UserBeta,double); 
    GETMEMBER(BinSize,double);
    GETMEMBER(ExtX,double); 
    GETMEMBER(ExtY,double); 
    GETMEMBER(PsiCut,double); 
    GETMEMBER(OSRadius,double); 
    GETMEMBER(AdaptFFT,bool); 
    GETMEMBER(AdaptCut_Alpha,bool); 
    GETMEMBER(ConstantArea,bool);
    GETMEMBER(ConstantThickness,bool);
    GETMEMBER(ConstantInnerRadius,bool);
    GETMEMBER(SmoothRings,bool);
    GETMEMBER(RingMinimalRadius,double);
    GETMEMBER(InnerRingMax,double);
    GETMEMBER(OuterRingMax,double);
    GETMEMBER(RingStep,double);
    GETMEMBER(RingParam_AreaOverPi,double); 
    GETMEMBER(RingParam_ExcFracMax,double); 
    GETMEMBER(RingParam_Thickness,double);
    GETMEMBER(RingParam_AlphaMax,double);
    GETMEMBER(StandardRingRadius,double);
    GETMEMBER(StandardRingThickness,double);
    GETMEMBER(AverageOff,bool);
    GETMEMBER(CorrectZenith,bool);
    GETMEMBER(EMin,double);
    GETMEMBER(EMax,double);
    GETMEMBER(FitAcceptance,bool);
    GETMEMBER(SelectAllRunsContributingToTheMap,bool);
    GETMEMBER(ApplySafeThreshold,bool);
    GETMEMBER(SafeThresholdFromAcceptance,bool);
    GETMEMBER(SafeThresholdFromPsiCut,bool);
    GETMEMBER(SafeThresholdFromPsiValue,double);
    GETMEMBER(SafeThresholdRatioParam,double);
    GETMEMBER(ProduceFluxProducts,bool);
    GETMEMBER(SpectralIndex,double);
    GETMEMBER(PointLikeFluxMaps,bool);
    GETMEMBER(ProduceSurfaceBrightnessMap,bool);
    GETMEMBER(ExposureMapsFromFits,bool);
    GETMEMBER(SaveResults,bool);
    GETMEMBER(Verbose,bool);
    GETMEMBER(ConfigureOutputsFlag, bool);
    GETMEMBER(ConfigureFluxProductsFlag, bool); 
    GETMEMBER(ConfigureRingMethodFlag, bool);
    GETMEMBER(ConfigureMapsFlag, bool); 
    GETMEMBER(ConfigureAcceptanceFlag, bool); 
    GETMEMBER(ConfigureExclusionsFlag, bool);
    GETMEMBER(StartConfigureFlag, bool);

    void Clear();

    void StartConfigure();

    bool ConfigureExclusions(bool ExclusionFromFits, 
			     bool ExclusionFromRegionFile);

    bool ConfigureAcceptance(bool EventsAndAcceptanceFromFile, 
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
			     double SafeThresholdRatioParam);
    
    bool ConfigureMaps(bool UseConfigTarget, 
		       bool UseConfigMapParams, 
		       double UserLambda, 
		       double UserBeta, 
		       double BinSize, 
		       double ExtX, 
		       double ExtY, 
		       double OSRadius, 
		       double EMin, 
		       double EMax);
    
    bool ConfigureRingMethod(bool AdaptFFT, 
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
			     bool AverageOff);

    bool ConfigureFluxProducts(bool ProduceFluxProducts,
			       double SpectralIndex, 
			       bool PointLikeFluxMaps, 
			       bool ProduceSurfaceBrightnessMap,
			       bool ExposureMapsFromFits,
			       const char *AnalysisConfig);

    bool ConfigureOutputs(bool SaveResults, bool Verbose);
    bool ConfigureRuns(const char *Resfile, int RunNumberMin, int RunNumberMax, bool UseRunListToMatch, bool UseRunsToForbid, const char *RunListToMatch, const char *RunsToForbid, std::string table_path);
    bool IsRunInFile(int RunNumber, std::vector<int> File, bool ReverseAnswer); 
    void RemoveBadTablesRuns(const char *Resfile, std::string table_path);
    void ComputeSafeThresholdPerRun(const char *RadialConfig, const char *Resfile, std::string table_path);
    bool EndConfigure(const char *RadialConfig, const char *Resfile, int RunNumberMin, int RunNumberMax, bool UseRunListToMatch, bool UseRunsToForbid, const char *RunListToMatch, const char *RunsToForbid, std::string table_path);
    void PrintConfigSummary();

    void RetrieveRunStat(const char *Resfile);    
    void DisplayEventsMap(const char *Resfile); // *MENU={Hierarchy="Clear All"}*
    void DrawRadialLookup(const char *Config, double ZenithAngle);
    void DisplayEventsAndAcceptanceMaps();
    void CreateEventsMap_FromFITS(const char *EventsMapFile);
    void CreateAcceptanceMap_FromFITS(const char *AccMapFile);
    void CreateExclusionMap_FromFITS(const char *ExclusionMapFile);
    void CreateEventsMap_FromResultFile(const char *Resfile);
    void CreateAcceptanceMap_FromResultFile(const char *Resfile);
    void CreateEventsAndAcceptanceMaps_FromRadialLookups(const char *RadialConfig, const char *Resfile, const char *Regionsfile, std::string table_path, const char *OutfilePrefix);
    void CreateExposureMaps_FromFITS(const char *ExposureMapFile);
    void CreateExposureMaps_FromFile(const char *ExposureMapFile);
    void CreateExtendedExposureMaps_FromFITS(const char *ExposureMapFile);
    void CreateExtendedExposureMaps_FromFile(const char *ExposureMapFile);
    void CreateExcessMaps_FromFITS(const char *ExcessMapFile);
    void CreateAlphaMaps_FromFITS(const char *AlphaMapFile);
    void CreateOffMaps_FromFITS(const char *OffMapFile);
    void CreateExclusionMask(const char *Resfile, const char *Regionsfile);
    void CreateEventsAndAcceptanceMaps(const char *RadialConfig, const char *Resfile, const char *Regionsfile, const char *Evtfile, const char *Accfile, std::string table_path, const char *OutfilePrefix);
    void CreateRunList(const char *Resfile, const char *Regionsfile, const char *OutfilePrefix, double SourceRadius, int sourceNumber);
    void CreateRunFileInfos(const char *Resfile, const char *Regionsfile, const char *OutfilePrefix);
    void CreateEventsInfosFile(const char *RadialConfig, const char *Resfile, const char *Regionsfile, std::string table_path, const char *OutfilePrefix);

    void CreateExposureMaps(const char *RadialConfig, const char *Resfile, const char *Regionsfile, std::string table_path, const char *OutfilePrefix);
    void CreatePSFMaps(const char *RadialConfig, const char *Resfile, const char *Regionsfile, std::string table_path, const char *OutfilePrefix);
    void CreateRingBgMaps(const char *AnalysisConfig, const char *Resfile, const char *Evtfile, const char *Accfile, const char *ExposureMapFile, const char *ExtendedExposureMapFile, const char *Regionsfile, std::string table_path, const char *OutfilePrefix);
    //void CreateFluxMaps_FromFITS(const char *AnalysisConfig, const char *ExcessMapFile, const char *AlphaMapFile, const char *OffMapFile, const char *ExposureMapFile, const char *OutfilePrefix);
    void CreateFluxMaps_FromResultFile(const char *AnalysisConfig, const char *ResultFile, const char *ExposureMapFile, const char *ExtendedExposureMapFile, const char *OutfilePrefix);

    void CreateEventsAndAcceptanceMaps_FromRadialLookups_PerRun(const char *RadialConfig, const char *Resfile, const char *Regionsfile, std::string table_path, const char *OutfilePrefix);    
    void CreateBgMaps_PerRun(const char *AnalysisConfig, const char *Resfile, const char *Evtfile, const char *Accfile, const char *ExposureMapFile, const char *ExtendedExposureMapFile, const char *Regionsfile, std::string table_path, const char *OutfilePrefix);

    void CreateONOFFTestMaps(const char *AnalysisConfig, const char *Resfile, const char *Evtfile, const char *Accfile, const char *ExposureMapFile, const char *ExtendedExposureMapFile, const char *Regionsfile, std::string table_path, const char *OutfilePrefix);

    bool CreateEventsAndAcceptanceMaps_FromRadialLookups_OneEntry(const char *RadialConfig, const char *Resfile, const char *Regionsfile, std::string table_path, const char *OutfilePrefix, int EntryNumber);
    bool CreateRingBgMaps_OneEntry(const char *AnalysisConfig, const char *Resfile, const char *Evtfile, const char *Accfile, const char *ExposureMapFile, const char *ExtendedExposureMapFile, const char *Regionsfile, std::string table_path, const char *OutfilePrefix, int EntryNumber);


  public:
    //SkyMaps
    Display::SkyHistogram2D *fEventsMap;
    Display::SkyHistogram2D *fAcceptanceMap;
    Display::SkyHistogram2D *fExclusionMap;
    Display::SkyHistogram2D *fExpectedCountsMap;
    Display::SkyHistogram2D *fExtendedExpectedCountsMap;
    Display::SkyHistogram2D *fExcessMap;
    Display::SkyHistogram2D *fUncorrelatedExcessMap;
    Display::SkyHistogram2D *fSignificanceMap;
    Display::SkyHistogram2D *fAlphaMap;
    Display::SkyHistogram2D *fOffMap;
    Display::SkyHistogram2D *fOffAxisMap;
    Display::SkyHistogram2D *fZenithAngleMap;
    Display::SkyHistogram2D *fMuonEfficiencyMap;
    Display::SkyHistogram2D *fMinSafeThresholdMap;
    Display::SkyHistogram2D *fAveragedSafeThresholdMap;
    Display::SkyHistogram2D *fAveragedLiveTimeMap;

    //values
    bool fEventsAndAcceptanceFromFile;
    bool fEventsAndAcceptanceFromRadial;
    bool fEventsAndAcceptanceFromRadialFile;
    bool fExclusionFromFits; 
    bool fExclusionFromRegionFile;
    bool fUseConfigTarget; 
    bool fUseConfigMapParams;
    double fUserLambda;
    double fUserBeta; 
    double fBinSize;
    double fExtX; 
    double fExtY; 
    double fPsiCut; 
    double fOSRadius; 
    bool fAdaptFFT;
    bool fAdaptCut_Alpha; 
    bool fConstantArea;
    bool fConstantThickness;
    bool fConstantInnerRadius;
    bool fSmoothRings;
    double fRingMinimalRadius;
    double fInnerRingMax;
    double fOuterRingMax;
    double fRingStep;
    double fRingParam_AreaOverPi; 
    double fRingParam_ExcFracMax; 
    double fRingParam_Thickness;
    double fRingParam_AlphaMax;
    double fStandardRingRadius;
    double fStandardRingThickness;
    bool fAverageOff;
    bool fCorrectZenith;
    double fEMin;
    double fEMax;
    bool fFitAcceptance;
    bool fSelectAllRunsContributingToTheMap;
    bool fApplySafeThreshold;
    bool fSafeThresholdFromAcceptance;
    bool fSafeThresholdFromPsiCut;
    double fSafeThresholdFromPsiValue;
    double fSafeThresholdRatioParam;
    bool fProduceFluxProducts;
    double fSpectralIndex;
    bool fPointLikeFluxMaps;
    bool fProduceSurfaceBrightnessMap;
    bool fExposureMapsFromFits;
    bool fSaveResults;
    bool fVerbose;

    //Configuration flags
    bool fConfigureOutputsFlag;
    bool fConfigureFluxProductsFlag;
    bool fConfigureRingMethodFlag;
    bool fConfigureMapsFlag;
    bool fConfigureAcceptanceFlag;
    bool fConfigureExclusionsFlag;
    bool fStartConfigureFlag;


    //Vectors, maps
    std::vector<int> fEntriesVector;
    std::map<int,double> fPerRunSafeThreshold;
    std::map<double,double> fScaleFactor;
    std::map<double,std::pair<double,double> > fRadialIntegrationEbin;

  protected:
    ClassDef(SurveySuite::MapMaker,1);  
  };  
}

