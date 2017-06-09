//void AnalyseRoundUpSource(char *Resfile, char *OutfilePrefix )
void AnalyseRoundUpSource(char *Resfile, char *OutfilePrefix,bool UseConfigTarget,double UserLambda,double UserBeta )
{
  gROOT->LoadMacro("RoundUpFancyMaps.C++");  
  SurveySuite::MapMaker *m = new SurveySuite::MapMaker();


  //Configure Run range and runs to forbid or to authorize
  int RunNumberMin = 0;
  int RunNumberMax = -1;
  bool UseRunListToMatch = 0;
  bool UseRunsToForbid = 0;
  char *RunListToMatch = ""; 
  char *RunsToForbid = "";
  
  double PsiCut = 2.0;
  bool SelectAllRunsContributingToTheMap = 1;

  //Configure ring background 
  bool AdaptFFT = 0; 
  bool AdaptCut_Alpha = 1;
  bool ConstantArea = 0;
  bool ConstantThickness = 1; 
  bool ConstantInnerRadius = 0; 
  bool SmoothRings = 1; 
  double RingMinimalRadius = 0.7;
  double InnerRingMax = 1.5;
  double OuterRingMax = 1.99;
  double RingStep = 0.08;
  double RingParam_AreaOverPi = 0.6;
  double RingParam_ExcFracMax = 0.5;
  double RingParam_Thickness = 0.44;
  double RingParam_AlphaMax = 0.25;
  double StandardRingRadius = 0.5;
  double StandardRingThickness = 0.1;
  bool AverageOff = 1;
  bool CorrectZenith = 0;
  double OSRadius = 0.1;

  bool FitAcceptance = 1;
  bool SaveResults = 1;
  bool Verbose = 1;
 
  char *AnalysisConfig = "Std";
  
  m->StartConfigure();
  bool ConfigExclusionSuccess = m->ConfigureExclusions(0, 0);
  bool ConfigAccSuccess = m->ConfigureAcceptance(0, 1, 0, PsiCut, 1, CorrectZenith, 1, 0, 1, 1, 0, 0.15);   
  bool ConfigMapSuccess = m->ConfigureMaps(UseConfigTarget, 1, UserLambda, UserBeta, 0, 0, 0, OSRadius, 0, -1);
  bool ConfigRingSuccess = m->ConfigureRingMethod(AdaptFFT, AdaptCut_Alpha, ConstantArea, ConstantThickness, ConstantInnerRadius, SmoothRings, RingMinimalRadius, InnerRingMax, OuterRingMax, RingStep, RingParam_AreaOverPi, RingParam_ExcFracMax, RingParam_Thickness, RingParam_AlphaMax, StandardRingRadius, StandardRingThickness, AverageOff);
  bool ConfigFluxSuccess = m->ConfigureFluxProducts(0, 0, 1, 0, 0, AnalysisConfig);
  bool ConfigOutSuccess = m->ConfigureOutputs(SaveResults, Verbose);
  bool ConfigEndSuccess = m->EndConfigure(AnalysisConfig, Resfile, RunNumberMin, RunNumberMax, UseRunListToMatch, UseRunsToForbid, RunListToMatch, RunsToForbid, "");
  bool ConfigSuccess = 0;
  if(ConfigExclusionSuccess && ConfigAccSuccess && ConfigMapSuccess && ConfigRingSuccess && ConfigFluxSuccess && ConfigOutSuccess && ConfigEndSuccess){
    ConfigSuccess = 1;
    m->PrintConfigSummary();
  }
  
  if(ConfigSuccess){
    m->CreateRingBgMaps(AnalysisConfig, Resfile, "", "", "", "", "", "", OutfilePrefix);
    m->CreateONOFFTestMaps(AnalysisConfig, Resfile, "", "", "", "", "", "", OutfilePrefix);
  }
  
}
