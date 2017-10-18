/**

\file  calibration_viewer.C

This macro displays histograms of pedestal and Nsb (estimated with HVI) timeline 
and flattens the broken pixels
Information are read from RunPedestal files (TotalJpT or PedNsb) and BrokenPixel files

It is possible to print some information using the function writeCAlibCoeffs() : this needs this macro to be compiled

*/

// #ifdef __MAKECINT__ // Part for dictionnary generation
// void calibration_viewer(int RunNo = -1, int MaxEvent=-1, bool UsePedNsb=false,const char *CalibVersion = "");
void writePedCoeffs(int telid);
void writeBpxCoeffs(int telid);
// #endif

// #if !defined (__CINT__) 
#include <sash/HESSArray.hh>
#include <sash/RunHeader.hh>
#include <sash/Telescope.hh>
#include <sash/Pixel.hh>
#include <sash/PixelConfig.hh>

#include <sash/DataSet.hh>
#include <sash/Folder.hh>
#include <sash/MakerChain.hh>
#include <sash/CINTMaker.hh>
#include <sash/Pointer.hh>
#include <sash/PointerSetIterator.hh>
#include <calibration/TelescopeCalibCoeffs.hh>
#include <calibration/CalibCoefficientsMaker.hh>
#include <pariscalibration/CalibrationFinder.hh>
#include <fromdb/FillRunHeader.hh>
#include <fromdb/CameraConfigurator.hh>
#include <fromdb/HeaderFixer.hh>
#include <pariscalibration/TelescopeLumi.hh>
#include <calibration/TelescopePedestalMonitor.hh>
#include <pariscalibration/LumiEstimator.hh>
#include <pariscalibrationmakers/CalibrationAccumulator.hh>
#include <iostream>
#include <valarray>
#include <TFile.h>
// #endif

// #ifndef  __MAKECINT__ 


Sash::HESSArray* hessarray=0;
Sash::RunHeader*  runheader=0;

int pixel=16; // pixel id to write the pedestal and NSB read at each pedestal event

void bpxviewer(int RunNo = -1,
			int MaxEvent=-1,
			bool UsePedNsb = false,
			const char *CalibVersion = "")
{
  if(RunNo == -1) {
    std::cout << "Please give run Number" << std::endl;
    return;
  }
    
  Sash::Folder * runf = Sash::Folder::GetFolder("run");
  hessarray=runf->GetHESSArray();
  runheader = hessarray->Handle((Sash::RunHeader*)0);
  runheader->GetRunNum() = RunNo;

  Sash::DataSet events("events");

  Sash::MakerChain *eventmanager = new Sash::MakerChain(true);
  eventmanager->UseMaker(new FromDB::FillRunHeader);
  eventmanager->UseMaker(new FromDB::CameraConfigurator);
  //eventmanager->UseMaker(new FromDB::HeaderFixer);
  eventmanager->UseMaker(new Calibration::CalibCoefficientsMaker("","","","","","","",1));
  ParisCalibration::CalibrationFinder *finder = new ParisCalibration::CalibrationFinder(ParisCalibration::CalibrationFinder::mAll,
											  ParisCalibration::CalibrationFinder::VerboseRun,
											  CalibVersion);
  if(UsePedNsb)
    finder->SetPedestalType(ParisCalibration::CalibrationFinder::PedNSB);
  finder->SetCalibType(ParisCalibration::CalibrationFinder::mPedestal|
		       ParisCalibration::CalibrationFinder::mBrokenPixel);
  finder->SetCalibOrigin(ParisCalibration::CalibrationFinder::mPedestal,
			 ParisCalibration::CalibrationFinder::fromFile);
  finder->SetCalibMode(ParisCalibration::CalibrationFinder::mPedestal,
		       ParisCalibration::CalibrationFinder::fromRun);
  finder->SetCalibOrigin(ParisCalibration::CalibrationFinder::mBrokenPixel,
			 ParisCalibration::CalibrationFinder::fromFile);
  finder->SetCalibMode(ParisCalibration::CalibrationFinder::mBrokenPixel,
		       ParisCalibration::CalibrationFinder::fromRun);   
  //finder->SetCalibType(ParisCalibration::CalibrationFinder::mPedestal);
  eventmanager->UseMaker(finder);
    
  // Estimate the NSB from HVI monitoring
  eventmanager->UseMaker(new ParisCalibration::LumiEstimator("^events$",ParisCalibration::LumiEstimator::mReadHV | ParisCalibration::LumiEstimator::mReadScaler | ParisCalibration::LumiEstimator::mReadPed, "NSB",120,false));

// #if !defined (__CINT__) 
  eventmanager->UseMaker(new Sash::CINTMaker("CT1_Pedestal","writePedCoeffs(1)"));
  eventmanager->UseMaker(new Sash::CINTMaker("CT2_Pedestal","writePedCoeffs(2)"));
  eventmanager->UseMaker(new Sash::CINTMaker("CT3_Pedestal","writePedCoeffs(3)"));
  eventmanager->UseMaker(new Sash::CINTMaker("CT4_Pedestal","writePedCoeffs(4)"));

  eventmanager->UseMaker(new Sash::CINTMaker("CT1_BrokenPixel","writeBpxCoeffs(1)"));
  eventmanager->UseMaker(new Sash::CINTMaker("CT2_BrokenPixel","writeBpxCoeffs(2)"));
  eventmanager->UseMaker(new Sash::CINTMaker("CT3_BrokenPixel","writeBpxCoeffs(3)"));
  eventmanager->UseMaker(new Sash::CINTMaker("CT4_BrokenPixel","writeBpxCoeffs(4)"));
// #endif    

  // We want to fill the calibration histogram
  eventmanager->UseMaker(new ParisCalibrationMakers::CalibrationAccumulator(true));
  
  // Now we loop over related datasets (Pedestal & BrokenPixel)
  eventmanager->Process(Sash::Folder::GetFolder("run"));
  TList *elist = events.GetListOfRelatedSets();
  TListIter iter(elist);
  TObject *o;
  while((o = iter()) != 0)
    {
      ((Sash::DataSet *)o)->EventLoop(eventmanager);
    }
  delete elist;
  eventmanager->Process(Sash::Folder::GetFolder("endrun"));
  eventmanager->Process(Sash::Folder::GetFolder("endanalysis"));
  delete eventmanager;

}

void writePedCoeffs(int telid)
{
  if(!hessarray) {
    std::cout << "writeCAlibCoeffs : NULL pointer HESSArray" << std::endl;
    return;
  }
  std::cout << "Received CT" << telid << "_Pedestal event" << std::endl; 
  const Sash::PointerSet<Sash::Telescope> & tels = runheader->GetTelsInRun();
  Sash::PointerSetIterator<Sash::Telescope> tel(tels.begin());
  for(;tel!=tels.end();++tel) {
    const Calibration::TelescopeCalibCoeffs * telcoeffs = (*tel)->Get((Calibration::TelescopeCalibCoeffs*)0);
    const ParisCalibration::TelescopeLumi *telnsb = (*tel)->Get((ParisCalibration::TelescopeLumi*)0);
    //    const Calibration::TelescopePedestalMonitor *telped = (*tel)->Get((Calibration::TelescopePedestalMonitor*)0);
    Sash::NonConstPointer<Sash::Pixel> begin = (*tel)->beginPixel();
    Sash::NonConstPointer<Sash::Pixel> end = (*tel)->endPixel();
    for(Sash::NonConstPointer<Sash::Pixel> pix = begin; pix!=end; ++pix) {
      /* Print one pixel (PixId=16) */
      if(pix->GetConfig()->GetPixelID()==pixel) {
	std::cout << "PIXEL = " << pix->GetConfig()->GetPixelID() << std::endl;
	if(telcoeffs)
	  {
	    const Calibration::PixelCalibCoeffs *pixcoeff = telcoeffs->GetMonitor(pix);
	    std::cout << pixcoeff->fPixelBroken[0] << "  " <<   pixcoeff->fPed[0]  << std::endl;
	    std::cout << pixcoeff->fPixelBroken[1] << "  " <<   pixcoeff->fPed[1]  << std::endl;
	  }
	
	if(telnsb)
	  {
	    const ParisCalibration::PixelLumi *pixnsb = telnsb->GetMonitor(pix);
	    std::cout << "NSB = " << pixnsb->fHviNsb << std::endl;
	  }
	std::cout << std::endl;
      }
    }
  }
}

void writeBpxCoeffs(int telid)
{
  if(!hessarray) {
    std::cout << "writeBpxCoeffs : NULL pointer HESSArray" << std::endl;
    return;
  }
  std::cout << "Received CT" << telid << "_BrokenPixel event" << std::endl; 
  const Sash::PointerSet<Sash::Telescope> & tels = runheader->GetTelsInRun();
  Sash::PointerSetIterator<Sash::Telescope> tel(tels.begin());
  for(;tel!=tels.end();++tel) {
    const Calibration::TelescopeCalibCoeffs * telcoeffs = (*tel)->Get((Calibration::TelescopeCalibCoeffs*)0);
    Sash::NonConstPointer<Sash::Pixel> begin = (*tel)->beginPixel();
    Sash::NonConstPointer<Sash::Pixel> end = (*tel)->endPixel();
    for(Sash::NonConstPointer<Sash::Pixel> pix = begin; pix!=end; ++pix) {
      if(telcoeffs) {
	const Calibration::PixelCalibCoeffs *pixcoeff = telcoeffs->GetMonitor(pix);
	/*Print the broken pixels */
	if(pixcoeff->fPixelBroken[0] || pixcoeff->fPixelBroken[1]) {
	  std::cout << "PIXEL = " << pix->GetConfig()->GetPixelID() << " is broken" << std::endl;
	  std::cout << pixcoeff->fPixelBroken[0] << "  " <<   pixcoeff->fPed[0]  << std::endl;
	  std::cout << pixcoeff->fPixelBroken[1] << "  " <<   pixcoeff->fPed[1]  << std::endl;	
	}
      }
    }
    std::cout << std::endl;
  }
}

