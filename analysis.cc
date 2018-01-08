#include <iostream> //Input and output
#include <string> //strings used to read in
#include <stdlib.h> 
#include <fstream> //Reading in - probably unecessary
#include <vector> //vectors used throughout
#include <cctype> //isdigit etc.
#include <algorithm>    // std::sort
#include <iostream>
#include <math.h>
#include "analysis.hh"


int main(int argc, char * argv[] ){
  if(argc == 3){
    std::string FileNameIn = argv[1];
    int FileEnergy =  strtol(argv[2], NULL, 10);
    Analysis Instance;//creates instance of class needed - calls constructor

    Instance.Initialise(FileNameIn,FileEnergy);//reads in values
    if(Instance.FileOpenFlag == false){ 
      Instance.Stats();//master analysis function (sorts data, calling other functions e.g. binning)
      Instance.Candidates();
      Instance.ReadOuts();
      Instance.OverMean();//Prints theta mean (actually calculated during file read)
      Instance.GRSIFile();
      Instance.ROOTfunc();

    } else {
      std::cout << "Wrong input parameters" << std::endl;
      std::cout << "Should be ./progname datafile.type int-energy-gated" << std::endl;
      std::cout << "E.G. ./kinan dir/1200+step.dat 1200" << std::endl;
    }
  }
  return 0;
}

Analysis::Analysis() {
  //initialising variables used - meanings listed when used for the first time 
  fDigitBool = false; FileOpenFlag = false; TCandidateFlag = false; PCandidateFlag = false;
  BananaFlag = false; SpotFlag = false; ScatterFlag = false;
  counter = 0; fCol = 0;
  pi = 3.14159265359;
  ThetaMax = pi; PhiMax = 2.*pi;
  MeanError = 0.5; StdDevError= 1.0;
  tHold = 0.; pHold = 0.; tCast = 0; pCast = 0;
  GatedThetaCount = 0; GatedPhiCount = 0; GatedTheta = 0.; GatedPhi = 0.;
  ThetaBins = 100, PhiBins=200;
  ThetaBinSize = ThetaMax/(double)ThetaBins, PhiBinSize = PhiMax/(double)PhiBins;
  
  //Allocating space
  for(int i=0;i<120;i++){//1 mode per segment
    BinCountMode.push_back(0);
    ThetaMode.push_back(0.);
    PhiMode.push_back(0.);
    SegMode.push_back(i+1);
    TSegMean.push_back(0);
    PSegMean.push_back(0);
    TSegStdDev.push_back(0);
    PSegStdDev.push_back(0);
    TSegOut.push_back(0);
    PSegOut.push_back(0);
    CandidatesArray[i]=0;
  }
  for(int i = 0; i < 20000; i++){//array not a vector
    BigArray[i] = 0;//space set-up and then used - no need to clear yet
  }
}

Analysis::~Analysis() {}

void Analysis::Initialise(std::string FileName, int FileEnergy){
    GateVal = (double) FileEnergy;
    std::ifstream kindata;//create input object
    kindata.open(FileName.c_str());//open data file
    
    if(kindata.is_open()) {//opens if gets in here
	std::cout << "\nFile stream open" << std::endl;
	std::cout<<"Reading the distribution from " << FileName << " ... "<<std::endl;
	
	while(fDigitBool == false){//skips a few lines if not data, counting as it goes
	  c = kindata.get();
	  if ( !std::isdigit(c) ) {
	      std::getline(kindata,fInfo);
	      counter++;
	  }else fDigitBool = true;
	}
	std::cout << "Number of text lines in data file: " << counter << std::endl;
	for(fCol = 0; std::getline(kindata, fInfo); ++fCol){}//counts data points when data found, rows to fCol
	std::cout << "Rows of data in file: " << fCol << std::endl;
	for(int i = 0; i < fCol; i++){ BigCount.push_back(0);}
	    
	kindata.clear();
	kindata.seekg(0,std::ios::beg);//reset istream to beginning of file
	

	for(fIter = 0; fIter<counter; fIter++){//skip text lines found by while loop above
	  std::getline(kindata,fInfo);
	}

	for(fIter=0; fIter<fCol; fIter++){
	  kindata >> Seg >> Theta >> Dep >> Phi;//Read in as int or double
	  SegIn.push_back(Seg); ThetaIn.push_back(Theta); DepIn.push_back(Dep); PhiIn.push_back(Phi);//Can then push back
	}
	std::cout << "File reading finished" << std::endl;
	kindata.close();
	std::cout << "File stream closed" << std::endl << std::endl;
    }else {
      std::cout << "File not opened" << std::endl << std::endl;//error condition
      FileOpenFlag = true;
    }
}

void Analysis::Stats(){//master analysis function (sorts data, calling other functions e.g. binning)

  fIter = 0;//segments start at 1 in data file so added on later
  for(fIter;fIter<10;fIter++){//iterate through rings then segs (30 total segs as reduced and mapped)
    for(fIter2 = 0; fIter2 < 3; fIter2++){
      SegT = 0.; SegP = 0.;//Holds mean values
      AngleCount = 0;//Holds counts (for the mean)
      
      for(int f=0; f<fCol; f++){//has to run through all data points
	if(SegIn[f] == fIter*12+fIter2+1){  //std::cout << "Good data " << f << " " <<fIter << std::endl;   
	  this->BinAngles(ThetaIn[f],PhiIn[f], f);//place data in bins (BigArray) if correct segment //counts if row contributes to a 1+ segment
	}//Big Array fully filled for the segment
      }//all necessary data from file for each individual seg is binned    

      for(int f=0; f<fCol; f++){//has to run through all rows
	  if(BigCount[f] == 1 && SegIn[f] == fIter*12+fIter2+1) {//if the value is part of a bin with more than one element
	    this->SegMean(ThetaIn[f], PhiIn[f]);//works out the mean per seg - running total, calculation in SegMean() if counts in bin are > 1
	  }
      }
      
      TSegMean[fIter*12+fIter2] = SegT/AngleCount;//works out the theta mean for the segment and stores it
      PSegMean[fIter*12+fIter2] = SegP/AngleCount;//works out the phi mean for the segment and stores it
      
      //std::cout << "MEANS SEGMENT: " << fIter*12+fIter2 << " T: " << TSegMean[fIter*12+fIter2] << " P: " << PSegMean[fIter*12+fIter2] << std::endl;
      
      this->Mode(fIter*12+fIter2);//work out the mode(s) per segment//ignores threshold
      
//       std::cout << "MODES SEGMENT: " << fIter*12+fIter2 << " T: " << ThetaMode[fIter*12+fIter2] << " P: " << PhiMode[fIter*12+fIter2] << std::endl;
      
      SegT = 0.; SegP=0.;//Resets values used lower down
      for(int f=0; f<fCol; f++){//has to run through all bins
	  if(BigCount[f] == 1 && SegIn[f] == fIter*12+fIter2+1) {//if the value is part of a bin with more than one element
	    this->StdDev(ThetaIn[f], PhiIn[f]);//works out the mean per seg - running total, calculation in SegMean() if counts in a bin are > 1
	  }
      }
      TSegStdDev[fIter*12+fIter2] = sqrt(SegT/(AngleCount-1));//SegT/P are the total squared differences for theta/phi
      PSegStdDev[fIter*12+fIter2] = sqrt(SegP/(AngleCount-1));//uses sample Std. Dev.
      
      std::cout << "SD's SEGMENT: " << fIter*12+fIter2 << " T: " << TSegStdDev[fIter*12+fIter2] << " P: " << PSegStdDev[fIter*12+fIter2] << std::endl;
      std::cout << "MMM SEGMENT: " << fIter*12+fIter2 << " T: " << TSegMean[fIter*12+fIter2]-ThetaMode[fIter*12+fIter2] << " P: " <<PSegMean[fIter*12+fIter2] - PhiMode[fIter*12+fIter2] << std::endl;
      this->ClearBinArray();//reset BigArray after analysis of each segment
    }
  }//all data in file is binned  
  
}

void Analysis::BinAngles(double ThetaToBin, double PhiToBin, int row){//row is the current column

    tHold = ThetaToBin/ThetaBinSize;//bin of theta element
    pHold = (PhiToBin+pi)/PhiBinSize;//bin of phi element
    tCast = (int) (tHold);//int bin after truncation allowances
    pCast = (int) (pHold);//int bin after truncation allowances
//       std::cout << "Bin request " << (pCast*100)+(tCast) << " against array of " << 100*200 << " elements, at segment" << fIter << std::endl;
    if(ThetaBins*PhiBins > ((pCast*ThetaBins)+tCast)){//check if bins within the limits (so no seg fault)
      BigArray[(pCast*ThetaBins)+tCast]+=1;//iterate bin of seg'd value
      if(BigArray[(pCast*ThetaBins)+tCast] > 1) BigCount[row] = 1;
    } else {
       std::cout << "Bin request " << (pCast*ThetaBins)+(tCast) << " larger than array with " << ThetaBins*PhiBins << " elements" << std::endl;
       std::cout << "With phi: " << pHold << " bin: " << pCast << " input phi: " << PhiToBin << std::endl;
       std::cout << "PhiMax" << PhiMax << " " << 2.*pi << " PhiBinSize: " << PhiBinSize << std::endl;
    }
}

void Analysis::Mode(int seg){
    for(int o=0; o<20000; o++){//iterate through master array for each seg bin
	if(BigArray[o] > BinCountMode[seg]){//if bin is greater than those in the mode segments
	    BinCountMode[seg] = BigArray[o];//set counts of mode to max
	    ThetaMode[seg] = (o%100)*ThetaBinSize + ThetaBinSize/2.;//set mode theta to current iterative theta value
	    PhiMode[seg] =(o/100)*PhiBinSize-pi + PhiBinSize/2.;
	 }
    }
}

void Analysis::ClearBinArray(){
  for(int i=0; i<ThetaBins*PhiBins; i++){//zeroes all elements on bin array
	BigArray[i] = 0; 
   }
}

void Analysis::ReadOuts(){
 /* std::string insert, spots = "Spotty segs: ", bananas = "Banana segs: ", bad = "Bad segs: ";
  
  for(int i = 0; i < 10; i++){
    for(int ii = 0; ii < 12; ii++){
      if(CandidatesArray[i*12+ii] == 1 || CandidatesArray[i*12+ii] == 2) insert = "a";
      else insert = "not a";
      std::cout << "R" << i <<"S" << ii << " is " << insert << " candidate" << std::endl;//could count candidates
      std::cout << "Theta mode:" << ThetaMode[i*12+ii] << " Phi mode:" << PhiMode[i*12+ii] << " With counts:" << BinCountMode[i*12+ii] << std::endl;
      std::cout << "Theta mean:" << TSegMean[i*12+ii] << " Phi mean:" << PSegMean[i*12+ii] << std::endl;
      std::cout <<"Theta std. dev.:" << TSegStdDev[i*12+ii] << " Phi std. dev.:" << PSegStdDev[i*12+ii] << std::endl;
      
      if(CandidatesArray[i*12+ii] == 1) spots += "R" +  to_string(i) + "S" + to_string(ii) + " ";
      if(CandidatesArray[i*12+ii] == 2) bananas += "R" +  to_string(i) + "S" + to_string(ii) + " ";
      if(CandidatesArray[i*12+ii] == 0) bad += "R" + to_string(i) + "S" + to_string(ii) + " ";
    }
  }
  std::cout << spots << std::endl;
  std::cout << bananas << std::endl;
  std::cout << bad << std::endl;*/
  std::cout << std::endl;
 
}

void Analysis::Candidates(){//if shape found, can tell us to use mode or mean
  double ModeMean = 0.;  
  bool StdDevCheckT = false, StdDevCheckP = false, ModeMeanCheckT = false, ModeMeanCheckP = false;
  
  for(fIter = 0; fIter<10; fIter++){//iterate through rings then segs (30 total segs as reduced and mapped)
    for(fIter2 = 0; fIter2 < 3; fIter2++){
      ModeMean = fabs(ThetaMode[fIter*12+fIter2] - TSegMean[fIter*12+fIter2]);
      ModeMeanCheckT = ModeMean < MeanError;
      StdDevCheckT = TSegStdDev[fIter*12+fIter2] < StdDevError;
      if (!StdDevCheckT) std::cout << "Theta Error!" << std::endl;
      if(ModeMeanCheckT && StdDevCheckT) TCandidateFlag = true;
      
      ModeMean = fabs(PhiMode[fIter*12+fIter2] - PSegMean[fIter*12+fIter2]);
      std::cout << fIter*12+fIter2 << " Phi MM: " << ModeMean << std::endl;
      ModeMeanCheckP = ModeMean < MeanError;
      StdDevCheckP = PSegStdDev[fIter*12+fIter2] < StdDevError;
      
      if(ModeMeanCheckP && StdDevCheckP) {
	SpotFlag = true;//MM < 0.5
      } else if(ModeMean < MeanError*3. && StdDevCheckP){
	BananaFlag = true;//MM > 0.5
      } else if (ModeMean < MeanError*3. && !StdDevCheckP) {
	ScatterFlag = true;
      }
      
      if(TCandidateFlag && SpotFlag) {//expand so which values (theta and phi) for output arrays 
	CandidatesArray[fIter*12+fIter2] = 1;
	CandidatesArray[fIter*12+fIter2+3] = 1;
	CandidatesArray[fIter*12+fIter2+6] = 1;
	CandidatesArray[fIter*12+fIter2+9] = 1;
      } else if (TCandidateFlag && BananaFlag) {
	  CandidatesArray[fIter*12+fIter2] = 2;
	  CandidatesArray[fIter*12+fIter2+3] = 2;
	  CandidatesArray[fIter*12+fIter2+6] = 2;
	  CandidatesArray[fIter*12+fIter2+9] = 2;
      } else if (TCandidateFlag && ScatterFlag) {
	  CandidatesArray[fIter*12+fIter2] = 3;
	  CandidatesArray[fIter*12+fIter2+3] = 3;
	  CandidatesArray[fIter*12+fIter2+6] = 3;
	  CandidatesArray[fIter*12+fIter2+9] = 3;
      }
      
      TCandidateFlag = false; SpotFlag = false; BananaFlag = false; ScatterFlag = false;//reset for next loop
    }
  }
}

void Analysis::OverMean(){
  std::cout << "Overall Theta Mean: " << GatedTheta / (double) GatedThetaCount << std::endl;//should move - global scope //also store variables
  std::cout << "Overall Phi Mean: " << GatedPhi / (double) GatedPhiCount << std::endl;//should move - global scope //also store variables
}

void Analysis::SegMean(double ThetaToUse, double PhiToUse){
  SegT += ThetaToUse;
  SegP += PhiToUse;
  AngleCount++;
}

void Analysis::StdDev(double TtoSD, double PtoSD){
  SegT+=pow(TtoSD-TSegMean[fIter*12+fIter2],2.);
  SegP+=pow(PtoSD-PSegMean[fIter*12+fIter2],2.);//no need to iterate as gated counts per seg unchanged
}

void Analysis::ROOTfunc(){
  std::ofstream ROOTdata;//create input object
  std::string ROOTout = "ROOTdata.dat";
  ROOTdata.open(ROOTout.c_str(), std::ofstream::out | std::ofstream::app);//open data file
    
  if(ROOTdata.is_open()) {//opens if gets in here
    for(int i = 0; i < 10; i++){
      for(int j = 0; j < 3; j++){ 
	if(TSegOut[i*12+j] != 0. || PSegOut[i*12+j] != 0.){
	  ROOTdata << GateVal << " " <<  i*12+j << " " << TSegOut[i*12+j] << " " << PSegOut[i*12+j] << std::endl;
	} else std::cout << "HELP" << std::endl;
      }
    }  
  } else std::cout << "File output stream not open" << std::endl;
  ROOTdata.close();
}

void Analysis::GRSIFile(){

  for(int i = 0; i < 10; i++){
    for(int j = 0; j < 3; j++){
      if(CandidatesArray[i*12+j] == 1) {
	  TSegOut[i*12+j]=ThetaMode[i*12+j];
	  TSegOut[i*12+j+3]=ThetaMode[i*12+j];
	  TSegOut[i*12+j+6]=ThetaMode[i*12+j];
	  TSegOut[i*12+j+9]=ThetaMode[i*12+j];
	  PSegOut[i*12+j]=PhiMode[i*12+j];
	  PSegOut[i*12+j+3]=PhiMode[i*12+j]+pi/2.;
	  if (PSegOut[i*12+j+3] < 0.) PSegOut[i*12+j+3] -= 2.*pi;
	  if (PSegOut[i*12+j+3] < -pi) PSegOut[i*12+j+3]+=2.*pi;
	  if (PSegOut[i*12+j+3] > pi) PSegOut[i*12+j+3]-=2.*pi;
	  PSegOut[i*12+j+6]=PhiMode[i*12+j]+pi;
	  if (PSegOut[i*12+j+6] < 0.) PSegOut[i*12+j+6] -= 2.*pi;
	  if (PSegOut[i*12+j+6] < -pi) PSegOut[i*12+j+6]+=2.*pi;
	  if (PSegOut[i*12+j+6] > pi) PSegOut[i*12+j+6]-=2.*pi;
	  PSegOut[i*12+j+9]=PhiMode[i*12+j]-pi/2.;
	  if (PSegOut[i*12+j+9] < -pi) PSegOut[i*12+j+9]+=2.*pi;
	  if (PSegOut[i*12+j+9] > pi) PSegOut[i*12+j+9]-=2.*pi;
      }else if(CandidatesArray[i*12+j] == 2) {
	  TSegOut[i*12+j]=TSegMean[i*12+j];
	  TSegOut[i*12+j+3]=TSegMean[i*12+j];
	  TSegOut[i*12+j+6]=TSegMean[i*12+j];
	  TSegOut[i*12+j+9]=TSegMean[i*12+j];
	  PSegOut[i*12+j]=PSegMean[i*12+j];
	  PSegOut[i*12+j+3]=PSegMean[i*12+j]+pi/2.;
	  if (PSegOut[i*12+j+3] < 0.) PSegOut[i*12+j+3] -= 2.*pi;
	  if (PSegOut[i*12+j+3] < -pi) PSegOut[i*12+j+3]+=2.*pi;
	  if (PSegOut[i*12+j+3] > pi) PSegOut[i*12+j+3]-=2.*pi;
	  PSegOut[i*12+j+6]=PSegMean[i*12+j]+pi;
	  if (PSegOut[i*12+j+6] < 0.) PSegOut[i*12+j+6] -= 2.*pi;
	  if (PSegOut[i*12+j+6] < -pi) PSegOut[i*12+j+6]+=2.*pi;
	  if (PSegOut[i*12+j+6] > pi) PSegOut[i*12+j+6]-=2.*pi;
	  PSegOut[i*12+j+9]=PSegMean[i*12+j]-pi/2.;
	  if (PSegOut[i*12+j+9] < -pi) PSegOut[i*12+j+9]+=2.*pi;
	  if (PSegOut[i*12+j+9] > pi) PSegOut[i*12+j+9]-=2.*pi;
      }else if(CandidatesArray[i*12+j] == 3) {
	  TSegOut[i*12+j]=ThetaMode[i*12+j];
	  TSegOut[i*12+j+3]=ThetaMode[i*12+j];
	  TSegOut[i*12+j+6]=ThetaMode[i*12+j];
	  TSegOut[i*12+j+9]=ThetaMode[i*12+j];
	  PSegOut[i*12+j]=PhiMode[i*12+j];
	  PSegOut[i*12+j+3]=PhiMode[i*12+j]+pi/2.;
	  if (PSegOut[i*12+j+3] < 0.) PSegOut[i*12+j+3] -= 2.*pi;
	  if (PSegOut[i*12+j+3] < -pi) PSegOut[i*12+j+3]+=2.*pi;
	  if (PSegOut[i*12+j+3] > pi) PSegOut[i*12+j+3]-=2.*pi;
	  PSegOut[i*12+j+6]=PhiMode[i*12+j]+pi;
	  if (PSegOut[i*12+j+6] < 0.) PSegOut[i*12+j+6] -= 2.*pi;
	  if (PSegOut[i*12+j+6] < -pi) PSegOut[i*12+j+6]+=2.*pi;
	  if (PSegOut[i*12+j+6] > pi) PSegOut[i*12+j+6]-=2.*pi;
	  PSegOut[i*12+j+9]=PhiMode[i*12+j]-pi/2.;
	  if (PSegOut[i*12+j+9] < -pi) PSegOut[i*12+j+9]+=2.*pi;
	  if (PSegOut[i*12+j+9] > pi) PSegOut[i*12+j+9]-=2.*pi;
      }else {	  
	  TSegOut[i*12+j] = 0.;
	  TSegOut[i*12+j+3] = 0.;
	  TSegOut[i*12+j+6] = 0.;
	  TSegOut[i*12+j+9] = 0.;
	  PSegOut[i*12+j] = 0.;
	  PSegOut[i*12+j+3] = 0.;
	  PSegOut[i*12+j+6] = 0.;
	  PSegOut[i*12+j+9] = 0.;
      }
    }
  }
   
  std::ofstream GRSIdata;//create input object
  std::string GRSIout = "GRSIdata.dat";
  GRSIdata.open(GRSIout.c_str(), std::ofstream::out | std::ofstream::app);//open data file
    
  if(GRSIdata.is_open()) {//opens if gets in here
    for(int i = 0; i < 120; i++){
      GRSIdata << GateVal << " " <<  SegMode[i]-1 << " " << TSegOut[i] << " " << PSegOut[i] << " " << CandidatesArray[i] << std::endl;
    }
  } else std::cout << "File output stream not open" << std::endl;
  GRSIdata.close();
}