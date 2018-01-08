#include <iostream> //Input and output
#include <string> //strings used to read in
#include <stdlib.h> 
#include <fstream> //Reading in - probably unecessary
#include <vector> //vectors used throughout
#include <cctype> //isdigit etc.

#include<sstream> //hack

int main(int, char**);

class Analysis
{
public:
    Analysis();
    ~Analysis();
  
    void Initialise(std::string, int);
    void Stats();
    void Mode(int);
    void BinAngles(double, double, int);
    void ClearBinArray();
    void ReadOuts();
    void Candidates();
    void OverMean();
    void StdDev(double, double);
    void SegMean(double,double);
    void ROOTfunc();
    void GRSIFile();
    
    //istream vars
    std::string fInfo;
    bool fDigitBool;
    char c, FileEnergy;
    double GateVal;//Energy of file

    //file input vars;
    int Seg;
    double Dep, Phi, Theta, MeanError, StdDevError; 
    
    //flags
    bool FileOpenFlag;
    
    //iterators
    int counter, fCol, fIter, fIter2;
    
    //binning vars (fixed)
    int ThetaBins, PhiBins;
    double pi, ThetaBinSize, PhiBinSize, ThetaMax, PhiMax;
    
    //binning vars (dynamic)
    double tHold, pHold;
    int tCast, pCast;
    
    int GatedThetaCount, GatedPhiCount;
    double GatedTheta, GatedPhi;
    
    //Stats vars (mean and std dev)
    double SegT, SegP;
    int AngleCount;
        
    //holds values of modes
    std::vector<int> SegMode, BinCountMode;
    std::vector<double> ThetaMode, PhiMode, TSegMean, PSegMean, TSegStdDev, PSegStdDev;
    std::vector<double> TSegOut, PSegOut;
    std::vector<bool> BigCount;
    
    //Shape Flags
    bool TCandidateFlag, PCandidateFlag, BananaFlag, SpotFlag, ScatterFlag;
    
    //vectors
    std::vector<int> SegIn;
    std::vector<double> ThetaIn, PhiIn, DepIn;
    int BigArray[20000];
    int CandidatesArray[120];
    
    
    template <typename T>
    std::string to_string(T value)
    {
      //create an output string stream
      std::ostringstream os ;

      //throw the value into the string stream
      os << value ;

      //convert the string stream into a string and return
      return os.str() ;
    }
    
    
    
};