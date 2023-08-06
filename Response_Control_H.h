#ifndef  Response_Control_H
#define  Response_Control_H 

#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>

#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TList.h"
#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"
#include "TChain.h"
#include "TF2.h"
#include "TH3F.h"
#include "TEnv.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCanvas.h"

//+I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+
class Response_vs_E_true {
public:
    Float_t eta_min,eta_max;
    TBranch* b_pt;
    TBranch* b_eta;
    TBranch* b_E;
    TBranch* b_NPV;
    TBranch* b_mu;
    TBranch* b_weight;
    TBranch* b_weight_tot;
    TBranch* b_pt_true;
    TBranch* b_jet_area;
    TBranch* b_jet_true_E;
    TBranch* b_rho;   
    TBranch* b_njet;   
    std::vector<float> *pt = 0;
    std::vector<float> *jet_eta; //jet eta
    std::vector<float>* pt_true;
    std::vector<float>* jet_area;
    std::vector<float> *jet_E;
    std::vector<float> *jet_true_E;
    Int_t           NPV;
    Int_t           njet;
    Float_t         mu;
    Double_t        weight_tot;
    // Float_t        weight_tot;
    TH2F * R_vs_E_true;
    Float_t         rho;
    Float_t pTmuNPVcorrection=0.0, calib4=0.0, term2=0.0;
    std::vector<double> NPVTerm =    {0.040, -0.002, 0.073, -0.065, -0.001, 0.103, 0.467, -0.212, 0.051, 0.345, 0.054, 0.015};
    std::vector<double> MuTerm =      {-0.012, -0.007, -0.054, -0.022, -0.074, -0.098, -0.345, 0.119, -0.056, -0.234, -0.019, 0.037};
    std::vector<double> ResidualAbsEtaBins = { 0, 0.9, 1.2, 1.5, 1.8, 2.4, 2.8, 3.2, 3.5, 4.0, 4.3, 6.0 };
    TChain* ch ;
    std::string TreeName = "IsolatedJet_tree";


    unsigned N_Jets, N_Jets_reco, N_Jets_true;
    Response_vs_E_true();
    ~Response_vs_E_true();
    TChain* Set_Up_TChain(std::string path_to_files, std::string TTree_name);

    void Response_Control(std::string path_to_files, std::string TTree_name , int increment , std::string Step, float eta_min , float eta_max, std::string LegendName , int Num_Leading_Jets  );



    void Create_Folder(std::string name);
    char * FSTCA(std::string & name);
    char * From_String_To_Char_Array( std::string & name);
    bool is_file_exist(const char *fileName);



    std::vector<double> ComputeOffsets(std::vector<double> etaBins, std::vector<double> term) ;
    double Interpolation (std::vector<double> etaBins, std::vector<double> term, std::vector<double> offset, Float_t eta) ;
    float correction1D(float pt, float eta , float mu, int NPV) ;
    float correctedPt1D(float pt, float eta, float area, float rho, float mu, int NPV );
    float correctionFactor(const float pt, const  float eta, const  float area, const  float rho, const  float mu, const  int NPV)  ;


    float Reco_Factor(float pt,float eta,float area , float rho , float mu, int NPV);
    float Area_Factor(float pt,float eta,float area , float rho , float mu, int NPV);
    template<class T>
    void displayVector(const std::vector<T> & l_vector);

    static const Int_t E_Bins_N = 28;
    const Double_t  E_Bins[E_Bins_N+1] = {10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 100, 120, 150, 200, 240, 300, 400, 500, 600, 800, 1000, 1200, 1500, 2000, 2500, 3000, 3500, 4000, 5500 };

    std::vector<std::function<void()>> functions;
    void displayProgress(const uint indexEntry, const uint N_Entries, const float percent);
    template<class T>    
    void fill2DHistograms(std::shared_ptr<std::vector<TH2D>> const &histogramVector, const std::vector<T> &bins, const T &threshold, const Double_t &x, const Double_t &y, const Double_t &weight);
    std::shared_ptr<std::vector<TH2D>> R_vs_E_true_vector = std::make_shared<std::vector<TH2D>>();
    std::vector<float> etaBins;
    uint numberOfBins;
    void makeBinning(float eta_min , float eta_max );
};
//+I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+
template<class T>
int findBinIndex (const std::vector<T> &binning, const T value) {
    typename std::vector<T>::const_iterator threshold = std::find_if(binning.begin(), binning.end(), [&value](const T &iter){
        return iter > value;
    });
    if(threshold == binning.end()){
        return -1;
    }
    else if (threshold == binning.begin()){
        return -2;
    }
    return std::distance(binning.begin(), threshold)-1;
}
//+I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+
float roundUp(float var){
    // we use array of chars to store number
    // as a string.
    char str[40];
    // Print in string the value of var
    // with two decimal point
    sprintf(str, "%.2f", var);
    // scan string value in var
    sscanf(str, "%f", &var);
    return var;
}
//+I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+
template<class T>
void Response_vs_E_true::displayVector(const std::vector<T> & l_vector){
  std::cout <<"{";
  for(const T & i : l_vector){
      std::cout << i<< " ";
  }
  std::cout <<"}";
  std::cout << " \n";
}
//+I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+
template<class T>
void Response_vs_E_true::fill2DHistograms(std::shared_ptr<std::vector<TH2D>> const &histogramVector, const std::vector<T> &bins, const T &threshold, const Double_t &x, const Double_t &y, const Double_t &weight){
    if (  (*histogramVector).size() != bins.size()){
        std::cout<<"ERRROR: wrong binning\n";
        exit(1);
    }
    const int binIndex = findBinIndex(bins, threshold);
    if(binIndex < 0){return;}
    (*histogramVector)[binIndex].Fill(x,y,weight);
    // std::shared_ptr<Base> const& s

    
    return;
}
//+I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+
void Response_vs_E_true::displayProgress(const uint indexEntry, const uint N_Entries, const float percent){
    const uint step = percent * N_Entries;
    static int prevStep;
    // std::cout <<indexEntry / step <<" "<<prevStep <<"\n";
    // indexEntry / step

    if(prevStep != indexEntry / step) std::cout<<"Done: "<< ((float) indexEntry / (float)N_Entries)* 100.0<<"%\n";
    prevStep = indexEntry / step;
    return;
}
//+I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+
void Response_vs_E_true::makeBinning(float eta_min, float eta_max ){
    if (eta_max < eta_min){
        float temp = eta_min;
        eta_min = eta_max;
        eta_max = temp;
    }
    float step = static_cast<float>( abs(eta_max - eta_min)/static_cast<Float_t>(numberOfBins));
    step = roundUp(step);
    std::cout<<step<<" "<<(eta_max - eta_min)<<"\n";
    for (  int i = 0 ; i < numberOfBins; i ++ )
        etaBins.push_back(roundUp(eta_min+i*step));
    std::cout<< etaBins.size()<<": ";
    displayVector(etaBins);
    TString nameR_vs_Etrue = "response vs E true_";
    // R_vs_E_true = new TH2F("response v e true","response v e true",  E_Bins_N,E_Bins, 180,-1.2,2.4);
    for ( uint i = 0; i < etaBins.size(); i++){
        nameR_vs_Etrue += std::to_string(i);
        TH2D histo(nameR_vs_Etrue, nameR_vs_Etrue, E_Bins_N, E_Bins, 180,-1.2,2.4);
        (*R_vs_E_true_vector).push_back(histo);
        nameR_vs_Etrue = "response vs E true";
    }
    std::cout<< etaBins.size()<<": "<< (*R_vs_E_true_vector).size()<<"\n";
    
    return;
}
//+I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+
Response_vs_E_true::Response_vs_E_true(){
    N_Jets = 0;
    N_Jets_reco = 0;
    N_Jets_true = 0;
    numberOfBins = 91;
}
//+I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+
char * Response_vs_E_true::From_String_To_Char_Array( std::string & name){
    char * char_name[500];
    for (int i =0; i < name.size();i++){
        //cout<<name.at(i);
        char_name[i] = & name.at(i);
    }
    
    return (*char_name);
}
//+I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+
char * Response_vs_E_true::FSTCA(std::string & name){
    return From_String_To_Char_Array(  name);
}
//+I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+
bool Response_vs_E_true::is_file_exist(const char *fileName)
{
    std::ifstream infile(fileName);
    bool exists = infile.good();
    return exists;
}
//+I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+
std::vector<double> Response_vs_E_true::ComputeOffsets(std::vector<double> etaBins, std::vector<double> term) {
    std::vector<double> offset;
    offset.push_back(term.at(0));  
    for(unsigned int i=1; i<etaBins.size();++i){
        double espace = etaBins.at(i)-etaBins.at(i-1);
        double ioffs = offset.at(i-1)+term.at(i)*espace;
        offset.push_back(ioffs);
    }
    return offset;
}
//+I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+
double Response_vs_E_true::Interpolation (std::vector<double> etaBins, std::vector<double> term, std::vector<double> offset, Float_t eta) {
    double correction=0;
    //  std::cout<<"eta : "<<eta<<std::endl;
    for(unsigned int i=0; i<(etaBins.size()-1);++i){
        if (eta >= etaBins.at(i) && eta < etaBins.at(i+1)){
            correction = offset.at(i) + (eta-etaBins.at(i))*term.at(i+1); // computing value of the linear function at bin i for a given value of eta
            //  std::cout<<"inside interpolation if"<<std::endl;
           break;
        }
    }
    return correction;      
}
//+I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+
float Response_vs_E_true::correction1D(float pt, float eta , float mu, int NPV)  {
    std::vector<double> NPVoffset;
    std::vector<double> Muoffset;           
    Double_t corr1D = 0.0;
    Double_t NPVtermAtEta, MutermAtEta;
    NPVoffset = ComputeOffsets(ResidualAbsEtaBins,NPVTerm);
    Muoffset = ComputeOffsets(ResidualAbsEtaBins,MuTerm);
    NPVtermAtEta = Interpolation(ResidualAbsEtaBins,NPVTerm,NPVoffset,abs(eta ));
    MutermAtEta = Interpolation(ResidualAbsEtaBins,MuTerm,Muoffset,abs(eta));
    corr1D = NPVtermAtEta*(NPV-1)+MutermAtEta*mu;
    return corr1D;
}
//+I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+
float Response_vs_E_true::correctedPt1D(float pt, float eta, float area, float rho, float mu, int NPV ) {
    float areaCorr; 
    areaCorr = area*rho*1.0;        
    float calibration1D = correction1D(pt - areaCorr , eta, mu, NPV);  
    return pt - areaCorr - calibration1D;      
}
//+I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+ +I+
float Response_vs_E_true::correctionFactor(const float pt, const  float eta, const  float area, const  float rho, const  float mu, const  int NPV)  {
    //std::cout<<"m_useRho = "<<m_useRho<<std::endl;
    // std::cout<<"1D\n";
    return correctedPt1D(pt,eta,area,rho,mu,NPV)/pt;
}
Response_vs_E_true::~Response_vs_E_true(){
    std::cout << "Finished the loop\n";
}
#endif
