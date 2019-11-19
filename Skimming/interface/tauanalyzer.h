#ifndef TAUANALYZER_H
#define TAUANALYZER_H

#include <ChargedAnalysis/Skimming/interface/baseanalyzer.h>

//Tau class to be saved in tree
struct Tau {
    TLorentzVector fourVec;
    Bool_t isMedium;
    Bool_t isTight;

    Bool_t isTriggerMatched;

    Float_t recoSF = 1.;
    Float_t mediumMvaSF = 1.;
    Float_t tightMvaSF = 1.;

    Float_t charge;
    Float_t isolation; 

    TLorentzVector genVec;
    Bool_t isgenMatched = false;
    Bool_t isFromHc = false;
};


class TauAnalyzer: public BaseAnalyzer {
    private:
        //Check if data
        bool isData;

        //Map for SF files
        std::map<int, std::string> tauIdSFfiles;
        std::map<int, std::string> antiMuSFfiles;
        std::map<int, std::string> antiEleSFfiles;

        //Hist with scale factors
        TF1* tauIdSFhist;
	TH1F* antiMuSFhist;
	TH1F* antiEleSFhist;

        //Kinematic cut criteria
        int era;
        float ptCut;
        float etaCut;

        //EDM Token for MINIAOD analysis
        tToken tauToken;
        trigObjToken triggerObjToken;
        genPartToken genParticleToken;

        //Vector with output varirables of the output tree
        std::vector<std::string> floatNames;
        std::vector<std::string> boolNames;
        std::vector<std::string> intNames;

        std::vector<std::vector<float>> floatVariables;
        std::vector<std::vector<bool>> boolVariables;
        std::vector<std::vector<int>> intVariables;

        //TTreeReader Values for NANO AOD analysis
        std::unique_ptr<TTreeReaderArray<float>> tauPt;
        std::unique_ptr<TTreeReaderArray<float>> tauEta;
        std::unique_ptr<TTreeReaderArray<float>> tauPhi;
        std::unique_ptr<TTreeReaderArray<float>> tauIso;
        std::unique_ptr<TTreeReaderArray<int>> tauCharge;
        std::unique_ptr<TTreeReaderArray<unsigned int>> nTau;
        std::unique_ptr<TTreeReaderArray<int>> tauDM;
        std::unique_ptr<TTreeReaderArray<unsigned char>> tauAntiEl;
        std::unique_ptr<TTreeReaderArray<unsigned char>> tauAntiMu;
        std::unique_ptr<TTreeReaderArray<unsigned char>> tauDM2017new;
        std::unique_ptr<TTreeReaderArray<unsigned char>> tauDM2017old;
        std::unique_ptr<TTreeReaderArray<bool>> tauIdDM;

      /*  std::unique_ptr<TTreeReaderArray<int>> nTaus;
        std::unique_ptr<TTreeReaderArray<int>> AntiEl;
        std::unique_ptr<TTreeReaderArray<int>> AntiMu;
        std::unique_ptr<TTreeReaderArray<int>> DMnew;
        std::unique_ptr<TTreeReaderArray<int>> DMold;*/

    public:
        TauAnalyzer(const int &era, const float &ptCut, const float &etaCut, tToken& tauToken, trigObjToken& triggerObjToken, genPartToken& genParticleToken);
        TauAnalyzer(const int &era, const float &ptCut, const float &etaCut, TTreeReader& reader);

	int SetGenParticles(const int &i, const int &pdgID);
        void BeginJob(std::vector<TTree*>& trees, bool &isData);
        void Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
