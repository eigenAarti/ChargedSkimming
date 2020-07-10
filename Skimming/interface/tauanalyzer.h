#ifndef TAUANALYZER_H
#define TAUANALYZER_H

#include <ChargedSkimming/Skimming/interface/baseanalyzer.h>
#include <DataFormats/PatCandidates/interface/Tau.h>

typedef edm::EDGetTokenT<std::vector<pat::Tau>> tToken;

//Tau class to be saved in tree		//NO STRUCT IN electronanalyzer.h?
struct Tau {
    Bool_t isMedium;
    Bool_t isTight;

    Bool_t isTriggerMatched;

    Float_t recoSF = 1.;
    Float_t mediumMvaSF = 1.;
    Float_t tightMvaSF = 1.;

    Float_t charge;
    Float_t isolation; 

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
        TF1* tauIdSFhistVL;
        TF1* tauIdSFhistL;
        TF1* tauIdSFhistM;
        TF1* tauIdSFhistT;
        TF1* tauIdSFhistVT;
        TF1* tauIdSFhistVVT;

	TH1F* antiMuSFhistL;
	TH1F* antiMuSFhistT;

	TH1F* antiEleSFhistVL;
	TH1F* antiEleSFhistL;
	TH1F* antiEleSFhistM;
	TH1F* antiEleSFhistT;
	TH1F* antiEleSFhistVT;

        //Kinematic cut criteria
        int era;
        float ptCut;
        float etaCut;

        //EDM Token for MINIAOD analysis
        tToken tauToken;
       // trigObjToken triggerObjToken;		//WHY IS THIS NOT NEEDED ANYMORE?
        genPartToken genParticleToken;

	//Vector with output varirables of the output tree
        std::map<std::string, std::vector<float>&> floatVar;
        std::map<std::string, std::vector<char>&> charVar;

        std::vector<float> Pt, Eta, Phi, ID_SF_VLoose, ID_SF_Loose, ID_SF_Medium, ID_SF_Tight, ID_SF_VTight, ID_SF_VVTight,
		 antiMuSF_Loose, antiMuSF_Tight, antiEleSF_VLoose, antiEleSF_Loose, antiEleSF_Medium, antiEleSF_Tight, antiEleSF_VTight;

        std::vector<char> Charge, DM, AntiEl, AntiMu, DMoldMVA2017, DMnewMVA2017, IdDM, isFrom_h;

        char nTaus;		

        //TTreeReader Values for NANO AOD analysis
        std::unique_ptr<TTreeReaderArray<float>> tauPt, tauEta, tauPhi, tauIso;
        std::unique_ptr<TTreeReaderArray<int>> tauCharge, tauDM;	
       // std::unique_ptr<TTreeReaderArray<unsigned int>> nTau;
        std::unique_ptr<TTreeReaderArray<unsigned char>> tauAntiEl, tauAntiMu, tauDM2017new, tauDM2017old;
        std::unique_ptr<TTreeReaderArray<bool>> tauIdDM;

    /*  std::unique_ptr<TTreeReaderArray<int>> nTaus;
        std::unique_ptr<TTreeReaderArray<int>> AntiEl;
        std::unique_ptr<TTreeReaderArray<int>> AntiMu;
        std::unique_ptr<TTreeReaderArray<int>> DMnew;
        std::unique_ptr<TTreeReaderArray<int>> DMold;*/

    public:
        TauAnalyzer(const int &era, const float &ptCut, const float &etaCut, tToken& tauToken, genPartToken& genParticleToken);
        TauAnalyzer(const int &era, const float &ptCut, const float &etaCut, TTreeReader& reader);

	int SetGenParticles(const int &i, const int &pdgID);
        void BeginJob(std::vector<TTree*>& trees, bool &isData, const bool& isSyst=false);
        void Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event);
        void EndJob(TFile* file);
};

#endif
