#ifndef BASEANALYZER_H
#define BASEANALYZER_H

#include <memory>
#include <map>
#include <vector>
#include <cmath>
#include <bitset>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <Rtypes.h>
#include <Math/Vector4Dfwd.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <FWCore/Framework/interface/Event.h>

#include <DataFormats/Common/interface/TriggerResults.h>
#include <FWCore/Common/interface/TriggerNames.h>

#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>

typedef edm::EDGetTokenT<edm::TriggerResults> trigToken;
typedef edm::EDGetTokenT<std::vector<reco::GenParticle>> genPartToken;

//Struct for cutflow
struct CutFlow {
    TH1F* hist;
    Float_t weight = 1.;

    unsigned int nMinEle=0;
    unsigned int nMinMu=0;
    unsigned int nMinJet=0;
    unsigned int nMinFatjet=0;
    
    bool passed = true;
};

class BaseAnalyzer {
    protected:
        //File path for SF etc.
        std::string filePath = std::string(std::getenv("CMSSW_BASE")) + "/src/ChargedSkimming/Skimming/data/";

        std::map<int, std::map<std::string, std::pair<int, int>>> runEras = {
              {2017, {
                        {"B", {297046, 299329}}, 
                        {"C", {299368, 302029}}, 
                        {"DE", {302030, 304797}},
                        {"F", {305040, 306462}}, 
                     }
              },
        };

        //Collection which are used in several analyzers if NANO AOD is analyzed
        TTreeReader* reader = NULL;
        bool isNANO;

        std::unique_ptr<TTreeReaderValue<unsigned int>> run;

        std::unique_ptr<TTreeReaderArray<float>> trigObjPt, trigObjPhi, trigObjEta;
        std::unique_ptr<TTreeReaderArray<int>> trigObjID, trigObjFilterBit;
    
        std::unique_ptr<TTreeReaderArray<float>> genPhi, genEta, genPt, genMass;
        std::unique_ptr<TTreeReaderArray<int>> genID, genMotherIdx, genStatus, eleGenIdx, muonGenIdx;

        //Set trihObj and Gen particle collection
        void SetCollection(bool &isData);
        bool isSyst=false;

        //Check for gen particle if it is last copy
        const reco::Candidate* FirstCopy(const reco::GenParticle& part, const int& pdgID);
        const reco::Candidate* FirstCopy(const reco::Candidate* part, const int& pdgID); //MINIAOD
        int FirstCopy(const int& index, const int& pdgID); //NANOAOD

        //Match Reco to gen particles
        bool SetGenParticles(ROOT::Math::PtEtaPhiMVector &lepton, const int &i, const int &pdgID, const std::vector<reco::GenParticle>& genParticle={});

    public:
        virtual ~BaseAnalyzer(){};
        BaseAnalyzer();
        BaseAnalyzer(TTreeReader* reader);
        virtual void BeginJob(std::vector<TTree*>& trees, bool &isData, const bool& isSyst=false) = 0;
        virtual void Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event = NULL) = 0;
        virtual void EndJob(TFile* file) = 0;
};
#endif
