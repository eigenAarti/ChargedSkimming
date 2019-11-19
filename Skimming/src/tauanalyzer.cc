#include <ChargedAnalysis/Skimming/interface/tauanalyzer.h>

TauAnalyzer::TauAnalyzer(const int &era, const float &ptCut, const float &etaCut, tToken& tauToken, trigObjToken& triggerObjToken, genPartToken& genParticleToken):	//for miniAOD
    BaseAnalyzer(),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut),
    tauToken(tauToken),
    triggerObjToken(triggerObjToken),
    genParticleToken(genParticleToken)
    {}

TauAnalyzer::TauAnalyzer(const int &era, const float &ptCut, const float &etaCut, TTreeReader& reader):	//for nanoAOD
    BaseAnalyzer(&reader),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut)
    {}

int TauAnalyzer::SetGenParticles(const int &i, const int &pdgID){

 bool isgenMatched = isNANO ? tauGenIdx->At(i) != -1 : false;
//std::cout << isgenMatched << std::endl;

  if(isgenMatched){
     int index = 0;
     index = FirstCopy(tauGenIdx->At(i), pdgID);	//gives the mother index

     int motherID = abs(genID->At(genMotherIdx->At(index)));		//getting motherID from the mother index

     if(motherID == 25){
         
         if(isNANO) index = FirstCopy(genMotherIdx->At(index), 25);	
         else std::cout<<"Not NanoAOD!";

         int motherID = abs(genID->At(genMotherIdx->At(index)));  
	//std::cout << motherID << std::endl;
         if(motherID == 37){
             return 1;
         }

         else{
             return 2;
         }
     }
  }
	
   return -1.;
}

void TauAnalyzer::BeginJob(std::vector<TTree*>& trees, bool &isData){		
    //SF files
    tauIdSFfiles = {
                    {2017, filePath + "/tauSF/TauID_SF_pt_MVAoldDM2017v2_2017ReReco.root"},
    };

    antiMuSFfiles = {
                    {2017, filePath + "tauSF/TauID_SF_eta_antiMu3_2017ReReco.root"},
    };

    antiEleSFfiles = {
                    {2017, filePath + "tauSF/TauID_SF_eta_antiEleMVA6_2017ReReco.root"},
    };

    //Set data bool
    this->isData = isData;

    //Hist with scale factors
    TFile* tauIdSFfile = TFile::Open(tauIdSFfiles[era].c_str());
    tauIdSFhist = (TF1*)tauIdSFfile->Get("VLoose_cent");

    TFile* antiMuSFfile = TFile::Open(antiMuSFfiles[era].c_str());
    antiMuSFhist = (TH1F*)antiMuSFfile->Get("Loose");

    TFile* antiEleSFfile = TFile::Open(antiEleSFfiles[era].c_str());
    antiEleSFhist = (TH1F*)antiEleSFfile->Get("VLoose");


    //Initiliaze TTreeReaderValues then using NANO AOD
    if(isNANO){
        tauPt = std::make_unique<TTreeReaderArray<float>>(*reader, "Tau_pt");
        tauEta = std::make_unique<TTreeReaderArray<float>>(*reader, "Tau_eta");
        tauPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "Tau_phi");
        tauCharge = std::make_unique<TTreeReaderArray<int>>(*reader, "Tau_charge");
        tauAntiEl = std::make_unique<TTreeReaderArray<unsigned char>>(*reader, "Tau_idAntiEle");
	tauAntiMu = std::make_unique<TTreeReaderArray<unsigned char>>(*reader, "Tau_idAntiMu");
	tauDM2017new = std::make_unique<TTreeReaderArray<unsigned char>>(*reader, "Tau_idMVAnewDM2017v2");
	tauDM2017old = std::make_unique<TTreeReaderArray<unsigned char>>(*reader, "Tau_idMVAoldDM2017v2");
	tauDM = std::make_unique<TTreeReaderArray<int>>(*reader, "Tau_decayMode");
	tauIdDM = std::make_unique<TTreeReaderArray<bool>>(*reader, "Tau_idDecayMode");

	//Set TTreeReader for genpart and trigger obj from baseanalyzer    
        SetCollection(this->isData);
    }

    //Set output names
    floatNames = {"E", "Px", "Py", "Pz", "Charge", "ID_SF_VLoose", "antiMuSF_Loose", "antiEleSF_VLoose", "isFrom_h"};
    boolNames = { "IdDM"};
    intNames = {"Charge_int", "DM", "AntiEl", "AntiMu", "DMoldMVA2017", "DMnewMVA2017"};

    floatVariables = std::vector<std::vector<float>>(floatNames.size(), std::vector<float>());
    boolVariables = std::vector<std::vector<bool>>(boolNames.size(), std::vector<bool>());
    intVariables = std::vector<std::vector<int>>(intNames.size(), std::vector<int>());

    //Set Branches of output tree
    for(TTree* tree: trees){
        for(unsigned int i=0; i<floatVariables.size(); i++){
            tree->Branch(("Tau_" + floatNames[i]).c_str(), &floatVariables[i]);
        }

        for(unsigned int i=0; i<boolVariables.size(); i++){
            tree->Branch(("Tau_" + boolNames[i]).c_str(), &boolVariables[i]);
        }

        for(unsigned int i=0; i<intVariables.size(); i++){
            tree->Branch(("Tau_" + intNames[i]).c_str(), &intVariables[i]);
        }
    }
}

void TauAnalyzer::Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event){
    //Clear variables vector
    for(std::vector<float>& variable: floatVariables){
        variable.clear();
    }

    for(std::vector<bool>& variable: boolVariables){
        variable.clear();
    }

    for(std::vector<int>& variable: intVariables){
        variable.clear();
    }

    edm::Handle<std::vector<pat::Tau>> taus;		
    edm::Handle<std::vector<pat::TriggerObjectStandAlone>> trigObjects;		//Not relevant for nanoAOD
    edm::Handle<std::vector<reco::GenParticle>> genParts;

    //Get Event info if using MINIAOD
    if(!isNANO){
        event->getByToken(tauToken, taus);
        event->getByToken(triggerObjToken, trigObjects);
    }

    float tauSize = isNANO ? tauPt->GetSize() : taus->size();
    
    //Loop over all taus
    for(unsigned int i = 0; i < tauSize; i++){
        float pt = isNANO ? tauPt->At(i) : 1.;
        float eta = isNANO ? tauEta->At(i) : 1.;
        float phi = isNANO ? tauPhi->At(i) : 1.;

        if(pt > ptCut && abs(eta) < etaCut){
            TLorentzVector lVec;
            lVec.SetPtEtaPhiM(pt, eta, phi, 1776.86*1e-3);	

	    //Tau four-momentum components
            floatVariables[0].push_back(lVec.E());   //Energy
            floatVariables[1].push_back(lVec.Px());  //Px
            floatVariables[2].push_back(lVec.Py());  //Py
            floatVariables[3].push_back(lVec.Pz());  //Pz     

            floatVariables[4].push_back(isNANO ? tauCharge->At(i) : 1.);  //charge  
		
            //Tau ID
            boolVariables[0].push_back(isNANO ? tauIdDM->At(i) : 1.);  //Medium MVA ID
 
	    if((int (tauAntiEl->At(i))&0) == 0){
		intVariables[2].push_back(isNANO ? int (tauAntiEl->At(i)) : 1.);	//with bitmask, just for demonstration
	    }
	    intVariables[0].push_back(isNANO ? tauCharge->At(i) : 1.);
	    intVariables[1].push_back(isNANO ? tauDM->At(i) : 1.);
	    //intVariables[2].push_back(isNANO ? int (tauAntiEl->At(i)) : 1.);
	    intVariables[3].push_back(isNANO ? int (tauAntiMu->At(i)) : 1.);		//without bitmask
	    intVariables[4].push_back(isNANO ? int (tauDM2017old->At(i)) : 1.);
	    intVariables[5].push_back(isNANO ? int (tauDM2017new->At(i)) : 1.); 

            if(!isData){
               //Fill scale factors
                floatVariables[5].push_back(tauIdSFhist->operator()(pt));
                floatVariables[6].push_back(antiMuSFhist->GetBinContent(antiMuSFhist->FindBin(abs(eta))));
                floatVariables[7].push_back(antiEleSFhist->GetBinContent(antiEleSFhist->FindBin(abs(eta))));

                //Save gen particle information
                floatVariables[8].push_back(SetGenParticles(i, 15));	
            } 

            //Fill tau in collection
        } 
    }

    for(CutFlow &cutflow: cutflows){			
        if(cutflow.nMinTau <= floatVariables[0].size()){
            if(cutflow.nMinTau!=0 and cutflow.passed){
                std::string cutName("N_{tau} >= " + std::to_string(cutflow.nMinTau) + " (no iso/ID req)");
                cutflow.hist->Fill(cutName.c_str(), cutflow.weight);
            }
        }

        else{
            cutflow.passed = false;
        }
    }
}


void TauAnalyzer::EndJob(TFile* file){
}
