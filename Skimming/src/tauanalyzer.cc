#include <ChargedSkimming/Skimming/interface/tauanalyzer.h>

TauAnalyzer::TauAnalyzer(const int &era, const float &ptCut, const float &etaCut, tToken& tauToken, genPartToken& genParticleToken):	//for miniAOD
    BaseAnalyzer(),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut),
    tauToken(tauToken),
    //triggerObjToken(triggerObjToken),
    genParticleToken(genParticleToken)
    {}

TauAnalyzer::TauAnalyzer(const int &era, const float &ptCut, const float &etaCut, TTreeReader& reader):	//for nanoAOD
    BaseAnalyzer(&reader),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut)
    {}

int TauAnalyzer::SetGenParticles(const int &i, const int &pdgID){		//IF NOT WORKING, JUST MOVE ON

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

void TauAnalyzer::BeginJob(std::vector<TTree*>& trees, bool &isData, const bool& isSyst){		
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

    //Hist with scale factors		//HOW DID I FIND THESE TFILE TYPES?							
    TFile* tauIdSFfile = TFile::Open(tauIdSFfiles[era].c_str());	//for tau ID
    tauIdSFhistVL = (TF1*)tauIdSFfile->Get("VLoose_cent");
    tauIdSFhistL = (TF1*)tauIdSFfile->Get("Loose_cent");
    tauIdSFhistM = (TF1*)tauIdSFfile->Get("Medium_cent");
    tauIdSFhistT = (TF1*)tauIdSFfile->Get("Tight_cent");
    tauIdSFhistVT = (TF1*)tauIdSFfile->Get("VTight_cent");
    tauIdSFhistVVT = (TF1*)tauIdSFfile->Get("VVTight_cent");

    TFile* antiMuSFfile = TFile::Open(antiMuSFfiles[era].c_str());	//for Antimu ID
    antiMuSFhistL = (TH1F*)antiMuSFfile->Get("Loose");
    antiMuSFhistT = (TH1F*)antiMuSFfile->Get("Tight");

    TFile* antiEleSFfile = TFile::Open(antiEleSFfiles[era].c_str());	//for Antielectron ID
    antiEleSFhistVL = (TH1F*)antiEleSFfile->Get("VLoose");
    antiEleSFhistL = (TH1F*)antiEleSFfile->Get("Loose");
    antiEleSFhistM = (TH1F*)antiEleSFfile->Get("Medium");
    antiEleSFhistT = (TH1F*)antiEleSFfile->Get("Tight");
    antiEleSFhistVT = (TH1F*)antiEleSFfile->Get("VTight");


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

//Variable name mapping to branch name
    floatVar = {				
            {"Pt", Pt},
            {"Eta", Eta},
            {"Phi", Phi},
	    {"ID_SF_VLoose", ID_SF_VLoose},
	    {"ID_SF_Loose", ID_SF_Loose},
	    {"ID_SF_Medium", ID_SF_Medium},
	    {"ID_SF_Tight", ID_SF_Tight},
	    {"ID_SF_VTight", ID_SF_VTight},
	    {"ID_SF_VVTight",ID_SF_VVTight},
	    {"antiMuSF_Loose", antiMuSF_Loose},
	    {"antiMuSF_Tight", antiMuSF_Tight},
	    {"antiEleSF_VLoose", antiEleSF_VLoose},
	    {"antiEleSF_Loose", antiEleSF_Loose},
	    {"antiEleSF_Medium", antiEleSF_Medium},
	    {"antiEleSF_Tight", antiEleSF_Tight},
	    {"antiEleSF_VTight", antiEleSF_VTight},
    };

    charVar = {
            {"Charge", Charge},
	    {"DM", DM},
	    {"AntiEl", AntiEl},
	    {"AntiMu", AntiMu},
	    {"DMoldMVA2017", DMoldMVA2017},
	    {"DMnewMVA2017", DMnewMVA2017},
	    {"IdDM", IdDM},			
            {"isFrom_h", isFrom_h},		
    };

 /* if(!isSyst){
        std::map<std::string, std::vector<float>&> SFvariations = {
            {"recoSFUp", recoSFUp},
            {"recoSFDown", recoSFDown},
            {"looseSFUp", looseSFUp},
            {"looseSFDown", looseSFDown},
            {"mediumSFUp", mediumSFUp},
            {"mediumSFDown", mediumSFDown},
            {"tightSFUp", tightSFUp},
            {"tightSFDown", tightSFDown},
        }; 

        floatVar.insert(SFvariations.begin(), SFvariations.end());   
    }

    //Set output names
    floatNames = {"E", "Px", "Py", "Pz", "Charge", "ID_SF_VLoose", "ID_SF_Loose", "ID_SF_Medium", "ID_SF_Tight", "ID_SF_VTight", "ID_SF_VVTight",
		 "antiMuSF_Loose", "antiMuSF_Tight", "antiEleSF_VLoose", "antiEleSF_Loose", "antiEleSF_Medium", "antiEleSF_Tight", "antiEleSF_VTight", "isFrom_h"};  
    //make both bool and int into char
    boolNames = { "IdDM"};
    intNames = {"Charge_int", "DM", "AntiEl", "AntiMu", "DMoldMVA2017", "DMnewMVA2017"};

    floatVariables = std::vector<std::vector<float>>(floatNames.size(), std::vector<float>());
    boolVariables = std::vector<std::vector<bool>>(boolNames.size(), std::vector<bool>());
    intVariables = std::vector<std::vector<int>>(intNames.size(), std::vector<int>()); */

    //Set Branches of output tree
    for(TTree* tree: trees){
        tree->Branch("Tau_Size", &nTaus);

        for(std::pair<const std::string, std::vector<float>&>& var: floatVar){
            tree->Branch(("Tau_" + var.first).c_str(), &var.second);
        }

        for(std::pair<const std::string, std::vector<char>&>& var: charVar){
            tree->Branch(("Tau_" + var.first).c_str(), &var.second);
        }
    }
}
   /* //Set Branches of output tree
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
}*/

void TauAnalyzer::Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event){
    //Clear variables vector
    for(std::pair<const std::string, std::vector<float>&>& var: floatVar){
        var.second.clear();
    }

    for(std::pair<const std::string, std::vector<char>&>& var: charVar){
        var.second.clear();
    }

    edm::Handle<std::vector<pat::Tau>> taus;		
  //  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> trigObjects;		//Not relevant for nanoAOD
    edm::Handle<std::vector<reco::GenParticle>> genParts;

    //Get Event info if using MINIAOD
    if(!isNANO){
        event->getByToken(tauToken, taus);
       // event->getByToken(triggerObjToken, trigObjects);
    }

    float tauSize = isNANO ? tauPt->GetSize() : taus->size();
    
    //Loop over all taus
    for(unsigned int i = 0; i < tauSize; i++){
        float pt = isNANO ? tauPt->At(i) : 1.;				
        float eta = isNANO ? tauEta->At(i) : 1.;
        float phi = isNANO ? tauPhi->At(i) : 1.;

        if(pt > ptCut && abs(eta) < etaCut){
            //TLorentzVector lVec;					
            //lVec.SetPtEtaPhiM(pt, eta, phi, 1776.86*1e-3);
	
	    //float variables
            Pt.push_back(pt);   				
            Eta.push_back(eta);  
            Phi.push_back(phi); 	    

            Charge.push_back(isNANO ? tauCharge->At(i) : 1.);  //charge  
		
            //Tau ID 
            IdDM.push_back(isNANO ? tauIdDM->At(i) : 1.);  //Medium MVA ID
 
	    if((int (tauAntiEl->At(i))&0) == 0){
		AntiEl.push_back(isNANO ? int (tauAntiEl->At(i)) : 1.);	//with bitmask, only for demonstration
	    }
	    Charge.push_back(isNANO ? tauCharge->At(i) : 1.);
	    DM.push_back(isNANO ? tauDM->At(i) : 1.);
	    //AntiEl.push_back(isNANO ? int (tauAntiEl->At(i)) : 1.);
	    AntiMu.push_back(isNANO ? int (tauAntiMu->At(i)) : 1.);		//without bitmask
	    DMoldMVA2017.push_back(isNANO ? int (tauDM2017old->At(i)) : 1.);
	    DMnewMVA2017.push_back(isNANO ? int (tauDM2017new->At(i)) : 1.); 

            if(!isData){
               //Fill scale factors											
                ID_SF_VLoose.push_back(tauIdSFhistVL->operator()(pt));			//WHAT IS THE PROBLEM HERE?			
                ID_SF_Loose.push_back(tauIdSFhistL->operator()(pt));						
                ID_SF_Medium.push_back(tauIdSFhistM->operator()(pt));						
                ID_SF_Tight.push_back(tauIdSFhistT->operator()(pt));						
                ID_SF_VTight.push_back(tauIdSFhistVT->operator()(pt));						
                ID_SF_VVTight.push_back(tauIdSFhistVVT->operator()(pt));
                antiMuSF_Loose.push_back(antiMuSFhistL->GetBinContent(antiMuSFhistL->FindBin(abs(eta))));
                antiMuSF_Tight.push_back(antiMuSFhistT->GetBinContent(antiMuSFhistT->FindBin(abs(eta))));
                antiEleSF_VLoose.push_back(antiEleSFhistVL->GetBinContent(antiEleSFhistVL->FindBin(abs(eta))));
                antiEleSF_Loose.push_back(antiEleSFhistL->GetBinContent(antiEleSFhistL->FindBin(abs(eta))));
                antiEleSF_Medium.push_back(antiEleSFhistM->GetBinContent(antiEleSFhistM->FindBin(abs(eta))));
                antiEleSF_Tight.push_back(antiEleSFhistT->GetBinContent(antiEleSFhistT->FindBin(abs(eta))));
                antiEleSF_VTight.push_back(antiEleSFhistVT->GetBinContent(antiEleSFhistVT->FindBin(abs(eta))));

                //Save gen particle information
                isFrom_h.push_back(SetGenParticles(i, 15));	
            } 

            //Fill tau in collection
        } 
    }

    //number of Taus
    nTaus = Pt.size();

    for(CutFlow &cutflow: cutflows){		//change this in nanoskimmer			
        if(cutflow.nMinTau <= nTaus){
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
