#include <ChargedSkimming/Skimming/interface/jetanalyzer.h>

JetAnalyzer::JetAnalyzer(const int &era, const float &ptCut, const float &etaCut, TTreeReader &reader):
    BaseAnalyzer(&reader),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut)
    {}

JetAnalyzer::JetAnalyzer(const int &era, const float &ptCut, const float &etaCut, std::vector<jToken>& jetTokens, std::vector<genjToken>& genjetTokens, mToken &metToken, edm::EDGetTokenT<double> &rhoToken, genPartToken& genParticleToken, secvtxToken& vertexToken, const std::string& systematic):
    BaseAnalyzer(),    
    era(era),
    ptCut(ptCut),
    etaCut(etaCut),
    jetTokens(jetTokens),
    metToken(metToken),
    rhoToken(rhoToken),
    genParticleToken(genParticleToken),
    vertexToken(vertexToken)
    {
        if(systematic.find("JEC") != std::string::npos){
            jecSyst = systematic;
            jecSyst.replace(jecSyst.find("JEC"), 3, "");

            if(jecSyst.find("Up") != std::string::npos){
                isUp = true;
                jecSyst.replace(jecSyst.find("Up"), 2, "");    
            }

            else{
                isUp = false;
                jecSyst.replace(jecSyst.find("Down"), 4, "");    
            }
        }

        if(systematic.find("JER") != std::string::npos){
            isJERsyst = true;

            if(systematic.find("Up") != std::string::npos) isUp = true;
            else isUp = false;
        }

        Deep = new TH1F("Deep", "Deep", 30, 0, 1);
    }


void JetAnalyzer::SetCorrector(const JetType &type, const int& runNumber){
    std::vector<JetCorrectorParameters> corrVec;

    for(std::string fileName: isData? JECDATA[era] : JECMC[era]){
        if(fileName.find("@") != std::string::npos){
            for(std::pair<std::string, std::pair<int, int>> eraNames: runEras[era]){
                if(eraNames.second.first <= runNumber and runNumber <= eraNames.second.second){
                    fileName.replace(fileName.find("@"), 1, eraNames.first);
                }
            }
        }

        fileName.replace(fileName.find("&"), 1, type == AK4 ? "AK4": "AK8");

        corrVec.push_back(JetCorrectorParameters(fileName));
    }
  
    jetCorrector[type] = new FactorizedJetCorrector(corrVec);
}

//https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorFWLite
float JetAnalyzer::CorrectEnergy(const ROOT::Math::PtEtaPhiMVector &jet, const float &rho, const float &area, const JetType &type){
    jetCorrector[type]->setJetPt(jet.Pt());
    jetCorrector[type]->setJetEta(jet.Eta());
    jetCorrector[type]->setRho(rho);
    jetCorrector[type]->setJetA(area);

    double correction = jetCorrector[type]->getCorrection();
    
    return correction;
}

//https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
//https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#L203-L263
float JetAnalyzer::SmearEnergy(const ROOT::Math::PtEtaPhiMVector &jet, const float &rho, const float &coneSize, const JetType &type, const std::vector<reco::GenJet> &genJets){
    //Configure jet SF reader
    jetParameter.setJetPt(jet.Pt()).setJetEta(jet.Eta()).setRho(rho);

    float reso = resolution[type].getResolution(jetParameter);
    float resoSF; 
    if(!isJERsyst) resoSF = resolution_sf[type].getScaleFactor(jetParameter);
    else resoSF = resolution_sf[type].getScaleFactor(jetParameter, isUp ? Variation::UP : Variation::DOWN);
        
    float smearFac = 1.; 
    float dR;
    float genPt, genPhi, genEta, genMass;
    unsigned int size;

    if(isNANO) size = (type == AK4) ? genJetPt->GetSize(): genFatJetPt->GetSize();
    else size = genJets.size();

    //Loop over all gen jets and find match
    for(unsigned int i = 0; i < size; i++){
        if(isNANO){
            genPt = (type == AK4) ? genJetPt->At(i): genFatJetPt->At(i);
            genPhi = (type == AK4) ? genJetPhi->At(i): genFatJetPhi->At(i);
            genEta = (type == AK4) ? genJetEta->At(i): genFatJetEta->At(i);
            genMass = (type == AK4) ? genJetMass->At(i): genFatJetMass->At(i);
        }

        else{
            genPt = genJets.at(i).pt();
            genPhi = genJets.at(i).phi();
            genEta = genJets.at(i).eta();
            genMass = genJets.at(i).mass();
        }

        dR = std::sqrt(std::pow(jet.Phi()-genPhi, 2) + std::pow(jet.Eta()-genEta, 2));

        //Check if jet and gen jet are matched
        if(dR < coneSize/2. and abs(jet.Pt() - genPt) < 3.*reso*jet.Pt()){
            ROOT::Math::PtEtaPhiMVector gJet(genPt, genEta, genPhi, genMass);
            genJet[type] = gJet;
            break;
        }
    }  

    //If you found gen matched 
    if(genJet[type] != ROOT::Math::PtEtaPhiMVector()){
        smearFac = 1.+(resoSF-1)*(jet.Pt() - genJet[type].Pt())/jet.Pt(); 
    }

    //If no match, smear with gaussian pdf
    else if(resoSF > 1.){
        std::default_random_engine generator;
        std::normal_distribution<> gaus(0, reso * std::sqrt(resoSF * resoSF - 1));

        smearFac = 1. + gaus(generator);
    }


    //Check if direction of jet not changed
    if(jet.Pt()*smearFac < 1e-2){
        smearFac = 1e-2/jet.Pt();
    }

    return smearFac;
}

int JetAnalyzer::SetGenParticles(ROOT::Math::PtEtaPhiMVector& validJet, const int &i, const int &pdgID, const JetType &type, const std::vector<reco::GenParticle>& genParticle){
    int nParton=0;
    bool isFromh1 = true, isFromh2 = true;

    //Check if gen matched particle exist
    if(genJet[type].Pt() != 0){
        float dR;
        
        //Find Gen particle to gen Jet
        int size=isNANO ? genPt->GetSize() : genParticle.size();

        for(int i=0; i < size; i++){
            const reco::Candidate* parton=NULL;
            int index=0;

            int ID = isNANO ? abs(genID->At(i)) : abs(genParticle.at(i).pdgId());
    
            if(ID == pdgID){
                if(isNANO) index = FirstCopy(i, pdgID);
                else parton = FirstCopy(genParticle.at(i), pdgID);
            }
    
            else continue;

            float eta, phi;
            phi = isNANO ? genPhi->At(index) : parton->phi();
            eta = isNANO ? genEta->At(index) : parton->eta();

            dR = std::sqrt(std::pow(phi - genJet[type].Phi(), 2) + std::pow(eta-genJet[type].Eta(), 2));
            float rMin = type == AK4 ? 0.3 : 0.4;

            if(dR <  rMin){
                int motherID = isNANO ? abs(genID->At(genMotherIdx->At(index))) : abs(parton->mother()->pdgId());

                if(motherID == 25){
                    const reco::Candidate* hBoson=NULL;

                    nParton++;
                                
                    if(isNANO) index = FirstCopy(genMotherIdx->At(index), 25);
                    else hBoson = FirstCopy(parton->mother(), 25);

                    int motherID = isNANO ? abs(genID->At(genMotherIdx->At(index))) : abs(hBoson->mother()->pdgId());  

                    if(motherID == 37){
                        isFromh1 = isFromh1 && true;
                        isFromh2 = isFromh2 && true;
                    }

                    else{
                        isFromh1 = isFromh1 && false;
                        isFromh2 = isFromh2 && true;
                    }
                }
            }

            if(nParton==1 and type==AK4){
                return !isFromh1 and !isFromh2 ? -1  : isFromh1 ? 1 : 2;
            }

            if(nParton==2 and type==AK8){
                return !isFromh1 and !isFromh2 ? -1  : isFromh1 ? 1 : 2;
            }
        }
    }

    return -1.;
}

void JetAnalyzer::BeginJob(std::vector<TTree*>& trees, bool &isData, const bool& isSyst){
    JECMC = {
                {2017, {filePath + "/JEC/Fall17_17Nov2017_V32_MC_L1FastJet_&PFchs.txt", 
                        filePath + "/JEC/Fall17_17Nov2017_V32_MC_L2Relative_&PFchs.txt",
                        filePath + "/JEC/Fall17_17Nov2017_V32_MC_L3Absolute_&PFchs.txt"}},
                            
    };

    JECDATA = {
                {2017, {filePath + "/JEC/Fall17_17Nov2017@_V32_DATA_L1FastJet_&PFchs.txt", 
                        filePath + "/JEC/Fall17_17Nov2017@_V32_DATA_L2Relative_&PFchs.txt",
                        filePath + "/JEC/Fall17_17Nov2017@_V32_DATA_L3Absolute_&PFchs.txt",
                        filePath + "/JEC/Fall17_17Nov2017@_V32_DATA_L2L3Residual_&PFchs.txt"}},
             
    };

    JECUNC = {
                {2017, filePath + "/JEC/Fall17_17Nov2017_V32_MC_UncertaintySources_&PFchs.txt"},
    };

    JMESF = {
                {2017, filePath + "/JME/Fall17_V3_MC_SF_&PFchs.txt"},     
    };

    JMEPtReso = {
                {2017, filePath + "/JME/Fall17_V3_MC_PtResolution_&PFchs.txt"},    
    };


    bTagSF = {
            {AK4, {
                {2017, filePath + "/btagSF/DeepFlavour_94XSF_V1_B_F.csv"},
                }
            },

            {AK8, {
                 {2017, filePath + "/btagSF/subjet_DeepCSV_94XSF_V4_B_F.csv"},
                }
            },
    };


    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94
    bTagCuts = {
            {AK4, {
                {2017, {0.051, 0.3033, 0.7489}}, //Loose, Medium, Tight
                }
            },

            {AK8, {
                {2017, {0.1522, 0.4941}}, //Loose, Medium, Tight
                }
            },
    };

    //Set data bool
    this->isData = isData;
    this->isSyst = isSyst;

    if(isNANO){
        //Initiliaze TTreeReaderValues
        fatJetPt = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_pt");
        fatJetEta = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_eta");
        fatJetPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_phi");
        fatJetMass = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_mass");
        fatJetArea = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_area");
        fatJetCSV = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_btagDeepB");
        fatJetTau1 = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_tau1");
        fatJetTau2 = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_tau2");
        fatJetTau3 = std::make_unique<TTreeReaderArray<float>>(*reader, "FatJet_tau3");

        jetPt = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_pt");
        jetEta = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_eta");
        jetPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_phi");
        jetMass = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_mass");
        jetArea = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_area");
        jetRho = std::make_unique<TTreeReaderValue<float>>(*reader, "fixedGridRhoFastjetAll");
        jetDeepBValue = std::make_unique<TTreeReaderArray<float>>(*reader, "Jet_btagDeepFlavB");
        
        metPhi = std::make_unique<TTreeReaderValue<float>>(*reader, "MET_phi");
        metPt = std::make_unique<TTreeReaderValue<float>>(*reader, "MET_pt");

        if(!this->isData){
            genJetPt = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJet_pt");
            genJetEta = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJet_eta");
            genJetPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJet_phi");
            genJetMass = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJet_mass");

            genFatJetPt = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJetAK8_pt");
            genFatJetEta = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJetAK8_eta");
            genFatJetPhi = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJetAK8_phi");
            genFatJetMass = std::make_unique<TTreeReaderArray<float>>(*reader, "GenJetAK8_mass");
            
            jetGenIdx = std::make_unique<TTreeReaderArray<int>>(*reader, "Jet_genJetIdx");
        }

        //Set TTreeReader for genpart and trigger obj from baseanalyzer
        SetCollection(this->isData);
    }
    
    for(JetType type: {AK4, AK8}){
        //Set configuration for bTagSF reader  ##https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
        calib[type] = BTagCalibration(std::to_string(type), bTagSF[type][era]);

        looseReader[type] = BTagCalibrationReader(BTagEntry::OP_LOOSE, "central", {"up", "down"});
        mediumReader[type] = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up", "down"});  
        tightReader[type] = BTagCalibrationReader(BTagEntry::OP_TIGHT, "central", {"up", "down"});  

        looseReader[type].load(calib[type],  BTagEntry::FLAV_B, (type==AK4) ? "comb": "lt");  
        mediumReader[type].load(calib[type],  BTagEntry::FLAV_B, (type==AK4) ? "comb": "lt");  
        
        if(type!=AK8){
            tightReader[type].load(calib[type],  BTagEntry::FLAV_B, (type==AK4) ? "comb": "lt");  
        }
    
        //Set configuration for JER tools
        std::string fileName = JMEPtReso[era];
        fileName.replace(fileName.find("&"), 1, type == AK4 ? "AK4": "AK8");
        resolution[type] = JME::JetResolution(fileName);

        fileName = JMESF[era];
        fileName.replace(fileName.find("&"), 1, type == AK4 ? "AK4": "AK8");
        resolution_sf[type] = JME::JetResolutionScaleFactor(fileName);

        //Set object to get JEC uncertainty
        if(jecSyst != ""){
            std::string fileName = JECUNC[era];
            fileName.replace(fileName.find("&"), 1, type == AK4 ? "AK4": "AK8");  

            jecUnc[type] = new JetCorrectionUncertainty(JetCorrectorParameters(fileName, "Total"));
        }

        else jecUnc[type] = NULL;
    }

    //Set output names
    variables = {
        {{"Jet", "E"}, E[AK4]},
        {{"Jet", "Px"}, Px[AK4]},
        {{"Jet", "Py"}, Py[AK4]},
        {{"Jet", "Pz"}, Pz[AK4]},
        {{"Jet", "loosebTagSF"}, loosebTagSF[AK4]},
        {{"Jet", "mediumbTagSF"}, mediumbTagSF[AK4]},
        {{"Jet", "tightbTagSF"}, tightbTagSF[AK4]},
        {{"Jet", "isFromh"}, isFromh[AK4]},
        {{"Jet", "FatJetIdx"}, FatJetIdx[AK4]},
        {{"FatJet", "E"}, Px[AK8]},
        {{"FatJet", "Px"}, Py[AK8]},
        {{"FatJet", "Py"}, Pz[AK8]},
        {{"FatJet", "Pz"}, E[AK8]},
        {{"FatJet", "Njettiness1"}, Njettiness1[AK8]},
        {{"FatJet", "Njettiness2"}, Njettiness2[AK8]},
        {{"FatJet", "Njettiness3"}, Njettiness3[AK8]},
        {{"FatJet", "topVsHiggs"}, topVsHiggs[AK8]},
        {{"FatJet", "QCDVsHiggs"}, QCDVsHiggs[AK8]},
        {{"FatJet", "isFromh"}, isFromh[AK8]},
        {{"JetParticle", "E"}, Px[PF]},
        {{"JetParticle", "Px"}, Py[PF]},
        {{"JetParticle", "Py"}, Pz[PF]},
        {{"JetParticle", "Pz"}, E[PF]},
        {{"JetParticle", "Vx"}, Vx[PF]},
        {{"JetParticle", "Vy"}, Vy[PF]},
        {{"JetParticle", "Vz"}, Vz[PF]},
        {{"JetParticle", "Charge"}, Charge[PF]},
        {{"JetParticle", "FatJetIdx"}, FatJetIdx[PF]},
        {{"SecondaryVertex", "E"}, Px[VTX]},
        {{"SecondaryVertex", "Px"}, Py[VTX]},
        {{"SecondaryVertex", "Py"}, Pz[VTX]},
        {{"SecondaryVertex", "Pz"}, E[VTX]},
        {{"SecondaryVertex", "Vx"}, Vx[VTX]},
        {{"SecondaryVertex", "Vy"}, Vy[VTX]},
        {{"SecondaryVertex", "Vz"}, Vz[VTX]},
        {{"SecondaryVertex", "FatJetIdx"}, FatJetIdx[VTX]},
    };

    bools = {
            {"isLooseB", isLooseB},
            {"isMediumB", isMediumB},
            {"isTightB", isTightB},
    };

    if(!isSyst){
        std::map<std::pair<std::string, std::string>, std::vector<float>&> SFvariations = {
            {{"Jet", "loosebTagSFSFUp"}, loosebTagSFUp[AK4]},
            {{"Jet", "loosebTagSFDown"}, loosebTagSFDown[AK4]},
            {{"Jet", "mediumbTagSFSFUp"}, mediumbTagSFUp[AK4]},
            {{"Jet", "mediumbTagSFDown"}, mediumbTagSFDown[AK4]},
            {{"Jet", "tightbTagSFSFUp"}, tightbTagSFUp[AK4]},
            {{"Jet", "tightbTagSFDown"}, tightbTagSFDown[AK4]},
        };

        variables.insert(SFvariations.begin(), SFvariations.end());   
    }

    //Set Branches of output tree
    for(TTree* tree: trees){
        for(std::pair<const std::pair<std::string, std::string>, std::vector<float>&>& variable: variables){
            tree->Branch((variable.first.first + "_" + variable.first.second).c_str(), &variable.second);
        }

        for(std::pair<const std::string, std::vector<bool>&>& variable: bools){
            tree->Branch(("Jet_" + variable.first).c_str(), &variable.second);
        }

        tree->Branch("MET_Px", &metPx);
        tree->Branch("MET_Py", &metPy);
        tree->Branch("HT", &HT);
    }
}

void JetAnalyzer::Analyze(std::vector<CutFlow> &cutflows, const edm::Event* event){
    //Clear variables vector
    for(std::pair<const std::pair<std::string, std::string>, std::vector<float>&>& variable: variables){
        variable.second.clear();
    }

    for(std::pair<const std::string, std::vector<bool>&>& variable: bools){
        variable.second.clear();
    }

    int nSubJets=0;
    HT=0;
    runNumber = isNANO ? *run->Get() : event->eventAuxiliary().id().run(); 

    if(jetCorrector.empty()){
        for(const JetType& type: {AK4, AK8}){
            SetCorrector(type, runNumber);
        }
    }

    //Get Event info is using MINIAOD
    edm::Handle<std::vector<pat::Jet>> jets, fatJets;
    edm::Handle<std::vector<reco::GenJet>> genJets, genfatJets;
    edm::Handle<std::vector<pat::MET>> MET;
    edm::Handle<double> rho;
    edm::Handle<std::vector<reco::GenParticle>> genParts;
    edm::Handle<std::vector<reco::VertexCompositePtrCandidate>> secVtx;

    if(!isNANO){
        event->getByToken(jetTokens[0], jets);
        event->getByToken(jetTokens[1], fatJets);
        event->getByToken(metToken, MET);
        event->getByToken(rhoToken, rho);
        event->getByToken(vertexToken, secVtx);

        if(!isData){
            event->getByLabel(edm::InputTag("slimmedGenJets"), genJets);
            event->getByLabel(edm::InputTag("slimmedGenJetsAK8"), genfatJets);
        }
    }

    //JER smearing
    float smearFac = 1., corrFac = 1.;

    //MET values not correct for JER yet
    metPx = isNANO ? *metPt->Get()*std::cos(*metPhi->Get()) : MET->at(0).uncorPx();
    metPy = isNANO ? *metPt->Get()*std::sin(*metPhi->Get()) : MET->at(0).uncorPy();
    
    float fatJetSize = isNANO ? fatJetPt->GetSize() : fatJets->size();
    float jetSize = isNANO ? jetPt->GetSize() : jets->size();

    //Check if initial jet collection has enough jets
    unsigned int nMin = 0;

    for(CutFlow& cutflow: cutflows){
        if(jetSize < cutflow.nMinJet){
            nMin++;
            cutflow.passed = false;
        }
    }

    if(nMin == cutflows.size()) return;
        
    //Loop over all fat jets
    for(unsigned int i = 0; i < fatJetSize; i++){
        float fatPt = isNANO ? fatJetPt->At(i) : fatJets->at(i).pt();
        float fatEta = isNANO ? fatJetEta->At(i) : fatJets->at(i).eta();
        float fatPhi = isNANO ? fatJetPhi->At(i) : fatJets->at(i).phi();
        float fatMass = isNANO ? fatJetMass->At(i) : fatJets->at(i).mass();

        ROOT::Math::PtEtaPhiMVector LV(fatPt, fatEta, fatPhi, fatMass);
    
        corrFac = isNANO ? CorrectEnergy(LV,  *jetRho->Get(), fatJetArea->At(i), AK8) : CorrectEnergy(LV, *rho, jets->at(i).jetArea(), AK8);

        //Get jet uncertainty
        if(jecUnc[AK8] !=  NULL){
            jecUnc[AK8]->setJetPt(corrFac*LV.Pt());
            jecUnc[AK8]->setJetEta(LV.Eta());
            float unc = jecUnc[AK8]->getUncertainty(isUp);
            corrFac *= isUp ? 1 + unc : 1 - unc;
        }

        //Smear pt if not data
        if(!isData){
            smearFac = isNANO ? SmearEnergy(LV*corrFac, *jetRho->Get(), 0.8, AK8) : SmearEnergy(LV*corrFac, *rho, 0.8, AK8, *genfatJets);
            LV *= smearFac*corrFac;
        }

        else LV *= corrFac;

        if(LV.Pt() > 170. and LV.M() > 40. and abs(LV.Eta()) < etaCut){
            float hScore = fatJets->at(i).bDiscriminator("pfDeepBoostedJetTags:probHbb");
            float topScore = 0.;

            for(const std::string& name : {"pfDeepBoostedJetTags:probTbcq", "pfDeepBoostedJetTags:probTbqq", "pfDeepBoostedJetTags:probTbc", "pfDeepBoostedJetTags:probTbq"}){
                topScore += fatJets->at(i).bDiscriminator(name);
            }

            topVsHiggs[AK8].push_back(hScore/(hScore + topScore));
            QCDVsHiggs[AK8].push_back(fatJets->at(i).bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:HbbvsQCD"));
    
            //Fatjet four momentum components
            E[AK8].push_back(LV.E());
            Px[AK8].push_back(LV.Px());  
            Py[AK8].push_back(LV.Py());  
            Pz[AK8].push_back(LV.Pz());  

            //Nsubjettiness
            Njettiness1[AK8].push_back(isNANO ? fatJetTau1->At(i) : fatJets->at(i).userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau1"));
            Njettiness2[AK8].push_back(isNANO ? fatJetTau2->At(i) : fatJets->at(i).userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau2"));
            Njettiness3[AK8].push_back(isNANO ? fatJetTau3->At(i) : fatJets->at(i).userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau3"));

            if(!isData){
                if(isNANO) isFromh[AK8].push_back(SetGenParticles(LV, i, 5, AK8));
                else{
                    event->getByToken(genParticleToken, genParts);
                    isFromh[AK8].push_back(SetGenParticles(LV, i, 5, AK8, *genParts));
                }
            }

            //Fill in particle flow candidates
            if(!isNANO){
                for(unsigned int j = 0; j < fatJets->at(i).numberOfDaughters(); j++){
                    reco::Candidate const * cand = fatJets->at(i).daughter(j);
                    
                    if(cand->numberOfDaughters() == 0){
                        //Jet particle four momentum components
                        E[PF].push_back(cand->energy());
                        Px[PF].push_back(cand->px());
                        Py[PF].push_back(cand->py());
                        Pz[PF].push_back(cand->pz());

                        //Jet particle vertex
                        Vx[PF].push_back(cand->vx());      
                        Vy[PF].push_back(cand->vy()); 
                        Vz[PF].push_back(cand->vz());
                        Charge[PF].push_back(cand->charge());
            
                        //Fat Jet Index
                        FatJetIdx[PF].push_back(E[AK8].size()-1);
                    }
                            
                    else{
                        for(unsigned int k = 0; k < cand->numberOfDaughters(); k++){
                            reco::Candidate const * cand2 = cand->daughter(k);

                            //Jet particle four momentum components
                            E[PF].push_back(cand2->energy());
                            Px[PF].push_back(cand2->px());
                            Py[PF].push_back(cand2->py());
                            Pz[PF].push_back(cand2->pz());

                            //Jet particle vertex
                            Vx[PF].push_back(cand2->vx());      
                            Vy[PF].push_back(cand2->vy()); 
                            Vz[PF].push_back(cand2->vz());
                            Charge[PF].push_back(cand2->charge());
                
                            //Fat Jet Index
                            FatJetIdx[PF].push_back(E[AK8].size()-1);
                        }
                    }
                }

                for(const reco::VertexCompositePtrCandidate &vtx: *secVtx){
                    ROOT::Math::PtEtaPhiMVector vtxLV(vtx.p4().Pt(), vtx.p4().Eta(), vtx.p4().Phi(), vtx.p4().M());

                    if(ROOT::Math::VectorUtil::DeltaR(LV, vtxLV) < 0.8){
                        //SV four momentum components
                        E[VTX].push_back(vtx.energy());
                        Px[VTX].push_back(vtx.px());
                        Py[VTX].push_back(vtx.py());
                        Pz[VTX].push_back(vtx.pz());

                        //SV vertex
                        Vx[VTX].push_back(vtx.vx());
                        Vy[VTX].push_back(vtx.vy()); 
                        Vz[VTX].push_back(vtx.vz());
            
                        //Fat Jet Index
                        FatJetIdx[VTX].push_back(E[AK8].size()-1);
                    }
                }
            }
        }
    }

    //Loop over all jets
    for(unsigned int i = 0; i < jetSize; i++){
        float pt = isNANO ? jetPt->At(i) : jets->at(i).pt();
        float eta = isNANO ? jetEta->At(i) : jets->at(i).eta();
        float phi = isNANO ? jetPhi->At(i) : jets->at(i).phi();
        float mass = isNANO ? jetMass->At(i) : jets->at(i).mass();

        //Define here already jet, because of smearing of 4-vec
        ROOT::Math::PtEtaPhiMVector LV(pt, eta, phi, mass);

        corrFac = isNANO ? CorrectEnergy(LV,  *jetRho->Get(), jetArea->At(i), AK4) : CorrectEnergy(LV, *rho, jets->at(i).jetArea(), AK4);

        //Get jet uncertainty
        if(jecUnc[AK8] !=  NULL){
            jecUnc[AK4]->setJetPt(corrFac*LV.Pt());
            jecUnc[AK4]->setJetEta(LV.Eta());
            float unc = jecUnc[AK4]->getUncertainty(isUp);
            corrFac *= isUp ? 1 + unc : 1 - unc;
        }

        //Smear pt if not data
        if(!isData){
            smearFac = isNANO ? SmearEnergy(LV*corrFac,  *jetRho->Get(), jetArea->At(i), AK4) : SmearEnergy(LV, *rho, 0.4, AK4, *genJets);

            LV*=smearFac*corrFac;
        }

        else LV*=corrFac;

        //Correct met
        metPx+= LV.Px()*(1-smearFac*corrFac);
        metPy+= LV.Py()*(1-smearFac*corrFac);

        LV *= smearFac*corrFac;

        if(LV.Pt() > ptCut and abs(LV.Eta()) < etaCut){
            HT+=LV.Pt();

            //Fatjet four momentum components
            E[AK4].push_back(LV.E());
            Px[AK4].push_back(LV.Px());
            Py[AK4].push_back(LV.Py());
            Pz[AK4].push_back(LV.Pz());

            //Check for btag
            float DeepBValue = 0;

            if(isNANO){ 
                DeepBValue = jetDeepBValue->At(i);
            } 

            else{
                for(std::string disc: {"pfDeepFlavourJetTags:probb", "pfDeepFlavourJetTags:probbb","pfDeepFlavourJetTags:problepb"}){
                    DeepBValue +=  jets->at(i).bDiscriminator(disc);
                }
            }

            isLooseB.push_back(bTagCuts[AK4][era][0] < DeepBValue);
            isMediumB.push_back(bTagCuts[AK4][era][1] < DeepBValue);
            isTightB.push_back(bTagCuts[AK4][era][2] < DeepBValue);

            if(!isData){
                //btag SF
                loosebTagSF[AK4].push_back(looseReader[AK4].eval_auto_bounds("central", BTagEntry::FLAV_B, abs(LV.Eta()), LV.Pt()));
                mediumbTagSF[AK4].push_back(mediumReader[AK4].eval_auto_bounds("central", BTagEntry::FLAV_B, abs(LV.Eta()), LV.Pt()));
                tightbTagSF[AK4].push_back(tightReader[AK4].eval_auto_bounds("central", BTagEntry::FLAV_B, abs(LV.Eta()), LV.Pt()));

                if(!isSyst){
                    loosebTagSFUp[AK4].push_back(looseReader[AK4].eval_auto_bounds("up", BTagEntry::FLAV_B, abs(LV.Eta()), LV.Pt()));
                    loosebTagSFDown[AK4].push_back(looseReader[AK4].eval_auto_bounds("down", BTagEntry::FLAV_B, abs(LV.Eta()), LV.Pt()));
                    mediumbTagSFUp[AK4].push_back(mediumReader[AK4].eval_auto_bounds("up", BTagEntry::FLAV_B, abs(LV.Eta()), LV.Pt()));
                    mediumbTagSFDown[AK4].push_back(mediumReader[AK4].eval_auto_bounds("down", BTagEntry::FLAV_B, abs(LV.Eta()), LV.Pt()));
                    tightbTagSFUp[AK4].push_back(tightReader[AK4].eval_auto_bounds("up", BTagEntry::FLAV_B, abs(LV.Eta()), LV.Pt()));
                    tightbTagSFDown[AK4].push_back(tightReader[AK4].eval_auto_bounds("down", BTagEntry::FLAV_B, abs(LV.Eta()), LV.Pt()));
                }


                if(isNANO) isFromh[AK4].push_back(SetGenParticles(LV, i, 6, AK4));
                else{
                    event->getByToken(genParticleToken, genParts);
                    isFromh[AK4].push_back(SetGenParticles(LV, i, 6, AK4, *genParts));
                }
            }

            //Check overlap with AK4 valid jets
            for(unsigned int j = 0; j < E[AK8].size(); j++){
                ROOT::Math::PtEtaPhiMVector FatLV(E[AK8][j], Px[AK8][j], Py[AK8][j], Pz[AK8][j]);

                if(ROOT::Math::VectorUtil::DeltaR(LV, FatLV) < 1.2){
                    FatJetIdx[AK4].push_back(j);
                    nSubJets++;
                }

                else{
                    FatJetIdx[AK4].push_back(-1.);
                }
            }
        } 
    }

    for(CutFlow& cutflow: cutflows){
        //Check if one combination of jet and fatjet number is fullfilled
        if(E[AK4].size() - nSubJets >= cutflow.nMinJet and E[AK8].size() == cutflow.nMinFatjet){
            if(cutflow.passed){
                std::string cutName("N^{AK4}_{jet} >= " + std::to_string(cutflow.nMinJet) + " && N^{AK8}_{jet} == " + std::to_string(cutflow.nMinFatjet));

                cutflow.hist->Fill(cutName.c_str(), cutflow.weight);
            }
        }

        else{
            cutflow.passed = false;
        }
    }
}

void JetAnalyzer::EndJob(TFile* file){
    for(const JetType& type: {AK4, AK8}){
        delete jetCorrector[type];
    }
}
