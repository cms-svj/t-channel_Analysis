#ifndef BASELINE_H
#define BASELINE_H

class Baseline
{
private:
    std::string myVarSuffix_;

    void baseline(NTupleReader& tr)
    {
        const auto& runtype                  = tr.getVar<std::string>("runtype");     
        const auto& filetag                  = tr.getVar<std::string>("filetag");
        const auto& runYear                  = tr.getVar<std::string>("runYear");
        const auto& RunNum                   = tr.getVar<unsigned int>("RunNum");
        const auto& blind                    = tr.getVar<bool>("blind");
        const auto& TriggerNames             = tr.getVec<std::string>("TriggerNames");
        const auto& TriggerPass              = tr.getVec<int>("TriggerPass");
        const auto& NNonIsoMuons             = tr.getVar<int>("NNonIsoMuons"+myVarSuffix_);
        const auto& NGoodLeptons             = tr.getVar<int>("NGoodLeptons"+myVarSuffix_);
        const auto& GoodLeptonsCharge        = tr.getVec<int>("GoodLeptonsCharge");
        const auto& NGoodMuons               = tr.getVar<int>("NGoodMuons"+myVarSuffix_);
        const auto& NGoodPlusMuons           = tr.getVar<int>("NGoodPlusMuons"+myVarSuffix_);
        const auto& NGoodMinusMuons          = tr.getVar<int>("NGoodMinusMuons"+myVarSuffix_);
        const auto& NGoodElectrons           = tr.getVar<int>("NGoodElectrons"+myVarSuffix_);
        const auto& NGoodPlusElectrons       = tr.getVar<int>("NGoodPlusElectrons"+myVarSuffix_);
        const auto& NGoodMinusElectrons      = tr.getVar<int>("NGoodMinusElectrons"+myVarSuffix_);
        const auto& NGoodLeptons_pt20        = tr.getVar<int>("NGoodLeptons_pt20"+myVarSuffix_);
        const auto& GoodLeptonsCharge_pt20   = tr.getVec<int>("GoodLeptonsCharge_pt20");
        const auto& HT_trigger               = tr.getVar<double>("HT_trigger"+myVarSuffix_);
        const auto& HT_trigger_pt30          = tr.getVar<double>("HT_trigger_pt30"+myVarSuffix_);
        const auto& HT_NonIsoMuon_pt30       = tr.getVar<double>("HT_NonIsoMuon_pt30"+myVarSuffix_);
        const auto& onZ                      = tr.getVar<bool>("onZ"+myVarSuffix_); 
        const auto& JetID                    = tr.getVar<bool>("JetID"+myVarSuffix_); 
        const auto& NGoodJets_pt45           = tr.getVar<int>("NGoodJets_pt45"+myVarSuffix_);
        const auto& NGoodBJets_pt45          = tr.getVar<int>("NGoodBJets_pt45"+myVarSuffix_);
        const auto& NGoodBJets_pt30          = tr.getVar<int>("NGoodBJets_pt30"+myVarSuffix_); 
        const auto& NGoodJets_pt30           = tr.getVar<int>("NGoodJets_pt30"+myVarSuffix_); 
        const auto& NNonIsoMuonJets_pt30     = tr.getVar<int>("NNonIsoMuonJets_pt30"+myVarSuffix_);         
        const auto& NGoodPhotons             = tr.getVar<int>("NGoodPhotons"+myVarSuffix_);
        const auto& Mbl                      = tr.getVar<double>("Mbl"+myVarSuffix_);
        const auto& passHEMVeto              = tr.getVar<bool>("passHEMVeto"+myVarSuffix_);
 
        //for(auto n : TriggerNames){std::cout<<n<<std::endl;}

        // ------------------------------
        // -- Data dependent stuff
        // ------------------------------
        bool passTriggerPhoton = PassTriggerPhoton(TriggerNames, TriggerPass);
        bool passTriggerAllHad = false, passTriggerMuon = false, passTriggerElectron = false, passTriggerNonIsoMuon = false;
        if (runYear == "2016")
        {
            passTriggerAllHad = PassTriggerAllHad2016(TriggerNames, TriggerPass);
            passTriggerMuon = PassTriggerMuon2016(TriggerNames, TriggerPass);
            passTriggerElectron = PassTriggerElectron2016(TriggerNames, TriggerPass);
            passTriggerNonIsoMuon = PassTriggerNonIsoMuon2016(TriggerNames, TriggerPass);
        }
        else if (runYear == "2017")
        {
            passTriggerAllHad = PassTriggerAllHad2017(TriggerNames, TriggerPass);
            passTriggerMuon = PassTriggerMuon2017(TriggerNames, TriggerPass);
            passTriggerElectron = PassTriggerElectron2017(TriggerNames, TriggerPass);
            passTriggerNonIsoMuon = PassTriggerNonIsoMuon2017(TriggerNames, TriggerPass);            
        }
        else if (runYear == "2018pre" || runYear == "2018post")
        {
            passTriggerAllHad = true;
            passTriggerMuon = PassTriggerMuon2018(TriggerNames, TriggerPass);
            passTriggerElectron = PassTriggerElectron2018(TriggerNames, TriggerPass);
            passTriggerNonIsoMuon = PassTriggerNonIsoMuon2018(TriggerNames, TriggerPass);
        }

        bool passTrigger   = true;
        bool passTriggerMC = true;
        bool passNonIsoTrigger = true;
        bool passNonIsoTriggerMC = true;
        bool passBlindHad_Good = true;
        bool passBlindLep_Good = true;        
        bool correct2018Split = true;
        if (runtype == "Data")
        {            
            if      (runYear == "2018pre"  && RunNum >= 319077) correct2018Split = false;
            else if (runYear == "2018post" && RunNum <  319077) correct2018Split = false;

            // Pass the right trigger
            if (filetag.find("Data_JetHT") != std::string::npos && !passTriggerAllHad) passTrigger = false;
            if (filetag.find("Data_SingleMuon") != std::string::npos && !passTriggerMuon) passTrigger = false;
            if (filetag.find("Data_SingleElectron") != std::string::npos && !passTriggerElectron) passTrigger = false;
            if (filetag.find("Data_SinglePhoton") != std::string::npos && !passTriggerPhoton) passTrigger = false;
            if (filetag.find("Data_SingleMuon") != std::string::npos && !passTriggerNonIsoMuon) passNonIsoTrigger = false;

            // Blinding data 
            if (NGoodJets_pt30 >= 9 && blind) passBlindHad_Good = false;
            if (NGoodJets_pt30 >= 9 && blind) passBlindLep_Good = false;
        }
        
        // ------------------------
        // -- MC dependent stuff - 
        // -----------------------
        bool passMadHT = true;
        if(runtype == "MC")
        {
            const auto& madHT  = tr.getVar<double>("madHT");

            // Exclude events with MadGraph HT > 100 from the DY & WJets inclusive samples
            if(filetag.find("DYJetsToLL_M-50_Incl") != std::string::npos && madHT > 100) passMadHT = false;
            if(filetag.find("WJetsToLNu_Incl") != std::string::npos      && madHT > 100) passMadHT = false;

            // Stitch TTbar samples together
            // remove HT overlap
            if((filetag.find("TTJets_Incl") != std::string::npos || 
                filetag.find("TTJets_SingleLeptFromT") != std::string::npos || 
                filetag.find("TTJets_SingleLeptFromTbar") != std::string::npos || 
                filetag.find("TTJets_DiLept") != std::string::npos) && madHT > 600) 
            {
                passMadHT = false;
            }
            // also remove lepton overlap from the inclusive sample
            const auto& GenElectrons        = tr.getVec<TLorentzVector>("GenElectrons");
            const auto& GenMuons            = tr.getVec<TLorentzVector>("GenMuons");
            const auto& GenTaus             = tr.getVec<TLorentzVector>("GenTaus");
            int NGenLeptons = GenElectrons.size() + GenMuons.size() + GenTaus.size();
            if (filetag.find("TTJets_Incl") != std::string::npos && NGenLeptons > 0) passMadHT = false;
            
            // MC modeling of the trigger
            if( !passTriggerMuon && !passTriggerElectron ) passTriggerMC = false;
            if( !passTriggerNonIsoMuon ) passNonIsoTriggerMC = false;
        }

        //------------------------
        // -- MET Filter dependent stuff
        // -----------------------
        const auto& globalSuperTightHalo2016Filter      = static_cast<bool>( tr.getVar<int>("globalSuperTightHalo2016Filter") );
        const auto& PrimaryVertexFilter                 = static_cast<bool>( tr.getVar<int>("PrimaryVertexFilter") );
        const auto& BadPFMuonFilter                     = static_cast<bool>( tr.getVar<bool>("BadPFMuonFilter") );
        const auto& EcalDeadCellTriggerPrimitiveFilter  = static_cast<bool>( tr.getVar<int>("EcalDeadCellTriggerPrimitiveFilter") );
        const auto& HBHEIsoNoiseFilter                  = static_cast<bool>( tr.getVar<int>("HBHEIsoNoiseFilter") );
        const auto& HBHENoiseFilter                     = static_cast<bool>( tr.getVar<int>("HBHENoiseFilter") );
        bool passMETFilters  = globalSuperTightHalo2016Filter && PrimaryVertexFilter && BadPFMuonFilter && EcalDeadCellTriggerPrimitiveFilter && HBHEIsoNoiseFilter && HBHENoiseFilter;
                
        // -----------------------------------
        // -- Define 1-lepton proto-baseline
        // -----------------------------------
        bool passBaselineGoodOffline1l = JetID                 &&
                                         passHEMVeto           &&
                                         correct2018Split      &&
                                         passMETFilters        &&
                                         HT_trigger_pt30 > 300 &&
                                         passMadHT             &&
                                         NGoodBJets_pt30 >= 1  &&
                                         (50 < Mbl && Mbl < 250);

        // -------------------------------
        // -- Define 0 Lepton Baseline
        // -------------------------------
        bool passBaseline0l_Good = JetID                &&
                                   passMETFilters       &&
                                   passMadHT            &&
                                   //passTrigger          &&
                                   //passTriggerMC        &&
                                   //passBlindHad_Good    &&
                                   NGoodLeptons == 0    &&
                                   HT_trigger > 500     && 
                                   NGoodBJets_pt45 >= 2 && 
                                   NGoodJets_pt45 >= 6;   
 
        // -------------------------------
        // -- Define 1 Lepton Baseline
        // -------------------------------
        bool passBaseline1mu_Good = passBaselineGoodOffline1l &&
                                    passTrigger               &&
                                    passTriggerMC             &&
                                    (runtype != "Data" || filetag.find("Data_SingleMuon") != std::string::npos) &&
                                    passBlindLep_Good         &&
                                    NGoodMuons == 1           && 
                                    NGoodElectrons == 0       &&
                                    NGoodJets_pt30 >= 7;

        bool passBaseline1el_Good = passBaselineGoodOffline1l &&
                                    passTrigger               &&
                                    passTriggerMC             &&
                                    (runtype != "Data" || filetag.find("Data_SingleElectron") != std::string::npos) &&
                                    passBlindLep_Good         &&
                                    NGoodElectrons == 1       &&
                                    NGoodMuons == 0           &&
                                    NGoodJets_pt30 >= 7;

        bool passBaseline1l_Good = passBaseline1mu_Good || passBaseline1el_Good;

        bool passBaseline1l_NonIsoMuon = HT_NonIsoMuon_pt30 > 300 &&
                                         correct2018Split         &&
                                         passHEMVeto              &&
                                         passMETFilters           &&
                                         passMadHT                &&
                                         passNonIsoTrigger        &&
                                         passNonIsoTriggerMC      &&
                                         (runtype != "Data" || filetag.find("Data_SingleMuon") != std::string::npos) &&
                                         NNonIsoMuons == 1        &&
                                         NGoodMuons == 0          &&
                                         NGoodElectrons == 0      &&
                                         JetID                    &&
                                         NNonIsoMuonJets_pt30 >= 7;

        // ----------------------------------
        // -- Define 2 Lepton onZ Baseline
        // ----------------------------------        
        bool passBaseline2lonZ_Good = JetID              &&
                                      passMETFilters     &&
                                      passMadHT          &&
                                      passTrigger        &&
                                      passTriggerMC      &&
                                      NGoodLeptons == 2 ? GoodLeptonsCharge[0]!=GoodLeptonsCharge[1] : false &&
                                      onZ                &&
                                      (runtype != "Data" || (NGoodMuons == 2 && filetag.find("Data_SingleMuon") != std::string::npos ) 
                                                         || (NGoodElectrons == 2 && filetag.find("Data_SingleElectron") != std::string::npos) ) &&
                                      NGoodJets_pt30 >= 7; 
        
        // -----------------------------------
        // -- Define 2 Lepton offZ Baseline
        // -----------------------------------
        bool passBaseline2l_Good = JetID                &&
                                   passMETFilters       &&
                                   passMadHT            &&
                                   passTrigger          &&
                                   passTriggerMC        &&
                                   passBlindLep_Good    &&
                                   NGoodLeptons == 2 ? GoodLeptonsCharge[0]!=GoodLeptonsCharge[1] : false && 
                                   !onZ                 &&
                                   (runtype != "Data"  || (NGoodMuons >= 1 && filetag.find("Data_SingleMuon") != std::string::npos ) 
                                                       || (NGoodElectrons == 2 && filetag.find("Data_SingleElectron") != std::string::npos) ) &&
                                   NGoodBJets_pt30 >= 1 &&
                                   NGoodJets_pt30 >= 4;

        // -----------------------------------
        // -- Define 2 Lepton pt20 Baseline
        // -----------------------------------

        bool passBaseline2l_pt20 =    JetID              &&
                                      passMETFilters     &&
                                      passMadHT          &&
            //passTrigger        &&
            //                          passTriggerMC      &&
                                      NGoodLeptons_pt20 == 2 ? GoodLeptonsCharge_pt20[0]!=GoodLeptonsCharge_pt20[1] : false &&
                                      (runtype != "Data" || (NGoodMuons == 2 && filetag.find("Data_SingleMuon") != std::string::npos ) 
                                                         || (NGoodElectrons == 2 && filetag.find("Data_SingleElectron") != std::string::npos) ) &&
                                      NGoodJets_pt30 >= 6 &&
                                      NGoodBJets_pt30 >= 1; 

        // -----------------------------------
        // -- Define 2 Lepton pt30 Baseline
        // -----------------------------------

        bool passBaseline2l_pt30 =    JetID              &&
                                      passMETFilters     &&
                                      passMadHT          &&
            //passTrigger        &&
            //                          passTriggerMC      &&
                                      NGoodLeptons == 2 ? GoodLeptonsCharge[0]!=GoodLeptonsCharge[1] : false &&
                                      (runtype != "Data" || (NGoodMuons == 2 && filetag.find("Data_SingleMuon") != std::string::npos ) 
                                                         || (NGoodElectrons == 2 && filetag.find("Data_SingleElectron") != std::string::npos) ) &&
                                      NGoodJets_pt30  >= 4 &&
                                      NGoodBJets_pt30 >= 1; 


        // -----------------------------------
        // -- Define 1 e and 1 m Baseline
        // -----------------------------------        
        bool passBaseline1e1m_Good = JetID                 &&
                                     passMETFilters        &&
                                     HT_trigger_pt30 > 300 &&
                                     passMadHT             &&
                                     passTriggerMuon       &&
                                     NGoodMuons == 1       &&
                                     NGoodElectrons == 1   &&
                                     (runtype != "Data" || filetag.find("Data_SingleMuon") != std::string::npos) &&
                                     (NGoodPlusMuons + NGoodPlusElectrons == 1)          &&
                                     (NGoodMinusMuons + NGoodMinusElectrons == 1)        &&
                                     NGoodJets_pt30 <= 6   && 
                                     NGoodBJets_pt30 >= 2;
        
        // -----------------------------------
        // -- Define 1 Photon Baseline
        // -----------------------------------
        bool passBaseline1photon_Good = passMadHT           &&
                                        passMETFilters      &&
                                        passTrigger         &&
                                        passTriggerMC       &&
                                        NGoodPhotons == 1   &&
                                        NGoodLeptons == 0   && 
                                        NGoodJets_pt30 >= 7; 
        
        tr.registerDerivedVar<bool>("passBaselineGoodOffline1l"+myVarSuffix_,   passBaselineGoodOffline1l);
        tr.registerDerivedVar<bool>("passBaseline0l_Good"+myVarSuffix_,         passBaseline0l_Good); 
        tr.registerDerivedVar<bool>("passBaseline1l_Good"+myVarSuffix_,         passBaseline1l_Good);
        tr.registerDerivedVar<bool>("passBaseline1l_NonIsoMuon"+myVarSuffix_,   passBaseline1l_NonIsoMuon);
        tr.registerDerivedVar<bool>("passBaseline2lonZ_Good"+myVarSuffix_,      passBaseline2lonZ_Good);
        tr.registerDerivedVar<bool>("passBaseline2l_Good"+myVarSuffix_,         passBaseline2l_Good);
        tr.registerDerivedVar<bool>("passBaseline2l_pt20"+myVarSuffix_,         passBaseline2l_pt20);
        tr.registerDerivedVar<bool>("passBaseline2l_pt30"+myVarSuffix_,         passBaseline2l_pt30);
        tr.registerDerivedVar<bool>("passBaseline1photon_Good"+myVarSuffix_,    passBaseline1photon_Good);
        tr.registerDerivedVar<bool>("passBaseline1e1m_Good"+myVarSuffix_,       passBaseline1e1m_Good);
        tr.registerDerivedVar<bool>("passBlindHad_Good"+myVarSuffix_,           passBlindHad_Good);
        tr.registerDerivedVar<bool>("passBlindLep_Good"+myVarSuffix_,           passBlindLep_Good);
        tr.registerDerivedVar<bool>("passTriggerAllHad"+myVarSuffix_,           passTriggerAllHad);
        tr.registerDerivedVar<bool>("passTriggerMuon"+myVarSuffix_,             passTriggerMuon);
        tr.registerDerivedVar<bool>("passTriggerElectron"+myVarSuffix_,         passTriggerElectron);
        tr.registerDerivedVar<bool>("passTrigger"+myVarSuffix_,                 passTrigger);
        tr.registerDerivedVar<bool>("passNonIsoTrigger"+myVarSuffix_,           passNonIsoTrigger);
        tr.registerDerivedVar<bool>("passTriggerMC"+myVarSuffix_,               passTriggerMC);
        tr.registerDerivedVar<bool>("passNonIsoTriggerMC"+myVarSuffix_,         passNonIsoTriggerMC);
        tr.registerDerivedVar<bool>("passMadHT"+myVarSuffix_,                   passMadHT);
        tr.registerDerivedVar<bool>("passMETFilters"+myVarSuffix_,              passMETFilters);
        tr.registerDerivedVar<bool>("correct2018Split"+myVarSuffix_,            correct2018Split);        
    }

    bool PassTriggerGeneral(std::vector<std::string>& mytriggers, const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        bool passTrigger = false;
        for(unsigned int i=0; i<TriggerNames.size(); ++i)
        {
            if(TriggerPass.at(i) != 1)
                continue;
            std::string trigname = TriggerNames.at(i);
            if( std::any_of(mytriggers.begin(), mytriggers.end(), [&] (std::string s) { return trigname.find(s)!=std::string::npos; }) )
            {
                passTrigger = true;
                break;
            }
        }
        return passTrigger;
    }

    bool PassTriggerAllHad2016(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {
            //"HLT_PFHT900",
            "HLT_PFHT450_SixJet40_BTagCSV",
            "HLT_PFHT400_SixJet30_DoubleBTagCSV",            
        };
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerAllHad2017(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {
            //"HLT_PFHT1050",
            "HLT_PFHT380_SixPFJet32_DoublePFBTagCSV",
            "HLT_PFHT430_SixPFJet40_PFBTagCSV",
        };
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerMuon2016(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {"HLT_IsoMu24","HLT_IsoTkMu24","HLT_Mu50","HLT_TkMu50"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerMuon2017(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        // IsoMu24 ended up prescaled for 3.6\fb, so also include the IsoMu27 here
        std::vector<std::string> mytriggers = {"HLT_IsoMu24","HLT_IsoMu27","HLT_Mu50"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerMuon2018(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        // Preliminary list for 2018 
        std::vector<std::string> mytriggers = {"HLT_IsoMu24","HLT_IsoMu27","HLT_Mu50"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerElectron2016(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {"HLT_Ele27_WPTight_Gsf","HLT_Ele115_CaloIdVT_GsfTrkIdT", "HLT_Photon175"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerElectron2017(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        // Ele35 fine for whole year, gets complicated for the other triggers
        std::vector<std::string> mytriggers = {"HLT_Ele35_WPTight","HLT_Ele115_CaloIdVT_GsfTrkIdT", "HLT_Photon200"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerElectron2018(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        // Preliminary list for 2018 
        std::vector<std::string> mytriggers = {"HLT_Ele35_WPTight","HLT_Ele115_CaloIdVT_GsfTrkIdT", "HLT_Photon200"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerPhoton(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {"HLT_Photon175_v"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }
    
    bool PassTriggerNonIsoMuon2016(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        std::vector<std::string> mytriggers = {"HLT_Mu50", "HLT_TkMu50"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }
    
    bool PassTriggerNonIsoMuon2017(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        //No track based muon 50 trigger in 2017 in our tuples
        std::vector<std::string> mytriggers = {"HLT_Mu50"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

    bool PassTriggerNonIsoMuon2018(const std::vector<std::string>& TriggerNames, const std::vector<int>& TriggerPass)
    {
        //No track based muon 50 trigger in 2018 in our tuples
        std::vector<std::string> mytriggers = {"HLT_Mu50"};
        return PassTriggerGeneral(mytriggers,TriggerNames,TriggerPass);
    }

public:
    Baseline(std::string myVarSuffix = "")
        : myVarSuffix_(myVarSuffix)
    {
        std::cout<<"Setting up Baseline"<<std::endl;
    }
    
    void operator()(NTupleReader& tr)
    {
        baseline(tr);
    }
};

#endif
