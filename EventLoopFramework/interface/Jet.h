#ifndef JET_H
#define JET_H

class Jet
{
private:
    std::string myVarSuffix_;

    void setJetVar(bool pass, std::vector<bool>& boolVec, int& n)
    {
        if( pass )
        {
            boolVec.push_back(true);                
            n++;
        }
        else
        {
            boolVec.push_back(false);
        }            
    }

    void jet(NTupleReader& tr)
    {
        const auto& Jets          = tr.getVec<TLorentzVector>(("Jets"+myVarSuffix_));
        const auto& etaCut        = tr.getVar<double>("etaCut");
        const auto& Muons         = tr.getVec<TLorentzVector>("Muons");
        const auto& GoodMuons     = tr.getVec<bool>("GoodMuons"+myVarSuffix_);
        const auto& NMuons        = tr.getVar<int>("NGoodMuons"+myVarSuffix_);
        const auto& Electrons     = tr.getVec<TLorentzVector>("Electrons");
        const auto& GoodElectrons = tr.getVec<bool>("GoodElectrons"+myVarSuffix_);
        const auto& NElectrons    = tr.getVar<int>("NGoodElectrons"+myVarSuffix_);
        const auto& NonIsoMuons   = tr.getVec<bool>("NonIsoMuons"+myVarSuffix_);

        //Adding code to create a vector of GoodJets -> defined as the jet collection that eliminates the closest jet to any good lepton (muon or electron) 
        //if that delta R is less than 0.4 and the pT of the jet and lepton is approximately the same
        auto tempGoodJets       = std::make_unique<std::vector<bool>>(Jets.size(), true);
        auto tempNonIsoMuonJets = std::make_unique<std::vector<bool>>(Jets.size(), true);

        if( NMuons > 0 ) 
        {
            for(unsigned int imu = 0; imu < Muons.size(); ++imu)
            {            
                TLorentzVector myMuon = Muons.at(imu);
                double         tempDeltaR = 10.0;
                int            tempJetIt  = -1;                
                for(unsigned int myJetIt = 0; myJetIt < Jets.size(); ++myJetIt ) 
                {                    
                    TLorentzVector myJet = Jets.at(myJetIt);
                    //First check pT matching between Jet and Muon
                    if( std::fabs( myJet.Pt() - myMuon.Pt() ) / myMuon.Pt() > 1.0 ) continue; 
                    double jetDeltaR = myMuon.DeltaR(myJet);                    
                    if( jetDeltaR < tempDeltaR ) 
                    {
                        tempDeltaR = jetDeltaR;
                        tempJetIt  = myJetIt;
                    }                    
                }//END of looping through jets

                if( tempDeltaR < 0.4 && tempJetIt != -1) 
                {
                    if(GoodMuons[imu])                     tempGoodJets->at(tempJetIt)       = false;
                    if(GoodMuons[imu] || NonIsoMuons[imu]) tempNonIsoMuonJets->at(tempJetIt) = false;
                }
            }//END of looping through muons
        }//END of NMuons if statement

        if( NElectrons > 0 ) 
        {
            for(unsigned int iel = 0; iel < Electrons.size(); ++iel)
            {
                if(!GoodElectrons[iel]) continue;
                TLorentzVector myElectron = Electrons.at(iel);
                double         tempDeltaR = 10.0;
                int            tempJetIt  = -1;                
                for(unsigned int myJetIt = 0; myJetIt < Jets.size(); ++myJetIt ) 
                {               
                    TLorentzVector myJet = Jets.at(myJetIt);
                    //Check pT matching between Jet and Electron
                    if( std::fabs( myJet.Pt() - myElectron.Pt() ) / myElectron.Pt() > 1.0 ) continue;
                    double jetDeltaR = myElectron.DeltaR(myJet);                    
                    if( jetDeltaR < tempDeltaR ) 
                    {
                        tempDeltaR = jetDeltaR;
                        tempJetIt  = myJetIt;
                    }
                }//END of looping through jets                
                if( tempDeltaR < 0.4 && tempJetIt != -1) 
                {
                    tempGoodJets->at(tempJetIt)       = false;
                    tempNonIsoMuonJets->at(tempJetIt) = false;
                }
            }//END of looping through electrons
        }//END of NElectrons if statement

        auto& jets_pt20_ = tr.createDerivedVec<bool>("Jets_pt20"+myVarSuffix_);
        auto& jets_pt30_ = tr.createDerivedVec<bool>("Jets_pt30"+myVarSuffix_);
        auto& jets_pt40_ = tr.createDerivedVec<bool>("Jets_pt40"+myVarSuffix_);
        auto& jets_pt45_ = tr.createDerivedVec<bool>("Jets_pt45"+myVarSuffix_);
        int NJets_pt20 = 0, NJets_pt30 = 0, NJets_pt40 = 0, NJets_pt45 = 0;

        auto& goodjets_      = tr.createDerivedVec<bool>("GoodJets"+myVarSuffix_);
        auto& goodjets_pt20_ = tr.createDerivedVec<bool>("GoodJets_pt20"+myVarSuffix_);
        auto& goodjets_pt30_ = tr.createDerivedVec<bool>("GoodJets_pt30"+myVarSuffix_);
        auto& goodjets_pt40_ = tr.createDerivedVec<bool>("GoodJets_pt40"+myVarSuffix_);
        auto& goodjets_pt45_ = tr.createDerivedVec<bool>("GoodJets_pt45"+myVarSuffix_);
        int NGoodJets = 0, NGoodJets_pt20 = 0, NGoodJets_pt30 = 0, NGoodJets_pt40 = 0, NGoodJets_pt45 = 0;

        auto& nonIsoMuonjets_      = tr.createDerivedVec<bool>("NonIsoMuonJets"+myVarSuffix_);
        auto& nonIsoMuonjets_pt30_ = tr.createDerivedVec<bool>("NonIsoMuonJets_pt30"+myVarSuffix_);
        auto& nonIsoMuonjets_pt40_ = tr.createDerivedVec<bool>("NonIsoMuonJets_pt40"+myVarSuffix_);
        auto& nonIsoMuonjets_pt45_ = tr.createDerivedVec<bool>("NonIsoMuonJets_pt45"+myVarSuffix_);
        int NNonIsoMuonJets = 0, NNonIsoMuonJets_pt30 = 0, NNonIsoMuonJets_pt40 = 0, NNonIsoMuonJets_pt45 = 0;

        for(unsigned int i = 0; i < Jets.size(); ++i ) 
        {               
            TLorentzVector lv = Jets.at(i);
           
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 20, jets_pt20_, NJets_pt20); 
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 30, jets_pt30_, NJets_pt30);
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 40, jets_pt40_, NJets_pt40);
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 45, jets_pt45_, NJets_pt45);

            setJetVar( abs(lv.Eta()) < etaCut &&             tempGoodJets->at(i), goodjets_,      NGoodJets     );
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 20 && goodjets_.at(i), goodjets_pt20_, NGoodJets_pt20);
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 30 && goodjets_.at(i), goodjets_pt30_, NGoodJets_pt30);
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 40 && goodjets_.at(i), goodjets_pt40_, NGoodJets_pt40);
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 45 && goodjets_.at(i), goodjets_pt45_, NGoodJets_pt45);

            setJetVar( abs(lv.Eta()) < etaCut &&             tempNonIsoMuonJets->at(i), nonIsoMuonjets_,      NNonIsoMuonJets     );
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 30 && nonIsoMuonjets_.at(i), nonIsoMuonjets_pt30_, NNonIsoMuonJets_pt30);
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 40 && nonIsoMuonjets_.at(i), nonIsoMuonjets_pt40_, NNonIsoMuonJets_pt40);
            setJetVar( abs(lv.Eta()) < etaCut && lv.Pt() > 45 && nonIsoMuonjets_.at(i), nonIsoMuonjets_pt45_, NNonIsoMuonJets_pt45);
        }

        tr.registerDerivedVar("NJets_pt20"+myVarSuffix_,  NJets_pt20);
        tr.registerDerivedVar("NJets_pt30"+myVarSuffix_,  NJets_pt30);
        tr.registerDerivedVar("NJets_pt40"+myVarSuffix_,  NJets_pt40);
        tr.registerDerivedVar("NJets_pt45"+myVarSuffix_,  NJets_pt45);

        tr.registerDerivedVar("NGoodJets"     +myVarSuffix_, NGoodJets);       
        tr.registerDerivedVar("NGoodJets_pt20"+myVarSuffix_, NGoodJets_pt20); 
        tr.registerDerivedVar("NGoodJets_pt30"+myVarSuffix_, NGoodJets_pt30);
        tr.registerDerivedVar("NGoodJets_pt30_double"+myVarSuffix_, static_cast<double>(NGoodJets_pt30));
        tr.registerDerivedVar("NGoodJets_pt40"+myVarSuffix_, NGoodJets_pt40);
        tr.registerDerivedVar("NGoodJets_pt45"+myVarSuffix_, NGoodJets_pt45);

        tr.registerDerivedVar("NNonIsoMuonJets"     +myVarSuffix_, NNonIsoMuonJets);        
        tr.registerDerivedVar("NNonIsoMuonJets_pt30"+myVarSuffix_, NNonIsoMuonJets_pt30);
        tr.registerDerivedVar("NNonIsoMuonJets_pt40"+myVarSuffix_, NNonIsoMuonJets_pt40);
        tr.registerDerivedVar("NNonIsoMuonJets_pt45"+myVarSuffix_, NNonIsoMuonJets_pt45);

        int NGoodJets_pt30_inclusive = NGoodJets_pt30 > 12 ? 12 : NGoodJets_pt30;
        tr.registerDerivedVar("NGoodJets_pt30_inclusive"+myVarSuffix_, NGoodJets_pt30_inclusive);
        tr.registerDerivedVar("NGoodJets_pt30_inclusive_shift"+myVarSuffix_, NGoodJets_pt30_inclusive - 7);

        int NNonIsoMuonJets_pt30_inclusive = NNonIsoMuonJets_pt30 > 12 ? 12 : NNonIsoMuonJets_pt30;
        tr.registerDerivedVar("NNonIsoMuonJets_pt30_inclusive"+myVarSuffix_, NNonIsoMuonJets_pt30_inclusive);
        tr.registerDerivedVar("NNonIsoMuonJets_pt30_inclusive_shift"+myVarSuffix_, NNonIsoMuonJets_pt30_inclusive - 7);
        tr.createDerivedVec<bool>("AllJets"+myVarSuffix_, Jets.size(), true);
    }

public:
    Jet(std::string myVarSuffix = "") 
        : myVarSuffix_(myVarSuffix)
    {
        std::cout<<"Setting up Jet"<<std::endl;
    }

    void operator()(NTupleReader& tr)
    {
        jet(tr);
    }
};

#endif
