#ifndef FATJETCOMBINE_H
#define FATJETCOMBINE_H

#include "t-channel_Analysis/EventLoopFramework/interface/Utility.h"
#include "TopTagger/TopTagger/interface/TopTaggerUtilities.h"
#include "TopTagger/TopTagger/interface/lester_mt2_bisect.h"

class FatJetCombine
{
private:
    std::string myVarSuffix_;

    void fatjetcombine(NTupleReader& tr)
    {
        const auto& JetsAK8               = tr.getVec<TLorentzVector>("JetsAK8");
        const auto& Tau1                  = tr.getVec<double>("JetsAK8_NsubjettinessTau1");
        const auto& Tau2                  = tr.getVec<double>("JetsAK8_NsubjettinessTau2");
        const auto& Tau3                  = tr.getVec<double>("JetsAK8_NsubjettinessTau3");
        const auto& softDropMass          = tr.getVec<double>("JetsAK8_softDropMass");
        const auto& prunedMass            = tr.getVec<double>("JetsAK8_prunedMass");
        const auto& Muons                 = tr.getVec<TLorentzVector>("Muons");
//        const auto& GoodMuons             = tr.getVec<bool>("GoodMuons"+myVarSuffix_);
        const auto& Electrons             = tr.getVec<TLorentzVector>("Electrons");
//        const auto& GoodElectrons         = tr.getVec<bool>("GoodElectrons"+myVarSuffix_);
        const auto& NMuons                = tr.getVar<int>("NGoodMuons"+myVarSuffix_);
        const auto& NElectrons            = tr.getVar<int>("NGoodElectrons"+myVarSuffix_);
        const auto& GoodBJets_pt30        = tr.getVec<bool>("GoodBJets_pt30"+myVarSuffix_);

        const auto& TwoLep_Mbl1_Idx      = tr.getVar<std::pair<int, int>>("TwoLep_Mbl1_Idx"+myVarSuffix_);
        const auto& TwoLep_Mbl2_Idx      = tr.getVar<std::pair<int, int>>("TwoLep_Mbl2_Idx"+myVarSuffix_);
        
        const auto& Jets                 = tr.getVec<TLorentzVector>("Jets"+myVarSuffix_);
//        const auto& GoodLeptons          = tr.getVec<std::pair<std::string,TLorentzVector>>("GoodLeptons"+myVarSuffix_);
        const auto& NGoodLeptons         = tr.getVar<int>("NGoodLeptons_pt20"+myVarSuffix_);

        const auto& subjets              = tr.getVec<std::vector<TLorentzVector>>("JetsAK8_subjets");

        //First clean out leptons from JetsAK8 collection

        auto& GoodJetsAK8         = tr.createDerivedVec<bool>("GoodJetsAK8"+myVarSuffix_);
        for (unsigned int j = 0; j < JetsAK8.size(); j++)
        {
            GoodJetsAK8.push_back(true);
        }                
        
        if(NMuons > 0)
        {
            for(unsigned int mu=0; mu < Muons.size(); mu++)
            {
                double minDeltaR = 0.8;
                TLorentzVector myMuon = Muons.at(mu);
                int muonCand = -1;
                for(unsigned int j=0; j < JetsAK8.size(); j++)
                {
                    TLorentzVector myJet = JetsAK8.at(j);
//                    if( std::fabs(myMuon.Pt() - myJet.Pt()) / myMuon.Pt() < 1 && myMuon.DeltaR(myJet) < minDeltaR)
                    if(myMuon.DeltaR(myJet) < minDeltaR)
                    {
                        minDeltaR = myMuon.DeltaR(myJet);
                        muonCand = j;
                    }
                }
                if(muonCand != -1) GoodJetsAK8.at(muonCand) = false;
            }
        }

        if(NElectrons > 0)
        {
            for(unsigned int el=0; el < Electrons.size(); el++)
            {
                double minDeltaR = 0.8;
                TLorentzVector myElec = Electrons.at(el);
                int elecCand = -1;
                for(unsigned int j=0; j < JetsAK8.size(); j++)
                {
                    TLorentzVector myJet = JetsAK8.at(j);
                    // if( std::fabs(myElec.Pt() - myJet.Pt()) / myElec.Pt() < 1 && myElec.DeltaR(myJet) < minDeltaR)
                    if(myElec.DeltaR(myJet) < minDeltaR)
                    {
                        minDeltaR = myElec.DeltaR(myJet);
                        elecCand = j;                  
                    }
                }
                if(elecCand != -1) GoodJetsAK8.at(elecCand) = false;
            }
        }
        int NGoodJetsAK8 = 0;

        //Now apply eta, pt, SDM cuts
        for(unsigned int j=0; j < JetsAK8.size(); j++)
        {
            if(abs(JetsAK8.at(j).Eta()) > 2.4) GoodJetsAK8.at(j) = false;
            if(JetsAK8.at(j).Pt() < 30) GoodJetsAK8.at(j) = false;
            if(softDropMass.at(j) < 20) GoodJetsAK8.at(j) = false;
            if(GoodJetsAK8.at(j)) NGoodJetsAK8 += 1;
        }

        auto& candidateLSP         = tr.createDerivedVec<bool>("candidateLSP"+myVarSuffix_);
        auto& candidateLSP_TLV     = tr.createDerivedVec<TLorentzVector>("candidateLSP_TLV"+myVarSuffix_);
        auto& candidateLSP_Idx     = tr.createDerivedVec<int>("candidateLSP_Idx"+myVarSuffix_);
        auto& candidateLSP_T21     = tr.createDerivedVec<double>("candidateLSP_T21"+myVarSuffix_);
        auto& candidateLSP_T32     = tr.createDerivedVec<double>("candidateLSP_T32"+myVarSuffix_);
        auto& candidateLSP_Pruned  = tr.createDerivedVec<double>("candidateLSP_Pruned"+myVarSuffix_);
        auto& candidateLSP_SDM     = tr.createDerivedVec<double>("candidateLSP_SDM"+myVarSuffix_);
        
            
        candidateLSP = GoodJetsAK8;

        TLorentzVector bottom1;
        TLorentzVector bottom2;

        if (NGoodLeptons == 2)
        {
            bottom1 = Jets[TwoLep_Mbl1_Idx.first];
            bottom2 = Jets[TwoLep_Mbl2_Idx.first];
        }
        else
        {
            int b1_idx = -1, b2_idx = -1;
            for (unsigned int j = 0; j < Jets.size(); j++)
            {
                if (GoodBJets_pt30.at(j) && b1_idx == -1)
                {
                    b1_idx = j;
                }
                else if (GoodBJets_pt30.at(j) && b1_idx != -1 && b2_idx == -1)
                {
                    b2_idx = j;
                }

            }
            if (b1_idx != -1) bottom1 = Jets.at(b1_idx);
            if (b2_idx != -1) bottom2 = Jets.at(b2_idx);
        }

        
        int NCandidateLSP = 0;
        bool bottom1Found = false;
        bool bottom2Found = false;

        for(unsigned int j=0; j < GoodJetsAK8.size(); j++)
        {
            if (candidateLSP.at(j))
            {
                if(bottom1.DeltaR(JetsAK8.at(j)) < 0.2 && abs(1 - bottom1.Pt()/JetsAK8.at(j).Pt()) < 0.2)
                {
                    bottom1Found = true;
                    candidateLSP.at(j) = false;
                }
                else if(bottom2.DeltaR(JetsAK8.at(j)) < 0.2 && abs(1 - bottom2.Pt()/JetsAK8.at(j).Pt()) < 0.2)
                {
                    bottom2Found = true;
                    candidateLSP.at(j) = false;
                }
            }
        }
        
        for(unsigned int j = 0; j < GoodJetsAK8.size(); j++)
        {
            if(candidateLSP.at(j))
            {
                NCandidateLSP++;
                TLorentzVector myJetAK8 = JetsAK8.at(j);
                candidateLSP_Idx.push_back(j);
                candidateLSP_T21.push_back(Tau2.at(j) / Tau1.at(j));
                candidateLSP_T32.push_back(Tau3.at(j) / Tau2.at(j));
                candidateLSP_Pruned.push_back(prunedMass.at(j) >= 0 ? prunedMass.at(j): 0.0);
                candidateLSP_SDM.push_back(softDropMass.at(j));

                int b1_sj_idx = -1, b2_sj_idx = -1;
                for(unsigned int s = 0; s < subjets.at(j).size(); s++)
                {
                    TLorentzVector mySJ = subjets.at(j).at(s);
                    if(!bottom1Found && bottom1.DeltaR(mySJ) < 0.2 && abs(1 - bottom1.Pt()/mySJ.Pt()) < 0.5)
                    {
                        b1_sj_idx = s;
                        bottom1Found = true;
                    }
                    else if(!bottom2Found && bottom2.DeltaR(mySJ) < 0.2 && abs(1 - bottom2.Pt()/mySJ.Pt()) < 0.5)
                    {
                        b2_sj_idx = s;
                        bottom2Found = true;
                    }
                }
                if (b1_sj_idx != -1) candidateLSP_TLV.push_back(myJetAK8 - subjets.at(j).at(b1_sj_idx));
                else if (b2_sj_idx != -1) candidateLSP_TLV.push_back(myJetAK8 - subjets.at(j).at(b2_sj_idx));
                else candidateLSP_TLV.push_back(myJetAK8);
            }
        } 

        asymm_mt2_lester_bisect::disableCopyrightMessage();
        tr.registerDerivedVar("NGoodJetsAK8"+myVarSuffix_, NGoodJetsAK8);
        tr.registerDerivedVar("NCandidateLSP"+myVarSuffix_, NCandidateLSP);
        
        }
    

public:
    FatJetCombine(std::string myVarSuffix = "") 
        : myVarSuffix_(myVarSuffix)
    {
        std::cout<<"Setting up FatJetCombine"<<std::endl;   
    }
    
    void operator()(NTupleReader& tr)
    {
        fatjetcombine(tr);
    }
};

#endif
