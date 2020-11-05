#ifndef FATJETCOMBINE_H
#define FATJETCOMBINE_H

#include "t-channel_Analysis/EventLoopFramework/interface/Utility.h"
#include "TMath.h"

class FatJetCombine
{
private:
    std::string myVarSuffix_;

    double mT(TLorentzVector v, TLorentzVector met)
    {
        double m = v.M();
        double et = met.Et();
        double MET = met.Pt();
        return TMath::Sqrt(m*m + 2*(TMath::Sqrt(m*m+MET*MET)*et - v.Vect().Dot(met.Vect())));
    }

    void fatjetcombine(NTupleReader& tr)
    {
        const auto& JetsAK8 = tr.getVec<TLorentzVector>("JetsAK8"+myVarSuffix_);
        const auto& MET     = tr.getVar<double>("MET");
        const auto& METPhi = tr.getVar<double>("METPhi");
        const auto& etaCut  = tr.getVar<double>("etaCut");
        const auto& HT_pt30 = tr.getVar<double>("HT_trigger_pt30");

        TLorentzVector lvMET;
        lvMET.SetPtEtaPhiM(MET, 0.0, METPhi, 0.0);

        auto& GoodJetsAK8  = tr.createDerivedVec<bool>("GoodJetsAK8_pt200"+myVarSuffix_, JetsAK8.size(), true);
        auto& NGoodJetsAK8 = tr.createDerivedVar<int>("NGoodJetsAK8_pt200"+myVarSuffix_, 0);
        auto& HT           = tr.createDerivedVar<double>("HT_pt200"+myVarSuffix_, 0.0);        

        //Now apply eta, pt cuts        
        std::vector<TLorentzVector> GoodJetsAK8_TLV;
        for(unsigned int j=0; j < JetsAK8.size(); j++)
        {
            double pt = JetsAK8.at(j).Pt();
            double eta = JetsAK8.at(j).Eta();
            bool cut = abs(eta) > etaCut || pt < 200;
            //bool cut = abs(eta) > 5.0 || pt < 170;
            if(cut)
            {
                GoodJetsAK8.at(j) = false;
            }
            else
            {
                NGoodJetsAK8 += 1;
                HT += pt;
                GoodJetsAK8_TLV.push_back(JetsAK8.at(j));
            }
        }
        std::sort( GoodJetsAK8_TLV.begin(), GoodJetsAK8_TLV.end(), [](const TLorentzVector& v1, const TLorentzVector& v2){return v1.Pt() > v2.Pt();} );

        //Register variables
        double dPhiMin = 10000.0;
        TLorentzVector leadingJetAK8, subleadingJetAK8;
        for(unsigned int i=0; i < GoodJetsAK8_TLV.size(); i++)
        {
            if(i==0) leadingJetAK8 = GoodJetsAK8_TLV.at(i);
            else if(i==1) subleadingJetAK8 = GoodJetsAK8_TLV.at(i);
            
            double dPhi = GoodJetsAK8_TLV.at(i).DeltaPhi(lvMET);
            dPhiMin = std::min(dPhiMin, dPhi);
        }

        tr.registerDerivedVar("leadingJetAK8_pt200"+myVarSuffix_, leadingJetAK8);
        tr.registerDerivedVar("subleadingJetAK8_pt200"+myVarSuffix_, subleadingJetAK8);               
        tr.registerDerivedVar("dRJ12"+myVarSuffix_, leadingJetAK8.DeltaR(subleadingJetAK8));
        tr.registerDerivedVar("dEtaJ12"+myVarSuffix_, leadingJetAK8.Eta() - subleadingJetAK8.Eta());
        tr.registerDerivedVar("dPhiJ1MET"+myVarSuffix_, leadingJetAK8.DeltaPhi(lvMET));
        tr.registerDerivedVar("dPhiJ2MET"+myVarSuffix_, subleadingJetAK8.DeltaPhi(lvMET));
        tr.registerDerivedVar("dPhiMinJMET"+myVarSuffix_, dPhiMin);
        tr.registerDerivedVar("ST_pt200"+myVarSuffix_, HT + MET);
        tr.registerDerivedVar("ST_pt30"+myVarSuffix_, HT_pt30 + MET);
        tr.registerDerivedVar("METrHT_pt30"+myVarSuffix_, MET / HT_pt30);
        tr.registerDerivedVar("METrST_pt30"+myVarSuffix_, MET /(HT_pt30 + MET));

        //Now calculate more intresting variables
        auto& mjjTLV  = tr.createDerivedVar<TLorentzVector>("mjjTLV"+myVarSuffix_, leadingJetAK8 + subleadingJetAK8);
        tr.registerDerivedVar("mT"+myVarSuffix_, mT( mjjTLV, lvMET) );
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
