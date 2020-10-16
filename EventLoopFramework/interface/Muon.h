#ifndef MUON_H
#define MUON_H

#include "t-channel_Analysis/EventLoopFramework/interface/Utility.h"

class Muon
{
private:
    std::string myVarSuffix_;

    void muon(NTupleReader& tr)
    {
        const auto& allMuons = tr.getVec<TLorentzVector>("Muons");
        const auto& allMuons_passIso = tr.getVec<bool>("Muons_passIso");
        const auto& allMuons_charge = tr.getVec<int>("Muons_charge");
        const auto& allMuons_medID = tr.getVec<bool>("Muons_mediumID");
        const auto& etaCut = tr.getVar<double>("etaCut");
        const auto& MET = tr.getVar<double>("MET"); 
        const auto& METPhi = tr.getVar<double>("METPhi");

        TLorentzVector lvMET;
        lvMET.SetPtEtaPhiM(MET, 0.0, METPhi, 0.0);

        auto* good_muons_ = new std::vector<bool>();
        auto* good_muons_pt20 = new std::vector<bool>();
        auto* muons_mtw_ = new std::vector<double>;
        int NGoodMuons = 0;
        int NGoodPlusMuons = 0;
        int NGoodMinusMuons = 0;
        int NGoodMuons_pt20 = 0;
        int NGoodPlusMuons_pt20 = 0;
        int NGoodMinusMuons_pt20 = 0;
        
        auto* noniso_muons_ = new std::vector<bool>();
        auto* GoodNonIsoMuons = new std::vector<std::pair<std::string, TLorentzVector>>();
        int NNonIsoMuons = 0;
        
        for(unsigned int imu = 0; imu < allMuons.size(); ++imu)
        {            
            TLorentzVector lvmu(allMuons.at(imu));
            double mtw = sqrt( 2*( lvMET.Pt()*lvmu.Pt() - (lvMET.Px()*lvmu.Px() + lvMET.Py()*lvmu.Py()) ) );
            muons_mtw_->push_back(mtw);
            if( abs(lvmu.Eta()) < etaCut && 
                lvmu.Pt() > 30 && 
                allMuons_passIso.at(imu) &&
                allMuons_medID.at(imu)
                )
            {
                good_muons_->push_back(true);
                NGoodMuons++;
                if( allMuons_charge.at(imu) ==  1 ) NGoodPlusMuons++;
                else if( allMuons_charge.at(imu) == -1 ) NGoodMinusMuons++;
                else std::cout<<utility::color("Charge values in nTuples are different", "red")<<std::endl;
            }
            else
            {
                good_muons_->push_back(false);
            }
            if( abs(lvmu.Eta()) < etaCut &&
                     lvmu.Pt() > 20 &&
                     allMuons_passIso.at(imu) &&
                     allMuons_medID.at(imu)
                )
            {
                good_muons_pt20->push_back(true);
                NGoodMuons_pt20 ++;
                if( allMuons_charge.at(imu) ==  1 ) NGoodPlusMuons_pt20++;
                else if( allMuons_charge.at(imu) == -1 ) NGoodMinusMuons_pt20++;
                else std::cout<<utility::color("Charge values in nTuples are different", "red")<<std::endl;
            }
            else
            {
                good_muons_pt20->push_back(false);
            }
            if( abs(lvmu.Eta()) < etaCut &&
                lvmu.Pt() > 55.0 &&
                !allMuons_passIso.at(imu) &&
                allMuons_medID.at(imu)
                )
            {
                TLorentzVector muon = allMuons.at(imu);
                noniso_muons_->push_back(true);
                GoodNonIsoMuons->push_back( std::make_pair("n", muon) );
                NNonIsoMuons++;
            }
            else
            {
                noniso_muons_->push_back(false);
            }
        }
        
        tr.registerDerivedVec("GoodMuons"+myVarSuffix_,  good_muons_);
        tr.registerDerivedVar("NGoodMuons"+myVarSuffix_, NGoodMuons );
        tr.registerDerivedVar("NGoodPlusMuons"+myVarSuffix_, NGoodPlusMuons );
        tr.registerDerivedVar("NGoodMinusMuons"+myVarSuffix_, NGoodMinusMuons );
        tr.registerDerivedVec("MuonsMTW"+myVarSuffix_,   muons_mtw_ );

        tr.registerDerivedVec("GoodMuons_pt20"+myVarSuffix_,  good_muons_pt20);
        tr.registerDerivedVar("NGoodMuons_pt20"+myVarSuffix_, NGoodMuons_pt20 );

        tr.registerDerivedVec("NonIsoMuons"+myVarSuffix_, noniso_muons_ );
        tr.registerDerivedVar("NNonIsoMuons"+myVarSuffix_, NNonIsoMuons );
        tr.registerDerivedVec("GoodNonIsoMuons"+myVarSuffix_, GoodNonIsoMuons );
    }

public:
    Muon(std::string myVarSuffix = "") 
        : myVarSuffix_(myVarSuffix)
    {
        std::cout<<"Setting up Muon"<<std::endl;   
    }
    
    void operator()(NTupleReader& tr)
    {
        muon(tr);
    }
};

#endif
