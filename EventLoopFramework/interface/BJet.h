#ifndef BJET_H
#define BJET_H

class BJet
{
private:
    std::string myVarSuffix_;

    void setVar(bool pass, std::vector<bool>& boolVec, int& n)
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

    void bjet(NTupleReader& tr)
    {
        const auto& Jets = tr.getVec<TLorentzVector>("Jets"+myVarSuffix_);
        const auto& Jets_bJetTagDeepCSVtotb = tr.getVec<double>("Jets"+myVarSuffix_+"_bJetTagDeepCSVtotb");
        const auto& etaCut = tr.getVar<double>("etaCut");
        const auto& JetsID = tr.getVec<bool>("Jets"+myVarSuffix_+"_ID");
        const auto& GoodJets = tr.getVec<bool>("GoodJets"+myVarSuffix_);
        const auto& loose = tr.getVar<double>("deepCSV_WP_loose");
        const auto& medium = tr.getVar<double>("deepCSV_WP_medium");
        const auto& tight = tr.getVar<double>("deepCSV_WP_tight");

        auto& bjets_loose_ = tr.createDerivedVec<bool>("BJets_loose"+myVarSuffix_);
        auto& bjets_pt30_loose_ = tr.createDerivedVec<bool>("BJets_pt30_loose"+myVarSuffix_);
        int NBJets_loose = 0, NBJets_pt30_loose = 0;

        auto& bjets_ = tr.createDerivedVec<bool>("BJets"+myVarSuffix_);
        auto& bjets_pt30_ = tr.createDerivedVec<bool>("BJets_pt30"+myVarSuffix_);
        auto& bjets_pt40_ = tr.createDerivedVec<bool>("BJets_pt40"+myVarSuffix_);
        auto& bjets_pt45_ = tr.createDerivedVec<bool>("BJets_pt45"+myVarSuffix_);
        int NBJets = 0, NBJets_pt30 = 0, NBJets_pt40 = 0, NBJets_pt45 = 0;

        auto& bjets_tight_ = tr.createDerivedVec<bool>("BJets_tight"+myVarSuffix_);
        auto& bjets_pt30_tight_ = tr.createDerivedVec<bool>("BJets_pt30_tight"+myVarSuffix_);
        auto& bjets_pt45_tight_ = tr.createDerivedVec<bool>("BJets_pt45_tight"+myVarSuffix_);
        int NBJets_tight = 0, NBJets_pt30_tight = 0, NBJets_pt45_tight = 0;

        auto& goodbjets_loose_ = tr.createDerivedVec<bool>("GoodBJets_loose"+myVarSuffix_);
        auto& goodbjets_pt30_loose_ = tr.createDerivedVec<bool>("GoodBJets_pt30_loose"+myVarSuffix_);
        int NGoodBJets_loose = 0, NGoodBJets_pt30_loose = 0;

        auto& goodbjets_ = tr.createDerivedVec<bool>("GoodBJets"+myVarSuffix_);
        auto& goodbjets_pt30_ = tr.createDerivedVec<bool>("GoodBJets_pt30"+myVarSuffix_);
        auto& goodbjets_pt40_ = tr.createDerivedVec<bool>("GoodBJets_pt40"+myVarSuffix_);
        auto& goodbjets_pt45_ = tr.createDerivedVec<bool>("GoodBJets_pt45"+myVarSuffix_);
        int NGoodBJets = 0, NGoodBJets_pt30 = 0, NGoodBJets_pt40 = 0, NGoodBJets_pt45 = 0;

        auto& goodbjets_tight_ = tr.createDerivedVec<bool>("GoodBJets_tight"+myVarSuffix_);
        auto& goodbjets_pt30_tight_ = tr.createDerivedVec<bool>("GoodBJets_pt30_tight"+myVarSuffix_);
        auto& goodbjets_pt45_tight_ = tr.createDerivedVec<bool>("GoodBJets_pt45_tight"+myVarSuffix_);
        int NGoodBJets_tight = 0, NGoodBJets_pt30_tight = 0, NGoodBJets_pt45_tight = 0;

        for (unsigned int ijet = 0; ijet < Jets.size(); ++ijet)
        {
            TLorentzVector lv = Jets.at(ijet);
            double bdisc = Jets_bJetTagDeepCSVtotb.at(ijet);

            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && bdisc > loose                , bjets_loose_,      NBJets_loose     );
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && bdisc > loose && lv.Pt() > 30, bjets_pt30_loose_, NBJets_pt30_loose);

            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && bdisc > loose                 && GoodJets.at(ijet), goodbjets_loose_,      NGoodBJets_loose     );
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && bdisc > loose && lv.Pt() > 30 && GoodJets.at(ijet), goodbjets_pt30_loose_, NGoodBJets_pt30_loose);

            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && bdisc > medium                , bjets_,      NBJets     );
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && bdisc > medium && lv.Pt() > 30, bjets_pt30_, NBJets_pt30);
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && bdisc > medium && lv.Pt() > 40, bjets_pt40_, NBJets_pt40);
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && bdisc > medium && lv.Pt() > 45, bjets_pt45_, NBJets_pt45);

            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && bdisc > medium                 && GoodJets.at(ijet), goodbjets_,      NGoodBJets     );
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && bdisc > medium && lv.Pt() > 30 && GoodJets.at(ijet), goodbjets_pt30_, NGoodBJets_pt30);
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && bdisc > medium && lv.Pt() > 40 && GoodJets.at(ijet), goodbjets_pt40_, NGoodBJets_pt40);
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && bdisc > medium && lv.Pt() > 45 && GoodJets.at(ijet), goodbjets_pt45_, NGoodBJets_pt45);

            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && bdisc > tight                , bjets_tight_,      NBJets_tight     );
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && bdisc > tight && lv.Pt() > 30, bjets_pt30_tight_, NBJets_pt30_tight);
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && bdisc > tight && lv.Pt() > 45, bjets_pt45_tight_, NBJets_pt45_tight);

            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && bdisc > tight                 && GoodJets.at(ijet), goodbjets_tight_,      NGoodBJets_tight     );
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && bdisc > tight && lv.Pt() > 30 && GoodJets.at(ijet), goodbjets_pt30_tight_, NGoodBJets_pt30_tight);
            setVar(JetsID.at(ijet) && abs(lv.Eta()) < etaCut && bdisc > tight && lv.Pt() > 45 && GoodJets.at(ijet), goodbjets_pt45_tight_, NGoodBJets_pt45_tight);            
        }

        tr.registerDerivedVar("NBJets_loose"+myVarSuffix_,       NBJets_loose);
        tr.registerDerivedVar("NBJets_pt30_loose"+myVarSuffix_,  NBJets_pt30_loose);

        tr.registerDerivedVar("NBJets"+myVarSuffix_,       NBJets);
        tr.registerDerivedVar("NBJets_pt30"+myVarSuffix_,  NBJets_pt30);
        tr.registerDerivedVar("NBJets_pt40"+myVarSuffix_,  NBJets_pt40);
        tr.registerDerivedVar("NBJets_pt45"+myVarSuffix_,  NBJets_pt45);

        tr.registerDerivedVar("NBJets_tight"+myVarSuffix_,       NBJets_tight);
        tr.registerDerivedVar("NBJets_pt30_tight"+myVarSuffix_,  NBJets_pt30_tight);
        tr.registerDerivedVar("NBJets_pt45_tight"+myVarSuffix_,  NBJets_pt45_tight);

        tr.registerDerivedVar("NGoodBJets_loose"+myVarSuffix_,       NGoodBJets_loose);
        tr.registerDerivedVar("NGoodBJets_pt30_loose"+myVarSuffix_,  NGoodBJets_pt30_loose);

        tr.registerDerivedVar("NGoodBJets"+myVarSuffix_,       NGoodBJets);
        tr.registerDerivedVar("NGoodBJets_pt30"+myVarSuffix_,  NGoodBJets_pt30);
        tr.registerDerivedVar("NGoodBJets_pt40"+myVarSuffix_,  NGoodBJets_pt40);
        tr.registerDerivedVar("NGoodBJets_pt45"+myVarSuffix_,  NGoodBJets_pt45);
        
        tr.registerDerivedVar("NGoodBJets_tight"+myVarSuffix_,       NGoodBJets_tight);
        tr.registerDerivedVar("NGoodBJets_pt30_tight"+myVarSuffix_,  NGoodBJets_pt30_tight);
        tr.registerDerivedVar("NGoodBJets_pt45_tight"+myVarSuffix_,  NGoodBJets_pt45_tight);
    }

public:
    BJet( std::string myVarSuffix = "" ) 
        : myVarSuffix_(myVarSuffix)
    {
        std::cout<<"Setting up BJet"<<std::endl;
    }

    void operator()(NTupleReader& tr)
    {
        bjet(tr);
    }
};

#endif
