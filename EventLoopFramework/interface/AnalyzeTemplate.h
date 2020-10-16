#ifndef AnalyzeTemplate_h
#define AnalyzeTemplate_h

#include "t-channel_Analysis/EventLoopFramework/interface/AnalyzeBase.h"
#include "t-channel_Analysis/EventLoopFramework/interface/Histo.h"
class NTupleReader;

class AnalyzeTemplate : public AnalyzeBase
{
public:
    void InitHistos()
    {
        TH1::SetDefaultSumw2();
        TH2::SetDefaultSumw2();

        const std::vector<std::string> weightVecAll = {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", 
                                                       "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"};
        const std::vector<std::string> weightVecAllQCD = {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", 
                                                          "totGoodMuonSF", "puWeightCorr", "prefiringScaleFactor"};

        my_Histos.emplace_back(new Histo1D("EventCounter", 2, -1.1, 1.1, "eventCounter", {}, {}));
        my_Histos.emplace_back(new Histo1D("h_njets_passBaseline1l",       20,   0.0,   20.0, "NGoodJets_pt30",             {"passBaseline1l_Good"},       weightVecAll));
        my_Histos.emplace_back(new Histo1D("h_njets_passQCDBaseline",      20,   0.0,   20.0, "NNonIsoMuonJets_pt30",       {"passBaseline1l_NonIsoMuon"}, weightVecAllQCD));
        my_Histos.emplace_back(new Histo1D("h_nb_passBaseline1l",          10,   0.0,   10.0, "NGoodBJets_pt30",            {"passBaseline1l_Good"},       weightVecAll));
        my_Histos.emplace_back(new Histo1D("h_ht_passBaseline1l",         500,   0.0, 5000.0, "HT_trigger_pt30",            {"passBaseline1l_Good"},       weightVecAll));
        my_Histos.emplace_back(new Histo1D("h_ht_passQCDBaseline",        500,   0.0, 5000.0, "HT_NonIsoMuon_pt30",         {"passBaseline1l_NonIsoMuon"}, weightVecAllQCD));
        my_Histos.emplace_back(new Histo1D("h_mbl_passBaseline1l",        300,   0.0,  300.0, "Mbl",                        {"passBaseline1l_Good"},       weightVecAll));
        my_Histos.emplace_back(new Histo1D("h_jPt_passBaseline1l",        200,   0.0, 2000.0, "Jets.Pt()", "GoodJets_pt30", {"passBaseline1l_Good"},       weightVecAll));
        my_Histos.emplace_back(new Histo1D("h_jEta_passBaseline1l",       200,  -6.0,    6.0, "Jets.Eta()","GoodJets_pt30", {"passBaseline1l_Good"},       weightVecAll));
        my_Histos.emplace_back(new Histo1D("h_jPhi_passBaseline1l",       200,  -4.0,    4.0, "Jets.Phi()","GoodJets_pt30", {"passBaseline1l_Good"},       weightVecAll));
        my_Histos.emplace_back(new Histo2D<int,double>("h_njets_ht_passBaseline1l", 20, 0.0, 20.0, "NGoodJets_pt30", 200, 0.0, 5000.0, "HT_trigger_pt30", {"passBaseline1l_Good"}, weightVecAll));
        my_Histos.emplace_back(new Histo3D<int,double,int>("h_njets_deepESM_nb_passBaseline1l", 20, 0.0, 20.0, "NGoodJets_pt30", 200, 0.0, 5000.0, "HT_trigger_pt30", 10, 0.0, 10.0, "NGoodBJets_pt30", {"passBaseline1l_Good"}, weightVecAll));
    }

    void Loop(NTupleReader& tr, double, int maxevents, bool)
    {
        //-------------------------------------
        //-- Initialize histograms to be filled
        //-------------------------------------
        InitHistos();

        while( tr.getNextEvent() )
        {
            //------------------------------------
            //-- Print Event Number
            //------------------------------------
            const bool breakLoop = printEventNum(maxevents, tr.getEvtNum());
            if(breakLoop) break;
               
            //-----------------------------------
            //-- Fill Histograms Below
            //-----------------------------------
            Fill(tr);

        }//END of while tr.getNextEvent loop   
    }//END of function
};

#endif
