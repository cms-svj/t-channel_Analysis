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

        const std::vector<std::string> weightVecAll = {"Lumi", "Weight"};

        my_Histos.emplace_back(new Histo1D("EventCounter", 2, -1.1, 1.1, "eventCounter", {}, {}));
        my_Histos.emplace_back(new Histo1D("h_njets",      20,    0.0,   20.0, "NGoodJetsAK8_pt200",                {"passBaselineSVJ"},       weightVecAll));
        my_Histos.emplace_back(new Histo1D("h_ht",         500,   0.0, 5000.0, "HT_trigger_pt30",                   {"passBaselineSVJ"},       weightVecAll));
        my_Histos.emplace_back(new Histo1D("h_htAK8",      500,   0.0, 5000.0, "HT_pt200",                          {"passBaselineSVJ"},       weightVecAll));
        my_Histos.emplace_back(new Histo1D("h_st",         500,   0.0, 5000.0, "ST_pt30",                           {"passBaselineSVJ"},       weightVecAll));
        my_Histos.emplace_back(new Histo1D("h_stAK8",      500,   0.0, 5000.0, "ST_pt200",                          {"passBaselineSVJ"},       weightVecAll));
        my_Histos.emplace_back(new Histo1D("h_met",        500,   0.0, 5000.0, "MET",                               {"passBaselineSVJ"},       weightVecAll));
        my_Histos.emplace_back(new Histo1D("h_jPt",        200,   0.0, 2000.0, "JetsAK8.Pt()", "GoodJetsAK8_pt200", {"passBaselineSVJ"},       weightVecAll));
        my_Histos.emplace_back(new Histo1D("h_jEta",       200,  -6.0,    6.0, "JetsAK8.Eta()","GoodJetsAK8_pt200", {"passBaselineSVJ"},       weightVecAll));
        my_Histos.emplace_back(new Histo1D("h_jPhi",       200,  -4.0,    4.0, "JetsAK8.Phi()","GoodJetsAK8_pt200", {"passBaselineSVJ"},       weightVecAll));
        my_Histos.emplace_back(new Histo1D("h_j1Pt",       200,   0.0, 2000.0, "leadingJetAK8_pt200.Pt()",          {"passBaselineSVJ"},       weightVecAll));
        my_Histos.emplace_back(new Histo1D("h_j2Pt",       200,   0.0, 2000.0, "subleadingJetAK8_pt200.Pt()",       {"passBaselineSVJ"},       weightVecAll));
        //my_Histos.emplace_back(new Histo1D("h_drJ12",      100,   0.0,   10.0, "drJ12",                             {"passBaselineSVJ"},       weightVecAll));
        //my_Histos.emplace_back(new Histo1D("h_drJ1MET",    100,   0.0,   10.0, "drJ1MET",                           {"passBaselineSVJ"},       weightVecAll));
        //my_Histos.emplace_back(new Histo1D("h_drJ2MET",    100,   0.0,   10.0, "drJ2MET",                           {"passBaselineSVJ"},       weightVecAll));
        my_Histos.emplace_back(new Histo1D("h_mjjPt",      200,   0.0, 2000.0, "mjjTLV.Pt()",                       {"passBaselineSVJ"},       weightVecAll));
        my_Histos.emplace_back(new Histo1D("h_mjjM",       200,   0.0, 2000.0, "mjjTLV.M()",                        {"passBaselineSVJ"},       weightVecAll));
        my_Histos.emplace_back(new Histo1D("h_mt",         200,   0.0, 2000.0, "mT",                                {"passBaselineSVJ"},       weightVecAll));
        //my_Histos.emplace_back(new Histo2D<double,double>("h_JetAK81_r_eta_phi", 200, -6.0, 6.0, "JetAK81_r_Eta", 200, -4.0, 4.0, "JetAK81_r_Phi", {"passBaselineSVJ"}, weightVecAll));
        //my_Histos.emplace_back(new Histo2D<double,double>("h_JetAK82_r_eta_phi", 200, -6.0, 6.0, "JetAK82_r_Eta", 200, -4.0, 4.0, "JetAK82_r_Phi", {"passBaselineSVJ"}, weightVecAll));
        //my_Histos.emplace_back(new Histo2D<double,double>("h_lvMET_r_eta_phi",   200, -6.0, 6.0, "lvMET_r_Eta",   200, -4.0, 4.0, "lvMET_r_Phi",   {"passBaselineSVJ"}, weightVecAll));
        //my_Histos.emplace_back(new Histo2D<double,double>("h_JetAK8_r_eta1_eta2", 200, -6.0, 6.0, "JetAK81_r_Eta", 200, -6.0, 6.0, "JetAK82_r_Eta", {"passBaselineSVJ"}, weightVecAll));
        //my_Histos.emplace_back(new Histo2D<double,double>("h_JetAK8_r_phi1_phi2", 200, -4.0, 4.0, "JetAK81_r_Phi", 200, -4.0, 4.0, "JetAK82_r_Phi", {"passBaselineSVJ"}, weightVecAll));
        //my_Histos.emplace_back(new Histo2D<double,double>("h_JetAK8MET_r_eta1_eta", 200, -6.0, 6.0, "JetAK81_r_Eta", 200, -6.0, 6.0, "lvMET_r_Eta", {"passBaselineSVJ"}, weightVecAll));
        //my_Histos.emplace_back(new Histo2D<double,double>("h_JetAK8MET_r_phi1_phi", 200, -4.0, 4.0, "JetAK81_r_Phi", 200, -4.0, 4.0, "lvMET_r_Phi", {"passBaselineSVJ"}, weightVecAll));
        //my_Histos.emplace_back(new Histo2D<double,double>("h_JetAK8MET_r_eta2_eta", 200, -6.0, 6.0, "JetAK82_r_Eta", 200, -6.0, 6.0, "lvMET_r_Eta", {"passBaselineSVJ"}, weightVecAll));
        //my_Histos.emplace_back(new Histo2D<double,double>("h_JetAK8MET_r_phi2_phi", 200, -4.0, 4.0, "JetAK82_r_Phi", 200, -4.0, 4.0, "lvMET_r_Phi", {"passBaselineSVJ"}, weightVecAll));

        //my_Histos.emplace_back(new Histo3D<int,double,int>("h_njets_deepESM_nb", 20, 0.0, 20.0, "NGoodJets_pt30", 200, 0.0, 5000.0, "HT_trigger_pt30", 10, 0.0, 10.0, "NGoodBJets_pt30", {}, weightVecAll));
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


            const auto& NGoodJetsAK8_pt200 = tr.getVar<int>("NGoodJetsAK8_pt200");
            tr.registerDerivedVar("passBaselineSVJ", NGoodJetsAK8_pt200 == 2);

               
            //-----------------------------------
            //-- Fill Histograms Below
            //-----------------------------------
            Fill(tr);

        }//END of while tr.getNextEvent loop   
    }//END of function
};

#endif
