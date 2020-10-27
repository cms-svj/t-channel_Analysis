#define Analyze_cxx
#include "t-channel_Analysis/EventLoopFramework/interface/Analyze.h"
#include "t-channel_Analysis/EventLoopFramework/interface/NTupleReader.h"
#include "t-channel_Analysis/EventLoopFramework/interface/Utility.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TProfile2D.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>

Analyze::Analyze() : initHistos(false)
{
}

void Analyze::InitHistos(const std::map<std::string, bool>& cutMap, const std::vector<TH1DInfo>& histInfos, 
                             const std::vector<TH2DInfo>& hist2DInfos,  const std::vector<TH2DProfileInfo>& hist2DProfileInfos)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    my_histos.emplace("EventCounter", std::make_shared<TH1D>("EventCounter","EventCounter", 2, -1.1, 1.1 ) );

    for(auto& mycut : cutMap)
    {
        for(const auto& hInfo : histInfos)
        { 
            my_histos.emplace(hInfo.name+mycut.first, 
                              std::make_shared<TH1D>((hInfo.name+mycut.first).c_str(),(hInfo.name+mycut.first).c_str(), hInfo.nBins, hInfo.low, hInfo.high));
        }

        for(const auto& h2dInfo : hist2DInfos)
        {
            my_2d_histos.emplace(h2dInfo.name+mycut.first, 
                                 std::make_shared<TH2D>((h2dInfo.name+mycut.first).c_str(),(h2dInfo.name+mycut.first).c_str(), 
                                                        h2dInfo.nBinsX, h2dInfo.lowX, h2dInfo.highX, h2dInfo.nBinsY, h2dInfo.lowY, h2dInfo.highY));
        }

        for(const auto& h2dProfile : hist2DProfileInfos)
        {
            my_2d_tp_histos.emplace(h2dProfile.name+mycut.first,
                                    std::make_shared<TProfile2D>((h2dProfile.name+mycut.first).c_str(),(h2dProfile.name+mycut.first).c_str(), 
                                                                 h2dProfile.nBinsX, h2dProfile.lowX, h2dProfile.highX, h2dProfile.nBinsY, h2dProfile.lowY, h2dProfile.highY, h2dProfile.lowZ, h2dProfile.highZ));
        }
    }

    my_histos.emplace( "h_cutFlow", std::make_shared<TH1D>("h_cutFlow", "h_cutFlow", 9,0,9));    
}

void Analyze::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        const auto& runtype                   = tr.getVar<std::string>("runtype");     
        const auto& filetag                   = tr.getVar<std::string>("filetag");
        const auto& eventCounter              = tr.getVar<int>("eventCounter");

        const auto& MET                       = tr.getVar<double>("MET");
        const auto& METPhi                    = tr.getVar<double>("METPhi");
        const auto& ntops                     = tr.getVar<int>("ntops");

        const auto& Jets                      = tr.getVec<TLorentzVector>("Jets");
        const auto& GoodJets_pt30             = tr.getVec<bool>("GoodJets_pt30");
        const auto& NGoodJets_pt30            = tr.getVar<int>("NGoodJets_pt30");
        const auto& NGoodBJets_pt30           = tr.getVar<int>("NGoodBJets_pt30");
        const auto& HT_trigger_pt30           = tr.getVar<double>("HT_trigger_pt30");
        const auto& ST_pt30                   = tr.getVar<double>("ST_pt30");
        const auto& METrHT_pt30               = tr.getVar<double>("METrHT_pt30");
        const auto& METrST_pt30               = tr.getVar<double>("METrST_pt30");
        const auto& JetID                     = tr.getVar<bool>("JetID");

        const auto& JetsAK8                   = tr.getVec<TLorentzVector>("JetsAK8");
        const auto& GoodJetsAK8_pt200         = tr.getVec<bool>("GoodJetsAK8_pt200");
        const auto& NGoodJetsAK8_pt200        = tr.getVar<int>("NGoodJetsAK8_pt200");
        const auto& JetIDAK8                  = tr.getVar<bool>("JetIDAK8");

        const auto& dEtaJ12                   = tr.getVar<double>("dEtaJ12");
        const auto& dRJ12                     = tr.getVar<double>("dRJ12");
        const auto& dPhiJ1MET                 = tr.getVar<double>("dPhiJ1MET");
        const auto& dPhiJ2MET                 = tr.getVar<double>("dPhiJ2MET");
        const auto& dPhiMinJMET               = tr.getVar<double>("dPhiMinJMET");
        const auto& mjjTLV                    = tr.getVar<TLorentzVector>("mjjTLV");
        const auto& mT                        = tr.getVar<double>("mT");

        // ------------------------
        // -- Print event number
        // ------------------------       
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if(tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() );

        // ------------------------
        // -- Define weight
        // ------------------------
        double weight=1.0, eventweight=1.0;
        if(runtype == "MC")
        {
            // Define Lumi weight
            const auto& Weight  = tr.getVar<double>("Weight");
            const auto& lumi = tr.getVar<double>("Lumi");
            eventweight = lumi*Weight;
            weight *= eventweight;
        }

        // -------------------------------
        // -- Define cuts
        // -------------------------------
        bool pass_general = true;
        
        // -------------------
        // --- Fill Histos ---
        // -------------------
        const std::map<std::string, bool> cut_map_1l 
        {
            { "",                    pass_general                                                                                       },
            { "_1AK8j",              pass_general && NGoodJetsAK8_pt200 == 1                                                            },
            { "_2AK8j",              pass_general && NGoodJetsAK8_pt200 == 2                                                            },
            { "_3AK8j",              pass_general && NGoodJetsAK8_pt200 == 3                                                            },
            { "_ge2AK8j_l2METrHT",   pass_general && NGoodJetsAK8_pt200 >= 2 && METrHT_pt30 < 2.0                                       },
            { "_ge2AK8j_lp6METrST",  pass_general && NGoodJetsAK8_pt200 >= 2 && METrST_pt30 < 0.6                                       },
            { "_ge2AK8j_l1p5dEta12", pass_general && NGoodJetsAK8_pt200 >= 2 && dEtaJ12 < 1.5                                           },
            { "_baseline",           pass_general && NGoodJetsAK8_pt200 >= 2 && METrHT_pt30 < 2.0 && METrST_pt30 < 0.6 && dEtaJ12 < 1.5 },
            { "_1j",                 pass_general && NGoodJets_pt30 == 1                                                                },
            { "_2j",                 pass_general && NGoodJets_pt30 == 2                                                                },
            { "_3j",                 pass_general && NGoodJets_pt30 == 3                                                                },
            { "_ge2j",               pass_general && NGoodJets_pt30 >= 2                                                                },
        };

        std::vector<TH1DInfo> histInfos = {
            {    "h_njets",             20,   0.0,   20.0},
            {    "h_njetsAK8",          20,   0.0,   20.0},
            {    "h_ntops",             10,   0.0,   10.0},
            {    "h_nb",                10,   0.0,   10.0},
            {    "h_ht",               500,   0.0, 5000.0},
            {    "h_st",               500,   0.0, 5000.0},
            {    "h_met",              500,   0.0, 5000.0},
            {    "h_jPt",              200,   0.0, 2000.0},
            {    "h_jEta",             200,  -6.0,    6.0},
            {    "h_jPhi",             200,  -4.0,    4.0},
            {    "h_jPtAK8",           200,   0.0, 2000.0},
            {    "h_jEtaAK8",          200,  -6.0,    6.0},
            {    "h_jPhiAK8",          200,  -4.0,    4.0},
            {    "h_weight",           200,  -5.0,    5.0},
            {    "h_dEtaJ12",          200,   0.0,   10.0},
            {    "h_dRJ12",            100,   0.0,   10.0},
            {    "h_dPhiJ1MET",        100,   0.0,    4.0},
            {    "h_dPhiJ2MET",        100,   0.0,    4.0},
            {    "h_dPhiMinJMET",      100,   0.0,    4.0},
            {    "h_mjjM",             500,   0.0, 5000.0}, 
            {    "h_mjjPt",            200,   0.0, 2000.0},
            {    "h_mjjEta",           200,  -6.0,    6.0},
            {    "h_mT",               500,   0.0, 5000.0}, 
            {    "h_METrHT_pt30",      100,   0.0,   20.0},
            {    "h_METrST_pt30",      100,   0.0,    1.0},
        };

        std::vector<TH2DInfo> hist2DInfos = {
            {    "h_jEta_jPhi",     100, -6.0,  6.0, 100,  -3.2,   3.2},
            {    "h_jEtaAK8_jPhiAK8",100, -6.0,  6.0, 100,  -3.2,   3.2},
            {    "h_njets_njetsAK8", 20, 0.0,  20.0, 20,  0.0,   20.0},
        };

        std::vector<TH2DProfileInfo> hist2DProfileInfos = {
        };

        // Initialize Histograms
        if(!initHistos)
        {
            InitHistos(cut_map_1l, histInfos, hist2DInfos, hist2DProfileInfos);
            initHistos = true;
        }

        my_histos["EventCounter"]->Fill(eventCounter);

        for(auto& kv : cut_map_1l)
        {
            if(kv.second)
            {
                double w = weight;
                my_histos["h_njets"               +kv.first]->Fill(NGoodJets_pt30, w);
                my_histos["h_njetsAK8"            +kv.first]->Fill(NGoodJetsAK8_pt200, w);
                my_histos["h_ntops"               +kv.first]->Fill(ntops, w);
                my_histos["h_nb"                  +kv.first]->Fill(NGoodBJets_pt30, w);
                my_histos["h_ht"                  +kv.first]->Fill(HT_trigger_pt30, w);
                my_histos["h_st"                  +kv.first]->Fill(ST_pt30, w);
                my_histos["h_met"                 +kv.first]->Fill(MET, w);
                my_histos["h_weight"              +kv.first]->Fill(weight, w);
                my_histos["h_dEtaJ12"             +kv.first]->Fill(abs(dEtaJ12), w);        
                my_histos["h_dRJ12"               +kv.first]->Fill(dRJ12, w);        
                my_histos["h_dPhiJ1MET"           +kv.first]->Fill(abs(dPhiJ1MET), w);    
                my_histos["h_dPhiJ2MET"           +kv.first]->Fill(abs(dPhiJ2MET), w);    
                my_histos["h_dPhiMinJMET"         +kv.first]->Fill(abs(dPhiMinJMET), w);    
                my_histos["h_mjjM"                +kv.first]->Fill(mjjTLV.M(), w);         
                my_histos["h_mjjPt"               +kv.first]->Fill(mjjTLV.Pt(), w);         
                my_histos["h_mjjEta"              +kv.first]->Fill(mjjTLV.Eta(), w);         
                my_histos["h_mT"                  +kv.first]->Fill(mT, w);           
                my_histos["h_METrHT_pt30"         +kv.first]->Fill(METrHT_pt30, w); 
                my_histos["h_METrST_pt30"         +kv.first]->Fill(METrST_pt30, w); 
                my_2d_histos["h_njets_njetsAK8"   +kv.first]->Fill(NGoodJets_pt30, NGoodJetsAK8_pt200, w);
                for(unsigned int j = 0; j < Jets.size(); j++)
                {
                    if(!GoodJets_pt30[j]) continue;
                    my_histos["h_jPt"+kv.first]->Fill(Jets.at(j).Pt(), w);
                    my_histos["h_jEta"+kv.first]->Fill(Jets.at(j).Eta(), w);
                    my_histos["h_jPhi"+kv.first]->Fill(Jets.at(j).Phi(), w);
                    my_2d_histos["h_jEta_jPhi"+kv.first]->Fill(Jets.at(j).Eta(), Jets.at(j).Phi(), w);
                }
                for(unsigned int j = 0; j < JetsAK8.size(); j++)
                {
                    if(!GoodJetsAK8_pt200[j]) continue;
                    my_histos["h_jPtAK8"+kv.first]->Fill(JetsAK8.at(j).Pt(), w);
                    my_histos["h_jEtaAK8"+kv.first]->Fill(JetsAK8.at(j).Eta(), w);
                    my_histos["h_jPhiAK8"+kv.first]->Fill(JetsAK8.at(j).Phi(), w);
                    my_2d_histos["h_jEtaAK8_jPhiAK8"+kv.first]->Fill(JetsAK8.at(j).Eta(), JetsAK8.at(j).Phi(), w);
                }
            }
        }

        // ------------
        // -- Cut flow
        // ------------
        if(true) my_histos["h_cutFlow"]->Fill(0.5, weight);

    } // end of event loop
}

void Analyze::WriteHistos(TFile* outfile)
{
    outfile->cd();
    
    for(const auto& p : my_histos) 
    {
        p.second->Write();
    }
    
    for(const auto& p : my_2d_histos) 
    {
        p.second->Write();
    }

    for(const auto& p : my_tp_histos) 
    {
        p.second->Write();
    }

    for(const auto& p : my_2d_tp_histos)
    {
        p.second->Write();
    }
    
    for(const auto& p : my_efficiencies) 
    {
        p.second->Write();
    }    
}
