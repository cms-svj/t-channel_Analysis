#define AnalyzeTest_cxx
#include "t-channel_Analysis/EventLoopFramework/interface/AnalyzeTest.h"
#include "t-channel_Analysis/EventLoopFramework/interface/Utility.h"
#include "t-channel_Analysis/EventLoopFramework/interface/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>
#include <TFile.h>

AnalyzeTest::AnalyzeTest()
{
    InitHistos();
}

//Define all your histograms here. 
void AnalyzeTest::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    //Define 1D histograms
    my_histos.emplace( "h_njets", std::make_shared<TH1D>( "h_njets", "h_njets", 8, 7, 15 ) ) ;

    //Define 2D histograms
    my_2d_histos.emplace( "h_njets_HT", std::make_shared<TH2D>( "h_njets_HT", "h_njets_HT", 8, 7, 15, 100, 0, 5000.0 ) );

    //Define TEfficiencies if you are doing trigger studies (for proper error bars) or cut flow charts.
    my_efficiencies.emplace("event_sel_weight", std::make_shared<TEfficiency>("event_sel_weight","event_sel_weight",9,0,9));
}

//Put everything you want to do per event here.
void AnalyzeTest::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        //This is added to count the number of events- do not change the next two lines.
        const auto& eventCounter        = tr.getVar<int>("eventCounter");
        my_histos["EventCounter"]->Fill( eventCounter );

        //--------------------------------------------------
        //-- Print Event Number 
        //--------------------------------------------------
        
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() & (10000 == 0) ) printf( " Event %i\n", tr.getEvtNum() );
        
        const auto& runtype             = tr.getVar<std::string>("runtype");     
        const auto& GoodLeptons         = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");

        const auto& JetID               = tr.getVar<bool>("JetID");
        const auto& NGoodLeptons        = tr.getVar<int>("NGoodLeptons");
        const auto& passTriggerMC       = tr.getVar<bool>("passTriggerMC");
        const auto& NGoodBJets_pt30     = tr.getVar<int>("NGoodBJets_pt30");
        const auto& Mbl                 = tr.getVar<double>("Mbl");
        const auto& HT_trigger_pt30     = tr.getVar<double>("HT_trigger_pt30");
        const auto& NGoodJets_pt30      = tr.getVar<int>("NGoodJets_pt30");
        
        const auto& passMadHT           = tr.getVar<bool>("passMadHT");
        const auto& passBaseline        = tr.getVar<bool>("passBaseline1l_Good");
       
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if ( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

        // ------------------------
        // -- Define weight
        // ------------------------
        double weight               = 1.0;
        double eventweight          = 1.0;
        double leptonScaleFactor    = 1.0;
        double bTagScaleFactor      = 1.0;
        double htDerivedScaleFactor = 1.0;
        double prefiringScaleFactor = 1.0;
        double puScaleFactor        = 1.0;
        
        if(runtype == "MC")
        {
            if( !passMadHT ) continue; //Make sure not to double count DY events
            // Define Lumi weight
            const auto& Weight  = tr.getVar<double>("Weight");
            const auto& lumi = tr.getVar<double>("Lumi");
            eventweight = lumi*Weight;
            
            // Define lepton weight
            if(NGoodLeptons == 1)
            {
                const auto& eleLepWeight = tr.getVar<double>("totGoodElectronSF");
                const auto& muLepWeight  = tr.getVar<double>("totGoodMuonSF");
                leptonScaleFactor = (GoodLeptons[0].first == "e") ? eleLepWeight : muLepWeight;
            }
            
            //PileupWeight = tr.getVar<double>("_PUweightFactor");
            bTagScaleFactor   = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            htDerivedScaleFactor = tr.getVar<double>("htDerivedweight");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor = tr.getVar<double>("puWeightCorr");
            
            weight *= eventweight*leptonScaleFactor*bTagScaleFactor*htDerivedScaleFactor*prefiringScaleFactor*puScaleFactor;
        }
        
        //Make cuts and fill histograms here
        if( passBaseline ) {
            my_histos["h_njets"]->Fill( NGoodJets_pt30, weight ); 
            my_2d_histos["h_njets_HT"]->Fill( NGoodJets_pt30, HT_trigger_pt30, weight );
        }

        // Example Fill event selection efficiencies
        my_efficiencies["event_sel_weight"]->SetUseWeightedEvents();
        my_efficiencies["event_sel_weight"]->FillWeighted(true,eventweight,0);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID,eventweight,1);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID && NGoodLeptons == 1,eventweight,2);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID && NGoodLeptons == 1 && passTriggerMC,eventweight,3);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID && NGoodLeptons == 1 && passTriggerMC && NGoodBJets_pt30 >= 1,eventweight,4);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID && NGoodLeptons == 1 && passTriggerMC && NGoodBJets_pt30 >= 1 && 50 < Mbl && Mbl < 250,eventweight,5);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID && NGoodLeptons == 1 && passTriggerMC && NGoodBJets_pt30 >= 1 && 50 < Mbl && Mbl < 250 && HT_trigger_pt30 > 300,eventweight,6);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID && NGoodLeptons == 1 && passTriggerMC && NGoodBJets_pt30 >= 1 && 50 < Mbl && Mbl < 250 && HT_trigger_pt30 > 300 && NGoodJets_pt30 >= 7,eventweight,7);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID && NGoodLeptons == 1 && passTriggerMC && NGoodBJets_pt30 >= 1 && 50 < Mbl && Mbl < 250 && HT_trigger_pt30 > 300 && NGoodJets_pt30 >= 7,weight,8);
    } 
}

void AnalyzeTest::WriteHistos(TFile* outfile)
{
    outfile->cd();

    for (const auto &p : my_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_2d_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_efficiencies) {
        p.second->Write();
    }    
}
