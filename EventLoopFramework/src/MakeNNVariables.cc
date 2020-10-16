#define MakeNNVariables_cxx
#include "t-channel_Analysis/EventLoopFramework/interface/MakeNNVariables.h"
#include "t-channel_Analysis/EventLoopFramework/interface/NTupleReader.h"
#include "t-channel_Analysis/EventLoopFramework/interface/MiniTupleMaker.h"
#include "t-channel_Analysis/EventLoopFramework/interface/Utility.h" 

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>
#include <stdio.h> 
//#include <fstream>
//#include <cstdio>

MakeNNVariables::MakeNNVariables()
{
    InitHistos();
}

void MakeNNVariables::InitHistos()
{
    my_histos.emplace( "EventCounterTrain", std::make_shared<TH1D>( "EventCounterTrain", "EventCounterTrain", 2, -1.1, 1.1 ) ); 
    my_histos.emplace( "EventCounterTest",  std::make_shared<TH1D>( "EventCounterTest",  "EventCounterTest",  2, -1.1, 1.1 ) ); 
    my_histos.emplace( "EventCounterVal",   std::make_shared<TH1D>( "EventCounterVal",   "EventCounterVal",   2, -1.1, 1.1 ) ); 
}//END of init histos

void MakeNNVariables::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    int count=0, numPassTrain=0, numPassTest=0, numPassVal=0;
    while( tr.getNextEvent() )
    {
        const auto& eventCounter        = tr.getVar<int>("eventCounter");
        //const auto& runtype             = tr.getVar<std::string>("runtype");
        const auto& isSignal            = tr.getVar<bool>("isSignal");
        const auto& passBaseline1l      = tr.getVar<bool>("passBaseline1l_Good");
        const auto& filetag             = tr.getVar<std::string>("filetag");

        auto& mass = tr.createDerivedVar<double>("mass", 0.0);
        if(!isSignal)
        {
            mass = 173.0;
        }
        else
        {
            for(unsigned int m = 300; m < 1500; m+=50)
            {
                mass = (filetag.find(std::to_string(m)) != std::string::npos) ? m : mass;
            }
        }
       
        //------------------------------------
        //-- Print Event Number
        //------------------------------------
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() % 1000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );        

        //-----------------------------------
        //  Initialize the tree
        //-----------------------------------       
        std::set<std::string> variables = {
            "Lumi",
            "isSignal",
            "deepESM_val",
            //"filetag",
            "NGoodJets_pt30_double",
            "NGoodBJets_pt30",
            "Mbl",
            "HT_trigger_pt30",
            //"HT",
            "NGoodElectrons",
            "NGoodMuons",
            //"NVtx",
            "Weight",
            "totalEventWeight",
            "fwm2_top6", 
            "fwm3_top6", 
            "fwm4_top6", 
            "fwm5_top6", 
            "fwm6_top6", 
            "fwm7_top6", 
            "fwm8_top6", 
            "fwm9_top6", 
            "fwm10_top6", 
            "jmt_ev0_top6", 
            "jmt_ev1_top6", 
            "jmt_ev2_top6",
            "Jet_pt_1",   "Jet_pt_2",   "Jet_pt_3",   "Jet_pt_4",   "Jet_pt_5",   "Jet_pt_6",   "Jet_pt_7",   "Jet_pt_8",   "Jet_pt_9",   "Jet_pt_10",   "Jet_pt_11",   "Jet_pt_12",  
            "Jet_eta_1",  "Jet_eta_2",  "Jet_eta_3",  "Jet_eta_4",  "Jet_eta_5",  "Jet_eta_6",  "Jet_eta_7",  "Jet_eta_8",  "Jet_eta_9",  "Jet_eta_10",  "Jet_eta_11",  "Jet_eta_12",
            "Jet_phi_1",  "Jet_phi_2",  "Jet_phi_3",  "Jet_phi_4",  "Jet_phi_5",  "Jet_phi_6",  "Jet_phi_7",  "Jet_phi_8",  "Jet_phi_9",  "Jet_phi_10",  "Jet_phi_11",  "Jet_phi_12",
            "Jet_m_1",    "Jet_m_2",    "Jet_m_3",    "Jet_m_4",    "Jet_m_5",    "Jet_m_6",    "Jet_m_7",    "Jet_m_8",    "Jet_m_9",    "Jet_m_10",    "Jet_m_11",    "Jet_m_12",
            "Jet_dcsv_1", "Jet_dcsv_2", "Jet_dcsv_3", "Jet_dcsv_4", "Jet_dcsv_5", "Jet_dcsv_6", "Jet_dcsv_7", "Jet_dcsv_8", "Jet_dcsv_9", "Jet_dcsv_10", "Jet_dcsv_11", "Jet_dcsv_12",
            "Jet_ptD_1",  "Jet_ptD_2",  "Jet_ptD_3",  "Jet_ptD_4",  "Jet_ptD_5",  "Jet_ptD_6",  "Jet_ptD_7",  "Jet_ptD_8",  "Jet_ptD_9",  "Jet_ptD_10",  "Jet_ptD_11",  "Jet_ptD_12",
            "Jet_axismajor_1", "Jet_axismajor_2", "Jet_axismajor_3", "Jet_axismajor_4", "Jet_axismajor_5", "Jet_axismajor_6", "Jet_axismajor_7", "Jet_axismajor_8", "Jet_axismajor_9", "Jet_axismajor_10", "Jet_axismajor_11", "Jet_axismajor_12",
            "Jet_axisminor_1", "Jet_axisminor_2", "Jet_axisminor_3", "Jet_axisminor_4", "Jet_axisminor_5", "Jet_axisminor_6", "Jet_axisminor_7", "Jet_axisminor_8", "Jet_axisminor_9", "Jet_axisminor_10", "Jet_axisminor_11", "Jet_axisminor_12",
            "Jet_multiplicity_1", "Jet_multiplicity_2", "Jet_multiplicity_3", "Jet_multiplicity_4", "Jet_multiplicity_5", "Jet_multiplicity_6", "Jet_multiplicity_7", "Jet_multiplicity_8", "Jet_multiplicity_9", "Jet_multiplicity_10", "Jet_multiplicity_11", "Jet_multiplicity_12",
            "GoodLeptons_pt_1",  "GoodLeptons_pt_2",
            "GoodLeptons_eta_1", "GoodLeptons_eta_2",
            "GoodLeptons_phi_1", "GoodLeptons_phi_2",
            "GoodLeptons_m_1",   "GoodLeptons_m_2",
            "JetsAK8Cands_pt_1",     "JetsAK8Cands_pt_2",     "JetsAK8Cands_pt_3",     "JetsAK8Cands_pt_4",     "JetsAK8Cands_pt_5",
            "JetsAK8Cands_eta_1",    "JetsAK8Cands_eta_2",    "JetsAK8Cands_eta_3",    "JetsAK8Cands_eta_4",    "JetsAK8Cands_eta_5",
            "JetsAK8Cands_phi_1",    "JetsAK8Cands_phi_2",    "JetsAK8Cands_phi_3",    "JetsAK8Cands_phi_4",    "JetsAK8Cands_phi_5",
            "JetsAK8Cands_m_1",      "JetsAK8Cands_m_2",      "JetsAK8Cands_m_3",      "JetsAK8Cands_m_4",      "JetsAK8Cands_m_5",
            "JetsAK8Cands_SDM_1",    "JetsAK8Cands_SDM_2",    "JetsAK8Cands_SDM_3",    "JetsAK8Cands_SDM_4",    "JetsAK8Cands_SDM_5",
            "JetsAK8Cands_Pruned_1", "JetsAK8Cands_Pruned_2", "JetsAK8Cands_Pruned_3", "JetsAK8Cands_Pruned_4", "JetsAK8Cands_Pruned_5",
            "JetsAK8Cands_T21_1",    "JetsAK8Cands_T21_2",    "JetsAK8Cands_T21_3",    "JetsAK8Cands_T21_4",    "JetsAK8Cands_T21_5",
            "lvMET_cm_pt",
            "lvMET_cm_eta",
            "lvMET_cm_phi",
            "lvMET_cm_m",
            "stop1_PtRank_1l_mass", "stop2_PtRank_1l_mass",
            "mass",
            "deepESM_valReg",
        };

        if( tr.isFirstEvent() ) 
        {
            std::string myTreeName = "myMiniTree";

            my_histos["EventCounterTrain"]->Fill( eventCounter );
            myTreeTrain = new TTree( (myTreeName).c_str() , (myTreeName).c_str() );
            myMiniTupleTrain = new MiniTupleMaker( myTreeTrain );
            myMiniTupleTrain->setTupleVars(variables);
            myMiniTupleTrain->initBranches(tr);

            my_histos["EventCounterTest"]->Fill( eventCounter );
            myTreeTest = new TTree( (myTreeName).c_str() , (myTreeName).c_str() );
            myMiniTupleTest = new MiniTupleMaker( myTreeTest );
            myMiniTupleTest->setTupleVars(variables);
            myMiniTupleTest->initBranches(tr);

            my_histos["EventCounterVal"]->Fill( eventCounter );
            myTreeVal = new TTree( (myTreeName).c_str() , (myTreeName).c_str() );
            myMiniTupleVal = new MiniTupleMaker( myTreeVal );
            myMiniTupleVal->setTupleVars(variables);
            myMiniTupleVal->initBranches(tr);
        }
        
        //-----------------------------------
        //-- Fill Histograms Below
        //-----------------------------------
        if( passBaseline1l ) 
        {
            int mod = count % 10;
            if(mod < 8)
            {
                myMiniTupleTrain->fill();
                numPassTrain++;
            }
            else if(mod == 8)
            {
                myMiniTupleTest->fill();
                numPassTest++;
            }
            else
            {
                myMiniTupleVal->fill();
                numPassVal++;
            }
            count++;
        }
    }//END of while tr.getNextEvent loop   
    std::cout<<"Total: "<<count<<" Train: "<<numPassTrain<<" Test: "<<numPassTest<<" Val: "<<numPassVal<<std::endl;
}//END of function
      
void MakeNNVariables::WriteHistos( TFile* outfile ) 
{
    const auto& outFileName = std::string(outfile->GetName());
    const auto& name = utility::split("first", outFileName, ".");

    TFile* outfileTrain = TFile::Open((name+"_Train.root").c_str(), "RECREATE");
    outfileTrain->cd();
    myTreeTrain->Write();
    my_histos["EventCounterTrain"]->Write();
    delete myTreeTrain;    
    delete myMiniTupleTrain;
    outfileTrain->Close();

    TFile* outfileTest = TFile::Open((name+"_Test.root").c_str(), "RECREATE");
    outfileTest->cd();
    myTreeTest->Write();
    my_histos["EventCounterTest"]->Write();
    delete myTreeTest;    
    delete myMiniTupleTest;
    outfileTest->Close();

    TFile* outfileVal = TFile::Open((name+"_Val.root").c_str(), "RECREATE");
    outfileVal->cd();
    myTreeVal->Write();
    my_histos["EventCounterVal"]->Write();
    delete myTreeVal;    
    delete myMiniTupleVal;
    outfileVal->Close();

    remove(outFileName.c_str());
}

