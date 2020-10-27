#ifndef Confg_h
#define Confg_h

#include "t-channel_Analysis/EventLoopFramework/interface/NTupleReader.h"
#include "t-channel_Analysis/EventLoopFramework/interface/PrepNTupleVars.h"
#include "t-channel_Analysis/EventLoopFramework/interface/RunTopTagger.h"
#include "t-channel_Analysis/EventLoopFramework/interface/Muon.h"
#include "t-channel_Analysis/EventLoopFramework/interface/Electron.h"
#include "t-channel_Analysis/EventLoopFramework/interface/Photon.h"
#include "t-channel_Analysis/EventLoopFramework/interface/Jet.h"
#include "t-channel_Analysis/EventLoopFramework/interface/BJet.h"
#include "t-channel_Analysis/EventLoopFramework/interface/CommonVariables.h"
#include "t-channel_Analysis/EventLoopFramework/interface/Baseline.h"
#include "t-channel_Analysis/EventLoopFramework/interface/BTagCorrector.h"
#include "t-channel_Analysis/EventLoopFramework/interface/ScaleFactors.h"
#include "t-channel_Analysis/EventLoopFramework/interface/PartialUnBlinding.h"
#include "t-channel_Analysis/EventLoopFramework/interface/FatJetCombine.h"

class Config
{
private:
    void registerModules(NTupleReader& tr, const std::vector<std::string>&& modules) const
    {
        const auto& runtype               = tr.getVar<std::string>("runtype");
        const auto& runYear               = tr.getVar<std::string>("runYear");
        const auto& leptonFileName        = tr.getVar<std::string>("leptonFileName");
        const auto& puFileName            = tr.getVar<std::string>("puFileName");
        const auto& meanFileName          = tr.getVar<std::string>("meanFileName");
        const auto& bjetFileName          = tr.getVar<std::string>("bjetFileName");
        const auto& bjetCSVFileName       = tr.getVar<std::string>("bjetCSVFileName");
        const auto& filetag               = tr.getVar<std::string>("filetag");
        const auto& TopTaggerCfg          = tr.getVar<std::string>("TopTaggerCfg");       
 
        for(const auto& module : modules)
        {
            if     (module=="PartialUnBlinding")            tr.emplaceModule<PartialUnBlinding>();
            else if(module=="PrepNTupleVars")               tr.emplaceModule<PrepNTupleVars>();
            else if(module=="RunTopTagger")                 tr.emplaceModule<RunTopTagger>(TopTaggerCfg);
            else if(module=="Muon")                         tr.emplaceModule<Muon>();
            else if(module=="Electron")                     tr.emplaceModule<Electron>();
            else if(module=="Photon")                       tr.emplaceModule<Photon>();
            else if(module=="Jet")                          tr.emplaceModule<Jet>();
            else if(module=="BJet")                         tr.emplaceModule<BJet>();
            else if(module=="CommonVariables")              tr.emplaceModule<CommonVariables>();
            else if(module=="Baseline")                     tr.emplaceModule<Baseline>();
            else if(module=="FatJetCombine")                tr.emplaceModule<FatJetCombine>();
            
            if(runtype == "MC")
            {
                if(module=="ScaleFactors")  tr.emplaceModule<ScaleFactors>(runYear, leptonFileName, puFileName, meanFileName);
                else if(module=="BTagCorrector")
                {
                    auto& bTagCorrector = tr.emplaceModule<BTagCorrector>(bjetFileName, "", bjetCSVFileName, filetag);
                    bTagCorrector.SetVarNames("GenParticles_PdgId", "Jets", "GoodJets_pt30", "Jets_bJetTagDeepCSVtotb", "Jets_partonFlavor");
                }
            }
        }
    }

public:
    Config() 
    {
    }

    void setUp(NTupleReader& tr) const
    {
        //Get and make needed info
        const auto& filetag = tr.getVar<std::string>("filetag");
        const auto& analyzer = tr.getVar<std::string>("analyzer");
        const bool isSignal = (filetag.find("_stop") != std::string::npos || filetag.find("_mStop") != std::string::npos || filetag.find("VLQ_2t4b") != std::string::npos || filetag.find("mMed") != std::string::npos || filetag.find("mZprime") != std::string::npos) ? true : false;

        std::string runYear, puFileName, leptonFileName, bjetFileName, bjetCSVFileName, meanFileName, TopTaggerCfg;
        
        double Lumi=0.0, deepCSV_WP_loose=0.0, deepCSV_WP_medium=0.0, deepCSV_WP_tight=0.0;
        bool blind = true;

        if(filetag.find("2016") != std::string::npos)
        {
            runYear = "2016";
            Lumi = 35900.0;
            deepCSV_WP_loose  = 0.2217;
            deepCSV_WP_medium = 0.6321;
            deepCSV_WP_tight  = 0.8953;            
            puFileName = "PileupHistograms_0121_69p2mb_pm4p6.root";
            leptonFileName = "allInOne_leptonSF_2016.root";
            bjetFileName = "allInOne_BTagEff.root";
            bjetCSVFileName = "DeepCSV_2016LegacySF_WP_V1.csv";
            meanFileName = "allInOne_SFMean.root";
            blind = false;
            TopTaggerCfg = "TopTaggerCfg_2016.cfg";
        }
        else if(filetag.find("2017") != std::string::npos)
        { 
            runYear = "2017";
            Lumi = 41525.0;
            deepCSV_WP_loose  = 0.1522;
            deepCSV_WP_medium = 0.4941;       
            deepCSV_WP_tight  = 0.8001;
            puFileName = "pu_ratio.root";
            leptonFileName = "allInOne_leptonSF_2017.root";
            bjetFileName = "allInOne_BTagEff.root";
            bjetCSVFileName = "DeepCSV_94XSF_WP_V4_B_F.csv";
            meanFileName = "allInOne_SFMean.root";
            blind = false;
            TopTaggerCfg = "TopTaggerCfg_2017.cfg";
        }
        else if(filetag.find("2018") != std::string::npos) 
        {
            runYear = "2018";
            Lumi = 21071.0+38654.0;
            deepCSV_WP_loose  = 0.1241;
            deepCSV_WP_medium = 0.4184;       
            deepCSV_WP_tight  = 0.7527;
            puFileName = "PileupHistograms_2018_69mb_pm5.root";
            leptonFileName = "allInOne_leptonSF_2018.root";
            bjetFileName = "allInOne_BTagEff.root";
            bjetCSVFileName = "DeepCSV_102XSF_WP_V1.csv";
            meanFileName = "allInOne_SFMean.root";
            blind = false;
            TopTaggerCfg = "TopTaggerCfg_2018.cfg";
        }

        tr.registerDerivedVar("runYear",runYear);
        tr.registerDerivedVar("Lumi",Lumi);
        tr.registerDerivedVar("deepCSV_WP_loose",deepCSV_WP_loose);
        tr.registerDerivedVar("deepCSV_WP_medium",deepCSV_WP_medium);
        tr.registerDerivedVar("deepCSV_WP_tight",deepCSV_WP_tight);
        tr.registerDerivedVar("isSignal",isSignal);
        tr.registerDerivedVar("puFileName",puFileName);
        tr.registerDerivedVar("leptonFileName",leptonFileName);        
        tr.registerDerivedVar("bjetFileName",bjetFileName);        
        tr.registerDerivedVar("bjetCSVFileName",bjetCSVFileName);        
        tr.registerDerivedVar("meanFileName",meanFileName);        
        tr.registerDerivedVar("etaCut",2.4); 
        tr.registerDerivedVar("blind",blind);
        tr.registerDerivedVar("TopTaggerCfg", TopTaggerCfg);

        //Register Modules that are needed for each Analyzer
        if(analyzer=="MakeNJetDists")
        {
            const std::vector<std::string> modulesList = {
                "PartialUnBlinding",
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "CommonVariables",
                "FatJetCombine",
                "Baseline",
                "BTagCorrector",
                "ScaleFactors"
            };
            registerModules(tr, std::move(modulesList));
        }
        else if (analyzer=="AnalyzeTest")
        {
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "PartialUnBlinding",
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "CommonVariables",
                "Baseline",
                "BTagCorrector",
                "ScaleFactors"
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="MakeNNVariables")
        {
            const std::vector<std::string> modulesList = {
                "PartialUnBlinding",
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "RunTopTagger",
                "CommonVariables",
                "Baseline",
                "FatJetCombine",
                "BTagCorrector",
                "ScaleFactors",
            };
            registerModules(tr, std::move(modulesList));
        }
        else
        {
            const std::vector<std::string> modulesList = {
                "PartialUnBlinding",
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "CommonVariables",
                "FatJetCombine",
                //"MakeMVAVariables",
                "Baseline",
                //"DeepEventShape",
                //"BTagCorrector",
                //"ScaleFactors",
                "RunTopTagger",
            };
            registerModules(tr, std::move(modulesList));
        }
    }
};

#endif
