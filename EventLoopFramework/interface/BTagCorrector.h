#ifndef BTAGCORRECTOR_H
#define BTAGCORRECTOR_H

//custom headers
#include "t-channel_Analysis/EventLoopFramework/interface/BTagCalibrationStandalone.h"
#include "t-channel_Analysis/EventLoopFramework/interface/NTupleReader.h"
#include "t-channel_Analysis/EventLoopFramework/interface/SATException.h"

//ROOT headers
#include <TROOT.h>
#include "TMath.h"
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include "TH2.h"

//STL headers
#include <string>
#include <vector>
#include <cmath>

class BTagCorrector
{
public:
    //constructor
    BTagCorrector(std::string file, std::string CSVFilePath, std::string CSVFile, std::string suffix = "") 
    : debug(false), btagSFunc(0), mistagSFunc(0), btagCFunc(0), ctagCFunc(0), mistagCFunc(0)
    , h_eff_b(nullptr), h_eff_c(nullptr), h_eff_udsg(nullptr), inFileName(file), MCBranch(""), JetsVec(""), JetMask(""), BJetsVec(""), JetsFlavor(""), myVarSuffix_("")
    {
        std::cout<<"Setting up BTagCorrector"<<std::endl;
        
        //Stops unwanted segfaults.
        TH1::AddDirectory(false);
        
        TFile inFile(inFileName.c_str());
        SetEffs(inFile, suffix);
        inFile.Close();
        
        if(CSVFilePath.size())
        {
            SetCalib((CSVFilePath + "/" + CSVFile).c_str());
        }
        else
        {
            SetCalib(CSVFile.c_str());
        }
    }

    //accessors
    void SetDebug(bool d) { debug = d; }

    void SetEffs(TFile& file, std::string suffix = "")
    {
        if(suffix.size())
        {
            std::string suffix2 = suffix;
            h_eff_b.reset( (TH2F*)file.Get(("n_eff_b_" + suffix2).c_str()) );
            
            if(!h_eff_b.get())
            {
                THROW_SATEXCEPTION("\033[1;31mError: Could not find \"n_eff_b_"+suffix2+"\" histogram in the Btag scale factor root file\033[0m");
            }
            h_eff_c.reset( (TH2F*)file.Get(("n_eff_c_" + suffix2).c_str()) );
            h_eff_udsg.reset( (TH2F*)file.Get(("n_eff_udsg_" + suffix2).c_str()) );
            std::unique_ptr<TH2F> d_eff_b( (TH2F*)file.Get(("d_eff_b_" + suffix2).c_str()) );
            std::unique_ptr<TH2F> d_eff_c( (TH2F*)file.Get(("d_eff_c_" + suffix2).c_str()) );
            std::unique_ptr<TH2F> d_eff_udsg( (TH2F*)file.Get(("d_eff_udsg_" + suffix2).c_str()) );
            
            h_eff_b->Divide(d_eff_b.get());
            h_eff_c->Divide(d_eff_c.get());
            h_eff_udsg->Divide(d_eff_udsg.get());            
        }
        else
        {
            h_eff_b.reset( (TH2F*)file.Get("n_eff_b") );
            h_eff_c.reset( (TH2F*)file.Get("n_eff_c") );
            h_eff_udsg.reset(  (TH2F*)file.Get("n_eff_udsg") );
	    std::unique_ptr<TH2F> d_eff_b( (TH2F*)file.Get("d_eff_b") );
            std::unique_ptr<TH2F> d_eff_c( (TH2F*)file.Get("d_eff_c") );
            std::unique_ptr<TH2F> d_eff_udsg( (TH2F*)file.Get("d_eff_udsg") );
            
            if(h_eff_b.get())
            {
                h_eff_b->Divide(d_eff_b.get());
                h_eff_c->Divide(d_eff_c.get());
                h_eff_udsg->Divide(d_eff_udsg.get());
            }
        }
    }

    void resetEffs(std::string suffix)
    {
        TFile inFile(inFileName.c_str());
        SetEffs(inFile, suffix);
        inFile.Close();
    }

    void SetCalib(std::string cfile)
    {        
        //initialize btag helper classes. Interface has been changed.
        calib = BTagCalibration("",cfile);
        reader = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central");
        reader.load(calib, BTagEntry::FLAV_B, "comb"); reader.load(calib, BTagEntry::FLAV_C, "comb");  reader.load(calib, BTagEntry::FLAV_UDSG, "incl");
        readerUp = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "up");
        readerUp.load(calib, BTagEntry::FLAV_B, "comb"); readerUp.load(calib, BTagEntry::FLAV_C, "comb");  readerUp.load(calib, BTagEntry::FLAV_UDSG, "incl");
        readerDown = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "down");
        readerDown.load(calib, BTagEntry::FLAV_B, "comb"); readerDown.load(calib, BTagEntry::FLAV_C, "comb");  readerDown.load(calib, BTagEntry::FLAV_UDSG, "incl");        
    }

    void SetVarNames(std::string MCBranchName, std::string JetsVecName, std::string JetMaskName, std::string BJetsVecName, std::string JetsFlavorName, std::string myVarSuffix = "")
    {
        MCBranch = MCBranchName;
        JetsVec = JetsVecName;
        JetMask = JetMaskName;
        BJetsVec = BJetsVecName;
        JetsFlavor = JetsFlavorName;
        myVarSuffix_ = myVarSuffix;
    }
    
    void SetBtagSFunc(int u) { btagSFunc = u; }
    //void SetCtagSFunc(int u) { ctagSFunc = u; }
    void SetCtagSFunc(int u) { btagSFunc = u; } //ctag and btag are correlated
    void SetMistagSFunc(int u) { mistagSFunc = u; }
    void SetBtagCFunc(int u) { btagCFunc = u; }
    void SetCtagCFunc(int u) { ctagCFunc = u; }
    void SetMistagCFunc(int u) { mistagCFunc = u; }

    /***********************************************************************************/
    // Method 1b) in twiki
    // https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
    /***********************************************************************************/
    std::vector<double>* GetCorrections(const std::vector<TLorentzVector>* Jets, const std::vector<bool>* jetMask,  const std::vector<int>* Jets_flavor)
    {
        //reset probabilities
        std::vector<double>* prob = new std::vector<double>(4,0.0);
        (*prob)[0] = 1.0;
        
        //first loop over jets
        std::vector<std::vector<double> > sfEffLists = std::vector<std::vector<double> >(Jets->size(),std::vector<double>());
        for(unsigned ja = 0; ja < Jets->size(); ++ja)
        {
            if(!jetMask->at(ja)) continue;
        
            //get sf and eff values (checks if already calculated)
            InitSFEff(Jets->at(ja).Pt(), Jets->at(ja).Eta(), Jets_flavor->at(ja), sfEffLists[ja]);
            double eps_a = sfEffLists[ja][0]*sfEffLists[ja][1]*sfEffLists[ja][2];
        
            //jet index, pt, eta, flavor, eff, sf, cf
            if(debug) std::cout << "Jet " << ja << ": " << Jets->at(ja).Pt() << ", " << fabs(Jets->at(ja).Eta()) << ", " << abs(Jets_flavor->at(ja)) 
                                << " sfEffLists[ja][0], " << sfEffLists[ja][0] << "  sfEffLists[ja][1], " << sfEffLists[ja][1] << "  sfEffLists[ja][2], " << sfEffLists[ja][2] << std::endl;
        
            //calculate prob(0 b-tags)
            (*prob)[0] *= (1-eps_a);
	 
            //sub-probabilities for following calculations
            double subprob1 = 1.0;
            double subprob2 = 0.0;
	 
            //second loop over jets
            for(unsigned jb = 0; jb < Jets->size(); ++jb)
            {
                //skip the same jet
                if(jb==ja || !jetMask->at(ja)) continue;	   
	   	   
                //get sf and eff values (checks if already calculated)
                InitSFEff(Jets->at(jb).Pt(), Jets->at(jb).Eta(), Jets_flavor->at(jb), sfEffLists[jb]);
	   
                double eps_b = sfEffLists[jb][0]*sfEffLists[jb][1]*sfEffLists[jb][2];
	   
                //jet index, pt, eta, flavor, eff, sf, cf
                if(debug) std::cout << "\tJet " << jb << ": " << Jets->at(jb).Pt() << ", " << fabs(Jets->at(jb).Eta()) << ", " << abs(Jets_flavor->at(jb)) 
                                    << ", " << sfEffLists[jb][0] << ", " << sfEffLists[jb][1] << ", " << sfEffLists[jb][2] << std::endl;
	   
                //calculate prob(1 b-tag)
                subprob1 *= (1-eps_b);
	   
                //sub-sub-probability for following calculations
                double subsubprob2 = 1.0;
	   
                //third loop over jets (only for jb>ja)
                if(jb<ja) continue;
                for(unsigned jc = 0; jc < Jets->size(); ++jc)
                {
                    //skip the same jet
                    if(jc==jb || jc==ja || !jetMask->at(ja)) continue;
	     
                    //get sf and eff values (checks if already calculated)
                    InitSFEff(Jets->at(jc).Pt(), Jets->at(jc).Eta(), Jets_flavor->at(jc), sfEffLists[jc]);
                    double eps_c = sfEffLists[jc][0]*sfEffLists[jc][1]*sfEffLists[jc][2];
	     
                    //jet index, pt, eta, flavor, eff, sf, cf		
                    if(debug) std::cout << "\t\tJet " << jc << ": " << Jets->at(jc).Pt() << ", " << fabs(Jets->at(jc).Eta()) << ", " << abs(Jets_flavor->at(jc)) 
                                        << ", " << sfEffLists[jc][0] << ", " << sfEffLists[jc][1] << ", " << sfEffLists[jc][2] << std::endl;
	     
                    //calculate prob(2 b-tags)
                    subsubprob2 *= (1-eps_c);
                }
	   
                //add up sub-sub-prob
                subprob2 += eps_b*subsubprob2;
            }
	 
            //add up sub-probs
            (*prob)[1] += eps_a*subprob1;
            (*prob)[2] += eps_a*subprob2;
        }
        
        //conserve probability       
        (*prob)[3] = 1 - (*prob)[0] -(*prob)[1] - (*prob)[2];
        if(debug) std::cout << (*prob)[0] << ", " << (*prob)[1] << ", " << (*prob)[2] << ", " << (*prob)[3] << std::endl;
        
        return prob;
    }

    /***********************************************************************************/
    //method 1a in twiki
    // https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods  
    /***********************************************************************************/
    double GetSimpleCorrection(const std::vector<TLorentzVector>* Jets, const std::vector<bool>* jetMask, const std::vector<int>* Jets_flavor, const std::vector<double>* Jets_bDiscriminatorCSV, const double wp)
    {
        double mcTag = 1.0, mcNoTag = 1.0, dataTag = 1.0, dataNoTag = 1.0;
        
        //loop over jets
        std::vector<std::vector<double> > sfEffLists = std::vector<std::vector<double> >(Jets->size(),std::vector<double>());
        for(unsigned ja = 0; ja < Jets->size(); ++ja)
        {
            if(!jetMask->at(ja)) continue;
            
            //get sf and eff values (checks if already calculated)
            InitSFEff(Jets->at(ja).Pt(), Jets->at(ja).Eta(), Jets_flavor->at(ja), sfEffLists[ja]);
            double eff_a = sfEffLists[ja][0]; //eff
            double cf_a = sfEffLists[ja][2]; //CF
            double sf_a = sfEffLists[ja][1];
            
            if( sfEffLists[ja][0] == 0.0 || sfEffLists[ja][1] == 0.0 || sfEffLists[ja][2] == 0.0 )
            {
                if(debug) std::cout<<"sfEffLists[ja][0] : "<<sfEffLists[ja][0]<<"  sfEffLists[ja][1] : "<<sfEffLists[ja][1]<<"  sfEffLists[ja][2] : "<<sfEffLists[ja][2]<<std::endl;
            }
            
            if(Jets_bDiscriminatorCSV->at(ja) > wp)
            {
                mcTag *= eff_a*cf_a;
                dataTag *= eff_a*cf_a*sf_a;
            } 
            else 
            {
                mcNoTag *= (1-eff_a*cf_a);
                dataNoTag *= (1-eff_a*cf_a*sf_a);
            }
        }
        
        double result = (mcNoTag * mcTag ==0) ? 1.0 : (dataNoTag * dataTag)/(mcNoTag * mcTag);
        return result;
    }

    void InitSFEff(double pt, double eta, int flav, std::vector<double>& sfEffList)
    {
        //avoid rerunning this
        sfEffList.clear();
        if(sfEffList.size()>0) return;
        
        //use abs(flav) always
        flav = abs(flav);
        
        sfEffList = std::vector<double>(3,1.0); //eff, sf (central, up, or down), cf (central, up, or down)
        
        if(flav==5)
        { //b-tag
            // data_t Uncertainty are now taken care automaticall with method eval_auto_bounds
            //in new interface.
            int pt_bin = h_eff_b->GetXaxis()->FindBin(pt); 
            if( pt_bin > h_eff_b->GetXaxis()->GetNbins() ) pt_bin = h_eff_b->GetXaxis()->GetNbins(); 
            int eta_bin = h_eff_b->GetYaxis()->FindBin(eta); 
            if ( eta_bin > h_eff_b->GetYaxis()->GetNbins() ) eta_bin = h_eff_b->GetYaxis()->GetNbins();
        
            sfEffList[0] = h_eff_b->GetBinContent(pt_bin, eta_bin);
        
            sfEffList[1] = (btagSFunc==0 ? reader.eval_auto_bounds("central",BTagEntry::FLAV_B,eta,pt) :
                            (btagSFunc==1 ? readerUp.eval_auto_bounds("up",BTagEntry::FLAV_B,eta,pt) :
                             readerDown.eval_auto_bounds("down",BTagEntry::FLAV_B,eta,pt) ) );       
        }
        else if(flav==4)
        { //charm mistag
            int pt_bin = h_eff_c->GetXaxis()->FindBin(pt); 
            if( pt_bin > h_eff_c->GetXaxis()->GetNbins() ) pt_bin = h_eff_c->GetXaxis()->GetNbins();
            int eta_bin = h_eff_c->GetYaxis()->FindBin(eta); 
            if ( eta_bin > h_eff_c->GetYaxis()->GetNbins() ) eta_bin = h_eff_c->GetYaxis()->GetNbins();
            sfEffList[0] =h_eff_c->GetBinContent(pt_bin, eta_bin);
        
            sfEffList[1] = (btagSFunc==0 ? reader.eval_auto_bounds("central",BTagEntry::FLAV_C,eta,pt) :
                            (btagSFunc==1 ? readerUp.eval_auto_bounds("up",BTagEntry::FLAV_C,eta,pt) :
                             readerDown.eval_auto_bounds("down", BTagEntry::FLAV_C,eta,pt) ) );
        }
        else if(flav<4 || flav==21)
        { //udsg mistag
            int pt_bin = h_eff_udsg->GetXaxis()->FindBin(pt); 
            if( pt_bin > h_eff_udsg->GetXaxis()->GetNbins() ) pt_bin = h_eff_udsg->GetXaxis()->GetNbins(); 
            int eta_bin = h_eff_udsg->GetYaxis()->FindBin(eta); 
            if ( eta_bin > h_eff_udsg->GetYaxis()->GetNbins() ) eta_bin = h_eff_udsg->GetYaxis()->GetNbins();
        
            sfEffList[0] = h_eff_udsg->GetBinContent( pt_bin, eta_bin);
        
            sfEffList[1] = (mistagSFunc==0 ? reader.eval_auto_bounds("central",BTagEntry::FLAV_UDSG,eta,pt) :
                            (mistagSFunc==1 ? readerUp.eval_auto_bounds("up",BTagEntry::FLAV_UDSG,eta,pt) :
                             readerDown.eval_auto_bounds("down", BTagEntry::FLAV_UDSG,eta,pt) ) );           
        }  
    }

    // To register Event weights/ Probabilities to FlatTuples
    void registerVarToNTuples(NTupleReader& tr)
    {
        //Check if this is data
        if( !tr.checkBranch(MCBranch) ) return;
        const auto& inputJets = tr.getVec<TLorentzVector>(JetsVec);
        const auto& jetMask = tr.getVec<bool>(JetMask);
        const auto& recoJetsBtag = tr.getVec<double>(BJetsVec);
        const auto& recoJetsFlavor = tr.getVec<int>(JetsFlavor);
        const auto& wp = tr.getVar<double>("deepCSV_WP_medium");
        
        /*************************************************/
        // Here we define which(up, down or central
        // to be stored central = 0, up = 1 down = else
        /*************************************************/
        
        /*************************************************/
        // Case 0: Central value;
        /*************************************************/        
        int switch_Unc = 0, switch_udsg_Unc = 0;
        SetBtagSFunc(switch_Unc); SetBtagCFunc(switch_Unc);
        SetCtagSFunc(switch_Unc); SetCtagCFunc(switch_Unc);
        SetMistagSFunc(switch_udsg_Unc); SetMistagCFunc(switch_udsg_Unc);
        //Method 1a) ignoring b-tag status 
        double evtWeightSimple_Central  = GetSimpleCorrection(&inputJets,&jetMask,&recoJetsFlavor,&recoJetsBtag,wp);
        if( std::isnan( evtWeightSimple_Central) || std::isinf(evtWeightSimple_Central) )
        {
            evtWeightSimple_Central = 1.0;
        } 
        //Register derived quantities to nTuples.
        tr.registerDerivedVar("bTagSF_EventWeightSimple_Central"+myVarSuffix_, evtWeightSimple_Central);        

        //// Method 1b) in different b-jet mullticipity bins.
        //std::vector<double> *evtWeightProb_Central = GetCorrections(&inputJets,&jetMask,&recoJetsFlavor);
        ////evtWeightProb[0] = probability of 0 Btags...... evtWeightProb[3] = probability of 3 Btags
        ////put event in each btag bin, weighted by evtWeightprob[0], evtWeightprob[1],
        //// evtWeightprob[2], evtWeightprob[3] for nb = 0, 1, 2, 3+
        //tr.registerDerivedVec("bTagSF_EventWeightProb_Central"+myVarSuffix_, evtWeightProb_Central);
        
        /*************************************************/
        // Case 1: Up  value;                            
        /*************************************************/
        switch_Unc = 1; switch_udsg_Unc = 0;
        SetBtagSFunc(switch_Unc); SetBtagCFunc(switch_Unc);
        SetCtagSFunc(switch_Unc); SetCtagCFunc(switch_Unc);
        SetMistagSFunc(switch_udsg_Unc); SetMistagCFunc(switch_udsg_Unc);
        double evtWeightSimple_Up  = GetSimpleCorrection(&inputJets,&jetMask,&recoJetsFlavor,&recoJetsBtag,wp);
        if( std::isnan( evtWeightSimple_Up) || std::isinf(evtWeightSimple_Up) )
        {
            evtWeightSimple_Up = 1.0;
        }
        tr.registerDerivedVar("bTagSF_EventWeightSimple_Up"+myVarSuffix_, evtWeightSimple_Up);

        //std::vector<double> *evtWeightProb_Up = GetCorrections(&inputJets,&jetMask,&recoJetsFlavor);
        //tr.registerDerivedVec("bTagSF_EventWeightProb_Up"+myVarSuffix_, evtWeightProb_Up);
        
        /*************************************************/
        // Case -1:Down  value;                            
        /*************************************************/
        switch_Unc = -1; switch_udsg_Unc = 0;
        SetBtagSFunc(switch_Unc); SetBtagCFunc(switch_Unc);
        SetCtagSFunc(switch_Unc); SetCtagCFunc(switch_Unc);
        SetMistagSFunc(switch_udsg_Unc); SetMistagCFunc(switch_udsg_Unc);
        double evtWeightSimple_Down  = GetSimpleCorrection(&inputJets,&jetMask,&recoJetsFlavor,&recoJetsBtag,wp);
        if( std::isnan( evtWeightSimple_Down) || std::isinf(evtWeightSimple_Down) )
        {
            evtWeightSimple_Down = 1.0;
        }
        tr.registerDerivedVar("bTagSF_EventWeightSimple_Down"+myVarSuffix_, evtWeightSimple_Down);

        //std::vector<double> *evtWeightProb_Down = GetCorrections(&inputJets,&jetMask,&recoJetsFlavor);
        //tr.registerDerivedVec("bTagSF_EventWeightProb_Down"+myVarSuffix_, evtWeightProb_Down);
        
        /*************************************************/
        // Mistag (udsg) Case 1: Up  value;                            
        /*************************************************/
        switch_Unc = 0; switch_udsg_Unc = 1;
        SetBtagSFunc(switch_Unc); SetBtagCFunc(switch_Unc);
        SetCtagSFunc(switch_Unc); SetCtagCFunc(switch_Unc);
        SetMistagSFunc(switch_udsg_Unc); SetMistagCFunc(switch_udsg_Unc);
        double evtWeightSimple_mistag_Up  = GetSimpleCorrection(&inputJets,&jetMask,&recoJetsFlavor,&recoJetsBtag,wp);
        if( std::isnan( evtWeightSimple_mistag_Up) || std::isinf(evtWeightSimple_mistag_Up) )
        {
            evtWeightSimple_mistag_Up = 1.0;
        }
        tr.registerDerivedVar("mistagSF_EventWeightSimple_Up"+myVarSuffix_, evtWeightSimple_mistag_Up);

        //std::vector<double> *evtWeightProb_mistag_Up =  GetCorrections(&inputJets,&jetMask,&recoJetsFlavor);
        //tr.registerDerivedVec("mistagSF_EventWeightProb_Up"+myVarSuffix_, evtWeightProb_mistag_Up);
        
        /*************************************************/
        // Case -1:Down  value;                            
        /*************************************************/
        switch_Unc = 0; switch_udsg_Unc = -1;
        SetBtagSFunc(switch_Unc); SetBtagCFunc(switch_Unc);
        SetCtagSFunc(switch_Unc); SetCtagCFunc(switch_Unc);
        SetMistagSFunc(switch_udsg_Unc); SetMistagCFunc(switch_udsg_Unc);
        double evtWeightSimple_mistag_Down  = GetSimpleCorrection(&inputJets,&jetMask,&recoJetsFlavor,&recoJetsBtag,wp);
        if( std::isnan( evtWeightSimple_mistag_Down) || std::isinf(evtWeightSimple_mistag_Down) )
        {
            evtWeightSimple_mistag_Down = 1.0;
        }
        tr.registerDerivedVar("mistagSF_EventWeightSimple_Down"+myVarSuffix_, evtWeightSimple_mistag_Down);

        //std::vector<double> *evtWeightProb_mistag_Down = GetCorrections(&inputJets,&jetMask,&recoJetsFlavor);
        //tr.registerDerivedVec("mistagSF_EventWeightProb_Down"+myVarSuffix_, evtWeightProb_mistag_Down);
    }

    //Operator
    void operator()(NTupleReader& tr)
    {
        registerVarToNTuples(tr);
    }

    //member variables
    bool debug;
    int btagSFunc, mistagSFunc;
    int btagCFunc, ctagCFunc, mistagCFunc;
    std::shared_ptr<TH2F> h_eff_b, h_eff_c, h_eff_udsg;
    std::string inFileName, MCBranch, JetsVec, JetMask, BJetsVec, JetsFlavor, myVarSuffix_;
    BTagCalibration calib;
    BTagCalibrationReader reader, readerUp, readerDown;
};

#endif
