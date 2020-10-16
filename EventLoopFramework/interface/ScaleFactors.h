#ifndef SCALEFACTORS_H
#define SCALEFACTORS_H

#include "TGraph.h"
#include "TFile.h"
#include "TKey.h"
#include "t-channel_Analysis/EventLoopFramework/interface/Utility.h"

class ScaleFactors
{
private:
    std::string myVarSuffix_;
    bool firstPrint_;
    bool printMeanError_;
    bool printGetBinError_;
    std::set<std::string> binError_message_;

    std::shared_ptr<TH2F> eleSFHistoTight_;        
    std::shared_ptr<TH2F> eleSFHistoIso_;
    std::shared_ptr<TH2F> eleSFHistoIP2D_;
    std::shared_ptr<TH2F> eleSFHistoReco_;
    std::shared_ptr<TH2F> eleSFHistoTrig_;
    std::shared_ptr<TH2F> muSFHistoMedium_;
    std::shared_ptr<TH2F> muSFHistoIso_;
    std::shared_ptr<TH2F> muSFHistoTrig_;
    std::shared_ptr<TH2F> nimuSFHistoTrig_;
    std::shared_ptr<TGraph> muSFHistoReco_;
    std::shared_ptr<TH1F> puSFHisto_;
    std::shared_ptr<TH1F> puSFUpHisto_;
    std::shared_ptr<TH1F> puSFDownHisto_;
    std::shared_ptr<TH2F> L1Prefireing_;
    std::map<std::string, double> sfMeanMap_;

    template<typename T> std::shared_ptr<T>& getHisto(TFile& f, std::shared_ptr<T>& h, const TString& name)
    {        
        if(name != "") h.reset( static_cast<T*>(f.Get(name)) );   
        else std::cerr<<utility::color("Warning: A needed scale factor histogram is set to nullptr and using 1.0 as the default", "red")<<std::endl;
        return h;
    }

    template<typename Th, typename Tb> int findBin(const std::shared_ptr<Th>& h, const Tb v, const std::string& axis, const std::string& message = "")
    {
        int bin = -1;
        if(h)
        {
            if(axis=="X")
            {
                bin = h->GetXaxis()->FindBin(v);
                if( v >= h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()) ) bin = h->GetNbinsX();
            }
            else if(axis=="Y")
            {
                bin = h->GetYaxis()->FindBin(v);
                if( v >= h->GetYaxis()->GetBinUpEdge(h->GetNbinsY()) ) bin = h->GetNbinsY();           
            }
        }
        if(bin == -1 && firstPrint_)
        {
            binError_message_.insert(message); 
            printGetBinError_ = true;
        }
        return bin;
    }

    double htScaleFactor(const int nJets, const double HT, const std::string& runYear) const
    {
        //All values updated for v1.0 on October 4, 2019. Commented out values are from the old setup (just 2016 and 2017).
        double norm = 1.0;
        double expo = 0.0;
        if( runYear == "2016" ) 
        {
            norm = 0.04176*nJets + 0.8915;
            expo = (-0.02358*nJets - 0.08940)/1000;
            //norm = 0.03422*nJets + 0.9367;
            //expo = (-0.02310*nJets - 0.0940 )/1000;
        }
        else if( runYear == "2017" ) 
        {
            norm = 0.03952*nJets + 0.9171;
            expo = (-0.03275*nJets - 0.03795)/1000;
            //norm = 0.02565*nJets + 0.9635;
            //expo = (-0.01418*nJets - 0.1101 )/1000;
        }
        else if( runYear == "2018pre" ) 
        {
            norm = 0.03953*nJets + 0.8934;
            expo = (-0.03832*nJets + 0.01108)/1000;
        }
        else if( runYear == "2018post" ) 
        {
            norm = 0.01953*nJets + 1.0150;
            expo = (-0.01409*nJets - 0.1325)/1000;
        }
        return norm*exp( expo*HT ); 
    }

    double htScaleFactorFlat2000(const int nJets, const double HT, const std::string& runYear) const
    {
        return ( HT > 2000.0 ) ? htScaleFactor(nJets,2000.0,runYear) : htScaleFactor(nJets,HT,runYear);
    }

    double htScaleFactorNJet8(const double HT, const std::string& runYear) const
    {
        double norm = 1.0;
        double expo = 0.0;      
        if( runYear == "2016" ) 
        {
            norm = 1.318;
            expo = -0.000349;
            //norm = 1.307;
            //expo = -0.0003416;
        }
        else if( runYear == "2017" )
        {
            norm = 1.216;
            expo = -0.0003109;
            //norm = 1.215;
            //expo = -0.0002613;
        }
        else if( runYear == "2018pre" )
        {
            norm = 1.229;
            expo = -0.0002752;
        }
        else if( runYear == "2018post" )
        {
            norm = 1.299;
            expo = -0.0003484;
        }
        return norm*exp( expo*HT );
    }

    double htScaleFactorMG(const int nJets, const double HT, const std::string& runYear) const
    {
        double norm = 1.0;
        double expo = 0.0;
        if( runYear == "2016" ) 
        {
            norm = -0.002429*nJets + 1.044;
            expo = (0.01214*nJets - 0.1198)/1000;
        }
        else if( runYear == "2017" ) 
        {
            norm = 0.1086*nJets + 0.5567;
            expo = (-0.8978*nJets + 0.2198)/1000;
        }
        else if( runYear == "2018pre" ) 
        {
            norm = 0.03349*nJets + 0.9464;
            expo = (-0.02058*nJets - 0.1488)/1000;
        }
        else if( runYear == "2018post" ) 
        {
            norm = 0.02729*nJets + 0.9632;
            expo = (-0.01304*nJets - 0.1637)/1000;
        }
        return norm*exp( expo*HT );
    }

    double htScaleFactorNJet8MG(const double HT, const std::string& runYear) const
    {
        double norm = 1.0;
        double expo = 0.0;
        if( runYear == "2016" ) 
        {
            norm = 1.073;
            expo = -0.00009169;
        }
        else if( runYear == "2017" )
        {
            norm = 1.236;
            expo = -0.0002949;
        }
        else if( runYear == "2018pre" )
        {
            norm = 1.2;
            expo = -0.0002516;
        }
        else if( runYear == "2018post" )
        {
            norm = 1.273;
            expo = -0.0003242;
        }
        return norm*exp( expo*HT );
    }

    double getMean(const std::string& name)
    {
        if( sfMeanMap_[name] == 0.0 )
        {
            printMeanError_ = true;
            return 1.0;
        }
        else
        {
            return sfMeanMap_[name];
        }
    }

    void scaleFactors(NTupleReader& tr)
    {
        // --------------------------------------------------------------------------------------
        // Calculate scale SF and variation
        // Following the example in SusyAnaTools PDFUncertainty.h
        // The scale weights are calculated using the envelope method and we ignore all anti-correlated variations (5 and 7)
        // --------------------------------------------------------------------------------------
        const auto& scaleWeights = tr.getVec<double>("ScaleWeights");
        const auto& filetag      = tr.getVar<std::string>("filetag");
        const auto& runYear      = tr.getVar<std::string>("runYear");
        
        double scaleWeightNominal = 1.0;
        std::vector<double>  myScaleWeights(6, 1.0);
        if( scaleWeights.size() == 9 ) 
        {
            //If there are not exactly 9 scale factors, then the vector  was filled incorrectly
            scaleWeightNominal = scaleWeights[0];
            myScaleWeights[0]  = scaleWeights[1];
            myScaleWeights[1]  = scaleWeights[2];
            myScaleWeights[2]  = scaleWeights[3];
            myScaleWeights[3]  = scaleWeights[4];
            myScaleWeights[4]  = scaleWeights[6];
            myScaleWeights[5]  = scaleWeights[8];
        }
        
        double scaleWeightUpperBound = *std::max_element( std::begin(myScaleWeights), std::end(myScaleWeights) );
        double scaleWeightLowerBound = *std::min_element( std::begin(myScaleWeights), std::end(myScaleWeights) );
        if( !std::isfinite(scaleWeightNominal) || !std::isfinite(scaleWeightUpperBound) || !std::isfinite(scaleWeightLowerBound) ) 
        {
            //TODO: There are some NaN/Inf values in the Diboson channel - still need to figure out why this is an issue.
            scaleWeightNominal = 1.0;
            scaleWeightUpperBound = 1.0;
            scaleWeightLowerBound = 1.0;
        }

        double scaleWeightUpperBound_corr = scaleWeightUpperBound;
        double scaleWeightLowerBound_corr = scaleWeightLowerBound;
        if(sfMeanMap_.find(filetag+"_sclUp") != sfMeanMap_.end() && sfMeanMap_.find(filetag+"_sclDown") != sfMeanMap_.end()) 
        {            
            const double meanUp = getMean(filetag+"_sclUp");
            const double meanDown = getMean(filetag+"_sclDown");
            scaleWeightUpperBound_corr = (1/meanUp)*scaleWeightUpperBound;
            scaleWeightLowerBound_corr = (1/meanDown)*scaleWeightLowerBound;
        }

        tr.registerDerivedVar("scaleWeightUp"+myVarSuffix_,   scaleWeightUpperBound_corr);
        tr.registerDerivedVar("scaleWeightDown"+myVarSuffix_, scaleWeightLowerBound_corr);
        tr.registerDerivedVar("scaleWeightUpUncor"+myVarSuffix_,   scaleWeightUpperBound);
        tr.registerDerivedVar("scaleWeightDownUncor"+myVarSuffix_, scaleWeightLowerBound);
        tr.registerDerivedVar("scaleWeightNom"+myVarSuffix_,  scaleWeightNominal);

        // --------------------------------------------------------------------------------------
        // Calculate parton shower variation
        // Note: not all samples have these weights stored, give them default value of 1. 
        // --------------------------------------------------------------------------------------        
        double PSweight_ISRUp = 1.0, PSweight_ISRDown = 1.0; 
        double PSweight_FSRUp = 1.0, PSweight_FSRDown = 1.0; 
        double PSweight_ISRUp_2 = 1.0, PSweight_ISRDown_2 = 1.0; 
        double PSweight_FSRUp_2 = 1.0, PSweight_FSRDown_2 = 1.0; 
        if(tr.hasVar("PSweights"))
        {
            const auto& PSweights = tr.getVec<double>("PSweights");
            if(PSweights.size() >= 12) // should have size of 14, but just put 12 or more to be able to use the sample with the bug
            {
                // Get nominal one so we can normalize it
                double MEweight = PSweights[0];
                // reduced variations, i.e. varying Pythia params isr:muRfac and fsr:muRfac with factor 1/sqrt(2) and sqrt(2)
                PSweight_ISRUp = PSweights[2]/MEweight;
                PSweight_FSRUp = (PSweights[3]/MEweight < 10.0) ? PSweights[3]/MEweight : 1.0;
                PSweight_ISRDown = PSweights[4]/MEweight;
                PSweight_FSRDown = (PSweights[5]/MEweight < 10.0) ? PSweights[5]/MEweight : 1.0;
                // nominal variations, i.e. varying Pythia params isr:muRfac and fsr:muRfac with factor 1/2 and 2
                PSweight_ISRUp_2 = PSweights[6]/MEweight;
                PSweight_FSRUp_2 = (PSweights[7]/MEweight < 10.0) ? PSweights[7]/MEweight : 1.0;
                PSweight_ISRDown_2 = PSweights[8]/MEweight;
                PSweight_FSRDown_2 = (PSweights[9]/MEweight < 10.0) ? PSweights[9]/MEweight : 1.0;
            }
        }
        tr.registerDerivedVar("PSweight_ISRUp"+myVarSuffix_,   PSweight_ISRUp);
        tr.registerDerivedVar("PSweight_ISRDown"+myVarSuffix_, PSweight_ISRDown);
        tr.registerDerivedVar("PSweight_FSRUp"+myVarSuffix_,   PSweight_FSRUp);
        tr.registerDerivedVar("PSweight_FSRDown"+myVarSuffix_, PSweight_FSRDown);
        tr.registerDerivedVar("PSweight_ISRUp_2"+myVarSuffix_,   PSweight_ISRUp_2);
        tr.registerDerivedVar("PSweight_ISRDown_2"+myVarSuffix_, PSweight_ISRDown_2);
        tr.registerDerivedVar("PSweight_FSRUp_2"+myVarSuffix_,   PSweight_FSRUp_2);
        tr.registerDerivedVar("PSweight_FSRDown_2"+myVarSuffix_, PSweight_FSRDown_2);

        // --------------------------------------------------------------------------------------
        // Now calculate the PDF scale factor and uncertainty 
        // Based on the 100 different replica values stored in PDFweights using envelope method and the median
        // --------------------------------------------------------------------------------------
        const auto& PDFweights = tr.getVec<double>("PDFweights");
        double central = 1.0, NNPDF_from_median_up = 1.0, NNPDF_from_median_down = 1.0;
        if(PDFweights.size() > 0)
        {
            const double reqCL = 0.68; //Choose a confidence level for the uncertainty
            std::vector<double> sortedPDFWeights = PDFweights; //Cannot sort a constant
            std::sort( sortedPDFWeights.begin() + 1, sortedPDFWeights.end() );
        
            const int upper = std::round( 0.5*(1 + reqCL)*100.0 );
            //const int lower = 1 + std::round( 0.5*(1 - reqCL)*100.0 );

            central  = 0.5*( sortedPDFWeights[50] + sortedPDFWeights[51] ); //Exactly 100 entries
            //const double errminus = abs(central - sortedPDFWeights[lower]);
            const double errplus  = abs(central - sortedPDFWeights[upper]);

            NNPDF_from_median_up   = central + errplus;
            NNPDF_from_median_down = central - errplus;
            NNPDF_from_median_up   = NNPDF_from_median_up/central   > 2.0 ? 1.0 : NNPDF_from_median_up/central   < -2.0 ? 1.0 : NNPDF_from_median_up/central;        
            NNPDF_from_median_down = NNPDF_from_median_down/central > 2.0 ? 1.0 : NNPDF_from_median_down/central < -2.0 ? 1.0 : NNPDF_from_median_down/central;
        }
        if( !std::isfinite(central) || !std::isfinite(NNPDF_from_median_up) || !std::isfinite(NNPDF_from_median_down) ) 
        {
            NNPDF_from_median_up    = 1.0;
            NNPDF_from_median_down  = 1.0;
            central                 = 1.0;
        }
        
        double NNPDF_from_median_up_corr = NNPDF_from_median_up;
        double NNPDF_from_median_down_corr = NNPDF_from_median_down;
        if(sfMeanMap_.find(filetag+"_pdf_Up") != sfMeanMap_.end() && sfMeanMap_.find(filetag+"_pdf_Down") != sfMeanMap_.end()) 
        {
            const double meanUp = getMean(filetag+"_pdf_Up");
            const double meanDown = getMean(filetag+"_pdf_Down");
            NNPDF_from_median_up_corr = (1/meanUp)*NNPDF_from_median_up;
            NNPDF_from_median_down_corr = (1/meanDown)*NNPDF_from_median_down;
        }

        tr.registerDerivedVar( "PDFweightUp"+myVarSuffix_,   NNPDF_from_median_up_corr );
        tr.registerDerivedVar( "PDFweightDown"+myVarSuffix_, NNPDF_from_median_down_corr );
        tr.registerDerivedVar( "PDFweightUpUncor"+myVarSuffix_,   NNPDF_from_median_up );
        tr.registerDerivedVar( "PDFweightDownUncor"+myVarSuffix_, NNPDF_from_median_down );
        tr.registerDerivedVar( "PDFweightNom"+myVarSuffix_,  central );
        
        // --------------------------------------------------------------------------------------
        // Now calculate the electron scale factor and uncertainty 
        // --------------------------------------------------------------------------------------
        const auto& electrons         = tr.getVec<TLorentzVector>("Electrons");
        const auto& goodElectrons     = tr.getVec<bool>("GoodElectrons"+myVarSuffix_);
        
        double totGoodElectronSF      = 1.0, totGoodElectronSF_Up   = 1.0, totGoodElectronSF_Down = 1.0;
        double totGoodElectronSFErr   = 0.0, totGoodElectronSFPErr2 = 0.0;
        double noTrigGoodElectronSF   = 1.0, noTrigGoodElectronSFErr = 0.0, noTrigGoodElectronSFPErr2 = 0.0;
        for( unsigned int iel = 0; iel < electrons.size(); iel++ ) 
        {
            if( !goodElectrons.at(iel) ) continue;
            
            //Get the scale factor from the rootfile
            const double elpt  = electrons.at(iel).Pt();
            const double eleta = electrons.at(iel).Eta();

            const int xbinElTight = findBin(eleSFHistoTight_, eleta, "X", "el id x");
            const int ybinElTight = findBin(eleSFHistoTight_, elpt,  "Y", "el id y");                
            const int xbinElIso   = findBin(eleSFHistoIso_,   eleta, "X", "el iso x");
            const int ybinElIso   = findBin(eleSFHistoIso_,   elpt,  "Y", "el iso y");
            const int xbinElReco  = findBin(eleSFHistoReco_,  eleta, "X", "el reco x");
            const int ybinElReco  = findBin(eleSFHistoReco_,  elpt,  "Y", "el reco y");
            const int xbinElTrig  = findBin(eleSFHistoTrig_,  elpt,  "X", "el trigger x");
            const int ybinElTrig  = findBin(eleSFHistoTrig_,  eleta, "Y", "el trigger y");
            if( xbinElTight != -1 && ybinElTight != -1 && xbinElIso != -1 && ybinElIso != -1 && xbinElReco != -1 && ybinElReco != -1 ) 
            {
                const double eleTightSF       = eleSFHistoTight_->GetBinContent( xbinElTight, ybinElTight );
                const double eleTightSFErr    = eleSFHistoTight_->GetBinError( xbinElTight, ybinElTight );
                const double eleTightPErr     = eleTightSFErr/eleTightSF;
                double eleIsoSF               = eleSFHistoIso_->GetBinContent( xbinElIso, ybinElIso );
                double eleIsoSFErr            = eleSFHistoIso_->GetBinError( xbinElIso, ybinElIso );
                double eleIsoPErr             = eleIsoSFErr/eleIsoSF;
                const double eleRecoSF        = eleSFHistoReco_->GetBinContent( xbinElReco, ybinElReco );
                const double eleRecoSFErr     = eleSFHistoReco_->GetBinError( xbinElReco, ybinElReco );
                const double eleRecoPErr      = eleRecoSFErr/eleRecoSF;
                const double eleTrigSF        = (eleSFHistoTrig_) ? eleSFHistoTrig_->GetBinContent( xbinElTrig, ybinElTrig ) : 1.0;
                const double eleTrigSFErr     = (eleSFHistoTrig_) ? eleSFHistoTrig_->GetBinError( xbinElTrig, ybinElTrig ) : 0.0;
                const double eleTrigPErr      = eleTrigSFErr/eleTrigSF;
                
                if( runYear == "2016" ) 
                { 
                    //The lepton scale factor is the multiplication of the three different scale factors. To get the proper error, you sum up the percentage errors in quadrature.
                    //If this is the year 2016, we need to add the IP2D histogram scale factors into the Iso scale factor
                    const double eleIP2DSF    = eleSFHistoIP2D_->GetBinContent( xbinElIso, ybinElIso );
                    const double eleIP2DSFErr = eleSFHistoIP2D_->GetBinError( xbinElIso, ybinElIso );
                    const double eleIP2DPErr  = eleIP2DSFErr/eleIP2DSF;

                    eleIsoSF    = eleIsoSF*eleIP2DSF;
                    eleIsoPErr  = utility::addInQuad( eleIsoPErr, eleIP2DPErr );
                    eleIsoSFErr = eleIsoPErr*eleIsoSF;
                }
                
                const double eleNoTrigSF      = eleTightSF*eleIsoSF*eleRecoSF;
                const double eleTotSF         = eleNoTrigSF*eleTrigSF; 
                const double eleNoTrigPErr    = utility::addInQuad( eleTightPErr, eleIsoPErr, eleRecoPErr );
                const double eleTotPErr       = utility::addInQuad( eleNoTrigPErr, eleTrigPErr );

                totGoodElectronSF       *= eleTotSF; 
                noTrigGoodElectronSF    *= eleNoTrigSF;
                totGoodElectronSFPErr2  += eleTotPErr*eleTotPErr;
                noTrigGoodElectronSFPErr2 += eleNoTrigPErr*eleNoTrigPErr;
            }            
        }

        totGoodElectronSFErr    = sqrt(totGoodElectronSFPErr2) * totGoodElectronSF;
        noTrigGoodElectronSFErr = sqrt(noTrigGoodElectronSFPErr2) * noTrigGoodElectronSF;
        totGoodElectronSF_Up    = totGoodElectronSF + totGoodElectronSFErr;
        totGoodElectronSF_Down  = totGoodElectronSF - totGoodElectronSFErr;

        tr.registerDerivedVar( "totGoodElectronSF"+myVarSuffix_,      totGoodElectronSF );
        tr.registerDerivedVar( "totGoodElectronSFErr"+myVarSuffix_,   totGoodElectronSFErr );
        tr.registerDerivedVar( "totGoodElectronSF_Up"+myVarSuffix_,   totGoodElectronSF_Up );
        tr.registerDerivedVar( "totGoodElectronSF_Down"+myVarSuffix_, totGoodElectronSF_Down );
        tr.registerDerivedVar( "noTrigGoodElectronSF"+myVarSuffix_,   noTrigGoodElectronSF );
        tr.registerDerivedVar( "noTrigGoodElectronSFErr"+myVarSuffix_, noTrigGoodElectronSFErr );

        // --------------------------------------------------------------------------------------
        // Adding code for implementing muon scale factors
        // --------------------------------------------------------------------------------------
        const auto& muons         = tr.getVec<TLorentzVector>("Muons");
        const auto& goodMuons     = tr.getVec<bool>("GoodMuons"+myVarSuffix_);
        const auto& nonisoMuons     = tr.getVec<bool>("NonIsoMuons"+myVarSuffix_);

        double totGoodMuonSF      = 1.0, totGoodMuonSF_Up     = 1.0, totGoodMuonSF_Down    = 1.0;
        double totGoodMuonSFErr   = 0.0, totGoodMuonSFPErr2   = 0.0;
        double noTrigGoodMuonSF   = 1.0, noTrigGoodMuonSFErr  = 0.0, noTrigGoodMuonSFPErr2 = 0.0;
        double totNonIsoMuonSF    = 1.0, totNonIsoMuonSF_Up   = 1.0, totNonIsoMuonSF_Down  = 1.0;
        double totNonIsoMuonSFErr = 0.0, totNonIsoMuonSFPErr2 = 0.0;        

        for( unsigned int imu = 0; imu < muons.size(); imu++ ) 
        {            
            //Get the scale factor from the rootfile
            const double mupt = muons.at(imu).Pt();
            const double mueta = muons.at(imu).Eta();

            const int xbinMuMedium     = findBin(muSFHistoMedium_, mupt,       "X", "mu id x");
            const int ybinMuMedium     = findBin(muSFHistoMedium_, abs(mueta), "Y", "mu id y");
            const int xbinMuIso        = findBin(muSFHistoIso_,    mupt,       "X", "mu iso x");
            const int ybinMuIso        = findBin(muSFHistoIso_,    abs(mueta), "Y", "mu iso y");            
            const int xbinMuTrig       = findBin(muSFHistoTrig_,   mupt,       "X", "mu trigger x");
            const int ybinMuTrig       = findBin(muSFHistoTrig_,   mueta,      "Y", "mu trigger y");
            const int xbinNonIsoMuTrig = findBin(nimuSFHistoTrig_, mupt,       "X", "mu iso trigger x");
            const int ybinNonIsoMuTrig = findBin(nimuSFHistoTrig_, mueta,      "Y", "mu iso trigger y");
            if( xbinMuMedium != -1 && ybinMuMedium != -1 && xbinMuIso != -1 && ybinMuIso != -1 ) 
            {
                //The SUSLepton Twiki claims that the errors in the histogrm are purely statistical and can be ignored and recommends a 3% error for each leg (ID+IP+ISO)
                const double muMediumSF         = muSFHistoMedium_->GetBinContent( xbinMuMedium, ybinMuMedium );
                const double muIsoSF            = muSFHistoIso_->GetBinContent( xbinMuIso, ybinMuIso );
                const double muTrigSF           = (muSFHistoTrig_) ? muSFHistoTrig_->GetBinContent( xbinMuTrig, ybinMuTrig ) : 1.0;
                const double muNonIsoTrigSF     = (nimuSFHistoTrig_) ? nimuSFHistoTrig_->GetBinContent( xbinNonIsoMuTrig, ybinNonIsoMuTrig ) : 1.0;
                const double muTrigSFErr        = (muSFHistoTrig_) ? muSFHistoTrig_->GetBinError( xbinMuTrig, ybinMuTrig ) : 0.0;
                const double muNonIsoTrigSFErr  = (nimuSFHistoTrig_) ? nimuSFHistoTrig_->GetBinError( xbinNonIsoMuTrig, ybinNonIsoMuTrig ) : 0.0;
                const double muTrigSFPErr       = muTrigSFErr/muTrigSF;
                const double muNonIsoTrigSFPErr = muNonIsoTrigSFErr/muNonIsoTrigSF;
                double muNoTrigSF               = muMediumSF*muIsoSF; 
                double muTotSF                  = muNoTrigSF*muTrigSF;
                double muNonIsoTotSF            = muNoTrigSF*muNonIsoTrigSF;
                const double muNoTrigSFPErr2    = 0.03*0.03;
                const double muTotSFPErr2       = muNoTrigSFPErr2 + muTrigSFPErr*muTrigSFPErr;
                const double muNonIsoTotSFPErr2 = muNoTrigSFPErr2 + muNonIsoTrigSFPErr*muNonIsoTrigSFPErr;
                
                if( runYear == "2016" ) 
                {
                    //For the general track reconstruction they claim that the errors for the systematic still need to be finalized - does not seem to have been finalized as of Dec 2018
                    //This reconstruction value only exists for 2016 - SUS SF people say the 3% will include the reco scale factor uncertainty for now
                    const double muRecoSF    = muSFHistoReco_->Eval( mueta );
                    muTotSF                 *= muRecoSF;
                    muNoTrigSF              *= muRecoSF;
                    muNonIsoTotSF           *= muRecoSF;
                }

                if( goodMuons.at(imu) )
                {                
                    totGoodMuonSF           *= muTotSF;
                    noTrigGoodMuonSF        *= muNoTrigSF;
                    totGoodMuonSFPErr2      += muTotSFPErr2;
                    noTrigGoodMuonSFPErr2   += muNoTrigSFPErr2;
                }
                if( nonisoMuons.at(imu) )
                {
                    totNonIsoMuonSF         *= muNonIsoTotSF;
                    totNonIsoMuonSFPErr2    += muNonIsoTotSFPErr2;
                }
            }
        }

        totGoodMuonSFErr     = sqrt(totGoodMuonSFPErr2)*totGoodMuonSF;
        totGoodMuonSF_Up     = totGoodMuonSF + totGoodMuonSFErr;
        totGoodMuonSF_Down   = totGoodMuonSF - totGoodMuonSFErr;
        noTrigGoodMuonSFErr  = sqrt(noTrigGoodMuonSFPErr2)*noTrigGoodMuonSF;
        totNonIsoMuonSFErr   = sqrt(totNonIsoMuonSFPErr2)*totNonIsoMuonSF;
        totNonIsoMuonSF_Up   = totNonIsoMuonSF + totNonIsoMuonSFErr;
        totNonIsoMuonSF_Down = totNonIsoMuonSF - totNonIsoMuonSFErr;

        tr.registerDerivedVar( "totGoodMuonSF"+myVarSuffix_,       totGoodMuonSF );
        tr.registerDerivedVar( "totGoodMuonSFErr"+myVarSuffix_,    totGoodMuonSFErr );
        tr.registerDerivedVar( "totGoodMuonSF_Up"+myVarSuffix_,    totGoodMuonSF_Up );
        tr.registerDerivedVar( "totGoodMuonSF_Down"+myVarSuffix_,  totGoodMuonSF_Down );
        tr.registerDerivedVar( "noTrigGoodMuonSF"+myVarSuffix_,    noTrigGoodMuonSF );
        tr.registerDerivedVar( "noTrigGoodMuonSFErr"+myVarSuffix_, noTrigGoodMuonSFErr );
        tr.registerDerivedVar( "totNonIsoMuonSF"+myVarSuffix_,     totNonIsoMuonSF );
        tr.registerDerivedVar( "totNonIsoMuonSFErr"+myVarSuffix_,  totNonIsoMuonSFErr );
        tr.registerDerivedVar( "totNonIsoMuonSF_Up"+myVarSuffix_,  totNonIsoMuonSF_Up );
        tr.registerDerivedVar( "totNonIsoMuonSF_Down"+myVarSuffix_,totNonIsoMuonSF_Down );

        // --------------------------------------------------------------------------------------
        // Adding a scale factor that corrects the disagreement between data and MC for Ht
        // --------------------------------------------------------------------------------------
        const auto& NGoodJets_pt30  = tr.getVar<int>("NGoodJets_pt30"+myVarSuffix_);
        const auto& HT_trigger_pt30 = tr.getVar<double>("HT_trigger_pt30"+myVarSuffix_);
        const auto& isSignal        = tr.getVar<bool>("isSignal");

        // Get the uncorrected ht scale factors
        const double htDerivedweightUncor         = htScaleFactor(NGoodJets_pt30, HT_trigger_pt30, runYear);
        const double htDerivedweightFlat2000Uncor = htScaleFactorFlat2000(NGoodJets_pt30, HT_trigger_pt30, runYear);
        const double htDerivedweightNJet7Uncor    = htScaleFactor(7, HT_trigger_pt30, runYear);
        const double htDerivedweightMGUncor       = htScaleFactorMG(NGoodJets_pt30, HT_trigger_pt30, runYear);
        const double fit2NJetBin8                 = htScaleFactorNJet8(HT_trigger_pt30, runYear);
        const double fit2NJetBin567               = htScaleFactor(8, HT_trigger_pt30, runYear);
        const double fit2NJetBin8MG               = htScaleFactorNJet8MG(HT_trigger_pt30, runYear);
        const double fit2NJetBin567MG             = htScaleFactor(8, HT_trigger_pt30, runYear);
                
        double htDerivedweight = 1.0, htDerivedweightFlat2000 = 1.0, htDerivedweightNJet7 = 1.0, htDerivedweightMG = 1.0;
        double htScaleUpUncor = 1.0, htScaleDownUncor = 1.0, htScaleUpMGUncor = 1.0, htScaleDownMGUncor = 1.0;
        double htScaleUp = 1.0, htScaleDown = 1.0, htScaleUpMG = 1.0, htScaleDownMG = 1.0;
        if( sfMeanMap_.find(filetag+"_ht") != sfMeanMap_.end() && !isSignal ) 
        {
            htDerivedweight         = (1.0/getMean(filetag+"_ht"         ))*htDerivedweightUncor;
            htDerivedweightFlat2000 = (1.0/getMean(filetag+"_ht_flat2000"))*htDerivedweightFlat2000Uncor;
            htDerivedweightNJet7    = (1.0/getMean(filetag+"_ht_njet7"   ))*htDerivedweightNJet7Uncor;
            htDerivedweightMG       = (1.0/getMean(filetag+"_ht_MG"      ))*htDerivedweightMGUncor;

            htScaleUpUncor          = htDerivedweight*(fit2NJetBin8/fit2NJetBin567);
            htScaleDownUncor        = htDerivedweight*(fit2NJetBin567/fit2NJetBin8);
            htScaleUpMGUncor        = htDerivedweightMG*(fit2NJetBin8MG/fit2NJetBin567MG);
            htScaleDownMGUncor      = htDerivedweightMG*(fit2NJetBin567MG/fit2NJetBin8MG);
            htScaleUp               = (1.0/getMean(filetag+"_htUp"     ))*htScaleUpUncor;
            htScaleDown             = (1.0/getMean(filetag+"_htDown"   ))*htScaleDownUncor;
            htScaleUpMG             = (1.0/getMean(filetag+"_ht_MGUp"  ))*htScaleUpMGUncor;
            htScaleDownMG           = (1.0/getMean(filetag+"_ht_MGDown"))*htScaleDownMGUncor;
        }

        tr.registerDerivedVar( "htDerivedweight"+myVarSuffix_, htDerivedweight);
        tr.registerDerivedVar( "htDerivedweightFlat2000"+myVarSuffix_, htDerivedweightFlat2000 );
        tr.registerDerivedVar( "htDerivedweightNJet7"+myVarSuffix_, htDerivedweightNJet7 );
        tr.registerDerivedVar( "htDerivedweightMG"+myVarSuffix_, htDerivedweightMG);
        tr.registerDerivedVar( "htDerivedweightUncor"+myVarSuffix_, htDerivedweightUncor);
        tr.registerDerivedVar( "htDerivedweightFlat2000Uncor"+myVarSuffix_, htDerivedweightFlat2000Uncor);
        tr.registerDerivedVar( "htDerivedweightNJet7Uncor"+myVarSuffix_, htDerivedweightNJet7Uncor );
        tr.registerDerivedVar( "htDerivedweightMGUncor"+myVarSuffix_, htDerivedweightMGUncor);
        tr.registerDerivedVar( "htScaleUpUncor"+myVarSuffix_, htScaleUpUncor);
        tr.registerDerivedVar( "htScaleDownUncor"+myVarSuffix_, htScaleDownUncor);
        tr.registerDerivedVar( "htScaleUpMGUncor"+myVarSuffix_, htScaleUpMGUncor);
        tr.registerDerivedVar( "htScaleDownMGUncor"+myVarSuffix_, htScaleDownMGUncor);
        tr.registerDerivedVar( "htScaleUp"+myVarSuffix_, htScaleUp);
        tr.registerDerivedVar( "htScaleDown"+myVarSuffix_, htScaleDown);
        tr.registerDerivedVar( "htScaleUpMG"+myVarSuffix_, htScaleUpMG);
        tr.registerDerivedVar( "htScaleDownMG"+myVarSuffix_, htScaleDownMG);

        //-----------------------------------------------------------------------------
        // Adding a scale factor for pileup 
        // For 2016 and 2018: Grab the individual pileup weight from the histogram found in PileupHistograms_*.root
        // For 2017: Using the puWeight stored in the nTuples
        // ----------------------------------------------------------------------------        
        const auto& puWeight  = tr.getVar<double>("puWeight");
        const auto& puSysUp   = tr.getVar<double>("puSysUp");
        const auto& puSysDown = tr.getVar<double>("puSysDown");
        const auto& tru_npv   = tr.getVar<double>("TrueNumInteractions");

        double puWeightUnCorr = 1.0, puSysUpUnCorr = 1.0, puSysDownUnCorr = 1.0;
        if( runYear == "2016" || runYear == "2018pre" || runYear == "2018post") 
        {
            puWeightUnCorr = puSFHisto_->GetBinContent( findBin(puSFHisto_, tru_npv, "X", "nom pu") );
            puSysUpUnCorr = puSFUpHisto_->GetBinContent( findBin(puSFUpHisto_, tru_npv, "X", "up pu") );
            puSysDownUnCorr = puSFDownHisto_->GetBinContent( findBin(puSFDownHisto_, tru_npv, "X", "down pu") );
        }
        else if( runYear == "2017") 
        {
            puWeightUnCorr = puWeight; 
            puSysUpUnCorr = puSysUp;
            puSysDownUnCorr = puSysDown;
        }        
        tr.registerDerivedVar( "puWeightUnCorr"+myVarSuffix_,  puWeightUnCorr);
        tr.registerDerivedVar( "puSysUpUnCorr"+myVarSuffix_,   puSysUpUnCorr);
        tr.registerDerivedVar( "puSysDownUnCorr"+myVarSuffix_, puSysDownUnCorr);

        // Adding correction to the pileup weight
        double puWeightCorr  = puWeight;
        double puSysUpCorr   = puSysUp;
        double puSysDownCorr = puSysDown;
        if( sfMeanMap_.find(filetag+"_pu") != sfMeanMap_.end() && !isSignal )
        {
            const double mean     = getMean(filetag+"_pu");
            const double meanUp   = getMean(filetag+"_pu_Up");
            const double meanDown = getMean(filetag+"_pu_Down");
            puWeightCorr = (1.0/mean)*puWeightUnCorr;
            puSysUpCorr  = (1.0/meanUp)*puSysUpUnCorr;
            puSysDownCorr= (1.0/meanDown)*puSysDownUnCorr;
        }
        tr.registerDerivedVar( "puWeightCorr"+myVarSuffix_, puWeightCorr);
        tr.registerDerivedVar( "puSysUpCorr"+myVarSuffix_, puSysUpCorr);
        tr.registerDerivedVar( "puSysDownCorr"+myVarSuffix_, puSysDownCorr);
        
        // --------------------------------------------------------------------------------------
        // Adding top pt reweighting for ttbar MC (Powheg)
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
        // 13 TeV all combined
        // --------------------------------------------------------------------------------------
        double topPtScaleFactor = 1.0;
        auto* topPtVec = new std::vector<double>();
        if(filetag == "TT" || filetag == "2016_TT" || filetag.find("TTTo") != std::string::npos)
        {
            const double a=0.0615, b=-0.0005;
            auto SF = [&](const double pt){return exp(a + b*pt);};
            
            const auto& GenParticles        = tr.getVec<TLorentzVector>("GenParticles");
            const auto& GenParticles_PdgId  = tr.getVec<int>("GenParticles_PdgId");

            for(unsigned int gpi=0; gpi < GenParticles.size(); gpi++)
            {
                if( abs(GenParticles_PdgId[gpi]) == 6 )
                {
                    topPtScaleFactor *= SF( GenParticles[gpi].Pt() );
                    topPtVec->push_back( GenParticles[gpi].Pt() );
                }
            }
            topPtScaleFactor = sqrt(topPtScaleFactor);
        }
        
        tr.registerDerivedVar( "topPtScaleFactor"+myVarSuffix_, topPtScaleFactor);
        tr.registerDerivedVec( "topPtVec"+myVarSuffix_, topPtVec);
        
        // --------------------------------------------------------------------------------------
        // Adding reweighting recipe to emulate Level 1 ECAL prefiring
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe
        // --------------------------------------------------------------------------------------
        double prefiringScaleFactor = 1.0, prefiringScaleFactorUp = 1.0, prefiringScaleFactorDown = 1.0;
        if( runYear == "2016" || runYear == "2017" )
        {
            const auto& Jets = tr.getVec<TLorentzVector>("Jets"+myVarSuffix_);            
            const auto& GoodJets_pt30 = tr.getVec<bool>("GoodJets_pt30"+myVarSuffix_);            
            for(unsigned int ijet = 0; ijet < Jets.size(); ++ijet)
            {            
                if(!GoodJets_pt30[ijet]) continue;
                const TLorentzVector& jet = Jets.at(ijet);
                int bin = L1Prefireing_->FindBin(jet.Eta(), jet.Pt());
                const double weight = L1Prefireing_->GetBinContent(bin);
                const double weightErr = std::max(0.2*weight, L1Prefireing_->GetBinError(bin));
                const double weightUp = weight + weightErr;
                const double weightDown = weight - weightErr;
                prefiringScaleFactor *= 1 - weight;
                prefiringScaleFactorUp *= 1 - weightDown;
                prefiringScaleFactorDown *= 1 - weightUp;
            }
        }
        tr.registerDerivedVar( "prefiringScaleFactor"+myVarSuffix_, prefiringScaleFactor);                    
        tr.registerDerivedVar( "prefiringScaleFactorUp"+myVarSuffix_, prefiringScaleFactorUp);                    
        tr.registerDerivedVar( "prefiringScaleFactorDown"+myVarSuffix_, prefiringScaleFactorDown);                    

        // --------------------------------------------------------------------------------------
        // Registering a variable that is the nominal total weight with lepton scale factor, btag scale factor, ht scale factor
        // --------------------------------------------------------------------------------------
        const auto& Weight         = tr.getVar<double>("Weight");
        const auto& bTagWeight     = tr.getVar<double>("bTagSF_EventWeightSimple_Central"+myVarSuffix_);
        const auto& NGoodElectrons = tr.getVar<int>("NGoodElectrons"+myVarSuffix_);
        const auto& NGoodMuons     = tr.getVar<int>("NGoodMuons"+myVarSuffix_);
        const auto& NNonIsoMuons   = tr.getVar<int>("NNonIsoMuons"+myVarSuffix_);
        
        double totalEventWeight         = -1.0;
        double totalEventWeightMG       = -1.0;
        double totalEventWeightNIM      = -1.0;
        double totalEventWeightNIM_ht   = -1.0;
        if( NGoodElectrons == 1 ) 
        {
            totalEventWeight = Weight*bTagWeight*totGoodElectronSF*htDerivedweight*prefiringScaleFactor*puWeightCorr;
            totalEventWeightMG = Weight*bTagWeight*totGoodElectronSF*htDerivedweightMG*prefiringScaleFactor*puWeightCorr;
        }
        else if ( NGoodMuons == 1 ) 
        {
            totalEventWeight = Weight*bTagWeight*totGoodMuonSF*htDerivedweight*prefiringScaleFactor*puWeightCorr;
            totalEventWeightMG = Weight*bTagWeight*totGoodMuonSF*htDerivedweightMG*prefiringScaleFactor*puWeightCorr;
        }
        if( NNonIsoMuons == 1 ) 
        {
            totalEventWeightNIM = Weight*totNonIsoMuonSF*prefiringScaleFactor*puWeightCorr;
            totalEventWeightNIM_ht = Weight*totNonIsoMuonSF*prefiringScaleFactor*puWeightCorr*htDerivedweight;
        }

        tr.registerDerivedVar( "totalEventWeight"+myVarSuffix_, totalEventWeight );
        tr.registerDerivedVar( "totalEventWeightMG"+myVarSuffix_, totalEventWeightMG );
        tr.registerDerivedVar( "totalEventWeightNIM"+myVarSuffix_, totalEventWeightNIM );
        tr.registerDerivedVar( "totalEventWeightNIM_ht"+myVarSuffix_, totalEventWeightNIM_ht );

        if(printGetBinError_) firstPrint_ = false;
    }

public:
    ScaleFactors( const std::string& runYear, const std::string& leptonFileName, const std::string& puFileName, const std::string& meanFileName, const std::string& myVarSuffix = "" )
        : myVarSuffix_(myVarSuffix), firstPrint_(true), printMeanError_(false), printGetBinError_(false)
    {
        std::cout<<"Setting up ScaleFactors"<<std::endl;
        TH1::AddDirectory(false); //According to Joe, this is a magic incantation that lets the root file close - if this is not here, there are segfaults?

        // Getting Lepton scale factor histograms
        TFile SFRootFile( leptonFileName.c_str() );
        TString eleSFHistoTightName, eleSFHistoIsoName, eleSFHistoRecoName, eleSFHistoTrigName;
        TString muSFHistoMediumName,  muSFHistoIsoName,  muSFHistoRecoName,  muSFHistoTrigName, nimuSFHistoTrigName;

        if( runYear == "2016")
        {
            eleSFHistoTightName = "Run2016_CutBasedTightNoIso94XV2";
            eleSFHistoIsoName = "Run2016_Mini";
            eleSFHistoRecoName = "EGamma_SF2D";
            eleSFHistoTrigName = "TrigEff_2016_num_el_pt40_trig_5jCut_htCut_DeepCSV";
            muSFHistoMediumName = "sf_mu_mediumID";
            muSFHistoIsoName = "sf_mu_mediumID_mini02";
            muSFHistoTrigName = "TrigEff_2016_num_mu_pt40_trig_5jCut_htCut_DeepCSV";
            //nimuSFHistoTrigName = "TrigEff_2016_num_nimu_pt40_trig_4jCut";
            nimuSFHistoTrigName = ""; //just for calculating non iso muon scale factors
        }
        else if( runYear == "2017")
        {
            eleSFHistoTightName = "Run2017_CutBasedTightNoIso94XV2";
            eleSFHistoIsoName = "Run2017_MVAVLooseTightIP2DMini";
            eleSFHistoRecoName = "EGamma_SF2D";
            eleSFHistoTrigName = "TrigEff_2017_num_el_pt40_trig_5jCut_htCut_DeepCSV";
            muSFHistoMediumName = "NUM_MediumID_DEN_genTracks_pt_abseta";
            muSFHistoIsoName = "TnP_MC_NUM_MiniIso02Cut_DEN_MediumID_PAR_pt_eta";
            muSFHistoTrigName = "TrigEff_2017_num_mu_pt40_trig_5jCut_htCut_DeepCSV";
            //nimuSFHistoTrigName = "TrigEff_2017_num_nimu_pt40_trig_4jCut";
            nimuSFHistoTrigName = ""; //just for calculating non iso muon scale factors
        }
        else if( runYear == "2018pre")
        {
            eleSFHistoTightName = "Run2018_CutBasedTightNoIso94XV2";
            eleSFHistoIsoName = "Run2018_Mini";
            eleSFHistoRecoName = "EGamma_SF2D";
            eleSFHistoTrigName = "TrigEff_2018pre_num_el_pt40_trig_5jCut_htCut_DeepCSV";
            muSFHistoMediumName = "NUM_MediumID_DEN_TrackerMuons_pt_abseta";
            muSFHistoIsoName = "TnP_MC_NUM_MiniIso02Cut_DEN_MediumID_PAR_pt_eta";
            muSFHistoTrigName = "TrigEff_2018pre_num_mu_pt40_trig_5jCut_htCut_DeepCSV";
            nimuSFHistoTrigName = "";
        }
        else if( runYear == "2018post")
        {
            eleSFHistoTightName = "Run2018_CutBasedTightNoIso94XV2";
            eleSFHistoIsoName = "Run2018_Mini";
            eleSFHistoRecoName = "EGamma_SF2D";
            eleSFHistoTrigName = "TrigEff_2018post_num_el_pt40_trig_5jCut_htCut_DeepCSV";
            muSFHistoMediumName = "NUM_MediumID_DEN_TrackerMuons_pt_abseta";
            muSFHistoIsoName = "TnP_MC_NUM_MiniIso02Cut_DEN_MediumID_PAR_pt_eta";
            muSFHistoTrigName = "TrigEff_2018post_num_mu_pt40_trig_5jCut_htCut_DeepCSV";
            nimuSFHistoTrigName = "";
        }

        getHisto(SFRootFile, eleSFHistoTight_, eleSFHistoTightName);
        getHisto(SFRootFile, eleSFHistoIso_,   eleSFHistoIsoName);
        getHisto(SFRootFile, eleSFHistoReco_,  eleSFHistoRecoName);
        getHisto(SFRootFile, eleSFHistoTrig_,  eleSFHistoTrigName);
        getHisto(SFRootFile, muSFHistoMedium_, muSFHistoMediumName);
        getHisto(SFRootFile, muSFHistoIso_,    muSFHistoIsoName);
        getHisto(SFRootFile, muSFHistoTrig_,   muSFHistoTrigName);
        getHisto(SFRootFile, nimuSFHistoTrig_, nimuSFHistoTrigName);
        if( runYear == "2016" ) 
        {
            getHisto(SFRootFile, eleSFHistoIP2D_, "Run2016_MVAVLooseIP2D");//In 2016, the isolation SF histogram is separate from the IP2D cut scale factor histogram.
            getHisto(SFRootFile, muSFHistoReco_,  "ratio_eff_aeta_dr030e030_corr");//Only 2016 requires the track reconstruction efficiency.
        }
        SFRootFile.Close();

        // Getting mean of some scale factors to keep total number of events the same after apply the scale factor
        TFile SFMeanRootFile( meanFileName.c_str() );
        TIter next(SFMeanRootFile.GetListOfKeys());
        TKey* key;
        while((key = static_cast<TKey*>(next())))
        {
            std::shared_ptr<TH1> h( static_cast<TH1*>(key->ReadObj()) );
            std::string name( h->GetTitle() );
            sfMeanMap_.insert(std::pair<std::string, double>(name, h->GetMean()));
        }
        SFMeanRootFile.Close();

        // Getting the pile up histograms
        TFile puRootFile( puFileName.c_str() );
        if( runYear == "2016" || runYear == "2018pre" || runYear == "2018post") 
        {
            getHisto(puRootFile, puSFHisto_,     "pu_weights_central");
            getHisto(puRootFile, puSFUpHisto_,   "pu_weights_up");
            getHisto(puRootFile, puSFDownHisto_, "pu_weights_down");
        }
        else if( runYear == "2017")
        {
            getHisto(puRootFile, puSFHisto_,     "pu_ratio_central");
            getHisto(puRootFile, puSFUpHisto_,   "pu_ratio_up");
            getHisto(puRootFile, puSFDownHisto_, "pu_ratio_down");
        }
        puRootFile.Close();

        // Get the L1prefiring scale factor histogram
        TFile L1PrefiringFile("L1prefiring_jetpt_2017BtoF.root");
        TString prefireHistoName = "L1prefiring_jetpt_2017BtoF";
        if(runYear == "2016")
        {
            prefireHistoName = "L1prefiring_jetpt_2016BtoH";
        }        
        getHisto(L1PrefiringFile, L1Prefireing_, prefireHistoName);
        L1PrefiringFile.Close();
    }

    ~ScaleFactors() 
    {
        if(printMeanError_)   std::cerr<<utility::color("Error: Scale Factor mean is 0.0 setting it to 1.0", "red")<<std::endl;

        std::string message = "";
        for(const auto& m : binError_message_) message += (message=="") ? m : ", "+m;
        if(printGetBinError_) std::cerr<<utility::color("Warning: There was an error in extracting the bin index for the following scale factor(s): \""+message+"\"", "red")<<std::endl;
    }

    void operator()(NTupleReader& tr)
    {
        scaleFactors(tr);
    }
};

#endif
