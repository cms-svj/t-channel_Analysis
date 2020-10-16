#ifndef AnalyzerBase_h
#define AnalyzerBase_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>
#include <TFile.h>
#include <map>
#include <string>

#include "t-channel_Analysis/EventLoopFramework/interface/Histo.h"
class NTupleReader;

class AnalyzeBase
{
protected:
    std::vector<std::unique_ptr<Histo_Base>> my_Histos;

public:
    virtual void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false) = 0;
    void WriteHistos(TFile* outfile)
    {
        outfile->cd();

        for(const auto& h : my_Histos)
        {
            h->Write();
        }
    }

    bool printEventNum(const int maxevents, const int evtNum, const int divisor = 1000)
    {
        bool b = false;
        if( maxevents != -1 && evtNum >= maxevents + 1 ) b = true;
        if( evtNum % divisor == 0 ) printf( " Event %i\n", evtNum );
        return b;
    }

    void Fill(NTupleReader& tr)
    {
        for(auto& h : my_Histos)
        {
            h->Fill(tr);
        }
    }
};

#endif
