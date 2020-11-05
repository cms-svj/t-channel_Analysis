#ifndef Analyze_h
#define Analyze_h

#include <TH1D.h>
#include <TH2D.h>
#include <TProfile2D.h>
#include <TEfficiency.h>
#include <TTree.h>
#include <TFile.h>

#include <map>
#include <string>

class NTupleReader;

class Analyze 
{
private:
    class TH1DInfo
    {
    public:
        std::string name;
        int nBins;
        double low;
        double high;        
    };

    class TH2DInfo
    {
    public:
        std::string name;
        int nBinsX;
        double lowX;
        double highX;        
        int nBinsY;
        double lowY;
        double highY;        
    };

    class TH2DProfileInfo
    {
    public:
        std::string name;
        int nBinsX;
        double lowX;
        double highX;        
        int nBinsY;
        double lowY;
        double highY;        
        double lowZ;
        double highZ;
    };

public:
    std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
    std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
    std::map<std::string, std::shared_ptr<TProfile>>  my_tp_histos;
    std::map<std::string, std::shared_ptr<TProfile2D>>  my_2d_tp_histos;
    std::map<std::string, std::shared_ptr<TEfficiency>>  my_efficiencies;
    bool initHistos;
    
    Analyze();
    ~Analyze(){};
    
    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
    void InitHistos(const std::map<std::string, bool>& cutMap, const std::vector<TH1DInfo>& histInfos, const std::vector<TH2DInfo>& hist2DInfos, const std::vector<TH2DProfileInfo>& hist2DProfileInfos);
    void WriteHistos(TFile* outfile);
};

#endif
