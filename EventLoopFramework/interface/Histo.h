#ifndef Histo_h
#define Histo_h

#include "t-channel_Analysis/EventLoopFramework/interface/NTupleReader.h"
#include "t-channel_Analysis/EventLoopFramework/interface/Utility.h"
#include <iostream>
#include <string>
#include <TH3.h>

class Histo_Base
{
public:
    virtual void Fill(const NTupleReader&) = 0;
    virtual void Write() const = 0;
};

template<typename Type>
class Histo_FirstChild : public Histo_Base
{
protected:
    std::string name_;
    int nBinsX_, nBinsY_, nBinsZ_;
    double lowX_, lowY_, lowZ_;
    double highX_, highY_, highZ_;
    std::string varX_, varY_, varZ_;
    std::string varXMask_;
    std::vector<std::string> cuts_;
    std::vector<std::string> weights_;
    std::unique_ptr<Type> histo_;    

    void Write() const
    {
        histo_->Write();
    }

    bool passCuts(const NTupleReader& tr) const
    {
        bool pass = true;
        for(const auto& cutName : cuts_)
        {
            const auto& cut = tr.getVar<bool>(cutName);
            pass = pass && cut;
            if(!pass) break;
        }
        return pass;
    }

    double getWeight(const NTupleReader& tr) const
    {
        double weight = 1.0;
        for(const auto& wName : weights_)
        {
            const auto& w = tr.getVar<double>(wName);
            weight *= w;
        }        
        return weight;
    }

    template<typename T> T tlvGetValue(const TLorentzVector& lv, const std::string& varType) const
    {
        T val = 0;
        if     (varType.find("P()")   != std::string::npos) val = lv.P();
        else if(varType.find("Pt()")  != std::string::npos) val = lv.Pt();
        else if(varType.find("Phi()") != std::string::npos) val = lv.Phi();
        else if(varType.find("Eta()") != std::string::npos) val = lv.Eta();
        else if(varType.find("M()")   != std::string::npos) val = lv.M();
        else if(varType.find("E()")   != std::string::npos) val = lv.E();
        else std::cout<<utility::color("No option for \""+varType+"\" found", "red")<<std::endl;
        return val;
    }
};

class Histo1D : public Histo_FirstChild<TH1>
{
private:
    template<typename T> void vecFill(const double weight, const NTupleReader& tr) const
    {
        const auto& vec = tr.getVec<T>(varX_);
        const auto& vecMask = (varXMask_ != "") ? tr.getVec<bool>(varXMask_) : std::vector<bool>();
        for(unsigned int j = 0; j < vec.size(); j++)
        {
            if(vecMask.size() > 0)
            {
                if(!vecMask[j]) continue;
            }
            histo_->Fill( vec[j], weight );
        }
    }

    void vecFilltlv(const double weight, const NTupleReader& tr) const
    {
        const std::string& varType = utility::split("last",  varX_, ".");
        const std::string&     var = utility::split("first", varX_, ".");
        const auto& vec = tr.getVec<TLorentzVector>(var);
        const auto& vecMask = (varXMask_ != "") ? tr.getVec<bool>(varXMask_) : std::vector<bool>();
        for(unsigned int j = 0; j < vec.size(); j++)
        {
            if(vecMask.size() > 0)
            {
                if(!vecMask[j]) continue;
            }
            histo_->Fill( tlvGetValue<double>(vec[j],varType), weight ); 
        }
    }

    void fillHisto(const NTupleReader& tr) const
    {
        std::string type;
        tr.getType(utility::split("first", utility::split("first", varX_, "."), "["), type);
        const double weight = getWeight(tr);

        if(type.find("std::vector<") != std::string::npos)
        {
            if     (type.find("double")         != std::string::npos) vecFill<double>      ( weight, tr );
            else if(type.find("unsigned int")   != std::string::npos) vecFill<unsigned int>( weight, tr );
            else if(type.find("int")            != std::string::npos) vecFill<int>         ( weight, tr );
            else if(type.find("bool")           != std::string::npos) vecFill<bool>        ( weight, tr );
            else if(type.find("float")          != std::string::npos) vecFill<float>       ( weight, tr );
            else if(type.find("char")           != std::string::npos) vecFill<char>        ( weight, tr );
            else if(type.find("short")          != std::string::npos) vecFill<short>       ( weight, tr );
            else if(type.find("long")           != std::string::npos) vecFill<long>        ( weight, tr );        
            else if(type.find("TLorentzVector") != std::string::npos) vecFilltlv           ( weight, tr );        
            else std::cout<<utility::color("Type \""+type+"\" is not an option for vector variable \""+varX_+"\"", "red")<<std::endl;
        }
        else
        {
            if     (type.find("double")         != std::string::npos) histo_->Fill( tr.getVar<double>(varX_),       weight );
            else if(type.find("unsigned int")   != std::string::npos) histo_->Fill( tr.getVar<unsigned int>(varX_), weight );
            else if(type.find("int")            != std::string::npos) histo_->Fill( tr.getVar<int>(varX_),          weight );
            else if(type.find("bool")           != std::string::npos) histo_->Fill( tr.getVar<bool>(varX_),         weight );
            else if(type.find("float")          != std::string::npos) histo_->Fill( tr.getVar<float>(varX_),        weight );
            else if(type.find("char")           != std::string::npos) histo_->Fill( tr.getVar<char>(varX_),         weight );
            else if(type.find("short")          != std::string::npos) histo_->Fill( tr.getVar<short>(varX_),        weight );
            else if(type.find("long")           != std::string::npos) histo_->Fill( tr.getVar<long>(varX_),         weight );        
            else if(type.find("TLorentzVector") != std::string::npos)
            {
                const std::string& varType = utility::split("last",  varX_, ".");
                const std::string&     var = utility::split("first", varX_, ".");
                histo_->Fill( tlvGetValue<double>(tr.getVar<TLorentzVector>(var), varType), weight );
            }
            else std::cout<<utility::color("Type \""+type+"\" is not an option for variable \""+varX_+"\"", "red")<<std::endl;
        }
    } 

public:
    void Fill(const NTupleReader& tr)
    {
        if(passCuts(tr))
        {
            fillHisto(tr);
        }
    }

    Histo1D(const std::string& name, const int nBinsX, const double lowX, const double highX, const std::string& varX, const std::vector<std::string>& cuts, const std::vector<std::string>& weights) 
        : Histo1D(name, nBinsX, lowX, highX, varX, "", cuts, weights)
    {
    }

    Histo1D(const std::string& name, const int nBinsX, const double lowX, const double highX, const std::string& varX, const std::string& varXMask, const std::vector<std::string>& cuts, const std::vector<std::string>& weights) 
    {
        name_ = name;
        nBinsX_ = nBinsX;
        lowX_ = lowX;
        highX_ = highX;
        varX_ = varX;
        varXMask_ = varXMask;
        cuts_ = cuts;
        weights_ = weights;
        histo_ = std::make_unique<TH1D>( name_.c_str(), name_.c_str(), nBinsX_, lowX_, highX_ );
    }
};

template<typename TX, typename TY>
class Histo2D : public Histo_FirstChild<TH2>
{
public:
    void Fill(const NTupleReader& tr)
    {
        if(passCuts(tr))
        {
            TX variableX = tr.getVar<TX>(varX_);
            TY variableY = tr.getVar<TY>(varY_);
            const double weight = getWeight(tr);
            histo_->Fill(variableX, variableY, weight);
        }
    }

    Histo2D(const std::string& name, 
            const int nBinsX, const double lowX, const double highX, const std::string& varX, 
            const int nBinsY, const double lowY, const double highY, const std::string& varY, 
            const std::vector<std::string>& cuts = {}, const std::vector<std::string>& weights = {})
    {
        name_ = name;
        nBinsX_ = nBinsX;
        lowX_ = lowX;
        highX_ = highX;
        varX_ = varX;
        nBinsY_ = nBinsY;
        lowY_ = lowY;
        highY_ = highY;
        varY_ = varY;
        cuts_ = cuts;
        weights_ = weights;
        histo_ = std::make_unique<TH2D>( name_.c_str(), name_.c_str(), nBinsX_, lowX_, highX_, nBinsY_, lowY_, highY_ );
    }
};

template<typename TX, typename TY, typename TZ>
class Histo3D : public Histo_FirstChild<TH3>
{
public:
    void Fill(const NTupleReader& tr)
    {
        if(passCuts(tr))
        {
            TX variableX = tr.getVar<TX>(varX_);
            TY variableY = tr.getVar<TY>(varY_);
            TZ variableZ = tr.getVar<TZ>(varZ_);
            const double weight = getWeight(tr);
            histo_->Fill(variableX, variableY, variableZ, weight);
        }
    }

    Histo3D(const std::string& name, 
            const int nBinsX, const double lowX, const double highX, const std::string& varX, 
            const int nBinsY, const double lowY, const double highY, const std::string& varY, 
            const int nBinsZ, const double lowZ, const double highZ, const std::string& varZ, 
            const std::vector<std::string>& cuts = {}, const std::vector<std::string>& weights = {})
    {
        name_ = name;
        nBinsX_ = nBinsX;
        lowX_ = lowX;
        highX_ = highX;
        varX_ = varX;
        nBinsY_ = nBinsY;
        lowY_ = lowY;
        highY_ = highY;
        varY_ = varY;
        nBinsZ_ = nBinsZ;
        lowZ_ = lowZ;
        highZ_ = highZ;
        varZ_ = varZ;
        cuts_ = cuts;
        weights_ = weights;
        histo_ = std::make_unique<TH3D>( name_.c_str(), name_.c_str(), nBinsX_, lowX_, highX_, nBinsY_, lowY_, highY_, nBinsZ_, lowZ_, highZ_ );
    }
};

#endif
