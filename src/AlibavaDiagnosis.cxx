// ***********************************************
// 
// Monitoring and diagnosis tool for the Alibava
// data-taking and (post-)processing 
//
// ***********************************************
// Author: Jordi Duarte-Campderros (IFCA-CERN), 
// jorge.duarte.campderros@cern.ch.
// 
// ***********************************************
#include "AlibavaDiagnosis.h"

// The alibava run and event headers
//#include "AuxiliaryStructures.h"
#include "ALIBAVA.h"

// ROOT 
#include "TFile.h"
#include "TCanvas.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TProfile.h"

#include "TROOT.h"

// System headers
//#include <fstream>
//#include <map>
//#include <algorithm>
//#include <numeric>
//#include <stdexcept>
//#include <utility>

AlibavaDiagnosis::AlibavaDiagnosis(const int & chipnumber):
    _chip_number(chipnumber)
{
}

AlibavaDiagnosis::~AlibavaDiagnosis()
{
    for(auto & h: _histos)
    {
        if(h.second != nullptr)
        {
            delete h.second;
            h.second = nullptr;
        }
    }
}

void AlibavaDiagnosis::book_plots()
{
    // The monitoring plot initialization
    // ----------------------------------
    // Note that most of the histograms are defined here
    // with a dummy number of bins and/or axis ranges, 
    // they are going to be fixed dynamically afterwards

    // 1. Calibration (needs a calibration plot)
    const std::string suffix_name("chip-"+std::to_string(_chip_number));
    _histos["calibration"] = new TH2F(std::string("histo_Calibration_"+suffix_name).c_str(),
            std::string("Calibration curves "+suffix_name+";Channel number;Injected pulse [e] x10^{3};Signal average [ADC]").c_str(),
            ALIBAVA::NOOFCHANNELS,0,ALIBAVA::NOOFCHANNELS-1,
            1,0,1);
    // 2. Pedestal/Noise (needs pedestal)
    _histos["pedestal"] = new TH2F(std::string("histo_Pedestals_"+suffix_name).c_str(),
            std::string("Pedestals "+suffix_name+";Channel number; Pedestal [ADC]").c_str(),
            ALIBAVA::NOOFCHANNELS,0,ALIBAVA::NOOFCHANNELS-1,
            1,0,1);
    _histos["noise"] = new TH2F(std::string("histo_Noise_"+suffix_name).c_str(),
            std::string("Noise "+suffix_name+";Channel number; Noise [ADC]").c_str(),
            ALIBAVA::NOOFCHANNELS,0,ALIBAVA::NOOFCHANNELS-1,
            1,0,1);
    // 3. Temperature
    _histos["temperature"] = new TH2F(std::string("histo_Temperature_"+suffix_name).c_str(),
            std::string("Temperature per event "+suffix_name+";Event number; Noise [ADC]").c_str(),
            1,0,1,
            1,0,1);
    // 4. TDC
    _histos["tdc"] = new TH2F(std::string("histo_TDC_"+suffix_name).c_str(),
            std::string("TDC per event "+suffix_name+";Event number; Time [ns]").c_str(),
            1,0,1,
            1,0,1);
    // 5. Signal (needs pedestal)
    _histos["signal"] = new TH1F(std::string("histo_Signal_"+suffix_name).c_str(),
            std::string("Signal "+suffix_name+";Signal [ADC]; Entries").c_str(),
            1,0,1);
    // 6. Hits   (needs pedestal)
    _histos["hits"] = new TH1F(std::string("histo_Hits_"+suffix_name).c_str(),
            std::string("Signal "+suffix_name+";Channel number; Entries").c_str(),
            ALIBAVA::NOOFCHANNELS,0,ALIBAVA::NOOFCHANNELS-1);
    // 7. Time profile (needs pedestal)
    _histos["timeprofile"] = new TProfile(std::string("histo_TimeProfile_"+suffix_name).c_str(),
            std::string("Time profile "+suffix_name+";Time [ns]; Signal average [ADC]").c_str(),
            50,0,100);
    // 8. Noise per event (needs pedestal)
    _histos["noiseevent"] = new TH2F(std::string("histo_NoiseEvent_"+suffix_name).c_str(),
            std::string("Noise per event "+suffix_name+";Event number [ns]; Noise [ADC]").c_str(),
            1,0,1,
            1,0,1);
    _histos["commonnoiseevent"] = new TH2F(std::string("histo_CommonNoiseEvent_"+suffix_name).c_str(),
            std::string("Common Noise per event "+suffix_name+";Event number [ns]; Common Noise [ADC]").c_str(),
            1,0,1,
            1,0,1);

    // Should it be assigned to no-where or belong to any file?
    // --> for all the histos, then AddDirectory(0)a
    for(auto & h: _histos)
    {
        dynamic_cast<TH1*>(h.second)->SetDirectory(0);
    }
}

void AlibavaDiagnosis::book_plot(const std::string & name, const TObject * theplot)
{
    // Consistency check
    if(_histos.find(name) != _histos.end())
    {
        std::cerr << "[AlibavaDiagnosis::book_plot] The plot object '"
            << name << "' already exists or it is already booked a different object "
            << " with the same name" << std::endl;
        // Throw exception?
        return;
    }

    // Get the object and clone it without appending it in the gDirectory of TROOT
    const std::string suffix_name("chip-"+std::to_string(this->_chip_number));
    _histos[name] = theplot->Clone(std::string("externalbook_histo_"+name+"_"+suffix_name).c_str());
    dynamic_cast<TH1*>(_histos[name])->SetDirectory(0);
}

template<class ROOTTYPE> 
    ROOTTYPE* AlibavaDiagnosis::get_diagnostic_plot(const std::string & plotname)
{
    if(_histos.find(plotname) == _histos.end())
    {
        std::cerr << "[AlibavaDiagnosis::get_diagnostic_plot ERROR] Invalid plotname"
            << " [" << plotname << "] " << std::endl;
        return nullptr;
    }
    return static_cast<ROOTTYPE*>(_histos[plotname]);
}

// Declaration of the used types
template TH2F*     AlibavaDiagnosis::get_diagnostic_plot(const std::string&);
template TProfile* AlibavaDiagnosis::get_diagnostic_plot(const std::string&);

const std::vector<TObject*> AlibavaDiagnosis::get_calibration_plots() const
{
    // Get the calibration plots booked externally
    // XXX FIXME - the name should be defined as a global (in the ALIBAVA.h probably)
    std::vector<TObject*> theobjects;
    std::string calname("calibration_profile_");
    for(int ichannel = 0; ichannel < ALIBAVA::NOOFCHANNELS; ++ichannel)
    {
        // -- first check whether they are here or not
        const std::string name(calname+std::to_string(ichannel));
        if(_histos.find(name) == _histos.end())
        {
            std::cerr << "[AlibavaDiagnosis::get_calibration_plots] The plot object '"
                << name << "' was not booked! Inconsistent use of this method...'" << std::endl;
            // XXX throw exception ?
            continue;
        }
        theobjects.push_back(_histos.at(name));
    }
    return theobjects;
}

void AlibavaDiagnosis::set_calibration_plot(const std::vector<TObject*> & curves)
{
    // First check at least there is something
    if(curves.size() == 0)
    {
        std::cerr << "[AlibavaDiagnosis::set_calibration_plots] No calibration plots "
            << "present in the argument. Incoherent use of this class... '" << std::endl;
        return;
    }
    // Get the info from the objects in order to re-bin the histogram
    TProfile * auxobj = static_cast<TProfile*>(curves[0]);
    const int bins_y = auxobj->GetNbinsX();
    const float ymin = auxobj->GetXaxis()->GetBinLowEdge(1);
    const float ymax = auxobj->GetXaxis()->GetBinUpEdge(bins_y);
    
    // a handler to work with the histogram
    TH2F * thehisto = static_cast<TH2F*>(_histos["calibration"]);
    
    // Rebin the histogram depending the profiles bins
    thehisto->SetBins(thehisto->GetNbinsX(),thehisto->GetXaxis()->GetBinLowEdge(1),thehisto->GetXaxis()->GetBinUpEdge(thehisto->GetNbinsX()),
            bins_y,ymin,ymax);
    // And fill the histogram (remember was it fill in channel order)
    for(unsigned int i=0; i < curves.size(); ++i)
    {
        // The vector elements orden corresponds to the channel number
        const int ich = static_cast<int>(i);
        const TProfile* prf = static_cast<TProfile*>(curves[i]);
        // for each channel, let's obtain the injected charge (x-profile) with
        // its measured ADC counts (y-profile)
        for(int ibin=1; ibin < prf->GetNbinsX()+1; ++ibin)
        {
            const int bin = thehisto->FindBin(ich,prf->GetBinCenter(ibin));
            thehisto->SetBinContent(bin,prf->GetBinContent(ibin));
        }
    }
}

void AlibavaDiagnosis::deliver_plots()
{
    // Store all external booked plots if any
    // XXX: a directory can be used to create a more ordered structure
    // in that case it is needed either a TFile argument or a gDirectory...
    for(auto & h: _histos)
    {
        if(h.second != nullptr && std::string(h.second->GetName()).find("externalbook") == 0)
        {
            h.second->Write();
        }
    }

    // Plot all the monitoring plots in a canvas and 
    // store that canvas in the file
    // ----------------------------------------------
    // XXX FIXME: Define the Style
    // Define the canvas
    const int wsize = 650;
    const int hsize = 850;
    TCanvas * canvas = new TCanvas(std::string("monitor_canvas_"+std::to_string(_chip_number)).c_str(),
            std::string("Monitor plots Beetle "+std::to_string(_chip_number)).c_str(),
            wsize,hsize);
    // and its layout: colummns, row, separation betwwen columns and rows
    canvas->Divide(2,4,0.01,0.01);
    // And plot (first define the order of plotting (except those pads
    // with more than 1 histo per pad)
    // Note that the map is defining the pad where the histogram should go
    // and the Draw option to be used
    const std::map<int,std::pair<std::string,std::string>> plotorder = 
    { 
        {1,{"calibration","SURF3Z"}},
        {3,{"signal",""}},
        {4,{"hits",""}},
        {5,{"timeprofile",""}},
        {6,{"temperature",""}},
        {7,{"tdc",""}}
    };
    // Define the pads and fill them
    for(const auto & padplot: plotorder)
    {
        canvas->cd(padplot.first);
        // Check if the plot is present
        if(_histos[padplot.second.first] == nullptr)
        {
            // This behaviour is also allowed (see calibration or pedestal
            // processing cases)
            delete canvas;
            canvas = nullptr;
            return;
        }
        _histos[padplot.second.first]->Draw(padplot.second.second.c_str());
    }
    // And now the composite pads
    // --------------------------
    // PAD 2: Pedestals - Noise
    auto pad_comp = canvas->cd(2);
    pad_comp->Divide(2,1,0.001,0.001);
    pad_comp->cd(1);
    _histos["pedestal"]->Draw();
    pad_comp->cd(2);
    _histos["noise"]->Draw();
    // PAD 8: Noise per Event 
    auto pad_comp8 = canvas->cd(8);
    pad_comp8->Divide(2,1,0.001,0.001);
    pad_comp8->cd(1);
    _histos["noiseevent"]->Draw();
    pad_comp8->cd(2);
    _histos["commonnoiseevent"]->Draw();
    // Save it to the TFile present in the gDirectory 
    // (managed by the client of this class: IOManager)
    canvas->Write();
    // and deallocate memory 
    delete canvas;
    canvas = nullptr;
    // Note that the pads are already deleted by the canvas
}
