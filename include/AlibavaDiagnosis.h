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
#ifndef ALIBAVADIAGNOSIS_H
#define ALIBAVADIAGNOSIS_H

//#include "IOManager.h"

// ROOT Classes
//

// System headers
#include<string>
#include<map>
#include<vector>

// forward declarations
// the canvas
class TObject;
class TTree;


class AlibavaDiagnosis
{
    public:
        AlibavaDiagnosis(const int & chipnumber);
        AlibavaDiagnosis() = delete;
        ~AlibavaDiagnosis();

        // Initialization of the histograms
        void book_plots();
        
        // Initialization of external plot objects, the object is cloned here
        void book_plot(const std::string & name, const TObject * plotobject);

        // Fill the set of predefined monitor plots (except the calibration
        // ones, which needs a specific function, see set_calibration_plot)
        void set_diagnostic_plots(const std::pair<std::vector<float>,std::vector<float> > & pednoise);

        // Get the objects needed to create the 2dim calibration plot (note
        // this function make sense only in the calibration file
        const std::vector<TObject*> get_calibration_plots() const;
        // Set 3dim calibration plot 
        void set_calibration_plot(const std::vector<TObject*> & curves);

        // Fill all the diagnostic plots
        void fill_diagnostic_plots(TTree * event_tree, TTree * header_tree);

        // Store the defined plots to a canvas
        void deliver_plots();

    private:
        int _chip_number;
        // The map of plots
        std::map<std::string,TObject*> _histos; 

        // Get the relation between the histogram name and the branched needed
        // to fill it (note that the branch change depending the processing or
        // not previously of the pedestal file
        std::map<std::string,std::string> get_branch_names(bool was_pedfile_proc);
        // Decide whether or not the current channel could be a seed cluster
        bool is_seed_cluster(const float & signal, const float & noise, const float & SoverN = 5.0) const;

        // the draw options
        inline std::map<std::string,std::string> get_draw_option()
        {
            return { {"temperature","L"} ,{"tdc","L"}, {"pedestal","L"},
                {"noise","L"},{"commonmode","L"},{"noiseevent","L"},
                {"signal",""},{"hits",""},{"timeprofile","prof"}    };
        }
};

#endif
