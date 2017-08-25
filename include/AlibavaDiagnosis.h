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

        // Get the objects needed to create the 3dim calibration plot (note
        // this function make sense only in the calibration file
        const std::vector<TObject*> get_calibration_plots() const;
        // Set 3dim calibration plot 
        void set_calibration_plot(const std::vector<TObject*> & curves);

        // Store the defined plots to a canvas
        void deliver_plots();

        //
        template<class ROOTTYPE>
            ROOTTYPE* get_diagnostic_plot(const std::string & plotname);

    private:
        int _chip_number;
        // The map of plots
        std::map<std::string,TObject*> _histos; 
};

#endif
