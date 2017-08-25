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
class TCanvas;
class TObject;


class AlibavaDiagnosis
{
    public:
        AlibavaDiagnosis(const int & chipnumber);
        AlibavaDiagnosis() = delete;
        ~AlibavaDiagnosis();

        // Initialization of the histograms
        void book_plots();
        
        // Initialization of external plot objects  which are cloned here
        void book_plot(const std::string & name, const TObject * plotobject);

        // Store the defined plots to a canvas
        void deliver_plots();

        //
        template<class ROOTTYPE>
            ROOTTYPE* get_diagnostic_plot(const std::string & plotname);

    private:
        int _chip_number;
        // The map of plots
        // The canvas to be stored
        TCanvas * _canvas;
        std::map<std::string,TObject*> _histos; 
};

#endif
