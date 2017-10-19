//  --- TESTING IOFortythieves

#include "IOFortythieves.h"

#include <string>
#include <iostream>

int main(int /*argc*/, char** /*argv[]*/) 
{
    const std::string fname("/home/duarte/Working/TB_CERN_2017/monitorplot_test/362_2017-05-20_06-15_gerva_MB2_iLGAD8533W1K05T_400V_146d0uA_20C_lat132_beam.root");
    IOFortythieves iof(fname,1);

    iof.initialize();
    std::cout << "Pedestal -----\n [";
    for(auto & ped: iof.pedestal())
    {
        std::cout << ped << " ";
    }
    std::cout << "]" <<  std::endl;
    std::cout << "Noise -----\n [";
    for(auto & ped: iof.noise())
    {
        std::cout << ped << " ";
    }
    std::cout << "]" <<  std::endl;
    std::cout << "Calibration -----\n [";
    for(auto & ped: iof.calibration())
    {
        std::cout << ped << " ";
    }
    std::cout << "]" <<  std::endl;

    std::cout << "N entries:" << iof.get_entries() << std::endl;
    std::cout << "ADC (corrected) for entry 1098" << std::endl;
    iof.process(34343);
    std::cout << "Event Number: " << iof.event_number() << " Time: " << iof.event_time() << " Temperature: " << iof.temperature() << std::endl;
    for(auto & ped: iof.adc_data())
    {
        std::cout << ped << " ";
    }
    std::cout << "]" <<  std::endl;

}
