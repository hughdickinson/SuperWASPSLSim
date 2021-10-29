//
//  main.cpp
//  SuperWASPSLSim
//
//  Created by Hugh Dickinson on 14/04/2021.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <algorithm>
#include <iterator>
#include <tuple>

#include <CCfits/FITS.h>
#include <CCfits/PHDU.h>
#include <CCfits/Table.h>
#include <CCfits/ExtHDU.h>
#include <CCfits/ColumnT.h>
#include <CCfits/FITSUtilT.h>

#include <boost/iterator/zip_iterator.hpp>

#include "EinsteinRadius.hpp"
#include "NonUniformCircularSource.h"

template<int NumSteps>
auto magnifications(double radiusUnitless, std::vector<double> range = {-1, 1}){
    
    auto magnifier  = NonUniformCircularSource(radiusUnitless);
    
    auto impactParam = 0.05;
    auto magnifications = std::array<double, NumSteps>{};
    
    auto elementCounter = 0;
    std::generate(magnifications.begin(), magnifications.end(), [&](){
        auto x = range[0] + elementCounter * (range[1] - range[0])/static_cast<float>(NumSteps);
        auto offset = std::hypot(x, impactParam);
        elementCounter++;
        return magnifier.magnification(offset);
    });
    
    return magnifications;
}


int testMagnificationSim(){
    auto solarRadius = 6.957e8*boost::units::si::meters;
    auto astronomicalUnit = 1.495978707e11*boost::units::si::meters;
    auto solarMass = 1.988409870698051e30*boost::units::si::kilograms;
    
    auto systemDistance = quantity<si::length>(2 * kiloparsec);
    auto systemSeparation = 1000.0*astronomicalUnit;
    auto lensMass = 20.0 * solarMass;
    auto sourceRadius = 5.0 * solarRadius;
    
    auto radiusUnitless = quantity<si::length>(sourceRadius) / *sourcePlaneEinsteinRadius<si::length>(systemDistance, systemDistance - systemSeparation, systemSeparation, lensMass);
    std::cout << radiusUnitless << "\t" << quantity<si::length>(sourceRadius)/radiusUnitless << "\t" << quantity<si::length>(sourceRadius) << std::endl;
    auto magValues = magnifications<200>(radiusUnitless);
    std::copy(magValues.begin(), magValues.end(), std::ostream_iterator<double>(std::cout, ", "));
    return 0;
}

auto testFitsLoad(std::string const & fileName){
    
    typedef std::vector<long> time_data_t;
    typedef std::vector<float> flux_data_t;
    typedef std::vector<float> hjd_data_t;
    
    auto fitsFile = std::unique_ptr<CCfits::FITS>(new CCfits::FITS(fileName, CCfits::Read));
    CCfits::Table & photometry = dynamic_cast<CCfits::Table &>(fitsFile->extension("PHOTOMETRY"));
    
    auto & timeColumn = photometry.column("TMID");
    auto & fluxColumn = photometry.column("TAMFLUX2");
    auto & fluxErrColumn = photometry.column("TAMFLUX2_ERR");
    
    auto timeData = std::unique_ptr<time_data_t>(new time_data_t(timeColumn.rows()));
    auto hjdData = std::unique_ptr<hjd_data_t>(new hjd_data_t(timeColumn.rows()));
    auto fluxData = std::unique_ptr<flux_data_t>(new flux_data_t(fluxColumn.rows()));
    auto fluxErrData = std::unique_ptr<flux_data_t>(new flux_data_t(fluxErrColumn.rows()));
    
    timeColumn.read<time_data_t::value_type>(*timeData, 0, timeData->size());
    fluxColumn.read<flux_data_t::value_type>(*fluxData, 0, fluxData->size());
    fluxErrColumn.read<flux_data_t::value_type>(*fluxErrData, 0, fluxErrData->size());
    
    float jdRef(0);
    float const secondsPerDay = 86400;
    fitsFile->pHDU().readKey("JD_REF", jdRef);
    
    std::transform(timeData->begin(), timeData->end(), hjdData->begin(), [secondsPerDay, jdRef](auto elem){
        return static_cast<double>(elem)/secondsPerDay + jdRef;
    });
    
    return std::make_tuple(std::move(timeData),
                           std::move(hjdData),
                           std::move(fluxData),
                           std::move(fluxErrData)
                           );
}

int main(int argc, const char * argv[]) {
    auto results = testFitsLoad("/Users/hugh.dickinson/Documents/Development/SuperWASPSLSim/SuperWASPSLSim/test_data/1SWASPJ000000.15320847.6.fits");
    
    auto zippedBegin = boost::make_zip_iterator(boost::make_tuple(std::get<0>(results)->begin(),
                                                                  std::get<1>(results)->begin(),
                                                                  std::get<2>(results)->begin(),
                                                                  std::get<3>(results)->begin()
                                                                  ));
    auto zippedEnd = boost::make_zip_iterator(boost::make_tuple(std::get<0>(results)->end(),
                                                                std::get<1>(results)->end(),
                                                                std::get<2>(results)->end(),
                                                                std::get<3>(results)->end()
                                                                ));
    std::for_each(zippedBegin, zippedEnd, [](auto elem){
        std::cout << elem.template get<0>() << "\t" << elem.template get<1>() <<  "\t" << elem.template get<2>() <<  "\t" << elem.template get<3>() << std::endl;
    });
    
    return 0;
}
