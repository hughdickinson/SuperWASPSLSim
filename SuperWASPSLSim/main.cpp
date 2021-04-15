//
//  main.cpp
//  SuperWASPSLSim
//
//  Created by Hugh Dickinson on 14/04/2021.
//

#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <array>
#include <algorithm>
#include <iterator>

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


int main(int argc, const char * argv[]) {
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
