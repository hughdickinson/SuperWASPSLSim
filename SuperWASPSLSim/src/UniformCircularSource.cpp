//
//  UniformCircularSource.cpp
//  SuperWASPSLSim
//
//  Created by Hugh Dickinson on 15/04/2021.
//

#include "UniformCircularSource.h"

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>

#include <cmath>

UniformCircularSource::UniformCircularSource(double radiusUnitless) :
radiusUnitless(radiusUnitless),
radiusUnitlessSquared(radiusUnitless * radiusUnitless){
}

UniformCircularSource::~UniformCircularSource(){
}

double UniformCircularSource::magnificationAtRadius(){
    
    auto kernel = (2 / this->radiusUnitless) + ((1 + this->radiusUnitlessSquared) / this->radiusUnitlessSquared) * ((M_PI / 2) + std::asin((this->radiusUnitlessSquared - 1) / (this->radiusUnitlessSquared + 1)));
    
    auto total = kernel / M_PI;
    
    return total;
}

double UniformCircularSource::magnification(double sourcePlaneCoordinate){
    
    auto kernel1 = (sourcePlaneCoordinate - radiusUnitless);
    kernel1 *= kernel1;
    auto kernel2 = std::sqrt(4.0 + kernel1);
    
    auto ellipticN = (
                      4.0
                      * this->radiusUnitless
                      * sourcePlaneCoordinate
                      / std::pow((sourcePlaneCoordinate + this->radiusUnitless), 2.0)
                      );
    
    auto ellipticK = std::sqrt(4.0 * ellipticN) / kernel2;
    
    auto firstTerm = ((sourcePlaneCoordinate + this->radiusUnitless) * kernel2) / (
                                                                                   2.0 * std::pow(this->radiusUnitless, 2.0));
    auto secondTerm = (
                       (sourcePlaneCoordinate - this->radiusUnitless)
                       * (4.0 + (0.5 * (std::pow(sourcePlaneCoordinate, 2.0) - this->radiusUnitlessSquared)))
                       / (kernel2 * this->radiusUnitlessSquared)
                       );
    auto thirdTerm = (
                      2.0
                      * kernel1
                      * (1.0 + this->radiusUnitlessSquared)
                      / (this->radiusUnitlessSquared * (sourcePlaneCoordinate + this->radiusUnitless) * kernel2)
                      );
    
    auto kernel3 = (
                    boost::math::ellint_2(ellipticK) * firstTerm
                    - boost::math::ellint_1(ellipticK) * secondTerm
                    + boost::math::ellint_3(ellipticK, ellipticN) * thirdTerm
                    );
    
    auto positiveSolutionReal=((kernel3 + M_PI) / (2.0 * M_PI));
    auto negativeSolutionReal=((kernel3 - M_PI) / (2.0 * M_PI));
    auto total = positiveSolutionReal + negativeSolutionReal;
    return total;
}

