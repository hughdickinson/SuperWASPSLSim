//
//  NonUniformCircularSource.cpp
//  SuperWASPSLSim
//
//  Created by Hugh Dickinson on 15/04/2021.
//

#include "NonUniformCircularSource.h"
#include "UniformCircularSource.h"

#include <cmath>

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/differentiation/finite_difference.hpp>

NonUniformCircularSource::NonUniformCircularSource(double radiusUnitless, double uLambda) :
UniformCircularSource(radiusUnitless),
uLambda(uLambda)
{
}

NonUniformCircularSource::~NonUniformCircularSource(){
}

double NonUniformCircularSource::profileIntegral(){
    double error(0);
    double result =  2 * M_PI * boost::math::quadrature::gauss_kronrod<double, 15>::integrate([this](double r){
        return r*this->radialProfile(r);
    }, 0.0, this->radiusUnitless, 15, 1e-4, &error);
    
    return result;
}

double NonUniformCircularSource::magnificationDerivative(double radialCoordinate, double sourcePlaneCoordinate){
    if(radialCoordinate < 0.0) {
        return 0;
    }
    return boost::math::differentiation::finite_difference_derivative([this, sourcePlaneCoordinate](double r){
        if(r < 0.0) {
            return 0.0;
        }
        return this->UniformCircularSource::magnification(sourcePlaneCoordinate);
    }, radialCoordinate);
}

double NonUniformCircularSource::magnificationProfileIntegral(double sourcePlaneCoordinate){
    double error(0);
    auto result =  2.0*M_PI*boost::math::quadrature::gauss_kronrod<double, 15>::integrate([this, sourcePlaneCoordinate](double r){
        return r * this->radialProfile(r) * ( this->UniformCircularSource::magnification(sourcePlaneCoordinate) + 0.5*r* this->magnificationDerivative(r, sourcePlaneCoordinate));
    }, 0.0, this->radiusUnitless, 15, 1e-4, &error);
    return result;
}

double NonUniformCircularSource::magnification(double sourcePlaneCoordinate){
    return this->magnificationProfileIntegral(sourcePlaneCoordinate) / this->profileIntegral();
}

double NonUniformCircularSource::radialProfile(double radialCoordinate){  // 0.6 => V-band
    return 1.0 - this->uLambda + this->uLambda * std::sqrt(1.0 - std::pow((radialCoordinate / this->radiusUnitless), 2));
}
