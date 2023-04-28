//
//  NonUniformCircularSource.cpp
//  SuperWASPSLSim
//
//  Created by Hugh Dickinson on 15/04/2021.
//

#include "NonUniformCircularSource.h"
#include "UniformCircularSource.h"

#include <cmath>
#include <limits>

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/differentiation/finite_difference.hpp>

NonUniformCircularSource::NonUniformCircularSource(double radiusUnitless, double u1=0.6, double u2=std::numeric_limits<double>::quiet_NaN()) :
UniformCircularSource(radiusUnitless)
{
    ldCoeffs = {u1};
    if (!std::isnan(u2)) {
        ldCoeffs.push_back(u2);
    }
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

std::vector<double> NonUniformCircularSource::magnifications(std::vector<double> const & sourcePlaneCoordinates){
    std::vector<double> res(sourcePlaneCoordinates.size());
    std::transform(sourcePlaneCoordinates.begin(), sourcePlaneCoordinates.end(), res.begin(), [this](double sourcePlaneCoordinate){
        return this->magnification(sourcePlaneCoordinate);
    });
    return res;
}


double NonUniformCircularSource::radialProfile(double radialCoordinate){
    double radiusRatio = (radialCoordinate / this->radiusUnitless);
    double mu = std::sqrt(1.0 - std::pow(radiusRatio, 2)); // mu=cos(theta), theta=arcsin(radiusRatio)

    switch(this->ldCoeffs.size()) {
        case 2:
            // Quadratic law
            // Zdenek Kopal, ‘Detailed Effects of Limb Darkening upon Light and Velocity Curves of Close Binary Systems’, Harvard College Observatory Circular 454 (1 January 1950): 1–12.            
            return 1.0 - this->ldCoeffs[0] * (1.0 - mu) - this->ldCoeffs[1] * std::pow((1.0 - mu), 2);

        case 1: 
            // Linear law
            // 0.6 => V-band
            // C. W. Allen, Astrophysical Quantities, London: University of London, 1973, https://ui.adsabs.harvard.edu/abs/1973asqu.book.....A.
            // via Hans J. Witt and Shude Mao, ‘Can Lensed Stars Be Regarded as Pointlike for Microlensing by MACHOs?’, The Astrophysical Journal 430 (1 August 1994): 505, https://doi.org/10.1086/174426.
            // Actually from K. Schwarzschild, ‘Ueber das Gleichgewicht der Sonnenatmosphäre’, Nachrichten von der Gesellschaft der Wissenschaften zu Göttingen, Mathematisch-Physikalische Klasse 1906 (1906): 41–53.
            return 1.0 - this->ldCoeffs[0] * (1.0 - mu);

        default:
            throw std::logic_error("Unsupported number of coefficients for limb-darkening law");
    }
}
