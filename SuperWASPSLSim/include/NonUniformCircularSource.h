//
//  NonUniformCircularSource.h
//  SuperWASPSLSim
//
//  Created by Hugh Dickinson on 15/04/2021.
//

#ifndef NonUniformCircularSource_h
#define NonUniformCircularSource_h

#include <vector>

#include "UniformCircularSource.h"

class NonUniformCircularSource : public UniformCircularSource {
    
private:
    std::vector<double> ldCoeffs;
    
    double radialProfile(double radialCoordinate);
    double profileIntegral();
    double magnificationProfileIntegral(double sourcePlaneCoordinate);
    double magnificationDerivative(double radialCoordinate, double sourcePlaneCoordinate);
    
public:
    NonUniformCircularSource(double radiusUnitless, double uLambda=0.6, std::vector<double> ldCoeffs={});
    ~NonUniformCircularSource();
    double magnification(double sourcePlaneCoordinate);
    std::vector<double> magnifications(std::vector<double> const & sourcePlaneCoordinates);
    
};

#endif /* NonUniformCircularSource_h */
