//
//  NonUniformCircularSource.h
//  SuperWASPSLSim
//
//  Created by Hugh Dickinson on 15/04/2021.
//

#ifndef NonUniformCircularSource_h
#define NonUniformCircularSource_h

#include "UniformCircularSource.h"

class NonUniformCircularSource : public UniformCircularSource {
    
private:
    
    double uLambda;
    
    double radialProfile(double radialCoordinate);
    double profileIntegral();
    double magnificationProfileIntegral(double sourcePlaneCoordinate);
    double magnificationDerivative(double radialCoordinate, double sourcePlaneCoordinate);
    
public:
    NonUniformCircularSource(double radiusUnitless, double uLambda=0.6);
    ~NonUniformCircularSource();
    double magnification(double sourcePlaneCoordinate);
    
};

#endif /* NonUniformCircularSource_h */
