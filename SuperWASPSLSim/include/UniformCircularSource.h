//
//  UniformCircularSource.h
//  SuperWASPSLSim
//
//  Created by Hugh Dickinson on 15/04/2021.
//

#ifndef UniformCircularSource_h
#define UniformCircularSource_h

class UniformCircularSource {
    
private:
    double magnificationAtRadius();
    
protected:
    double radiusUnitless;
    double radiusUnitlessSquared;
    
public:
    
    UniformCircularSource(double unitlessRadius);
    virtual ~UniformCircularSource();
    
    double magnification(double sourcePlaneCoordinate);
};

#endif /* UniformCircularSource_h */
