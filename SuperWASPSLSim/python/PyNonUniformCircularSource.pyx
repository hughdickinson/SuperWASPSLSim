# distutils: language = c++

from NonUniformCircularSource cimport NonUniformCircularSource

cdef class PyNonUniformCircularSource:
    cdef NonUniformCircularSource * c_non_uniform_circular_source

    def __cinit__(self, double radiusUnitless, double uLambda=0.6, double u2=float('nan')):
        self.c_non_uniform_circular_source = new NonUniformCircularSource(radiusUnitless, uLambda, u2)
        
    def __dealloc__(self):
        del self.c_non_uniform_circular_source

    def magnification(self, sourcePlaneCoordinate):
        return self.c_non_uniform_circular_source.magnification(sourcePlaneCoordinate)
        
    def magnifications(self, sourcePlaneCoordinates):
        return self.c_non_uniform_circular_source.magnifications(sourcePlaneCoordinates)
