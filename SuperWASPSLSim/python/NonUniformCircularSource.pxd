from libcpp.vector cimport vector

cdef extern from "../src/NonUniformCircularSource.cpp":
    pass

cdef extern from "../src/UniformCircularSource.cpp":
    pass

# Declare the class with cdef
cdef extern from "../include/NonUniformCircularSource.h":
    cdef cppclass NonUniformCircularSource:
        NonUniformCircularSource(double, double, double) except +
        double magnification(double)
        vector[double] magnifications(vector[double])
        
cdef extern from "../include/UniformCircularSource.h":
    cdef cppclass UniformCircularSource:
        UniformCircularSource(double, double) except +
        double magnification(double)

