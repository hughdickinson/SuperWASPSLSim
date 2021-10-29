//
//  EinsteinRadius.h
//  SuperWASPSLSim
//
//  Created by Hugh Dickinson on 15/04/2021.
//

#ifndef EinsteinRadius_h
#define EinsteinRadius_h

#include <memory>

#include <boost/units/base_units/astronomical/parsec.hpp>
#include <boost/units/static_rational.hpp>
#include <boost/units/make_scaled_unit.hpp>
#include <boost/units/systems/si/codata_constants.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/mass.hpp>
#include <boost/units/cmath.hpp>
#include <boost/units/io.hpp>

using namespace boost::units;

typedef unit<length_dimension, make_system<astronomical::parsec_base_unit>::type> parsec_unit;
BOOST_UNITS_STATIC_CONSTANT(parsec,parsec_unit);

typedef make_scaled_unit<parsec_unit, scale<10, static_rational<3> > >::type kiloparsec_unit;
BOOST_UNITS_STATIC_CONSTANT(kiloparsec,kiloparsec_unit);

template <typename EinsteinRadiusUnit, typename SourceToObsUnit, typename LensToObsUnit, typename SourceToLensUnit, typename CompactObjectMassUnit>
auto sourcePlaneEinsteinRadius(
                               quantity<SourceToObsUnit> const & sourceToObsDist,
                               quantity<LensToObsUnit> const & lensToObsDist,
                               quantity<SourceToLensUnit> const & sourceToLensDist,
                               quantity<CompactObjectMassUnit> const & lensMass){
    auto G = boost::units::si::constants::codata::G;
    auto C = boost::units::si::constants::codata::c;
    auto einsteinRadius = boost::units::sqrt((4.0 * G * lensMass * sourceToLensDist * lensToObsDist)
                                             / (sourceToObsDist * C * C));
    
    return std::unique_ptr<quantity<EinsteinRadiusUnit> >(new quantity< EinsteinRadiusUnit>(einsteinRadius));
    
}

template <typename EinsteinRadiusUnit, typename SourceToObsUnit, typename LensToObsUnit, typename SourceToLensUnit, typename CompactObjectMassUnit>
auto lensPlaneEinsteinRadius(
                             quantity<SourceToObsUnit> const & sourceToObsDist,
                             quantity<LensToObsUnit> const & lensToObsDist,
                             quantity<SourceToLensUnit> const & sourceToLensDist,
                             quantity<CompactObjectMassUnit> const & lensMass){
    
    auto einsteinRadius = (sourceToObsDist / lensToObsDist) * *sourcePlaneEinsteinRadius<EinsteinRadiusUnit>(sourceToObsDist,lensToObsDist, sourceToLensDist, lensMass);
    
    return std::unique_ptr<quantity<EinsteinRadiusUnit> >(new quantity< EinsteinRadiusUnit>(einsteinRadius));
    
}

#endif /* EinsteinRadius_h */
