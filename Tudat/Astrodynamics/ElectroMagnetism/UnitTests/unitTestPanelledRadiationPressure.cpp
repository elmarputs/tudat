#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"

#include "Tudat/Astrodynamics/ElectroMagnetism/panelledRadiationPressure.h"
#include "Tudat/Astrodynamics/Ephemerides/keplerEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/constantRotationalEphemeris.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/PropagationSetup/createAccelerationModels.h"

namespace tudat
{

using namespace tudat::basic_astrodynamics;
using namespace tudat::simulation_setup;
using namespace tudat::ephemerides;
using namespace tudat::basic_astrodynamics;
using namespace tudat::electro_magnetism;

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_PanelledRadiationPressure )

BOOST_AUTO_TEST_CASE( testSimpleGeometryPanelledRadiationPressure )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create bodies needed in simulation
    double initialEphemerisTime = -86400.0;
    double finalEphemerisTime = 1.1 * 365.25 * 86400.0;
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Sun" );
    NamedBodyMap bodyMap = createBodies(
                getDefaultBodySettings( bodyNames,initialEphemerisTime, finalEphemerisTime ) );

    // Create vehicle
    double vehicleMass = 2500.0;
    bodyMap[ "Vehicle" ] = std::make_shared< Body >( );
    bodyMap[ "Vehicle" ]->setConstantBodyMass( vehicleMass );

    // Put vehicle on circular orbit around Sun
    Eigen::Vector6d initialStateInKeplerianElements = Eigen::Vector6d::Zero( );
    initialStateInKeplerianElements[ 0 ] = physical_constants::ASTRONOMICAL_UNIT;
    bodyMap[ "Vehicle" ]->setEphemeris(
                std::make_shared< KeplerEphemeris >(
                    initialStateInKeplerianElements, 0.0, spice_interface::getBodyGravitationalParameter( "Sun" ),
                    "Sun", "ECLIPJ2000", 1 ) );
    Eigen::Vector7d rotationalStateVehicle;
    rotationalStateVehicle.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond( Eigen::Matrix3d::Identity() ));
    rotationalStateVehicle.segment( 4, 3 ) = Eigen::Vector3d::Zero();
//    bodyMap[ "Vehicle" ]->setRotationalEphemeris(
//                std::make_shared< ConstantRotationalEphemeris >(
//                    Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ), "ECLIPJ2000", "VehicleFixed" ) );
    bodyMap[ "Vehicle" ]->setRotationalEphemeris(
                std::make_shared< ConstantRotationalEphemeris >(
                    rotationalStateVehicle, "ECLIPJ2000", "VehicleFixed" ) );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    for( unsigned int test = 0; test <= 2; test++ )
    {
        // Create radiation pressure properties
        std::vector< double > areas;
        areas.push_back( 2.532 );
        areas.push_back( 3.254 );
        areas.push_back( 8.654 );
        areas.push_back( 1.346 );
        areas.push_back( 2.454 );
        areas.push_back( 5.345 );

        std::vector< double > emissivities;
        if( test == 0 )
        {
            emissivities.push_back( 0.0 );
            emissivities.push_back( 1.0 );
            emissivities.push_back( 0.0 );
            emissivities.push_back( 1.0 );
        }
        else if( test == 1 )
        {
            emissivities.push_back( 0.2 );
            emissivities.push_back( 0.3 );
            emissivities.push_back( 0.4 );
            emissivities.push_back( 0.5 );
        }
        emissivities.push_back( 0.5 );
        emissivities.push_back( 0.2 );

        std::vector< double > diffuseReflectionCoefficients;
        if( test == 0 )
        {
            diffuseReflectionCoefficients.push_back( 0.0 );
            diffuseReflectionCoefficients.push_back( 0.0 );
            diffuseReflectionCoefficients.push_back( 0.0 );
            diffuseReflectionCoefficients.push_back( 0.0 );
            diffuseReflectionCoefficients.push_back( 0.0 );
            diffuseReflectionCoefficients.push_back( 0.0 );
        }
        else
        {
            diffuseReflectionCoefficients.push_back( 0.6 );
            diffuseReflectionCoefficients.push_back( 0.5 );
            diffuseReflectionCoefficients.push_back( 0.4 );
            diffuseReflectionCoefficients.push_back( 0.3 );
            diffuseReflectionCoefficients.push_back( 0.2 );
            diffuseReflectionCoefficients.push_back( 0.1 );
        }

        if( test == 2 )
        {
            bodyMap[ "Vehicle" ]->setRotationalEphemeris(
                        std::make_shared< tudat::ephemerides::SimpleRotationalEphemeris >(
                            0.2, 0.4, -0.2, 1.0E-5, 0.0, "ECLIPJ2000", "VehicleFixed" ) );
        }

        std::vector< Eigen::Vector3d > panelSurfaceNormals;
        panelSurfaceNormals.push_back( -Eigen::Vector3d::UnitX( ) );
        panelSurfaceNormals.push_back( -Eigen::Vector3d::UnitY( ) );
        panelSurfaceNormals.push_back( Eigen::Vector3d::UnitX( ) );
        panelSurfaceNormals.push_back( Eigen::Vector3d::UnitY( ) );
        panelSurfaceNormals.push_back( Eigen::Vector3d::UnitZ( ) );
        panelSurfaceNormals.push_back( -Eigen::Vector3d::UnitZ( ) );

        std::shared_ptr< PanelledRadiationPressureInterfaceSettings > radiationPressureInterfaceSettings =
                std::make_shared< PanelledRadiationPressureInterfaceSettings >(
                    "Sun", emissivities, areas, diffuseReflectionCoefficients, panelSurfaceNormals );
        std::shared_ptr< PanelledRadiationPressureInterface > radiationPressureInterface =
                std::dynamic_pointer_cast< PanelledRadiationPressureInterface >(
                    createRadiationPressureInterface( radiationPressureInterfaceSettings, "Vehicle", bodyMap ) );
        bodyMap[ "Vehicle" ]->setRadiationPressureInterface( "Sun", radiationPressureInterface );

        // Define accelerations
        SelectedAccelerationMap accelerationMap;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOnVehicle;
        accelerationsOnVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                       panelled_radiation_pressure_acceleration ) );
        accelerationMap[ "Vehicle" ] = accelerationsOnVehicle;
        std::map< std::string, std::string > centralBodyMap;
        centralBodyMap[ "Vehicle" ] = "Sun";
        AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodyMap, accelerationMap, centralBodyMap );
        std::shared_ptr< AccelerationModel< Eigen::Vector3d > > accelerationModel =
                accelerationModelMap.at( "Vehicle" ).at( "Sun" ).at( 0 );

        // Compute spacecraft orbital period, and compute test times
        double orbitalPeriod =
                2.0 * mathematical_constants::PI * std::sqrt(
                    std::pow( physical_constants::ASTRONOMICAL_UNIT, 3.0 ) /
                    spice_interface::getBodyGravitationalParameter( "Sun" ) );
        std::vector< double > testTimes = { 0.0, orbitalPeriod / 4.0, orbitalPeriod / 2.0, 3.0 * orbitalPeriod / 4.0 };

        // Compute panelled radiation pressure for various relative Sun positions
        Eigen::Vector3d calculatedAcceleration, expectedAcceleration;
        Eigen::Vector3d sunCenteredVehiclePosition;
        std::shared_ptr< Ephemeris > vehicleEphemeris = bodyMap[ "Vehicle" ]->getEphemeris( );
        for( unsigned int i = 0; i < testTimes.size( ); i++ )
        {
            // Update environment and acceleration
            bodyMap[ "Sun" ]->setStateFromEphemeris( testTimes[ i ] );
            bodyMap[ "Vehicle" ]->setStateFromEphemeris( testTimes[ i ] );
            bodyMap[ "Vehicle" ]->setCurrentRotationToLocalFrameFromEphemeris( testTimes[ i ] );
            radiationPressureInterface->updateInterface( testTimes[ i ] );
            accelerationModel->updateMembers( testTimes[ i ] );

            // Retrieve acceleration
            calculatedAcceleration = accelerationModel->getAcceleration( );

            // Manually compute acceleration
            double radiationPressure = radiationPressureInterface->getCurrentRadiationPressure( );
            sunCenteredVehiclePosition = vehicleEphemeris->getCartesianState( testTimes[ i ] ).segment( 0, 3 );
            if( test == 0 )
            {
                if( i == 0 || i == 2 )
                {
                    expectedAcceleration = radiationPressure / vehicleMass  *
                            radiationPressureInterface->getArea( i ) * sunCenteredVehiclePosition.normalized( );
                }
                else if( i == 1 || i == 3 )
                {
                    expectedAcceleration = -radiationPressure / vehicleMass * (
                                2.0 * radiationPressureInterface->getArea( i ) *
                                radiationPressureInterface->getCurrentSurfaceNormal( i ) );
                }
            }
            else
            {
                expectedAcceleration = radiationPressure / vehicleMass *
                        ( 1.0 + radiationPressureInterface->getEmissivity( i ) + 2.0 * diffuseReflectionCoefficients.at( i ) / 3.0 ) *
                        radiationPressureInterface->getArea( i ) * sunCenteredVehiclePosition.normalized( );
            }

            if( test == 0 || test == 1 )
            {
                // Check computed acceleration
                for( unsigned int j = 0; j < 3; j++ )
                {
                    BOOST_CHECK_SMALL( std::fabs( calculatedAcceleration( j ) - expectedAcceleration( j ) ), 2.0E-23 );
                }
            }
            else
            {
                std::shared_ptr< PanelledRadiationPressureAcceleration > panelledRadiationPressureAcceleration =
                        std::dynamic_pointer_cast< PanelledRadiationPressureAcceleration >( accelerationModel );
                Eigen::Quaterniond currentRotationToInertialFrame =
                        bodyMap[ "Vehicle" ]->getRotationalEphemeris( )->getRotationToBaseFrame( testTimes.at( i ) );
                for( unsigned int j = 0; j < panelSurfaceNormals.size( ); j++ )
                {
                   Eigen::Vector3d currentPanelNormal =
                           panelledRadiationPressureAcceleration->getCurrentPanelSurfaceNormalInPropagationFrame( j );
                   Eigen::Vector3d expectedPanelNormal =
                           currentRotationToInertialFrame * panelSurfaceNormals.at( j );
                   for( unsigned int k = 0; k < 3; k++ )
                   {
                       BOOST_CHECK_SMALL( std::fabs( currentPanelNormal( k ) - expectedPanelNormal( k ) ), 1.0E-15 );
                   }
                }
            }
        }
    }
}


BOOST_AUTO_TEST_CASE( testPanelledRadiationPressureMontenbruckData )
{

    // Box-and-wings model is obtained from Montenbruck, O., 2017.
    // Semi-analytical solar radiation pressure modeling for QZS-1 orbit-normal and yaw-steering attitude.
    // Advances in Space Research, 59(8), pp.2088-2100..

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create bodies needed in simulation
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = 1.1 * 365.25 * 86400.0;
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Sun" );
    NamedBodyMap bodyMap = createBodies(
                getDefaultBodySettings( bodyNames,initialEphemerisTime, finalEphemerisTime ) );

    // Create vehicle
    double vehicleMass = 2000.0;
    bodyMap[ "Vehicle" ] = std::make_shared< Body >( );
    bodyMap[ "Vehicle" ]->setConstantBodyMass( vehicleMass );


    for ( int testCase = 0 ; testCase < 4 ; testCase++){

        // Put vehicle on circular orbit around Sun

        if ( testCase == 0 || testCase == 3 ){
            Eigen::Vector6d initialStateInKeplerianElements = Eigen::Vector6d::Zero( );
            initialStateInKeplerianElements[ 0 ] = physical_constants::ASTRONOMICAL_UNIT;
            bodyMap[ "Vehicle" ]->setEphemeris( std::make_shared< KeplerEphemeris >( initialStateInKeplerianElements, 0.0,
                                                         spice_interface::getBodyGravitationalParameter( "Sun" ), "Sun", "ECLIPJ2000", 1 ) );
        }
        else if ( testCase == 1 ){
            Eigen::Vector6d initialStateInKeplerianElements = Eigen::Vector6d::Zero( );
            initialStateInKeplerianElements[ 0 ] = physical_constants::ASTRONOMICAL_UNIT;
            initialStateInKeplerianElements[ orbital_element_conversions::inclinationIndex ] = unit_conversions::convertDegreesToRadians( 90.0 );
            bodyMap[ "Vehicle" ]->setEphemeris( std::make_shared< KeplerEphemeris >( initialStateInKeplerianElements, 0.0,
                                                         spice_interface::getBodyGravitationalParameter( "Sun" ), "Sun", "ECLIPJ2000", 1 ) );
        }
        else if ( testCase == 2 ){
            Eigen::Vector6d initialStateInKeplerianElements = Eigen::Vector6d::Zero( );
            initialStateInKeplerianElements[ 0 ] = physical_constants::ASTRONOMICAL_UNIT;
            initialStateInKeplerianElements[ orbital_element_conversions::inclinationIndex ] = unit_conversions::convertDegreesToRadians( 20.0 );
            bodyMap[ "Vehicle" ]->setEphemeris( std::make_shared< KeplerEphemeris >( initialStateInKeplerianElements, 0.0,
                                                         spice_interface::getBodyGravitationalParameter( "Sun" ), "Sun", "ECLIPJ2000", 1 ) );
        }

        // Setup rotational ephemeris for vehicle.
        if ( testCase < 3 ){
            Eigen::Vector7d rotationalStateVehicle;
            rotationalStateVehicle.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond( Eigen::Matrix3d::Identity() ));
            rotationalStateVehicle.segment( 4, 3 ) = Eigen::Vector3d::Zero();
            bodyMap[ "Vehicle" ]->setRotationalEphemeris( std::make_shared< ConstantRotationalEphemeris >( rotationalStateVehicle, "ECLIPJ2000",
                                                                                                           "VehicleFixed" ) );
        }
        else if ( testCase == 3 ){
            bodyMap[ "Vehicle" ]->setRotationalEphemeris( std::make_shared< tudat::ephemerides::SimpleRotationalEphemeris >(
                            0.2, 0.4, -0.2, 1.0E-5, 0.0, "ECLIPJ2000", "VehicleFixed" ) );
        }

        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );


        // Create radiation pressure properties.
        std::vector< double > areas;
        if ( testCase < 3 ){
            areas.push_back( 2.0 );
            areas.push_back( 4.0 );
            areas.push_back( 6.0 );
            areas.push_back( 9.9 );
            areas.push_back( 2.3 );
            areas.push_back( 9.9 );
            areas.push_back( 2.3 );
            areas.push_back( 4.6 );
            areas.push_back( 2.7 );
            areas.push_back( 5.8 );
            areas.push_back( 2.7 );
        }
        else if ( testCase == 3 ){
            areas.push_back( 9.9 );
            areas.push_back( 2.3 );
            areas.push_back( 9.9 );
            areas.push_back( 2.3 );
        }

        std::vector< double > emissivities;
        if ( testCase < 3 ){
            emissivities.push_back( 0.0 );
            emissivities.push_back( 0.1 );
            emissivities.push_back( 0.0 );
            emissivities.push_back( 0.0 );
            emissivities.push_back( 0.1 );
            emissivities.push_back( 0.0 );
            emissivities.push_back( 0.1 );
            emissivities.push_back( 0.0 );
            emissivities.push_back( 0.1 );
            emissivities.push_back( 0.0 );
            emissivities.push_back( 0.1 );
        }
        else if ( testCase == 3 ){
            emissivities.push_back( 0.0 );
            emissivities.push_back( 0.1 );
            emissivities.push_back( 0.0 );
            emissivities.push_back( 0.1 );
        }

        std::vector< double > diffuseReflectionCoefficients;
        if ( testCase < 3 ){
            diffuseReflectionCoefficients.push_back( 0.06 );
            diffuseReflectionCoefficients.push_back( 0.46 );
            diffuseReflectionCoefficients.push_back( 0.06 );
            diffuseReflectionCoefficients.push_back( 0.06 );
            diffuseReflectionCoefficients.push_back( 0.46 );
            diffuseReflectionCoefficients.push_back( 0.06 );
            diffuseReflectionCoefficients.push_back( 0.46 );
            diffuseReflectionCoefficients.push_back( 0.06 );
            diffuseReflectionCoefficients.push_back( 0.46 );
            diffuseReflectionCoefficients.push_back( 0.06 );
            diffuseReflectionCoefficients.push_back( 0.46 );
        }
        else if ( testCase == 3 ){
            diffuseReflectionCoefficients.push_back( 0.06 );
            diffuseReflectionCoefficients.push_back( 0.46 );
            diffuseReflectionCoefficients.push_back( 0.06 );
            diffuseReflectionCoefficients.push_back( 0.46 );
        }

        std::vector< Eigen::Vector3d > panelSurfaceNormals;
        if ( testCase < 3 ){
            panelSurfaceNormals.push_back( Eigen::Vector3d::UnitZ( ) );
            panelSurfaceNormals.push_back( Eigen::Vector3d::UnitZ( ) );
            panelSurfaceNormals.push_back( - Eigen::Vector3d::UnitZ( ) );
            panelSurfaceNormals.push_back( Eigen::Vector3d::UnitX( ) );
            panelSurfaceNormals.push_back( Eigen::Vector3d::UnitX( ) );
            panelSurfaceNormals.push_back( - Eigen::Vector3d::UnitX( ) );
            panelSurfaceNormals.push_back( - Eigen::Vector3d::UnitX( ) );
            panelSurfaceNormals.push_back( Eigen::Vector3d::UnitY( ) );
            panelSurfaceNormals.push_back( Eigen::Vector3d::UnitY( ) );
            panelSurfaceNormals.push_back( - Eigen::Vector3d::UnitY( ) );
            panelSurfaceNormals.push_back( - Eigen::Vector3d::UnitY( ) );
        }
        else if ( testCase == 3 ){
            panelSurfaceNormals.push_back( Eigen::Vector3d::UnitX( ) );
            panelSurfaceNormals.push_back( Eigen::Vector3d::UnitX( ) );
            panelSurfaceNormals.push_back( - Eigen::Vector3d::UnitX( ) );
            panelSurfaceNormals.push_back( - Eigen::Vector3d::UnitX( ) );
        }


        // Create panelled radiation pressure interface.
        std::shared_ptr< PanelledRadiationPressureInterfaceSettings > radiationPressureInterfaceSettings =
                std::make_shared< PanelledRadiationPressureInterfaceSettings >(
                    "Sun", emissivities, areas, diffuseReflectionCoefficients, panelSurfaceNormals );

        std::shared_ptr< PanelledRadiationPressureInterface > radiationPressureInterface =
                std::dynamic_pointer_cast< PanelledRadiationPressureInterface >(
                    createRadiationPressureInterface( radiationPressureInterfaceSettings, "Vehicle", bodyMap ) );

        bodyMap[ "Vehicle" ]->setRadiationPressureInterface( "Sun", radiationPressureInterface );



        // Define accelerations.
        SelectedAccelerationMap accelerationMap;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOnVehicle;
        accelerationsOnVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                       panelled_radiation_pressure_acceleration ) );

        accelerationMap[ "Vehicle" ] = accelerationsOnVehicle;

        std::map< std::string, std::string > centralBodyMap; centralBodyMap[ "Vehicle" ] = "Sun";
        AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodyMap, accelerationMap, centralBodyMap );
        std::shared_ptr< AccelerationModel< Eigen::Vector3d > > accelerationModel =
                accelerationModelMap.at( "Vehicle" ).at( "Sun" ).at( 0 );


        // Compute spacecraft orbital period, and compute test times
        double orbitalPeriod = 2.0 * mathematical_constants::PI * std::sqrt( std::pow( physical_constants::ASTRONOMICAL_UNIT, 3.0 ) /
                                                                             spice_interface::getBodyGravitationalParameter( "Sun" ) );
        std::vector< double > testTimes = { 0.0, orbitalPeriod / 4.0, orbitalPeriod / 2.0, 3.0 * orbitalPeriod / 4.0 };


        // Compute panelled radiation pressure for various relative Sun positions.
        Eigen::Vector3d calculatedAcceleration, expectedAcceleration;

        Eigen::Vector3d sunCenteredVehiclePosition;
        std::shared_ptr< Ephemeris > vehicleEphemeris = bodyMap[ "Vehicle" ]->getEphemeris( );


        for( unsigned int i = 0; i < testTimes.size( ) ; i++ )
        {
            // Update environment and acceleration
            bodyMap[ "Sun" ]->setStateFromEphemeris( testTimes[ i ] );
            bodyMap[ "Vehicle" ]->setStateFromEphemeris( testTimes[ i ] );
            bodyMap[ "Vehicle" ]->setCurrentRotationToLocalFrameFromEphemeris( testTimes[ i ] );
            radiationPressureInterface->updateInterface( testTimes[ i ] );
            accelerationModel->updateMembers( testTimes[ i ] );

            // Retrieve acceleration.
            calculatedAcceleration = accelerationModel->getAcceleration( );

            double radiationPressure = radiationPressureInterface->getCurrentRadiationPressure( );
            Eigen::Vector3d expectedVehicleToSunNormalisedVector = ( bodyMap[ "Sun" ]->getState( ) - bodyMap[ "Vehicle" ]->getState( ) ).segment( 0, 3 )
                    .normalized();


            if ( testCase == 0 ){
                if ( i == 0 ){
                    // Sun-Vehicle vector along +X-axis.
                    Eigen::Vector3d expectedAccelerationDirection = Eigen::Vector3d::UnitX( );
                    double expectedAccelerationMagnitude = radiationPressure / bodyMap[ "Vehicle" ]->getBodyMass()
                            * ( areas[5] * ( ( 1 - emissivities[5] ) + 2.0 / 3.0 * diffuseReflectionCoefficients[5] )
                            + areas[6] * ( ( 1 - emissivities[6] ) + 2.0 / 3.0 * diffuseReflectionCoefficients[6] + 2 * emissivities[6] ) );
                    expectedAcceleration = expectedAccelerationMagnitude * expectedAccelerationDirection;
                }
                else if ( i == 1 ){
                    // Sun-Vehicle vector along +Y-axis.
                    Eigen::Vector3d expectedAccelerationDirection = Eigen::Vector3d::UnitY( );
                    double expectedAccelerationMagnitude = radiationPressure / bodyMap[ "Vehicle" ]->getBodyMass()
                            * ( areas[9] * ( ( 1 - emissivities[9] ) + 2.0 / 3.0 * diffuseReflectionCoefficients[9] )
                            + areas[10] * ( ( 1 - emissivities[10] ) + 2.0 / 3.0 * diffuseReflectionCoefficients[10] + 2 * emissivities[10] ) );
                    expectedAcceleration = expectedAccelerationMagnitude * expectedAccelerationDirection;
                }
                else if ( i == 2 ){
                    // Sun-Vehicle vector along -X-axis.
                    Eigen::Vector3d expectedAccelerationDirection = - Eigen::Vector3d::UnitX( );
                    double expectedAccelerationMagnitude = radiationPressure / bodyMap[ "Vehicle" ]->getBodyMass()
                            * ( areas[3] * ( ( 1 - emissivities[3] ) + 2.0 / 3.0 * diffuseReflectionCoefficients[3] )
                            + areas[4] * ( ( 1 - emissivities[4] ) + 2.0 / 3.0 * diffuseReflectionCoefficients[4] + 2 * emissivities[4] ) );
                    expectedAcceleration = expectedAccelerationMagnitude * expectedAccelerationDirection;
                }
                else if ( i == 3 ){
                    // Sun-Vehicle vector along -Y-axis.
                    Eigen::Vector3d expectedAccelerationDirection = - Eigen::Vector3d::UnitY( );
                    double expectedAccelerationMagnitude = radiationPressure / bodyMap[ "Vehicle" ]->getBodyMass()
                            * ( areas[7] * ( ( 1 - emissivities[7] ) + 2.0 / 3.0 * diffuseReflectionCoefficients[7] )
                            + areas[8] * ( ( 1 - emissivities[8] ) + 2.0 / 3.0 * diffuseReflectionCoefficients[8] + 2 * emissivities[8] ) );
                    expectedAcceleration = expectedAccelerationMagnitude * expectedAccelerationDirection;
                }
            }


            if ( testCase == 1 ){

                if ( i == 0 ){
                    // Sun-Vehicle vector along +X-axis.
                    Eigen::Vector3d expectedAccelerationDirection = Eigen::Vector3d::UnitX( );
                    double expectedAccelerationMagnitude = radiationPressure / bodyMap[ "Vehicle" ]->getBodyMass()
                            * ( areas[5] * ( ( 1 - emissivities[5] ) + 2.0 / 3.0 * diffuseReflectionCoefficients[5] )
                            + areas[6] * ( ( 1 - emissivities[6] ) + 2.0 / 3.0 * diffuseReflectionCoefficients[6] + 2 * emissivities[6] ) );
                    expectedAcceleration = expectedAccelerationMagnitude * expectedAccelerationDirection;
                }
                else if ( i == 1 ){
                    // Sun-Vehicle vector along +Z-axis.
                    Eigen::Vector3d expectedAccelerationDirection = Eigen::Vector3d::UnitZ( );
                    double expectedAccelerationMagnitude = radiationPressure / bodyMap[ "Vehicle" ]->getBodyMass()
                            * ( areas[2] * ( ( 1 - emissivities[2] ) + 2.0 / 3.0 * diffuseReflectionCoefficients[2] ) );
                    expectedAcceleration = expectedAccelerationMagnitude * expectedAccelerationDirection;
                }
                else if ( i == 2 ){
                    // Sun-Vehicle vector along -X-axis.
                    Eigen::Vector3d expectedAccelerationDirection = - Eigen::Vector3d::UnitX( );
                    double expectedAccelerationMagnitude = radiationPressure / bodyMap[ "Vehicle" ]->getBodyMass()
                            * ( areas[3] * ( ( 1 - emissivities[3] ) + 2.0 / 3.0 * diffuseReflectionCoefficients[3] )
                            + areas[4] * ( ( 1 - emissivities[4] ) + 2.0 / 3.0 * diffuseReflectionCoefficients[4] + 2 * emissivities[4] ) );
                    expectedAcceleration = expectedAccelerationMagnitude * expectedAccelerationDirection;
                }
                else if ( i == 3 ){
                    // Sun-Vehicle vector along -Z-axis.
                    Eigen::Vector3d expectedAccelerationDirection = - Eigen::Vector3d::UnitZ( );
                    double expectedAccelerationMagnitude = radiationPressure / bodyMap[ "Vehicle" ]->getBodyMass()
                            * ( areas[0] * ( ( 1 - emissivities[0] ) + 2.0 / 3.0 * diffuseReflectionCoefficients[0] )
                            + areas[1] * ( ( 1 - emissivities[1] ) + 2.0 / 3.0 * diffuseReflectionCoefficients[1] + 2 * emissivities[1] ) );
                    expectedAcceleration = expectedAccelerationMagnitude * expectedAccelerationDirection;
                }
            }

            if ( testCase == 2 ){

                double cosinusPanelInclinationPositiveYaxis = expectedVehicleToSunNormalisedVector.dot( Eigen::Vector3d::UnitY() );
                double cosinusPanelInclinationNegativeYaxis = expectedVehicleToSunNormalisedVector.dot( - Eigen::Vector3d::UnitY() );
                double cosinusPanelInclinationPositiveZaxis = expectedVehicleToSunNormalisedVector.dot( Eigen::Vector3d::UnitZ() );
                double cosinusPanelInclinationNegativeZaxis = expectedVehicleToSunNormalisedVector.dot( - Eigen::Vector3d::UnitZ() );


                if ( i == 0 ){
                    // Sun-Vehicle vector along +X-axis.
                    Eigen::Vector3d expectedAccelerationDirection = Eigen::Vector3d::UnitX( );
                    double expectedAccelerationMagnitude = radiationPressure / bodyMap[ "Vehicle" ]->getBodyMass()
                            * ( areas[5] * ( ( 1 - emissivities[5] ) + 2.0 / 3.0 * diffuseReflectionCoefficients[5] )
                            + areas[6] * ( ( 1 - emissivities[6] ) + 2.0 / 3.0 * diffuseReflectionCoefficients[6] + 2 * emissivities[6] ) );
                    expectedAcceleration = expectedAccelerationMagnitude * expectedAccelerationDirection;
                }

                else if ( i == 1 ){
                    expectedAcceleration = - radiationPressure / bodyMap[ "Vehicle" ]->getBodyMass()
                            * ( areas[2] * cosinusPanelInclinationNegativeZaxis * ( ( 1 - emissivities[2] ) * expectedVehicleToSunNormalisedVector
                                        + 2.0 / 3.0 * diffuseReflectionCoefficients[2] * panelSurfaceNormals[2] )
                            + areas[9] * cosinusPanelInclinationNegativeYaxis * ( ( 1 - emissivities[9] ) * expectedVehicleToSunNormalisedVector
                                        + 2.0 / 3.0 * diffuseReflectionCoefficients[9] * panelSurfaceNormals[9])
                            + areas[10] * cosinusPanelInclinationNegativeYaxis * ( ( 1 - emissivities[10] ) * expectedVehicleToSunNormalisedVector
                                        + ( 2.0 / 3.0 * diffuseReflectionCoefficients[10]
                                        + 2 * emissivities[10] * cosinusPanelInclinationNegativeYaxis ) * panelSurfaceNormals[10] ) );
                }
                else if ( i == 2 ){
                    // Sun-Vehicle vector along -X-axis.
                    Eigen::Vector3d expectedAccelerationDirection = - Eigen::Vector3d::UnitX( );
                    double expectedAccelerationMagnitude = radiationPressure / bodyMap[ "Vehicle" ]->getBodyMass()
                            * ( areas[3] * ( ( 1 - emissivities[3] ) + 2.0 / 3.0 * diffuseReflectionCoefficients[3] )
                            + areas[4] * ( ( 1 - emissivities[4] ) + 2.0 / 3.0 * diffuseReflectionCoefficients[4] + 2 * emissivities[4] ) );
                    expectedAcceleration = expectedAccelerationMagnitude * expectedAccelerationDirection;
                }

                else if ( i == 3 ){
                    expectedAcceleration = - radiationPressure / bodyMap[ "Vehicle" ]->getBodyMass()
                            * ( areas[0] * cosinusPanelInclinationPositiveZaxis * ( ( 1 - emissivities[0] ) * expectedVehicleToSunNormalisedVector
                                        + 2.0 / 3.0 * diffuseReflectionCoefficients[0] * panelSurfaceNormals[0] )
                            + areas[1] * cosinusPanelInclinationPositiveZaxis * ( ( 1 - emissivities[1] ) * expectedVehicleToSunNormalisedVector
                                        + ( 2.0 / 3.0 * diffuseReflectionCoefficients[1]
                                        + 2 * emissivities[1] * cosinusPanelInclinationPositiveZaxis) * panelSurfaceNormals[1] )
                            + areas[7] * cosinusPanelInclinationPositiveYaxis * ( ( 1 - emissivities[7] ) * expectedVehicleToSunNormalisedVector
                                        + 2.0 / 3.0 * diffuseReflectionCoefficients[7] * panelSurfaceNormals[7])
                            + areas[8] * cosinusPanelInclinationPositiveYaxis * ( ( 1 - emissivities[8] ) * expectedVehicleToSunNormalisedVector
                                        + ( 2.0 / 3.0 * diffuseReflectionCoefficients[8]
                                        + 2 * emissivities[8] * cosinusPanelInclinationPositiveYaxis) * panelSurfaceNormals[8] ) );
                }
            }

            if ( testCase == 3 ){

               Eigen::Quaterniond currentRotationToInertialFrame = bodyMap[ "Vehicle" ]->getRotationalEphemeris( )
                        ->getRotationToBaseFrame( testTimes.at( i ) );
               Eigen::Vector3d expectedPanelNormalPositiveXaxis = currentRotationToInertialFrame * Eigen::Vector3d::UnitX();
               Eigen::Vector3d expectedPanelNormalNegativeXaxis = currentRotationToInertialFrame * - Eigen::Vector3d::UnitX();


               double cosinusPanelInclinationPositiveXaxis = expectedVehicleToSunNormalisedVector.dot( expectedPanelNormalPositiveXaxis );

               Eigen::Vector3d expectedPanelNormal;
               if ( cosinusPanelInclinationPositiveXaxis >= 0 ){
                   expectedPanelNormal = expectedPanelNormalPositiveXaxis;
               }
               else{
                   expectedPanelNormal = expectedPanelNormalNegativeXaxis;
               }

               expectedAcceleration = - radiationPressure / bodyMap[ "Vehicle" ]->getBodyMass()
                       * ( areas[0] * abs(cosinusPanelInclinationPositiveXaxis) * ( ( 1 - emissivities[0] ) * expectedVehicleToSunNormalisedVector
                       + 2.0 / 3.0 * diffuseReflectionCoefficients[0] * expectedPanelNormal )
                       + areas[1] * abs(cosinusPanelInclinationPositiveXaxis) * ( ( 1 - emissivities[1] ) * expectedVehicleToSunNormalisedVector
                       + ( 2.0 / 3.0 * diffuseReflectionCoefficients[1] + 2.0 * abs(cosinusPanelInclinationPositiveXaxis) * emissivities[1] )
                       * expectedPanelNormal ) );
            }


            std::cout << "calculated acceleration: " << calculatedAcceleration << "\n\n";
            std::cout << "expected acceleration: " << expectedAcceleration << "\n\n";


            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( std::fabs( calculatedAcceleration( j ) - expectedAcceleration( j ) ), 3.0E-23 );
            }

            }
    }

}



BOOST_AUTO_TEST_CASE( testPanelledRadiationPressureTimeVaryingPanelOrientation )
{

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create bodies needed in simulation
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = 1.1 * 365.25 * 86400.0;
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Sun" );
    NamedBodyMap bodyMap = createBodies(
                getDefaultBodySettings( bodyNames,initialEphemerisTime, finalEphemerisTime ) );

    // Create vehicle
    double vehicleMass = 2000.0;
    bodyMap[ "Vehicle" ] = std::make_shared< Body >( );
    bodyMap[ "Vehicle" ]->setConstantBodyMass( vehicleMass );

    std::vector< double > rightAscensionPole;
    rightAscensionPole.push_back( 0.0 ); rightAscensionPole.push_back( 0.2 );
    std::vector< double > declinationPole;
    declinationPole.push_back( mathematical_constants::PI / 2.0 ); declinationPole.push_back( 0.4 );
    std::vector< double > primeMeridianLongitude;
    primeMeridianLongitude.push_back( - mathematical_constants::PI / 2.0 ); primeMeridianLongitude.push_back( - 0.2 );
    std::vector< double > rotationalRate;
    rotationalRate.push_back( 1.0E-5 ); rotationalRate.push_back( 1.0E-5 );
    std::vector< double > numberSecondsSinceEpoch;
    numberSecondsSinceEpoch.push_back( 0.0 ); numberSecondsSinceEpoch.push_back( 0.0 );



    for ( int testCase = 0 ; testCase < 2 ; testCase++){

        // Compute spacecraft orbital period, and compute test times
        double orbitalPeriod = 2.0 * mathematical_constants::PI * std::sqrt( std::pow( physical_constants::ASTRONOMICAL_UNIT, 3.0 ) /
                                                                             spice_interface::getBodyGravitationalParameter( "Sun" ) );
        std::vector< double > testTimes = { 0.0, orbitalPeriod / 4.0, orbitalPeriod / 2.0, 3.0 * orbitalPeriod / 4.0 };

        // Put vehicle on circular orbit around Sun
        Eigen::Vector6d initialStateInKeplerianElements = Eigen::Vector6d::Zero( );
        initialStateInKeplerianElements[ 0 ] = physical_constants::ASTRONOMICAL_UNIT;
        bodyMap[ "Vehicle" ]->setEphemeris( std::make_shared< KeplerEphemeris >( initialStateInKeplerianElements, 0.0,
                                                     spice_interface::getBodyGravitationalParameter( "Sun" ), "Sun", "ECLIPJ2000", 1 ) );


        /// First calculation with simple rotational ephemeris and constant panel orientation

        // Define simple rotational ephemeris.
        bodyMap[ "Vehicle" ]->setRotationalEphemeris( std::make_shared< tudat::ephemerides::SimpleRotationalEphemeris >(
                        rightAscensionPole[ testCase ], declinationPole[ testCase ],  primeMeridianLongitude[ testCase ],
                        rotationalRate[ testCase ], numberSecondsSinceEpoch[ testCase ], "ECLIPJ2000", "VehicleFixed" ) );

        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );


        // Create radiation pressure properties.
        std::vector< double > areas;
        areas.push_back( 2.0 );
        areas.push_back( 4.0 );

        std::vector< double > emissivities;
        emissivities.push_back( 0.0 );
        emissivities.push_back( 0.1 );

        std::vector< double > diffuseReflectionCoefficients;
        diffuseReflectionCoefficients.push_back( 0.06 );
        diffuseReflectionCoefficients.push_back( 0.46 );

        std::vector< Eigen::Vector3d > panelSurfaceNormals;
        panelSurfaceNormals.push_back( Eigen::Vector3d::UnitX( ) );
        panelSurfaceNormals.push_back( - Eigen::Vector3d::UnitX( ) );


        // Create panelled radiation pressure interface.
        std::shared_ptr< PanelledRadiationPressureInterfaceSettings > radiationPressureInterfaceSettings =
                std::make_shared< PanelledRadiationPressureInterfaceSettings >(
                    "Sun", emissivities, areas, diffuseReflectionCoefficients, panelSurfaceNormals );

        std::shared_ptr< PanelledRadiationPressureInterface > radiationPressureInterface =
                std::dynamic_pointer_cast< PanelledRadiationPressureInterface >(
                    createRadiationPressureInterface( radiationPressureInterfaceSettings, "Vehicle", bodyMap ) );

        bodyMap[ "Vehicle" ]->setRadiationPressureInterface( "Sun", radiationPressureInterface );


        // Define accelerations.
        SelectedAccelerationMap accelerationMap;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOnVehicle;
        accelerationsOnVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                       panelled_radiation_pressure_acceleration ) );

        accelerationMap[ "Vehicle" ] = accelerationsOnVehicle;

        std::map< std::string, std::string > centralBodyMap; centralBodyMap[ "Vehicle" ] = "Sun";
        AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodyMap, accelerationMap, centralBodyMap );
        std::shared_ptr< AccelerationModel< Eigen::Vector3d > > accelerationModel =
                accelerationModelMap.at( "Vehicle" ).at( "Sun" ).at( 0 );


        // Compute radiation pressure acceleration for different Sun positions.
        std::vector< Eigen::Vector3d > calculatedAcceleration;

        Eigen::Vector3d sunCenteredVehiclePosition;
        std::shared_ptr< Ephemeris > vehicleEphemeris = bodyMap[ "Vehicle" ]->getEphemeris( );


        for( unsigned int i = 0; i < testTimes.size( ) ; i++ )
        {
            // Update environment and acceleration
            bodyMap[ "Sun" ]->setStateFromEphemeris( testTimes[ i ] );
            bodyMap[ "Vehicle" ]->setStateFromEphemeris( testTimes[ i ] );
            bodyMap[ "Vehicle" ]->setCurrentRotationToLocalFrameFromEphemeris( testTimes[ i ] );
            radiationPressureInterface->updateInterface( testTimes[ i ] );
            accelerationModel->updateMembers( testTimes[ i ] );

            // Retrieve acceleration.
            calculatedAcceleration.push_back( accelerationModel->getAcceleration( ) );

        }



        /// Second calculation with constant rotational ephemeris and time-varying panel orientation

        // Define constant rotational ephemeris
        Eigen::Vector7d rotationalStateVehicle;
        rotationalStateVehicle.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond( Eigen::Matrix3d::Identity() ));
        rotationalStateVehicle.segment( 4, 3 ) = Eigen::Vector3d::Zero();
        bodyMap[ "Vehicle" ]->setRotationalEphemeris( std::make_shared< ConstantRotationalEphemeris >( rotationalStateVehicle, "ECLIPJ2000",
                                                                                                       "VehicleFixed" ) );

        // Define time-varying panel orientation.
        std::vector< std::function< Eigen::Vector3d ( const double ) > > timeVaryingPanelSurfaceNormals;
        for ( int j = 0 ; j < panelSurfaceNormals.size() ; j++ ){

            timeVaryingPanelSurfaceNormals.push_back( [ = ]( const double currentTime ){

                Eigen::Vector3d currentPanelSurfaceOrientation = (reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion( basic_mathematics::computeModulo(
                          ( currentTime - numberSecondsSinceEpoch[ testCase ] ) * rotationalRate[ testCase ], 2.0 * mathematical_constants::PI ) )
                          * reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion(
                          declinationPole[ testCase ], rightAscensionPole[ testCase ], primeMeridianLongitude[ testCase ] ) )
                          .toRotationMatrix().inverse() * panelSurfaceNormals[ j ];

                return currentPanelSurfaceOrientation; } );

        }

        // Create panelled radiation pressure interface.
        std::shared_ptr< PanelledRadiationPressureInterfaceSettings > radiationPressureInterfaceSettingsWithTimeVaryingPanelSurfaceNormal =
                std::make_shared< PanelledRadiationPressureInterfaceSettings >(
                    "Sun", emissivities, areas, diffuseReflectionCoefficients, timeVaryingPanelSurfaceNormals );

        std::shared_ptr< PanelledRadiationPressureInterface > radiationPressureInterfaceTimeVaryingSurfaceNormal =
                std::dynamic_pointer_cast< PanelledRadiationPressureInterface >(
                    createRadiationPressureInterface( radiationPressureInterfaceSettingsWithTimeVaryingPanelSurfaceNormal, "Vehicle", bodyMap ) );

        bodyMap[ "Vehicle" ]->setRadiationPressureInterface( "Sun", radiationPressureInterfaceTimeVaryingSurfaceNormal );

        accelerationModelMap = createAccelerationModelsMap( bodyMap, accelerationMap, centralBodyMap );

        std::shared_ptr< AccelerationModel< Eigen::Vector3d > > accelerationModelTimeVaryingPanelSurfaceNormal =
                accelerationModelMap.at( "Vehicle" ).at( "Sun" ).at( 0 );


        // Compute radiation pressure acceleration for different Sun positions.
        std::vector< Eigen::Vector3d > calculatedAccelerationTimeVaryingPanelOrientation;

        for( unsigned int i = 0; i < testTimes.size( ) ; i++ )
        {
            // Update environment and acceleration
            bodyMap[ "Sun" ]->setStateFromEphemeris( testTimes[ i ] );
            bodyMap[ "Vehicle" ]->setStateFromEphemeris( testTimes[ i ] );
            bodyMap[ "Vehicle" ]->setCurrentRotationToLocalFrameFromEphemeris( testTimes[ i ] );
            radiationPressureInterfaceTimeVaryingSurfaceNormal->updateInterface( testTimes[ i ] );
            accelerationModelTimeVaryingPanelSurfaceNormal->updateMembers( testTimes[ i ] );

            // Retrieve acceleration.
            calculatedAccelerationTimeVaryingPanelOrientation.push_back( accelerationModelTimeVaryingPanelSurfaceNormal->getAcceleration( ) );

        }


        for( unsigned int j = 0; j < testTimes.size() ; j++ )
        {
            for ( unsigned int i = 0 ; i < 3 ; i++ ){
                BOOST_CHECK_SMALL( std::fabs(
                                       calculatedAcceleration[ j ][ i ] - calculatedAccelerationTimeVaryingPanelOrientation[ j ][ i ] ), 3.0E-23 );
            }
        }
    }

}




BOOST_AUTO_TEST_SUITE_END( )


}

}
