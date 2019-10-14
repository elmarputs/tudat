/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef HYBRIDMETHOD_H
#define HYBRIDMETHOD_H

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include <map>
#include "Tudat/Astrodynamics/LowThrustDirectMethods/hybridMethodLeg.h"
#include "Tudat/Astrodynamics/LowThrustDirectMethods/optimisationSettings.h"

namespace tudat
{
namespace low_thrust_direct_methods
{

//! Transform thrust model as a function of time into hybrid method thrust model.
Eigen::Matrix< double, 10 , 1 > convertToHybridMethodThrustModel( std::function< Eigen::Vector3d( const double ) > thrustModelWrtTime );

class HybridMethod : public transfer_trajectories::LowThrustLeg
{
public:

    //! Constructor.
    HybridMethod(
            const Eigen::Vector6d& stateAtDeparture,
            const Eigen::Vector6d& stateAtArrival,
            const double maximumThrust,
            const double specificImpulse,
            const double timeOfFlight,
            simulation_setup::NamedBodyMap& bodyMap,
            const std::string bodyToPropagate,
            const std::string centralBody,
            std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
            std::shared_ptr< transfer_trajectories::OptimisationSettings > optimisationSettings ) :
        LowThrustLeg( stateAtDeparture, stateAtArrival, timeOfFlight, bodyMap, bodyToPropagate, centralBody ),
        maximumThrust_( maximumThrust ),
        specificImpulse_( specificImpulse ),
        integratorSettings_( integratorSettings ),
        optimisationSettings_( optimisationSettings )
    {

        // Store initial spacecraft mass.
        initialSpacecraftMass_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass();

        // Convert the thrust model proposed as initial guess into simplified thrust model adapted to the hybrid method.
        if ( optimisationSettings_->initialGuessThrustModel_.first != nullptr )
        {
            initialGuessThrustModel_.first = convertToHybridMethodThrustModel( optimisationSettings_->initialGuessThrustModel_.first );
        }
        else
        {
            Eigen::VectorXd emptyVector;
            initialGuessThrustModel_.first = emptyVector; // emptyVector;
        }
        initialGuessThrustModel_.second = optimisationSettings_->initialGuessThrustModel_.second;
        // Perform optimisation
        std::pair< std::vector< double >, std::vector< double > > bestIndividual = performOptimisation( );
        championFitness_ = bestIndividual.first;
        championDesignVariables_ = bestIndividual.second;


        // Transform vector of design variables into 3D vector of throttles.
        Eigen::VectorXd initialCostates; initialCostates.resize( 5 );
        Eigen::VectorXd finalCostates; finalCostates.resize( 5 );
        for ( unsigned int i = 0 ; i < 5 ; i++ )
        {
            initialCostates[ i ] = championDesignVariables_[ i ];
            finalCostates[ i ] = championDesignVariables_[ i + 5 ];
        }

        bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );

        // Create Sims-Flanagan leg from the best optimisation individual.
        hybridMethodLeg_ = std::make_shared< HybridMethodLeg >( stateAtDeparture_, stateAtArrival_, initialCostates, finalCostates, maximumThrust_,
                                                                specificImpulse_, timeOfFlight_, bodyMap_, bodyToPropagate_, centralBody_, integratorSettings );

    }

    //! Default destructor.
    ~HybridMethod( ) { }

    //! Convert time to independent variable.
    double convertTimeToIndependentVariable( const double time )
    {
        return time;
    }

    //! Convert independent variable to time.
    double convertIndependentVariableToTime( const double independentVariable )
    {
        return independentVariable;
    }

    //! Perform optimisation.
    std::pair< std::vector< double >, std::vector< double > > performOptimisation( );

    //! Compute DeltaV.
    double computeDeltaV( )
    {
        return hybridMethodLeg_->getTotalDeltaV( );
    }

    //! Compute current cartesian state.
    Eigen::Vector6d computeCurrentStateVector( const double currentTime );

    //! Compute state history.
    void getTrajectory(
            std::vector< double >& epochsVector,
            std::map< double, Eigen::Vector6d >& propagatedTrajectory )
    {
        propagatedTrajectory = hybridMethodLeg_->propagateTrajectory( epochsVector, propagatedTrajectory );
    }


    Eigen::Vector3d computeCurrentThrust( double time,
                                          std::function< double ( const double ) > specificImpulseFunction,
                                          std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Return thrust profile.
    void getThrustProfile( std::vector< double >& epochsVector,
                           std::map< double, Eigen::VectorXd >& thrustProfile,
                           std::function< double ( const double ) > specificImpulseFunction,
                           std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Compute direction thrust acceleration in cartesian coordinates.
    Eigen::Vector3d computeCurrentThrustAccelerationDirection(
            double currentTime, std::function< double ( const double ) > specificImpulseFunction,
            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Compute magnitude thrust acceleration.
    double computeCurrentThrustAccelerationMagnitude(
            double currentTime, std::function< double ( const double ) > specificImpulseFunction,
            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Return thrust acceleration profile.
    void getThrustAccelerationProfile(
            std::vector< double >& epochsVector,
            std::map< double, Eigen::VectorXd >& thrustAccelerationProfile,
            std::function< double ( const double ) > specificImpulseFunction,
            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );


    //! Return best individual.
    std::vector< double > getBestIndividual( )
    {
        return championDesignVariables_;
    }

    //! Return fitness of best individual.
    std::vector< double > getBestIndividualFitness( )
    {
        return championFitness_;
    }

    //! Return best hybrid method leg after optimisation.
    std::shared_ptr< HybridMethodLeg > getOptimalHybridMethodLeg( )
    {
        return hybridMethodLeg_;
    }


    //! Retrieve acceleration map (thrust and central gravity accelerations).
    basic_astrodynamics::AccelerationMap retrieveLowThrustAccelerationMap(
            std::function< double ( const double ) > specificImpulseFunction,
            std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Define appropriate translational state propagator settings for the full propagation.
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > createLowThrustTranslationalStatePropagatorSettings(
            basic_astrodynamics::AccelerationMap accelerationModelMap,
            std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave );


private:

    //! Maximum allowed thrust.
    double maximumThrust_;

    //! Specific impulse.
    double specificImpulse_;

    //! Integrator settings.
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings_;

    //! Optimisation settings.
    std::shared_ptr< transfer_trajectories::OptimisationSettings > optimisationSettings_;

    //! Initial guess for the optimisation.
    //! The first element contains the thrust throttles corresponding to the initial guess for the thrust model.
    //! The second element defines the bounds around the initial time (in percentage).
    std::pair< Eigen::VectorXd, double > initialGuessThrustModel_;

    //! Fitness vector of the optimisation best individual.
    std::vector< double > championFitness_;

    //! Design variables vector corresponding to the optimisation best individual.
    std::vector< double > championDesignVariables_;

    //! Initial mass of the spacecraft.
    double initialSpacecraftMass_;

    //! Hybrid method leg corresponding to the best optimisation output.
    std::shared_ptr< HybridMethodLeg > hybridMethodLeg_;

};


} // namespace low_thrust_direct_methods
} // namespace tudat

#endif // HYBRIDMETHOD_H