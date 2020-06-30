/*    Copyright (c) 2010-2020, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATBUNDLE_YARKOVSKYACCELERATION_H
#define TUDATBUNDLE_YARKOVSKYACCELERATION_H

#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

namespace tudat
{
namespace electro_magnetism
{

	class YarkovskyAcceleration: public basic_astrodynamics::AccelerationModel3d
	{
	private:
		//! Typedef for double-returning function.
		typedef std::function< double( ) > DoubleReturningFunction;

		//! Typedef for Eigen::Vector3d-returning function.
		typedef std::function< Eigen::Vector3d( ) > Vector3dReturningFunction;

	public:

		YarkovskyAcceleration( Vector3dReturningFunction sourcePositionFunction,
				Vector3dReturningFunction acceleratedBodyPositionFunction,
				Eigen::Vector3d spinAxis, const double stellarLuminosity,
				const double efficiencyFactor, const double diameter,
				const double density, const double phaseLag ):
			sourcePositionFunction_( sourcePositionFunction ),
			acceleratedBodyPositionFunction_( acceleratedBodyPositionFunction ),
			stellarLuminosity_( stellarLuminosity ), efficiencyFactor_(efficiencyFactor), phaseLag_( phaseLag ),
			diameter_( diameter ), density_( density ), spinAxis_( spinAxis )
		{
			this->updateMembers( );
		}

		Eigen::Vector3d getAcceleration( )
		{
			double positionNorm = currentVectorToSource_.norm( );

			return ( efficiencyFactor_ * 3.0 * tudat::mathematical_constants::PI * stellarLuminosity_ /
					 ( 8.0 * diameter_ * density_ * tudat::physical_constants::SPEED_OF_LIGHT *
					   std::pow( positionNorm, 3 ))) * ( rotationMatrix_ * currentVectorToSource_ );
		}

		void updateMembers( const double currentTime = TUDAT_NAN )
		{
			// Calculate the current position of the accelerated body w.r.t. the source body
			currentVectorToSource_ = sourcePositionFunction_( ) - acceleratedBodyPositionFunction_( );

			rotationMatrix_ = Eigen::AngleAxisd( phaseLag_, spinAxis_ );
		}

	private:

		//! Function pointer returning position of source.
		/*!
		 * Function pointer returning position of source (3D vector).
		 */
		const Vector3dReturningFunction sourcePositionFunction_;

		//! Function pointer returning position of accelerated body.
		/*!
		 * Function pointer returning position of accelerated body (3D vector).
		 */
		const Vector3dReturningFunction acceleratedBodyPositionFunction_;

		// Initialized by updateMembers()
		Eigen::Vector3d currentVectorToSource_;

		double stellarLuminosity_;

		double efficiencyFactor_;

		double phaseLag_;

		// Initialized by updateMembers()
		Eigen::AngleAxisd rotationMatrix_;

		double diameter_;

		double density_;

		Eigen::Vector3d spinAxis_;
	};
}
}

#endif //TUDATBUNDLE_YARKOVSKYACCELERATION_H
