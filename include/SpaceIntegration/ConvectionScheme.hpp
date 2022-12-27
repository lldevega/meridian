/*
* Convection schemes implementation
*/

#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cassert>
#include <memory>

/// Base class for convection schemes
/*
 *
 */
class ConvectionSchemeBase
{
public:
	/// Default constructor.
	ConvectionSchemeBase() {};

	/// Numerical flux computation from left and right states.
	/*
	 * Requires the left and right states.
	 */
	virtual FieldVector ComputeNumericalFluxVectorDotNormal(const FieldVector &leftState,
	    const FieldVector &rightState, const double nx, const double ny) const = 0 ;

	/// Physical flux computation from left and right states.
	/*
	 * Requires the left and right states.
	 */
	FluxVector ComputePhysicalFluxVector(const FieldVector &state) const
	{
		const FluxVector physicalFluxVector = FluxVector(state);
		return physicalFluxVector;
	}

protected:
	/// Compute the delta of the state between left and right element.
	/*
	 * Requires the left and right states.
	 */
	FieldVector ComputeDeltaState(const FieldVector &leftState,  const FieldVector &rightState) const
	{
		FieldVector deltaState = rightState;
		deltaState -= leftState;
		return deltaState;
	}

	/// Compute the flux vector dot product with the normal.
	/*
	 * Requires the flux vector and the normal components
	 */
	FieldVector ComputeFluxVectorDotNormal(const FluxVector &fluxVector, double nx, double ny) const
	{
		// initialize the return
		FieldVector fluxVectorDotNormal;

		// compute the flux vector dot product with the normal components
		fluxVectorDotNormal.SetDensity(fluxVector.GetDensityFlux()[0] * nx +
			fluxVector.GetDensityFlux()[1] * ny);

		fluxVectorDotNormal.SetMomentumX(fluxVector.GetMomentumXFlux()[0] * nx +
			fluxVector.GetMomentumXFlux()[1] * ny);

		fluxVectorDotNormal.SetMomentumY(fluxVector.GetMomentumYFlux()[0] * nx +
			fluxVector.GetMomentumYFlux()[1] * ny);

		fluxVectorDotNormal.SetStagnationEnergy(fluxVector.GetStagnationEnergyFlux()[0] * nx +
			fluxVector.GetStagnationEnergyFlux()[1] * ny);

		return fluxVectorDotNormal;
	}
};

/// Implementation of a basic central scheme.
/*
 * The computation of the numerical flux is the average of the left and right ones.
 */
class CentralScheme: public ConvectionSchemeBase
{
public:
	using BaseClass = ConvectionSchemeBase;
	using BaseClass::BaseClass;

	/// Numerical flux computation.
	/*
	 * Requires the left and right states.
	 */
	FieldVector ComputeNumericalFluxVectorDotNormal(const FieldVector &leftState,
		const FieldVector &rightState, const double nx, const double ny) const override
	{
		// averaged flux vector
		FluxVector numericalFluxVector = ComputePhysicalFluxVector(leftState);
		numericalFluxVector += ComputePhysicalFluxVector(rightState);
		numericalFluxVector *= 0.5;

		// compute the dot product
		FieldVector numericalFluxVectorDotNormal = this->ComputeFluxVectorDotNormal(numericalFluxVector,
			nx, ny);

		return numericalFluxVectorDotNormal;
	}
};

/// Implementation of the HLL scheme.
/*
 * Bibliography: TODO
 */
class HLLScheme: public ConvectionSchemeBase
{
public:
	using BaseClass = ConvectionSchemeBase;
	using BaseClass::BaseClass;

	/// Numerical flux computation
	/*
	 * Requires the left and right states.
	 */
	FieldVector ComputeNumericalFluxVectorDotNormal(const FieldVector &leftState, const FieldVector &rightState,
		const double nx, const double ny) const override
	{
		/// left and right fluxes
		FluxVector numericalFluxVectorLeft = ComputePhysicalFluxVector(leftState);
	    FluxVector numericalFluxVectorRight = ComputePhysicalFluxVector(rightState);

	    // dot product with normal
	    FieldVector numericalFluxVectorDotNormalLeft = this->ComputeFluxVectorDotNormal(numericalFluxVectorLeft,
			nx, ny);
	    FieldVector numericalFluxVectorDotNormalRight = this->ComputeFluxVectorDotNormal(numericalFluxVectorRight,
			nx, ny);

	    // wave speeds
	    std::vector<double> waveSpeeds = EstimateWaveSpeeds(leftState, rightState, nx, ny);
	    double sL = waveSpeeds[0];
	    double sR = waveSpeeds[1];

	    // case 1
	    if (0 <= sL)
	    	return numericalFluxVectorDotNormalLeft;
	    // case 2
	    if (0 >= sR)
	    	return numericalFluxVectorDotNormalRight;
	    // case 3
	    else
	    {
	    	// sanity check
	    	assert(sL <= 0 && 0 <= sR &&
	    		"HLLScheme::ComputeNumericalFluxVector Something went wrong with the wave speeds...");

	    	// compute the inverse of the denominator
	    	double invDenominator = 1. / (sR - sL);

	    	// compute delta state contribution
	    	FieldVector deltaState = this->ComputeDeltaState(leftState, rightState);
            deltaState *= sL * sR;

	    	// compute numerator
	    	FieldVector numerator = numericalFluxVectorDotNormalLeft * sR -
	    	    numericalFluxVectorDotNormalRight * sL;
	    	numerator += deltaState;

	    	FieldVector numericalFluxVectorDotNormal = numerator * invDenominator;

	    	return numericalFluxVectorDotNormal;
	    }
	}

protected:
	std::vector<double> EstimateWaveSpeeds(const FieldVector &leftState, const FieldVector &rightState,
		const double nx, const double ny) const
	{
		std::vector<double> averagedState = ComputeRoeAveragedFaceState(leftState, rightState, nx, ny);
		double sL = std::abs(averagedState[4]) - averagedState[3]; // v~ - c~
		double sR = std::abs(averagedState[4]) + averagedState[3]; // v~ + c~

		std::vector<double> waveSpeeds{sL, sR};

		return waveSpeeds;
	}

	std::vector<double> ComputeRoeAveragedFaceState(const FieldVector &leftState, const FieldVector &rightState,
		const double nx, const double ny) const
	{
		// averaged density: ro^~ = sqrt(roL * roR)
		double averagedDensity = std::pow(leftState.GetDensity() * rightState.GetDensity(), 0.5);

		// left and right density square root
		double sqrtRootLeftDensity = std::pow(leftState.GetDensity(), 0.5);
		double sqrtRootRightDensity = std::pow(rightState.GetDensity(), 0.5);

		// left and right velocity x and y
		double leftVelocityX = leftState.ConservativeToPrimitive()[1];
		double leftVelocityY = leftState.ConservativeToPrimitive()[2];
		double rightVelocityX = rightState.ConservativeToPrimitive()[1];
		double rightVelocityY = rightState.ConservativeToPrimitive()[2];

		// left and right entalpy
		double leftEntalpy = leftState.ComputeStagnationEntalpy();
		double rightEntalpy = rightState.ComputeStagnationEntalpy();

		// common averaging denominator
		double denominator = sqrtRootLeftDensity + sqrtRootRightDensity;

		// averaged velocity x
		double averagedVelocityX = (sqrtRootLeftDensity * leftVelocityX +
			sqrtRootRightDensity * rightVelocityX) / denominator;

		// averaged velocity y
		double averagedVelocityY = (sqrtRootLeftDensity * leftVelocityY +
			sqrtRootRightDensity * rightVelocityY) / denominator;

		// averaged kinetic energy
		double averagedK = 0.5 * (averagedVelocityX * averagedVelocityX +
			averagedVelocityY * averagedVelocityY);

		// averaged entalpy
		double averagedEntalpy = (sqrtRootLeftDensity *  leftEntalpy +
			sqrtRootRightDensity * rightEntalpy) / denominator;

		// averaged speed of sound
		double averagedSpeedOfSound = std::pow(0.4 * (averagedEntalpy - averagedK), 0.5);

		// averaged velocity as seen by the left cell
		double averagedVelocity = nx * averagedVelocityX + ny * averagedVelocityY;

		return std::vector<double>{averagedDensity, averagedVelocityX, averagedVelocityY,
			averagedSpeedOfSound, averagedVelocity};
	}
};
