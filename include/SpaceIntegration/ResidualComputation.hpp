/*
* Residual computation implementation.
*/

#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cassert>
#include <memory>

template <class ConvectionScheme>
class ResidualComputation
{
public:
	/// Forward types from mesh class.
	using CellContainer = typename Mesh::CellContainer;
	using FaceContainer = typename Mesh::FaceContainer;

	/// Constructor
	/*
	 * Requires the inner and boundary face containers and their respective stencils.
	 */
	ResidualComputation(const CellContainer &cellContainer, const FaceContainer &innerFaceContainer,
		Stencil &innerFaceStencil, const BoundaryConditionContainer &bcContainer,
		const ConvectionScheme &convectionScheme):
			_cellContainer(cellContainer), _innerFaceContainer(innerFaceContainer),
	        _innerFaceStencil(innerFaceStencil), _bcContainer(bcContainer), _convectionScheme(convectionScheme)
	{ }

	/// Compute residual method.
	void ComputeResidual(const State &solution, State &residual)
	{
		assert(solution.size() == residual.size() && solution.size() > 0 &&
			"Solution and residual sizes do not match or are not initialized");

		// inner face loop
		for (int i = 0; i < _innerFaceContainer.size(); ++i)
		{
			// get the face and its normal
			const Face &innerFace = _innerFaceContainer[i];
			const auto &normal = innerFace.GetNormal();
			double length = innerFace.GetLength();
			const auto faceIdx = innerFace.GetCenter().GetNodeID();

			// get the left and right element indices
			auto leftElemIdx = _innerFaceStencil.GetStencil()[faceIdx][0];
			auto rightElemIdx = _innerFaceStencil.GetStencil()[faceIdx][1];

			// compute the center-to-face vector for the left element. TODO: this must be precomputed data
			double leftRadiusXComp = innerFace.GetCenter().GetXCoord()
				- this->_cellContainer[leftElemIdx].GetCenter().GetXCoord();
			double leftRadiusYComp = innerFace.GetCenter().GetYCoord()
				- this->_cellContainer[leftElemIdx].GetCenter().GetYCoord();

			// get the sign of the face normal as seen from the left element (convention: pointing outwards)
			double leftRadiusXCompSign = 1.0 - 2.0 * static_cast<double> (std::signbit(leftRadiusXComp));
			double leftRadiusYCompSign = 1.0 - 2.0 * static_cast<double> (std::signbit(leftRadiusYComp));

			// scale the normal components
			double nx = std::abs(normal[0]) * leftRadiusXCompSign;
			double ny = std::abs(normal[1]) * leftRadiusYCompSign;

			// get the left and right states
			const FieldVector &leftState = solution.GetFieldVectorAt(leftElemIdx);
			const FieldVector &rightState = solution.GetFieldVectorAt(rightElemIdx);

            // compute the numerical flux vector
			FieldVector numericalFluxVectorDotNormal = this->_convectionScheme.
				ComputeNumericalFluxVectorDotNormal(leftState, rightState, nx, ny);

			// the residual contribution of the current face for the cells in its stencil
			// is the dot product of the numerical flux vector and the face normal
			double res_ro = numericalFluxVectorDotNormal.GetDensity() * length;
			double res_rou = numericalFluxVectorDotNormal.GetMomentumX()* length;
			double res_rov = numericalFluxVectorDotNormal.GetMomentumY() * length;
			double res_roE = numericalFluxVectorDotNormal.GetStagnationEnergy() * length;

			// write the face contribution to the residual using the field vector data structure
			// for the left element, the negative sign is used as the face normal points outward wrt this element
			const FieldVector faceResidualContributionToLeft(-res_ro, -res_rou, -res_rov, -res_roE);
			// for the right element, the positive sign is used as the face normal points inwards wrt this element
			const FieldVector faceResidualContributionToRight(res_ro, res_rou, res_rov, res_roE);

			// write the contribution to the element residual
			residual.AddFieldVectorAt(leftElemIdx, faceResidualContributionToLeft);
			residual.AddFieldVectorAt(rightElemIdx, faceResidualContributionToRight);
		}

		// boundary condition loop
		for (auto bc : this->_bcContainer)
		{
			//bc->Print();
			BoundaryConditionLoop(bc, solution, residual);
		}
	}

protected:
	/// Members.
	const CellContainer &_cellContainer;
	const FaceContainer &_innerFaceContainer;
	Stencil &_innerFaceStencil;
	const BoundaryConditionContainer &_bcContainer;
	const ConvectionScheme &_convectionScheme;

	template <class BoundaryConditionType>
	void BoundaryConditionLoop(BoundaryConditionType bc_, const State &solution, State &residual)
	{
		// dereference bc
		auto &bc = *bc_;

		// retrieve the face container for this boundary codntion
		const FaceContainer &bcFaceContainer = bc.GetBoundaryFaceContainer();
		Stencil bcFaceStencil = bc.GetBoundaryFaceStencil();

		// loop over boundary faces of the current boundary condition
		for (int i = 0; i < bcFaceContainer.size(); ++i)
		{
			const Face &bcFace = bcFaceContainer[i];
			const auto &normal = bcFace.GetNormal();
			double length = bcFace.GetLength();

			// get the element index
			auto elemIdx = bcFaceStencil.GetStencil()[bcFace.GetCenter().GetNodeID()][0];

			// compute the center-to-face vector for the left element. TODO: this must be precomputed data
			double radiusXComp = bcFace.GetCenter().GetXCoord()
				- this->_cellContainer[elemIdx].GetCenter().GetXCoord();
			double radiusYComp = bcFace.GetCenter().GetYCoord()
				- this->_cellContainer[elemIdx].GetCenter().GetYCoord();

			// get the sign of the face normal as seen from the left element (convention: pointing outwards)
			double radiusXCompSign = 1.0 - 2.0 * static_cast<double> (std::signbit(radiusXComp));
			double radiusYCompSign = 1.0 - 2.0 * static_cast<double> (std::signbit(radiusYComp));

			// scale the normal components
			double nx = std::abs(normal[0]) * radiusXCompSign;
			double ny = std::abs(normal[1]) * radiusYCompSign;

			// get the face state
			const FieldVector &faceState = bc.ApplyTreatmentAt(i, solution);

			// get the attached element state
			const FieldVector &cellCenteredState = solution.GetFieldVectorAt(elemIdx);

			// compute the numerical flux vector
			FieldVector numericalFluxVectorDotNormal = this->_convectionScheme.
				ComputeNumericalFluxVectorDotNormal(cellCenteredState, faceState, nx, ny);

			// the residual contribution of the current face for the cells in its stencil
			// is the dot product of the numerical flux vector and the face normal
			double res_ro = numericalFluxVectorDotNormal.GetDensity() * length;
			double res_rou = numericalFluxVectorDotNormal.GetMomentumX()* length;
			double res_rov = numericalFluxVectorDotNormal.GetMomentumY() * length;
			double res_roE = numericalFluxVectorDotNormal.GetStagnationEnergy() * length;

			// write the face contribution to the residual using the field vector data structure
			// for the element, the negative sign is used as the face normal points outward wrt this element
			const FieldVector faceResidualContribution(-res_ro, -res_rou, -res_rov, -res_roE);

			// write the contribution to the element residual
			residual.AddFieldVectorAt(elemIdx, faceResidualContribution);
		}
	}
};
