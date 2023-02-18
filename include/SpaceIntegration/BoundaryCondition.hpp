/*
* Boundary condition treatment implementation
*/

#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cassert>
#include <memory>
#include <map>

/// Boundary condition class
/*
 * This class provides the data structure hosting a list of boundary faces that share a treatment.
 * It also provides the functionality to retrieve the face state after the BC treatment is applied.
 */
class BoundaryConditionBase
{
public:
	/// Forward types from mesh class.
	using CellContainer = typename Mesh::CellContainer;
	using FaceContainer = typename Mesh::FaceContainer;

	/// Constructor.
	/*
	 *
	 */
	BoundaryConditionBase(const CellContainer &cellContainer, const FaceContainer &bdryFaceContainer,
	    Stencil &bdryFaceStencil): _cellContainer(cellContainer), _bdryFaceContainer(bdryFaceContainer),
	    _bdryFaceStencil(bdryFaceStencil) { }

	/// Default constructor.
	BoundaryConditionBase(): _cellContainer(), _bdryFaceContainer(), _bdryFaceStencil(){ }

//	/// Asignment operator.
//	void operator=(const BoundaryConditionBase &other)
//	{
//        this->_cellContainer = other.GetCellContainer();
//        this->_bdryFaceContainer = other.GetBoundaryFaceContainer();
//        this->_bdryFaceStencil = other.GetBoundaryFaceStencil();
//	}

	const CellContainer &GetCellContainer() const
	{
		return this->_cellContainer;
	}

	const FaceContainer &GetBoundaryFaceContainer() const
	{
		return this->_bdryFaceContainer;
	}

	Stencil GetBoundaryFaceStencil() const
	{
		return this->_bdryFaceStencil;
	}

	/// Applies the treatment at a face idx
	/*
	 *
	 */
	FieldVector ApplyTreatmentAt(int idx, const State &solution)
	{
		// get the associated face object
		const Face &bdryFace = _bdryFaceContainer[idx];

		// get the element the BC face is attached to and the state
		int elemIdx = _bdryFaceStencil.GetStencil()[idx][0];
        const Cell &cell = _cellContainer[elemIdx];
        const FieldVector &cellCenteredState = solution.GetFieldVectorAt(elemIdx);

	    const FieldVector faceState = ComputeFaceState(bdryFace, cell, cellCenteredState);

		return faceState;
	}

	/// Purely virtual method to compute the face state.
	/*
	 * Requires the boundary face object, the cell object to which the face is attached and its cell-centered state
	 * TODO: make it a void and throw expcetion if called.
	 */
	virtual FieldVector ComputeFaceState (const Face &bdryFace, const Cell &cell, const FieldVector &cellCenterdState)
	    const
	{
		return FieldVector();
	}

protected:
	/// Members.
	const CellContainer _cellContainer; // gets copied, but otherwise default ctor not possible..
	const FaceContainer _bdryFaceContainer; // gets copied, but otherwise default ctor not possible..
	Stencil _bdryFaceStencil; // gets copied, but otherwise default ctor not possible..
};

/// TODO
class BCWallInviscid: public BoundaryConditionBase
{
public:
	/// Forward constructor.
	using BaseClass = BoundaryConditionBase;
	using BaseClass::BaseClass;

    /// Implementation of `BaseClass::ComputeFaceState`
	FieldVector ComputeFaceState(const Face &bdryFace, const Cell &cell, const FieldVector &cellCenteredState)
	    const override
	{
		// retrieve the normal
		const auto &normal = bdryFace.GetNormal();

		// compute the face to center vector for the attached element. TODO: this must be precomputed data
		double radiusXComp = cell.GetCenter().GetXCoord() - bdryFace.GetCenter().GetXCoord();
		double radiusYComp = cell.GetCenter().GetYCoord() - bdryFace.GetCenter().GetYCoord();

		// get the sign of the face normal as seen from the element (convention: pointing inwards)
		double radiusXCompSign = 1.0 - 2.0 * static_cast<double> (std::signbit(radiusXComp));
		double radiusYCompSign = 1.0 - 2.0 * static_cast<double> (std::signbit(radiusYComp));

		// scale the normal components
		double nx = std::abs(normal[0]) * radiusXCompSign;
		double ny = std::abs(normal[1]) * radiusYCompSign;

		// get the cell centered state extrapolated to face (constant reconstruction)
		double density =  cellCenteredState.GetDensity();
		double velocityX = cellCenteredState.GetMomentumX() / cellCenteredState.GetDensity();
		double velocityY = cellCenteredState.GetMomentumY() / cellCenteredState.GetDensity();
		double pressure = cellCenteredState.ComputePressure();

		// retrieve the momentum projection on the normal
		double vN = velocityX * nx + velocityY * ny;

		// remove the previous projection from the velocity, projected on the normal pointing inwards
		// this yield the velocity to be imposed at the face
		double faceVelocityX = velocityX - 2. * vN * nx;
		double faceVelocityY = velocityY - 2. * vN * ny;

		/// TODO
		const auto gas = GasThermodynamicProperties();

		// compute the face state: density and pressure are extrapolated
		double faceDensity = density;
		double facePressure = pressure;
		double faceTemperature = facePressure / (faceDensity * gas._R);
		double faceMomentumX = faceDensity * faceVelocityX;
		double faceMomentumY = faceDensity * faceVelocityY;
		double faceVelocitySquared = faceVelocityX * faceVelocityX + faceVelocityY * faceVelocityY;
		double faceStagnationEnergy = faceDensity * (gas._Cv * faceTemperature + 0.5 * faceVelocitySquared);

		// and set them in the face state
		FieldVector bdryFaceState(faceDensity, faceMomentumX, faceMomentumY, faceStagnationEnergy);

		// check whether the boundary is working as expected
		double computedMomentumDotNormal = bdryFaceState.GetMomentumX() * nx +
			bdryFaceState.GetMomentumY() * ny;
		double delta = std::abs(computedMomentumDotNormal);
		//assert(delta < 0.0001 && "BCWallInviscid::ComputeFaceState() failed to computed face state");
		double delta2 = std::abs(bdryFaceState.ComputePressure() - cellCenteredState.ComputePressure());
		//assert(delta2 < 0.0001 && "BCWallInviscid::ComputeFaceState() failed to computed face state");

		return bdryFaceState;
	}
};

/// Boundary condition for imposed exterior state
/*
 * This boundary condition imposes a reference exterior state that is enforced. This state is defined using
 * the temperature, pressure, mach and flow direction.
 */
class BCExteriorState: public BoundaryConditionBase
{
public:
	/// Forward constructor.
	using BaseClass = BoundaryConditionBase;

	/// Constructor allowing to specify a reference state to be imposed (overrides BaseClass' default one).
	BCExteriorState(const CellContainer &cellContainer, const FaceContainer &bdryFaceContainer,
        Stencil &bdryFaceStencil, double Mach, double pressure, double temperature, double AoA /*in degrees*/):
        BaseClass(cellContainer, bdryFaceContainer, bdryFaceStencil),
		_refState(ComputeReferenceState(Mach, pressure, temperature, AoA)) { }

    /// Implementation of `BaseClass::ComputeFaceState`
	FieldVector ComputeFaceState(const Face &bdryFace, const Cell &/*cell*/, const FieldVector &/*cellCenteredState*/)
	    const override
	{
		return GetReferenceState();
	}

	/// Get private member
	const FieldVector &GetReferenceState() const
	{
		return this->_refState;
	}

protected:
	/// Members.
	const FieldVector _refState;

	/// Reference state computation
	/*
	 * Retrieves the conservative state from the Mach, pressure, temperature and angle of attack information.
	 */
	const FieldVector ComputeReferenceState(double Mach, double pressure, double temperature,
		double AoA /*in degrees*/)
	{
		// get the gas properties. TODO: fix dangling reference problem
		const auto gas = GasThermodynamicProperties(); //this->_refState.GetGas();

		// compute the associated state
		double velocity = Mach * std::pow(gas._gamma * gas._R * temperature, 0.5);
		double sinAoA = std::sin(AoA * M_PI / 180.0);
		double cosAoA = std::cos(AoA * M_PI / 180.0);
		double velocityX = velocity * cosAoA;
		double velocityY = velocity * sinAoA;
		double density = pressure / (gas._R * temperature);
		double stagnationEnergy = density * gas._Cv * temperature + 0.5 * density * velocity * velocity;

		FieldVector refState(density, density * velocityX, density * velocityY, stagnationEnergy);

		return refState;
	}
};

/// TODO
class BCStagnationInflow: public BoundaryConditionBase
{
public:
	/// Forward base class.
	using BaseClass = BoundaryConditionBase;
	using BaseClass::BaseClass;

	/// Constructor allowing to specify an outflow pressure (overrides BaseClass' default one).
	BCStagnationInflow(const CellContainer &cellContainer, const FaceContainer &bdryFaceContainer,
	    Stencil &bdryFaceStencil,
		double stagnationPressure, double stagnationTemperature, std::vector<double> direction):
	    BaseClass(cellContainer, bdryFaceContainer, bdryFaceStencil),
		_Pi(stagnationPressure), _Ti(stagnationTemperature), _dir(direction) { }

    /// Implementation of `BaseClass::ComputeFaceState`
	FieldVector ComputeFaceState(const Face &bdryFace, const Cell &cell, const FieldVector &cellCenteredState)
	    const override
	{
		/// TODO
		const auto gas = GasThermodynamicProperties(); //this->_refState.GetGas();

		// compute the normal velocity
		double velocityX = cellCenteredState.ConservativeToPrimitive()[1];
		double velocityY = cellCenteredState.ConservativeToPrimitive()[2];

	    double velocity = std::sqrt(velocityX * velocityX + velocityY * velocityY);
	    double speedOfSound = cellCenteredState.ComputeSpeedOfSound();

	    double faceHt = gas._Cp * _Ti;
	    double jMinus = velocity - 2 * speedOfSound / (gas._gamma - 1.);
	    double jPlus = (-(3. - gas._gamma) * jMinus + 4 * std::sqrt((gas._gamma + 1.) * faceHt - 0.5 * (gas._gamma - 1.) * jMinus * jMinus))
	    		/ (gas._gamma + 1.);

	    double faceSpeedOfSound = (gas._gamma - 1.) * 0.25 * (jPlus - jMinus);
	    double faceMach = 0.5 * (jPlus + jMinus) / faceSpeedOfSound;
	    double facePressure = _Pi / std::pow(1.0 + 0.5 * (gas._gamma - 1.) * faceMach * faceMach, 3.5 /*CHANGE THIS!*/);
	    double faceDensity = gas._gamma * facePressure / (faceSpeedOfSound * faceSpeedOfSound);
	    double faceTemperature = facePressure / (faceDensity * gas._R);
	    double faceVelocity = faceMach * faceSpeedOfSound;
	    double faceVelocityX = faceVelocity * _dir[0];
	    double faceVelocityY = faceVelocity * _dir[1];
	    double faceStagnationEnergy = faceDensity * (gas._Cv * faceTemperature +
	    	0.5 * faceVelocity * faceVelocity);

	    FieldVector faceState(faceDensity, faceDensity * faceVelocityX, faceDensity * faceVelocityY,
	    	faceStagnationEnergy);

	    // sanity check: the face state has the right Pi, Ti
	    double deltaPi = std::abs((faceState.ComputeStagnationPressure() -  _Pi) / _Pi);
	    double deltaTi = std::abs((faceState.ComputeStagnationTemperature() - _Ti) / _Ti);

	    assert(deltaPi < 0.0001 && deltaTi < 0.0001
	    	&& "BCStagnationInflow::ComputeFaceState() failed to computed face state");

        return faceState;
	}

protected:
	/// Members.
	const double _Pi;
	const double _Ti;
    const std::vector<double> _dir;
};

/// TODO
class BCOutflowPressure: public BoundaryConditionBase
{
public:
	/// Forward base class.
	using BaseClass = BoundaryConditionBase;
	using BaseClass::BaseClass;

	/// Constructor allowing to specify an outflow pressure (overrides BaseClass' default one).
	BCOutflowPressure(const CellContainer &cellContainer, const FaceContainer &bdryFaceContainer,
	    Stencil &bdryFaceStencil, double backPressure):
	    BaseClass(cellContainer, bdryFaceContainer, bdryFaceStencil), _backPressure(backPressure) { }

    /// Implementation of `BaseClass::ComputeFaceState`
	FieldVector ComputeFaceState(const Face &bdryFace, const Cell &cell, const FieldVector &cellCenteredState)
	    const override
	{
		/// TODO
		const auto gas = GasThermodynamicProperties(); //this->_refState.GetGas();

		// retrieve the normal
		const auto &normal = bdryFace.GetNormal();

		// compute the face to center vector for the attached element. TODO: this must be precomputed data
		double radiusXComp = bdryFace.GetCenter().GetXCoord() - cell.GetCenter().GetXCoord();
		double radiusYComp = bdryFace.GetCenter().GetYCoord() - cell.GetCenter().GetYCoord();

		// get the sign of the face normal as seen from the element (convention: pointing outwards)
		double radiusXCompSign = 1.0 - 2.0 * static_cast<double> (std::signbit(radiusXComp));
		double radiusYCompSign = 1.0 - 2.0 * static_cast<double> (std::signbit(radiusYComp));

		// scale the normal components
		double nx = std::abs(normal[0]) * radiusXCompSign;
		double ny = std::abs(normal[1]) * radiusYCompSign;

	    // get the interior pressure, density, speed of sound and entropy
		double pressure = cellCenteredState.ComputePressure();
		double density = cellCenteredState.GetDensity();
		double speedOfSound = cellCenteredState.ComputeSpeedOfSound();
		double entropy = pressure * std::pow(density, -gas._gamma);

		// compute the normal velocity
		double velocityX = cellCenteredState.ConservativeToPrimitive()[1];
		double velocityY = cellCenteredState.ConservativeToPrimitive()[2];
		double vN = velocityX * nx + velocityY * ny;

	    // outcoming invariant
	    double jPlus  = vN + 2.0 / (gas._gamma - 1.) * speedOfSound;

	    // face state
	    double faceDensity  = std::pow(this->_backPressure / entropy, 1. / gas._gamma);
	    double faceSpeedOfSound = std::sqrt(gas._gamma * this->_backPressure / faceDensity);
	    double faceVN = jPlus - 2. / (gas._gamma - 1.) * faceSpeedOfSound;
	    double faceVelocityX = velocityX + (faceVN - vN) * nx;
	    double faceVelocityY = velocityY + (faceVN - vN) * ny;
	    double faceTemperature = this->_backPressure / (faceDensity * gas._R);

		// mach number
		double mach = std::abs(vN) / speedOfSound;
	    if (mach > 1.0)
	    	return cellCenteredState;

	    double faceStagnationEnergy = faceDensity * (gas._Cv * faceTemperature + 0.5 * (
	    	faceVelocityX * faceVelocityX + faceVelocityY * faceVelocityY));

        FieldVector faceState(faceDensity, faceVelocityX * faceDensity, faceVelocityY * faceDensity,
        	faceStagnationEnergy);

	    // sanity check: the face state has the right pressure
	    double deltaP = std::abs((faceState.ComputePressure() -  _backPressure) / _backPressure);
	    assert(deltaP < 0.0001
	    	&& "BCOutflowPressure::ComputeFaceState() failed to computed face state");

        return faceState;
	}

protected:
	/// Members.
	const double _backPressure;
};


class BoundaryTreatmentMap
{
public:
	using StringVector = typename std::vector<std::string>;
	using ValueDict = typename std::map<std::string, double>;

	/// Constructor from values.
	BoundaryTreatmentMap(const std::string treatmentType, const StringVector markerTags,
	    const ValueDict valueDict): _treatmentType(treatmentType), _markerTags(markerTags),
	    _valueDict(valueDict) { }

	/// Copy constructor.
	BoundaryTreatmentMap(const BoundaryTreatmentMap &other)
	{
		this->_treatmentType = other.GetTreatmentType();
		this->_markerTags = other.GetMarkerTags();
		this->_valueDict = other.GetValueDict();
	}

	const std::string GetTreatmentType() const
	{
		return this->_treatmentType;
	}

	const StringVector GetMarkerTags() const
	{
		return this->_markerTags;
	}

	const ValueDict GetValueDict() const
	{
		return this->_valueDict;
	}

protected:
    /// Members (get copied).
	std::string _treatmentType;
	StringVector _markerTags;
	ValueDict _valueDict;
};

///// Boundary condition container class.
/*
 *
 */
class BoundaryConditionContainer: public std::vector<std::shared_ptr<BoundaryConditionBase> >
{
public:
	/// Define BC map type.
	/// Example of use: `BCMapping bcMapping {"BCWallInviscid", {"upper side", "lower side"}};`.
	using BoundaryConditionMap = typename std::vector<BoundaryTreatmentMap>;
	using StringVector = typename BoundaryTreatmentMap::StringVector;
	using ValueDict = typename BoundaryTreatmentMap::ValueDict;

	// Forward type from mesh class
	using FaceContainer = typename Mesh::FaceContainer;
	using CellContainer = typename Mesh::CellContainer;

	/// Default constructor, as the base class has a default one.
	/*
	 * Requires the boundary condition mapping.
	 */
	BoundaryConditionContainer(const BoundaryConditionMap &bcMapping) : _bcMapping(bcMapping)
	{
		// expect a non-empty bcMapping
		assert(!_bcMapping.empty() && "BC map is empty");
		// reserve the adequate size of this
		this->reserve(_bcMapping.size());
	}

	/// Compute method.
	/*
	 *
	 */
	void Compute(const CellContainer &cellContainer, const FaceContainer &bdryFaceContainer,
		Stencil &bdryFaceStencil)
	{
        // retrieve the keys of the BC map
		std::vector<std::string> requestedBCTreatments;
		for (auto it = this->_bcMapping.begin(); it != this->_bcMapping.end(); ++it)
			requestedBCTreatments.push_back(it->GetTreatmentType());

		// loop over all boundary condition treatments that have been asked
		for (int j = 0; j < requestedBCTreatments.size(); ++j)
		{
			// get the requested BC treatment
			std::string requestedBCTreatment = requestedBCTreatments[j];
			// get the marker tags for which the current BC treatment must be applied
			const StringVector &markerTags = this->_bcMapping[j].GetMarkerTags();
			// retrieve the value dict
			const ValueDict &valueDict = this->_bcMapping[j].GetValueDict();

			// initialize the message displaying information
			std::cout << std::endl;
			std::cout << "Boundary condition treatment " << requestedBCTreatment << " will be applied"
				" on boundary faces:";

			// start our face container
			FaceContainer bcFacesContainer;
			Stencil bcFaceStencil(bdryFaceContainer.size());

			for (int i = 0; i < bdryFaceContainer.size(); ++i)
			{
				const std::string &tag = bdryFaceContainer[i].GetTag();
				if (std::find(markerTags.begin(), markerTags.end(), tag) != markerTags.end())
				{
					const Face &face = bdryFaceContainer[i];
					const std::vector<int> &faceStencilEntry = bdryFaceStencil.GetStencil()[i];


					bcFacesContainer.push_back(face);
					bcFaceStencil.FillAt(face.GetCenter().GetNodeID(), faceStencilEntry);

					// display the boundary face ID
					std::cout << " " << face.GetCenter().GetNodeID();
				}
			}
			std::cout << std::endl;

			if (requestedBCTreatment == "BCExteriorState")
			{
				// check value dict
				meridian::CheckValueDict(valueDict, {"Mach", "Pressure", "Temperature", "AoA"},
					"Wrong value dict given to BCExteriorState");

				// retrieve data
				double mach = valueDict.at("Mach");
				double pressure = valueDict.at("Pressure");
				double temperature = valueDict.at("Temperature");
				double AoA = valueDict.at("AoA");

				// instantiate the requested boundary treatment
			    auto bc = std::make_shared<BCExteriorState>(cellContainer, bcFacesContainer,
			    	bdryFaceStencil, mach, pressure, temperature, AoA);

	            // finally, place the boundary condition that has just been created into this
			    this->push_back(bc);
			}

			else if (requestedBCTreatment == "BCWallInviscid")
			{
				// instantiate the requested boundary treatment
				 auto bc = std::make_shared<BCWallInviscid>(
					cellContainer, bcFacesContainer, bdryFaceStencil);

	            // finally, place the boundary condition that has just been created into this
			    this->push_back(bc);
			}

			else if (requestedBCTreatment == "BCStagnationInflow")
			{
				// check value dict
				meridian::CheckValueDict(valueDict,
					{"StagnationPressure", "StagnationTemperature", "AoA"},
					"Wrong value dict given to BCStagnationInflow");

				// retrieve data
				double Pi = valueDict.at("StagnationPressure");
				double Ti = valueDict.at("StagnationTemperature");
				double AoA = valueDict.at("AoA");

				// transform the AoA to vector
				std::vector<double> dir{std::cos(AoA * M_PI / 180.0), std::sin(AoA * M_PI / 180.0)};

				// instantiate the requested boundary treatment
				auto bc = std::make_shared<BCStagnationInflow>(
					cellContainer, bcFacesContainer, bdryFaceStencil,
					Pi, Ti, dir);

	            // finally, place the boundary condition that has just been created into this
		        this->push_back(bc);
			}

			else if (requestedBCTreatment == "BCOutflowPressure")
			{
				// check value dict
				meridian::CheckValueDict(valueDict, {"Pressure"},
					"Wrong value dict given to BCOutflowPressure");

				// retrieve data
				double pressure = valueDict.at("Pressure");

				// instantiate the requested boundary treatment
				auto bc = std::make_shared<BCOutflowPressure>(cellContainer, bcFacesContainer,
					bdryFaceStencil, pressure);

	            // finally, place the boundary condition that has just been created into this
		        this->push_back(bc);
			}

			// throw an exception otherwise
			else
				throw std::invalid_argument(
					"Requested boundary condition treatment misspelled or not implemented yet");
		}
	}

protected:
	/// Members.
	const BoundaryConditionMap & _bcMapping;
};
