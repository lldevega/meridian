/*
* Field vector and state class.
*/

#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cassert>
#include <memory>

/// Gas thermodynamic properties structure.
/*
 * This struct holds the specific heat ratio and the gas constant
 */
struct GasThermodynamicProperties
{
	/// Constructor.
	GasThermodynamicProperties(double R, double gamma) : _R(R), _gamma(gamma),
	    _Cp(ComputeCp(gamma, R)), _Cv(ComputeCv(gamma, R)) { };

	/// Default constructor with air values.
	GasThermodynamicProperties() : _R(287.06), _gamma(1.4),
		_Cp(ComputeCp(1.4, 287.06)),  _Cv(ComputeCv(1.4, 287.06)) { };

	/// Get the constant pressure heat capacity.
	double ComputeCp(double gamma, double R)
	{
		return R * gamma / (gamma - 1.0);
	}

	/// Get the constant volume heat capacity.
	double ComputeCv(double gamma, double R)
	{
		return R / (gamma - 1.0);
	}

	/// The gas constant.
	double _R;
	/// The specific heat ratio.
	double _gamma;
	/// The constant pressure heat capacity.
	double _Cp;
	/// The constant volume heat capacity.
	double _Cv;
};

/// Field vector class.
/*
 * A field vector is the data structure holding the
 * state variables of a single point.
 */
class FieldVector
{
public:
	/// Constructor
	/*
	 * Requires the ID of the point and the gas thermodynamic properties object.
	 */
	FieldVector(int ID, double ro, double rou, double rov, double roE, const GasThermodynamicProperties &gas) :
		_ID(ID), _ro(ro), _rou(rou), _rov(rov), _roE(roE), _gas(gas) { }

	/// Constructor with default gas properties
	FieldVector(double ro, double rou, double rov, double roE) : _ID(-1), _ro(ro), _rou(rou), _rov(rov),
	    _roE(roE), _gas(GasThermodynamicProperties()) { }

	/// Default constructor
	FieldVector() : _gas(GasThermodynamicProperties()) { }

	/// Assignment operator from other instance of this class.
	void operator=(const FieldVector &other)
	{
        this->_ro = other.GetDensity();
        this->_rou = other.GetMomentumX();
        this->_rov = other.GetMomentumY();
        this->_roE = other.GetStagnationEnergy();
        this->_gas = other.GetGas();
	}

	/// Addition operator from other instance of this class.
	void operator+=(const FieldVector &other)
	{
        this->_ro += other.GetDensity();
        this->_rou += other.GetMomentumX();
        this->_rov += other.GetMomentumY();
        this->_roE += other.GetStagnationEnergy();
	}

	/// Subtraction operator from other instance of this class.
	void operator-=(const FieldVector &other)
	{
        this->_ro -= other.GetDensity();
        this->_rou -= other.GetMomentumX();
        this->_rov -= other.GetMomentumY();
        this->_roE -= other.GetStagnationEnergy();
	}

	/// Scaling this by a float value.
	void operator*=(double scalar)
	{
        this->_ro *= scalar;
        this->_rou *= scalar;
        this->_rov *= scalar;
        this->_roE *= scalar;
	}

	/// Scaling this by a float value.
	FieldVector operator*(double scalar)
	{
		FieldVector result = *this;
        result.SetDensity(result.GetDensity() * scalar);
        result.SetMomentumX(result.GetMomentumX() * scalar);
        result.SetMomentumY(result.GetMomentumY() * scalar);
        result.SetStagnationEnergy(result.GetStagnationEnergy() * scalar);

        return result;
	}

	/// Addition operator from other instance of this class.
	FieldVector operator+(const FieldVector &other)
	{
		FieldVector result = *this;
        result += other;

        return result;
	}

	/// Subtraction operator from other instance of this class.
	FieldVector operator-(const FieldVector &other)
	{
		FieldVector result = *this;
		result -= other;

        return result;
	}

	/// Square operator (useful for the residual computation).
	void Square()
	{
		this->_ro *= this->_ro;
		this->_rou *= this->_rou;
		this->_rov *= this->_rov;
		this->_roE *= this->_roE;
	}

	/// Square root operator (useful for the residual computation).
	void SquareRoot()
	{
		this->_ro = std::sqrt(this->_ro);
		this->_rou = std::sqrt(this->_rou);
		this->_rov = std::sqrt(this->_rov);
		this->_roE = std::sqrt(this->_roE);
	}

	/// Get and set methods of private members
	int GetID() const
	{
		return this->_ID;
	}

	const GasThermodynamicProperties &GetGas() const
	{
		return this->_gas;
	}

	double GetDensity() const
	{
		return this->_ro;
	}

	void SetDensity(double value)
	{
		this->_ro = value;
	}

	double GetMomentumX() const
	{
		return this->_rou;
	}

	void SetMomentumX(double value)
	{
		this->_rou = value;
	}

	double GetMomentumY() const
	{
		return this->_rov;
	}

	void SetMomentumY(double value)
	{
		this->_rov = value;
	}

	double GetStagnationEnergy() const
	{
		return this->_roE;
	}

	void SetStagnationEnergy(double value)
	{
        this->_roE = value;
	}

	/// Computation of thermodynamic quantities
	double ComputeVelocityX() const
	{
		return this->_rou / this->_ro;
	}

	double ComputeVelocityY() const
	{
		return this->_rov / this->_ro;
	}

	double ComputeTemperature() const
	{
		double sqrVelocity = ComputeSquaredVelocity();
		double temperature = (this->_roE / this->_ro - 0.5 * sqrVelocity) / this->_gas._Cv;

		return temperature;
	}

	double ComputePressure() const
	{
	    double temperature = ComputeTemperature();
		return this->_ro * this->_gas._R * temperature;
	}

	double ComputeSpeedOfSound() const
	{
	    double temperature = ComputeTemperature();
		return std::pow(this->_gas._gamma * this->_gas._R * temperature, 0.5);
	}

	double ComputeMach() const
	{
	    double temperature = ComputeTemperature();
	    double sqrVelocity = ComputeSquaredVelocity();
		return std::pow(sqrVelocity / (this->_gas._gamma * this->_gas._R * temperature), 0.5);
	}

	double ComputeStagnationTemperature() const
	{
		double temperature = ComputeTemperature();
		double Mach = ComputeMach();
		return temperature * (1.0 + 0.5 * (this->_gas._gamma - 1.0) * Mach * Mach);
	}

	double ComputeStagnationPressure() const
	{
		double pressure = ComputePressure();
		double Mach = ComputeMach();
		double exponent = this->_gas._gamma / (this->_gas._gamma - 1.0);
		return pressure * std::pow((1.0 + 0.5 * (this->_gas._gamma - 1.0) * Mach * Mach), exponent);
	}

	double ComputeStagnationEntalpy() const
	{
		return ComputeStagnationTemperature() * this->_gas._Cp;
	}

	/// Returns a vector holding the primitive state corresponding to this conservative state.
	/*
	 * Arbitrary ordered as follows; density, VelocityX, VelocitY, pressure.
	 */
	std::vector<double> ConservativeToPrimitive() const
	{
		double density = GetDensity();
		double VelocityX = ComputeVelocityX();
		double VelocityY = ComputeVelocityY();
		double pressure = ComputePressure();
		return std::vector<double> {density, VelocityX, VelocityY, pressure};
	}

protected:
	///Members
	int _ID;
	double _ro;
	double _rou;
	double _rov;
	double _roE;
	GasThermodynamicProperties _gas;

	/// Helper function returning the squared modulus of the velocity.
	double ComputeSquaredVelocity() const
	{
		return ComputeVelocityX() * ComputeVelocityX() +
			ComputeVelocityY() * ComputeVelocityY();
	}
};

/// Flux vector class
/*
 * This class provides the data structure to store the 2D-components of the flux vector.
 */
class FluxVector
{
public:
	/// Constructor.
	FluxVector(const FieldVector &fieldVector): _fluxro(ComputeDensityFluxVector(fieldVector)),
	    _fluxrou(ComputeMomentumXFluxVector(fieldVector)), _fluxrov(ComputeMomentumYFluxVector(fieldVector)),
		_fluxroE(ComputeStagnationEnergyFluxVector(fieldVector))
	{ }

	/// Assignment operator from other instance of this class.
	void operator=(const FluxVector &other)
	{
		for (int i = 0; i < 2; ++i)
		{
	        this->_fluxro[i] = other.GetDensityFlux()[i];
	        this->_fluxrou[i] = other.GetMomentumXFlux()[i];
	        this->_fluxrov[i] = other.GetMomentumYFlux()[i];
	        this->_fluxroE[i] = other.GetStagnationEnergyFlux()[i];
		}
	}

	/// Addition operator from other instance of this class.
	void operator+=(const FluxVector &other)
	{
		for (int i = 0; i < 2; ++i)
		{
	        this->_fluxro[i] += other.GetDensityFlux()[i];
	        this->_fluxrou[i] += other.GetMomentumXFlux()[i];
	        this->_fluxrov[i] += other.GetMomentumYFlux()[i];
	        this->_fluxroE[i] += other.GetStagnationEnergyFlux()[i];
		};
	}

	/// Substraction operator from other instance of this class.
	void operator-=(const FluxVector &other)
	{
		for (int i = 0; i < 2; ++i)
		{
	        this->_fluxro[i] -= other.GetDensityFlux()[i];
	        this->_fluxrou[i] -= other.GetMomentumXFlux()[i];
	        this->_fluxrov[i] -= other.GetMomentumYFlux()[i];
	        this->_fluxroE[i] -= other.GetStagnationEnergyFlux()[i];
		}
	}

	/// Multiplication by a scalar.
	void operator*=(double scalar)
	{
		for (int i = 0; i < 2; ++i)
		{
	        this->_fluxro[i] *= scalar;
	        this->_fluxrou[i] *= scalar;
	        this->_fluxrov[i] *= scalar;
	        this->_fluxroE[i] *= scalar;
		}
	}

	/// Sum operator.
	FluxVector operator+(const FluxVector &other)
	{
		FluxVector result = *this;

		for (int i = 0; i < 2; ++i)
		{
	        result.GetDensityFlux()[i] += other.GetDensityFlux()[i];
	        result.GetMomentumXFlux()[i] += other.GetMomentumXFlux()[i];
	        result.GetMomentumYFlux()[i] += other.GetMomentumYFlux()[i];
	        result.GetStagnationEnergyFlux()[i] += other.GetStagnationEnergyFlux()[i];
		}
		return result;
	}

	/// Subtraction operator.
	FluxVector operator-(const FluxVector &other)
	{
		FluxVector result = *this;

		for (int i = 0; i < 2; ++i)
		{
	        result.GetDensityFlux()[i] -= other.GetDensityFlux()[i];
	        result.GetMomentumXFlux()[i] -= other.GetMomentumXFlux()[i];
	        result.GetMomentumYFlux()[i] -= other.GetMomentumYFlux()[i];
	        result.GetStagnationEnergyFlux()[i] -= other.GetStagnationEnergyFlux()[i];
		}
		return result;
	}

	/// Scale operator.
	/// TODO: wrong!
	FluxVector operator*(double scalar)
	{
		FluxVector result = *this;

		for (int i = 0; i < 2; ++i)
		{
	        result.GetDensityFlux()[i] *= scalar;
	        result.GetMomentumXFlux()[i] *= scalar;
	        result.GetMomentumYFlux()[i] *= scalar;
	        result.GetStagnationEnergyFlux()[i] *= scalar;
		}
		return result;
	}

	/// Get the private members.
	std::vector<double> GetDensityFlux() const
	{
		return this->_fluxro;
	}

	/// Get the private members.
	std::vector<double> GetMomentumXFlux() const
	{
		return this->_fluxrou;
	}

	/// Get the private members.
	std::vector<double> GetMomentumYFlux() const
	{
		return this->_fluxrov;
	}

	/// Get the private members.
	std::vector<double> GetStagnationEnergyFlux() const
	{
		return this->_fluxroE;
	}

protected:
	/// Members.
	std::vector<double> _fluxro;
	std::vector<double> _fluxrou;
	std::vector<double> _fluxrov;
	std::vector<double> _fluxroE;

	/// Computation of members.
	std::vector <double> ComputeDensityFluxVector(const FieldVector &fieldVector)
	{
		double fluxro_i = fieldVector.GetMomentumX();
		double fluxro_j = fieldVector.GetMomentumY();

		return std::vector<double> {fluxro_i, fluxro_j};
	}

	std::vector <double> ComputeMomentumXFluxVector(const FieldVector &fieldVector)
	{
		// ro * u² + p
		double fluxrou_i = fieldVector.GetMomentumX() * fieldVector.GetMomentumX() / fieldVector.GetDensity()
			+ fieldVector.ComputePressure();
		// ro * u * v
		double fluxrou_j = fieldVector.GetMomentumX() * fieldVector.GetMomentumY() / fieldVector.GetDensity();

		return std::vector<double> {fluxrou_i, fluxrou_j};
	}

	std::vector <double> ComputeMomentumYFluxVector(const FieldVector &fieldVector)
	{
		// ro * u * v
		double fluxrov_i = fieldVector.GetMomentumY() * fieldVector.GetMomentumX() / fieldVector.GetDensity();

		// ro * v² + p
		double fluxrov_j = fieldVector.GetMomentumY() * fieldVector.GetMomentumY() / fieldVector.GetDensity()
			+ fieldVector.ComputePressure();

		return std::vector<double> {fluxrov_i, fluxrov_j};
	}

	std::vector <double> ComputeStagnationEnergyFluxVector(const FieldVector &fieldVector)
	{
		double fluxroE_i = fieldVector.GetMomentumX() * fieldVector.ComputeStagnationEntalpy();
		double fluxroE_j = fieldVector.GetMomentumY() * fieldVector.ComputeStagnationEntalpy();

		return std::vector<double> {fluxroE_i, fluxroE_j};
	}

};

/// State class.
/*
 * The state class is the data structure holding all the Field Vector instances
 * of our mesh.
 */
class State: public std::vector<FieldVector>
{
public:
	/// Default constructor, as the base class has a default one.
	/*
	 * Requires the number of entries and an instance of the gas thermodynamic properties.
	 */
	State(int nEntries, const GasThermodynamicProperties &gas) : _nEntries(nEntries), _gas(gas)
	{
		this->resize(nEntries);
	}

	/// Sets the conservative variables to zero for all entries.
	void SetEntries()
	{
		for (int i = 0; i < this->size(); ++i)
		{
			FieldVector entry(i,
				0. /*Density*/, 0./*MomentumX*/, 0./*MomentumY*/, 0./*StagnationEnergy*/,
				this->_gas);
			this->at(i) = entry;
		}
	}

	/// Uniformly sets the conservative states according to a stagnation state.
	/*
	 * The conservative state is deduced from the Mach number, the pressure,
	 * the temperature and the angle of attack.
	 */
	void SetState(double Mach, double pressure, double temperature, double AoA /*in degrees*/)
	{
		// check first if we need to initialize the state
		if (!IsInitialized())
			Initialize(_nEntries);

		// instantiate the field vectors that will be the entries of this
		SetEntries();

		// compute the associated state
		double velocity = Mach * std::pow(this->_gas._gamma * this->_gas._R * temperature, 0.5);
		double sinAoA = std::sin(AoA * M_PI / 180.0);
		double cosAoA = std::cos(AoA * M_PI / 180.0);
		double velocityX = velocity * cosAoA;
		double velocityY = velocity * sinAoA;
		double density = pressure / (this->_gas._R * temperature);
		double stagnationEnergy = density * this->_gas._Cv * temperature + 0.5 * density * velocity * velocity;

		// loop over all entries and set their conservative state to the compute one
		for (FieldVector &fieldVector : *this)
		{
			fieldVector.SetDensity(density);
			fieldVector.SetMomentumX(density * velocityX);
			fieldVector.SetMomentumY(density * velocityY);
			fieldVector.SetStagnationEnergy(stagnationEnergy);
		}
	}

	/// Get field vector at specific point.
	const FieldVector &GetFieldVectorAt(int idx) const
	{
		return this->at(idx);
	}

	/// Set field vector at specific point.
	void SetFieldVectorAt(int idx, const FieldVector other)
	{
		this->at(idx) = other;
	}

	/// Set field vector at specific point.
	void AddFieldVectorAt(int idx, const FieldVector &other)
	{
		this->at(idx) += other;
	}

	/// Get private members.
	int GetNumEntries() const
	{
		return this->_nEntries;
	}

	/// Get private members.
	const GasThermodynamicProperties &GetGas() const
	{
		return this->_gas;
	}

protected:
	/// Members.
	int _nEntries;
	const GasThermodynamicProperties &_gas;

	/// Returns true if this has at least one entry allocated.
	bool IsInitialized()
	{
		return this->size() > 0;
	}

	/// Resizes this to the number of entries required.
	void Initialize(int nEntries)
	{
		this->resize(nEntries);
	}
};
