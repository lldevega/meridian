/*
* Time integration classes.
*/

#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cassert>
#include <memory>
#include <iomanip>

/// Time step controller class.
/*
 * This class provides the time step to be applied at each cell, based on the CFL number,
 * a reference mach number and the user criteria (either local or global time step).
 */
class TimeStepController
{
public:
	/// Forward type from mesh.
	using CellContainer = typename Mesh::CellContainer;

    /// Constructor.
	/*
	 * Requires the CFL and the number of entries (equal to the number of cells in the mesh).
	 */
	TimeStepController(double CFL, const CellContainer &cellContainer) : _CFL(CFL), _cellContainer(cellContainer),
		_timeStep()

	{
		InitTimeStep(cellContainer.size());
	}

	/// Compute the local time step
	/*
	 * The time step is computed for each cell based on its reference length.
	 * Local time step must only be used for reaching steady states.
	 */
	void ComputeLocalTimeStep() /*,const FieldVector &refState)*/
	{
		//double refMach = refState.ComputeMach();
        double refSoundSpeed = 340.0;

		for (int i = 0; i < this->_cellContainer.size(); ++i)
		{
			double refLength = std::pow(this->_cellContainer[i].GetSurface(), 0.5);
			double dt = this->_CFL * refLength / refSoundSpeed;
			this->_timeStep[i] = dt;
		}
	}

	/// Compute the global time step
	/*
	 * The time step is set globally for all the mesh cells, using the most restrictive (smallest)
	 * value.
	 */
	void ComputeGlobalTimeStep()/*, const FieldVector &refState)*/
	{
		//double refMach = refState.ComputeMach();
        double refSoundSpeed = 340.0;
        double dt = 1.0e10;

		for (int i = 0; i < this->_cellContainer.size(); ++i)
		{
			double refLength = std::pow(this->_cellContainer[i].GetSurface(), 0.5);
			double dtLocal = this->_CFL * refLength / refSoundSpeed;

			if (dtLocal < dt)
				dt = dtLocal;
		}

		for (int i = 0; i < this->_cellContainer.size(); ++i)
			this->_timeStep[i] = dt;
	}

	/// Set physical time step
	/*
	 * Sets an user defined time step
	 */
	void SetPhysicalTimeStep(double dt)
	{
		for (int i = 0; i < this->_cellContainer.size(); ++i)
			this->_timeStep[i] = dt;
	}

    /// Get time step.
	const std::vector<double> &GetTimeStep() const
	{
		return this->_timeStep;
	}


protected:
	/// Members.
	double _CFL;
	const CellContainer &_cellContainer;
	std::vector<double> _timeStep;

	/// Initialize the time step vector
	void InitTimeStep(int nEntries)
	{
		this->_timeStep.reserve(nEntries);
	}
};

/// Base class for time integration methods.
/*
 * Purely virtual class.
 */
class TimeIntegrationSchemeBase
{
public:
	/// Forward type from mesh class.
	using CellContainer = typename Mesh::CellContainer;
	using ConvergenceMonitor = typename std::vector<FieldVector>;

	/// Constructor.
	/*
	 * Requires the cell container and the time step controller.
	 */
	TimeIntegrationSchemeBase(const CellContainer &cellContainer, const TimeStepController &timeStepController):
		_cellContainer(cellContainer), _timeStepController(timeStepController), _convergenceMonitor (ConvergenceMonitor())
	{ }

	/// Purely virtual method to update the solution based on the compute residual.
	virtual void Update(int iteration, const State &residual, State &solution) = 0;

	/// Get convergence monitor
	const ConvergenceMonitor &GetConvergenceMonitor()
	{
		return this->_convergenceMonitor;
	}

	void PrintHeader()
	{
		// print header
		std::cout<< std::endl;
		std::cout<< "It   " << "ro           " << "rou          " << "rov          " << "roE " << std::endl;
		std::cout<< "-------------------------------------------------------" << std::endl;
	}

	/// Write the convergence monitor.
	/*
	 * Requires the filename to write the residual evolution.
	 */
	void ExportConvergenceMonitor(const char *filename)
	{
		// open the output file
		std::fstream outfile;
		outfile.open(filename, std::ios::out);

		// header
		outfile << "VARIABLES=" << "\t" << "Iteration" << "\t" << "DensityResidual" << "\t"
			    << "MomentumXResidual" << "\t" << "MomentumYResidual" << "\t"
				<< "StagnationEnergyResidual" << "\n";

		// write lines
		for(int i = 0; i < this->_convergenceMonitor.size(); i++)
		{
			outfile << i + 1 /*iteration number starting at 1 */ << "\t"
			    << this->_convergenceMonitor[i].GetDensity() / this->_convergenceMonitor[0].GetDensity() << "\t"
			    << this->_convergenceMonitor[i].GetMomentumX() / this->_convergenceMonitor[0].GetMomentumX()<< "\t"
				<< this->_convergenceMonitor[i].GetMomentumY() / this->_convergenceMonitor[0].GetMomentumY()<< "\t"
				<< this->_convergenceMonitor[i].GetStagnationEnergy() / this->_convergenceMonitor[0].GetStagnationEnergy()
				<< "\n";

		}
		outfile.close();
	}

protected:
	/// Members.
	const TimeStepController &_timeStepController;
	const CellContainer &_cellContainer;
	ConvergenceMonitor _convergenceMonitor;

	/// Print current residual norm reduction.
	void PrintResidualReduction(int iteration)
	{
		assert(this->_convergenceMonitor.size() > 0 && "Convergence monitor is empty");

		// comopute norm reduction
        double roResidualNormReduction = this->_convergenceMonitor[iteration - 1].GetDensity()
		    / this->_convergenceMonitor[0].GetDensity();
		double rouResidualNormReduction = this->_convergenceMonitor[iteration - 1].GetMomentumX()
			/ this->_convergenceMonitor[0].GetMomentumX();
		double rovResidualNormReduction = this->_convergenceMonitor[iteration - 1].GetMomentumY()
			/ this->_convergenceMonitor[0].GetMomentumY();
		double roEResidualNormReduction = this->_convergenceMonitor[iteration - 1].GetStagnationEnergy()
			/ this->_convergenceMonitor[0].GetStagnationEnergy();

		// print to screen
		std::cout << iteration << "    ";
		std::cout << std::setprecision(5) << std::scientific;
		std::cout << roResidualNormReduction << "  " << rouResidualNormReduction
			<< "  " << rovResidualNormReduction << "  " << roEResidualNormReduction << std::endl;
	}
};

/// Implementation of the first order backward euler time integration method
/*
 * surf * (Q^n+1 - Q^n) / dt = -R(Q^n) ==> Q^n+1 = Q^n -R(Q^n) * dt / surf
 */
class BackwardEuler: public TimeIntegrationSchemeBase
{
public:
	/// Forward constructor.
	using BaseClass = TimeIntegrationSchemeBase;
	using BaseClass::BaseClass;

	/// Implementation of the time integration update.
	void Update(int iteration, const State &residual, State &solution) override
	{
		// init the norm of the residuals
		FieldVector residualNorm(0., 0., 0., 0.);

		// loop over all entries
		for (int i = 0; i < solution.size(); ++i)
		{
			// compute update
			FieldVector update = residual[i];
			update *= (_timeStepController.GetTimeStep()[i] / _cellContainer[i].GetSurface());

			// add to solution
		    solution[i] += update;

		    // compute the residual norm
		    FieldVector squaredResidual = residual[i];
		    squaredResidual.Square();
		    residualNorm += squaredResidual;
		}

		// now the residual norm can be computed and stored
		residualNorm.SquareRoot();
	    this->_convergenceMonitor.push_back(residualNorm);

	    // display residual reduction
	    PrintResidualReduction(iteration);
	}
};
