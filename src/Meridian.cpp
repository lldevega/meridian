/// Meridian functions
/*
 *
 */

#include "../include/Utils/Utils.hpp"
#include "../include/SpaceIntegration/Mesh.hpp"
#include "../include/SpaceIntegration/State.hpp"
#include "../include/SpaceIntegration/BoundaryCondition.hpp"
#include "../include/SpaceIntegration/ResidualComputation.hpp"
#include "../include/SpaceIntegration/ConvectionScheme.hpp"
#include "../include/TimeIntegration/TimeIntegrationScheme.hpp"
#include "../include/SolverIO/SolverIO.hpp"

/// Generate the mesh representation.
/*
 * Requires the file name of an existing SU2 format mesh file.
 */
Mesh CreateMesh(const std::string &meshFilename, const bool displayInfo = false,
	const bool exportToTecplot = true)
{
	// init the SU2 mesh format reader
	SU2Reader reader(meshFilename);

	// read the mesh file
	reader.ReadMesh();

	// get containers
	auto &nodes = reader.GetNodeList();
	auto &elements = reader.GetElementList();
	auto &boundaries = reader.GetBoundaryList();

	// generate a mesh object
	Mesh mesh(nodes, elements, boundaries);
	mesh.Compute();

	// display into to screen
	if (displayInfo)
	    mesh.DisplayInfo();

	// dump to tecplot format (coordinates only)
	if (exportToTecplot)
	    mesh.ExportToTecplot("grid.dat");

	return mesh;
}

/// Create the gas model.
/*
 * Requires the gas constant R and the specific heat ratio gamma.
 */
GasThermodynamicProperties CreateGasModel(double R = 287.06, double gamma = 1.4)
{
	// instantiate class and return
	GasThermodynamicProperties gas(R, gamma);
	return gas;
}

/// Generate the state, holding the flow variables at each cell.
/*
 *
 */
State CreateInitialState(const Mesh &mesh, const GasThermodynamicProperties &gas,
	const meridian::ValueDict params)
{
	// number of entries
	int nEntries = mesh.GetCellContainer().size();

	// generate the state
	State solution(nEntries, gas);

	// get the reference state
	meridian::CheckValueDict(params, {"Mach", "Pressure", "Temperature", "AoA"},
		"Wrong value dict passed for creation of the initial state");

	double mach = params.at("Mach");
	double pressure = params.at("Pressure");
	double temperature = params.at("Temperature");
	double AoA = params.at("AoA");

	// initialize all entries to a constant state
	solution.SetState(mach, pressure, temperature, AoA);

	return solution;
}

/// Create the inital residual (evals to zero for all entries
/*
 *
 */
State CreateInitialResidual(const Mesh &mesh, const GasThermodynamicProperties &gas)
{
	State residual(mesh.GetCellContainer().size(), gas);
	residual.SetEntries();
	return residual;
}

///Create the boundary condition container
/*
 * Note that the user inputs (bcMaps) are checked at
 * BoundaryConditionContainer::Compute() level.
 */
BoundaryConditionContainer CreateBoundaryConditionContainer(Mesh &mesh,
	const BoundaryConditionContainer::BoundaryConditionMap bcMap)
{
	// instantiate the boundary condition container
    BoundaryConditionContainer bcContainer(bcMap);

    // compute the boundaries
    bcContainer.Compute(mesh.GetCellContainer(), mesh.GetBoundaryFaceContainer(),
        mesh.GetBoundaryFaceStencil());

    return bcContainer;
}

/// Create the time step controller.
/*
 *
 */
TimeStepController CreateTimeStepController(const Mesh &mesh, const meridian::ValueDict params)
{
	// check user inputs
	meridian::CheckValueDict(params, {"CFL", "Activate local time step"},
		"Wrong value dict passed to the time step controller");

	// get inputs
	double CFL = params.at("CFL");
	bool useLocalTimeStep = static_cast<bool>(params.at("Activate local time step"));

	// instantiate our class
	TimeStepController timeStepController(CFL, mesh.GetCellContainer());

	// compute the time step size per cell
	if (useLocalTimeStep)
	    timeStepController.ComputeLocalTimeStep();
	else
		timeStepController.ComputeGlobalTimeStep();

	return timeStepController;
}

/// Create the time integration.
/*
 *
 */
BackwardEuler CreateTimeIntegration(const Mesh &mesh, const TimeStepController &timeStepController)
{
	// instantiate and return
	BackwardEuler timeIntegration(mesh.GetCellContainer(), timeStepController);
	return timeIntegration;
}


/// Run the solver using the space and time integration methods.
/*
 * Requires a residual computation, a time integration, a solution field,
 * a residual field, a mesh object and the number of iterations to run.
 */
template <class ResidualComputationType, class TimeIntegrationType>
void RunCFD(ResidualComputationType &residualComputation,
	TimeIntegrationType &timeIntegration, State &solution, State &residual,
	Mesh &mesh, const meridian::ValueDict params)
{
	// check user inputs
	meridian::CheckValueDict(params, {"Number of iterations"},
		"Wrong value dict passed run the CFD loop");

	// get inputs
	int nIter = params.at("Number of iterations");

	// loop for the given number of iterations
	for (int iteration = 1; iteration < nIter + 1; ++ iteration)
	{
		// compute the residual field using the previous solution
		residualComputation.ComputeResidual(solution, residual);
		// update the solution using the integrated residual field
		timeIntegration.Update(iteration, residual, solution);
		// set the entries of the residual to zero
		residual.SetEntries();
	}

	// compute the residual using the final solution
    residualComputation.ComputeResidual(solution, residual);
    // export to tecplot format
	ExportSolutionToTecplot("solution.dat", mesh, solution, residual);
	timeIntegration.ExportConvergenceMonitor("convergence.dat");
}
