/*
 * ******************  MERIDIAN CFD SOLVER ***************************
 *
 * Meridian is a 2D finite volume CFD solver of the Euler equations,
 * intended to solve easy-to-medium complexity problems with a basic
 * accuracy level.
 *
 * This code features a series of classes providing:
 *     1. Unstructured mesh representation, generated from SU2 format
 *        files.
 *     2. Space integration, using first order discretization
 *        (constant reconstruction) and both centered / upwinding
 *        schemes. Boundary conditions include solid walls, farfield
 *        state, stagnation inflow and backpressure outflow.
 *     3. Time integration, using the explicit backward Euler scheme.
 *     4. Input / Output functionalities, allowing to dump the mesh
 *        and the convergence history in tecplot format.
 *
 * Author: Luis Lopez de Vega (l.lopezdevega@gmail.com), 2022
 *
 * ********************************************************************
 */

#include <iostream>
#include "src/Meridian.cpp"

int main(int argc, char *argv[])
{
	// check that the right number of command line args has been passed
	if (argc < 2)
		throw std::invalid_argument("Please include a configuration file to the call of meridian!");

	// print logo for fun
	PrintLogo();

	// create an object that reads configuration file and has access to its information
	ConfigFileReader params = CreateConfigFileReader(argv[1]);

	// create the computational mesh
	Mesh mesh = CreateMesh(params.GetMeshFilename());

	// instantiate the gas model
	GasThermodynamicProperties gas = CreateGasModel(params.GetGasConstant(), params.GetGasSpecificHeatRatio());

	// create the initial state
	State solution = CreateInitialState(mesh, gas,
	    {{"Mach", params.GetFreeStreamMach()}, {"Pressure", params.GetFreeStreamPressure()},
	        {"Temperature", params.GetFreeStreamTemperature()}, {"AoA", params.GetFreeStreamAoA()}});

	// create the boundary container
	BoundaryConditionContainer::BoundaryConditionMap bcMap(params.GetBCMap());

    // create the boundary condition container
	BoundaryConditionContainer bcContainer = CreateBoundaryConditionContainer(mesh, bcMap);

	// create the convection scheme- and the residual computation
	HLLScheme hllScheme;
	ResidualComputation<HLLScheme> residualComputation(mesh.GetCellContainer(), mesh.GetInternalFaceContainer(),
		mesh.GetInternalFaceStencil(), bcContainer, hllScheme);

	// create the residual
	State residual = CreateInitialResidual(mesh, gas);

	// create the time step controller.
	TimeStepController timeStepController = CreateTimeStepController(
		mesh, {{"CFL", params.GetCFLNumber()}}, params.GetTimeStep());

	// create the time integration
	BackwardEuler timeIntegration = CreateTimeIntegration(mesh, timeStepController);

	// run the CFD loop
	RunCFD(residualComputation, timeIntegration, solution, residual, mesh,
		{{"Number of iterations", params.GetNbIterations()}});

	// successfull run
	std::cout << std::endl;
	std::cout << "End of meridian" << std::endl;

	return 0;
}
