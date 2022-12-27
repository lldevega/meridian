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

int main()
{
	// print logo for fun
	PrintLogo();

	// create the computational mesh
	//Mesh mesh = CreateMesh("examples/pipe/pipe.su2");
	Mesh mesh = CreateMesh("examples/cylinder/cylinder.su2");
	//Mesh mesh = CreateMesh("examples/naca0012/naca0012.su2");

	// instantiate the gas model
	GasThermodynamicProperties gas = CreateGasModel(287.06, 1.4);

	// create the initial state
	State solution = CreateInitialState(mesh, gas,
	    {{"Mach", 0.2}, {"Pressure", 101325.0}, {"Temperature", 288.15}, {"AoA", 0.0}});

	// create the boundary container
//	BoundaryConditionContainer::BoundaryConditionMap bcMap{
//	   {"BCWallInviscid", {"upper", "lower"}, {} },
//	   {"BCStagnationInflow", {"inlet"}, {{"StagnationPressure", 101325.0 * 1.02828}, {"StagnationTemperature", 288.15*1.008}, {"AoA", 0.0}}},
//	   {"BCOutflowPressure", {"outlet"}, {{"Pressure", 101325.0}}}
//	};
	BoundaryConditionContainer::BoundaryConditionMap bcMap{
	   {"BCWallInviscid", {"cylinder"}, {} },
	   {"BCExteriorState", {"farfield"}, {{"Mach", 0.2}, {"Pressure", 101325.0}, {"Temperature", 288.15}, {"AoA", 0.0}} },
	};

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
		mesh, {{"CFL", 0.4}, {"Activate local time step", 1.}}
	);

	// create the time integration
	BackwardEuler timeIntegration = CreateTimeIntegration(mesh, timeStepController);

	// run the CFD loop
	RunCFD(residualComputation, timeIntegration, solution, residual, mesh,
		{{"Number of iterations", 5000}});

	return 0;
}
