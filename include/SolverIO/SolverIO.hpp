/*
* Solver Input/Output functionalities
*/

#include <string.h>

/// Base class for reading text files and keeping information as class members
/*
 *
 */
class FileReaderBase
{
public:
	/// Constructor.
	FileReaderBase(const std::string name) : _inputFilename(name) { }

	/// Virtual method for reading the file.
	/*
	 * This method implements the necessary actions to read the contents of the file and store the relevant
	 * information as class members.
	 */
	virtual void Read() = 0;

protected:
	/// Members.
	const std::string _inputFilename;

    /// Method reading and storing numeric entries from a file.
	/*
	 *
	 */
    template <class ContainerType>
    void String2VectorOfNumerics(const std::string& str, ContainerType& cont, char delimiter = ' ')
    {
    	std::stringstream ss(str);
    	std::string token;
    	while(std::getline(ss, token, delimiter))
    	{
    		cont.push_back(std::stod(token.c_str()));
    	}
    }

    /// Method reading and storing string entries from a file.
    /*
     *
     */
    template <class ContainerType>
    void String2VectorOfString(const std::string& str, ContainerType& cont, char delimiter = ' ')
    {
    	std::stringstream ss(str);
    	std::string token;
    	while(std::getline(ss, token, delimiter))
    	{
    		cont.push_back(token.c_str());
    	}
    }
};

/// Reads the configuration file.
class ConfigFileReader: public FileReaderBase
{
public:
	// Type definition.
	using StringVector = typename std::vector<std::string>;
	using ValueDict = typename std::map<std::string, double>;

	/// Forward constructor.
	using BaseClass = FileReaderBase;
	using BaseClass::BaseClass;

	/// Reads the configuration file and stores all the information in members.
	void Read() override
	{
    	std::ifstream file;
    	std::string line;

    	file.open(this->_inputFilename.c_str());

    	if (file.is_open())
    	{
    		// set booleans to false
    		bool readingMeshParams = false;
    		bool readingGasParams = false;
    		bool readingFreeStreamParams = false;
    		bool readingBoundaryConditionParams = false;
    		bool readingSpatialIntegrationParams = false;
    		bool readingTimeIntegrationParams = false;

    		bool readingBC = false;
			std::string bcName;
			StringVector markers;
			ValueDict valueDict;

    		// read lines of the input file
    		while (getline(file, line))
    		{
    			std::vector<std::string> vect;
    			this->String2VectorOfString(line, vect);

    	       	// comment mark is % for the configuration file
    	        if (vect[0] == "%")
    	        	continue;

    	        // select the current reading format
    			SelectReader(vect, readingMeshParams, readingGasParams, readingFreeStreamParams,
    				readingBoundaryConditionParams, readingSpatialIntegrationParams, readingTimeIntegrationParams);

    			// case reading the mesh params
    			if (readingMeshParams && vect[0] != "[MESH]")
    			{
    				// we only have one parameter, which is the MeshFile
    				if (vect[0] == "MeshFile=")
    					_meshFilename = vect[1];
    				// otherwise we skip the line
    				else
    					SetAndSkip(vect[0], "[MESH]");
    			}
    			// case reading the gas properties
    			else if (readingGasParams && vect[0] != "[GAS_PROPERTIES]")
    			{
    				// gas constant, R and gamma
    				if (vect[0] == "R=")
    					_gasConstant = std::stod(vect[1]);
    				else if (vect[0] == "Gamma=")
    					_gasSpecificHeatRatio = std::stod(vect[1]);
    				// otherwise we skip the line
    				else
    					SetAndSkip(vect[0], "[GAS_PROPERTIES]");
    			}
    			// case reading the free stream params
    			else if (readingFreeStreamParams && vect[0] != "[FREESTREAM_STATE]")
    			{
    				// pressure, temperature, mach and AoA
    				if (vect[0] == "Pressure=")
    					_freeStreamPressure = std::stod(vect[1]);
    				else if (vect[0] == "Temperature=")
    					_freeStreamTemperature = std::stod(vect[1]);
    				else if (vect[0] == "Mach=")
    					_freeStreamMach = std::stod(vect[1]);
    				else if (vect[0] == "AoA=")
    					_freeStreamAoA = std::stod(vect[1]);
    				// otherwise we skip the line
    				else
    					SetAndSkip(vect[0], "[FREESTREAM_STATE]");
    			}
    			// case reading the boundary conditions
    			else if (readingBoundaryConditionParams && vect[0] != "[BOUNDARY_CONDITIONS]")
    			{

    				// set the boundary condition we are currently reading
    				if (IsBoundaryTreatmentAvailable(vect[0]) && !readingBC)
    				{
    					bcName = SplitBoundaryTreatmentName(vect[0]);
    					readingBC = true;
    				}

    				// general case
    				else if (readingBC)
    				{
    					// read the markers
    					if (vect[0] == "Markers=")
    					{
    						markers.reserve(10);
    						for (int i = 1; i < vect.size(); ++i)
    							markers.push_back(vect[i]);
    					}
    					// read the variables
    					else if (vect[0] == "Variables=")
    					{
    						for (int i = 1; i < vect.size(); ++i)
    							valueDict.insert({vect[i], -1.0});
    					}
    					else if (vect[0] == "Values=")
    					{
    						int i = 1;
    						for (auto &variable : valueDict)
    						{
    							variable.second = std::stod(vect[i]);
    							++i;
    						}
    					}
    					// if the new line has another BC name, then create it
    				    else if (IsBoundaryTreatmentAvailable(vect[0])
    				    	|| vect[0] == "[END_BOUNDARY_CONDITIONS]")
    					{
    						// check that at least there is a marker
    						assert(markers.size() > 0 && "No markers specified for boundary condition...");
    						// fill in the boundary condition map with this entry
    						_bcMap.push_back(BoundaryTreatmentMap(bcName, markers, valueDict));
    						// set to false to allow reading another bc
    						readingBC = true;
    						bcName = SplitBoundaryTreatmentName(vect[0]);
    						// reset markers
    						markers.resize(0);
    						// reset value dict
    						valueDict = ValueDict();
    					}
    				}
    				// otherwise we skip the line
    				else
    				   SetAndSkip(vect[0], "[BOUNDARY_CONDITIONS]");
    			}
    			// case reading the spatial scheme
    			else if (readingSpatialIntegrationParams && vect[0] != "[SPATIAL_INTEGRATION]")
    			{
    				// pressure, temperature, mach and AoA
    				if (vect[0] == "ConvectionScheme=")
    					_convectionSchemeName = vect[1];
    				// otherwise we skip the line
    				else
    					SetAndSkip(vect[0], "[SPATIAL_INTEGRATION]");
    			}
    			// case reading the time integration scheme
    			else if (readingTimeIntegrationParams && vect[0] != "[TIME_INTEGRATION]")
    			{
    				// pressure, temperature, mach and AoA
    				if (vect[0] == "CFL=")
    					_CFLNumber = std::stod(vect[1]);
    				else if(vect[0] == "NbIterations=")
    					_nbIterations = std::stod(vect[1]);
    				else if(vect[0] == "TimeStep=")
    					_timeStep = vect[1];
    				// otherwise we skip the line
    				else
    					SetAndSkip(vect[0], "[TIME_INTEGRATION]");
    			}
    		}
    	}
	}

	/// Get private members.
	std::string GetMeshFilename()
	{
		return this->_meshFilename;
	}

	double GetGasConstant()
	{
		return this->_gasConstant;
	}

	double GetGasSpecificHeatRatio()
	{
		return this->_gasSpecificHeatRatio;
	}

	double GetFreeStreamPressure()
	{
		return this->_freeStreamPressure;
	}

	double GetFreeStreamTemperature()
	{
		return this->_freeStreamTemperature;
	}

	double GetFreeStreamMach()
	{
		return this->_freeStreamMach;
	}

	double GetFreeStreamAoA()
	{
		return this->_freeStreamAoA;
	}

	BoundaryConditionContainer::BoundaryConditionMap GetBCMap()
	{
		return this->_bcMap;
	}

	std::string GetConvectionSchemeName()
	{
		return this->_convectionSchemeName;
	}

	double GetCFLNumber()
	{
		return this->_CFLNumber;
	}

	int GetNbIterations()
	{
		return this->_nbIterations;
	}

	std::string GetTimeStep()
	{
		return this->_timeStep;
	}

	void PrettyPrint()
	{
		std::cout << std::endl;
		std::cout << "    Mesh:" << std::endl;
		std::cout << "        File: " << GetMeshFilename() << std::endl;
		std::cout << "    Gas properties: " << std::endl;
		std::cout << "        Gas constant [J/(kg K)]: " << GetGasConstant() <<std::endl;
		std::cout << "        Gas specific heat ratio [-]: " << GetGasSpecificHeatRatio() <<std::endl;
		std::cout << "    Free stream state: " << std::endl;
		std::cout << "        Pressure [Pa]: " << GetFreeStreamPressure() <<std::endl;
		std::cout << "        Temperature [K]: " << GetFreeStreamTemperature() <<std::endl;
		std::cout << "        Mach [-]: " << GetFreeStreamMach() <<std::endl;
		std::cout << "        AoA [deg]: " << GetFreeStreamAoA() <<std::endl;
		std::cout << "    Spatial integration: " << std::endl;
		std::cout << "        Convection scheme: " << GetConvectionSchemeName() <<std::endl;
		std::cout << "    Time integration: " << std::endl;
		std::cout << "        CFL number: " << GetCFLNumber() <<std::endl;
		std::cout << "        Number of iterations: " << GetNbIterations() <<std::endl;
		std::cout << "        Time step: " << GetTimeStep() <<std::endl;
		std::cout << "    Boundary conditions: " << std::endl;
		for (int i = 0; i < _bcMap.size(); ++i)
		{
			std::cout << "        Treatment type: " << _bcMap[i].GetTreatmentType() <<std::endl;
		    std::cout << "            Markers: ";
		    for (int j = 0; j < _bcMap[i].GetMarkerTags().size(); ++j)
		    	std::cout << _bcMap[i].GetMarkerTags()[j] << " ";
		    std::cout << std::endl;
		    std::cout << "            Boundary state: " << std::endl;
		    for (auto state : _bcMap[i].GetValueDict())
		    	std::cout << "                 " << state.first << ": " << state.second << std::endl;
		}
		std::cout << std::endl;

	}

protected:
	/// Members.
	// name of the mesh file (SU2 format)
	std::string _meshFilename;

	// gas properties
	double _gasConstant;
	double _gasSpecificHeatRatio;

	// free stream values
	double _freeStreamPressure;
	double _freeStreamTemperature;
	double _freeStreamMach;
	double _freeStreamAoA;

	// boundary condition params
	BoundaryConditionContainer::BoundaryConditionMap _bcMap;
	const std::vector<std::string> _availableBoundaryTreatment {"[[BCWallInviscid]]", "[[BCExteriorState]]",
		"[[BCStagnationInflow]]", "[[BCOutflowPressure]]"};

	// convection scheme name
	std::string _convectionSchemeName;

	// time integration params
	double _CFLNumber;
	int _nbIterations;
	std::string _timeStep;

    void SelectReader(const std::vector<std::string> &vect, bool &readingMeshParams, bool &readingGasParams,
    	bool &readingFreeStreamParams, bool &readingBoundaryConditionParams, bool &readingSpatialIntegrationParams,
		bool &readingTimeIntegrationParams)
    {
    	if (vect[0] == "[MESH]")
    	{
    		readingMeshParams = true;
    		readingGasParams = false;
    		readingFreeStreamParams = false;
    		readingBoundaryConditionParams = false;
    		readingSpatialIntegrationParams = false;
    		readingTimeIntegrationParams = false;
    	}

    	else if (vect[0] == "[GAS_PROPERTIES]")
    	{
    		readingMeshParams = false;
    		readingGasParams = true;
    		readingFreeStreamParams = false;
    		readingBoundaryConditionParams = false;
    		readingSpatialIntegrationParams = false;
    		readingTimeIntegrationParams = false;
    	}

    	else if (vect[0] == "[FREESTREAM_STATE]")
    	{
    		readingMeshParams = false;
    		readingGasParams = false;
    		readingFreeStreamParams = true;
    		readingBoundaryConditionParams = false;
    		readingSpatialIntegrationParams = false;
    		readingTimeIntegrationParams = false;
    	}

    	else if (vect[0] == "[BOUNDARY_CONDITIONS]")
    	{
    		readingMeshParams = false;
    		readingGasParams = false;
    		readingFreeStreamParams = false;
    		readingBoundaryConditionParams = true;
    		readingSpatialIntegrationParams = false;
    		readingTimeIntegrationParams = false;
    	}

    	else if (vect[0] == "[SPATIAL_INTEGRATION]")
    	{
    		readingMeshParams = false;
    		readingGasParams = false;
    		readingFreeStreamParams = false;
    		readingBoundaryConditionParams = false;
    		readingSpatialIntegrationParams = true;
    		readingTimeIntegrationParams = false;
    	}

    	else if (vect[0] == "[TIME_INTEGRATION]")
    	{
    		readingMeshParams = false;
    		readingGasParams = false;
    		readingFreeStreamParams = false;
    		readingBoundaryConditionParams = false;
    		readingSpatialIntegrationParams = false;
    		readingTimeIntegrationParams = true;
    	}
    }

    /// Writes a message to screen when a key name is skipped
    /*
     *
     */
    void SetAndSkip(const std::string &keyName, const std::string &sectionName)
    {
    	std::cout << "    Entry " << keyName << " of " << sectionName << " not understood and ignored..." << std::endl;
    }

    /// Returns true if a given boundary treatment is available
    bool IsBoundaryTreatmentAvailable(const std::string keyName)
    {
    	bool isAvailable = false;
    	for (int i=0; i < _availableBoundaryTreatment.size(); ++i)
    		isAvailable = isAvailable || keyName == _availableBoundaryTreatment[i];

    	return isAvailable;
    }

    /// Retrieve the boundary treatment name from file
    const std::string SplitBoundaryTreatmentName(const std::string &keyName)
    {
    	// rear delimiter and split
        std::string delimiter1 = "]]";
        std::string token = keyName.substr(0, keyName.find(delimiter1));
        // front delimiter and split
        std::string delimiter2 = "[[";
        token = token.substr(token.find(delimiter2) + 2, token.find(delimiter2) + token.length());

        return token;
    }
};

///
/*
 *
 */
class SU2Reader: public FileReaderBase
{
public:
    /// Forward constructor.
    using BaseClass = FileReaderBase;
    using BaseClass::BaseClass;

    /// Reads a SU2 mesh and stores the elements, nodes and boundary faces in vectors.
    void Read() override
    {
    	std::ifstream file;
    	std::string line;

    	file.open(this->_inputFilename.c_str());
    	if (file.is_open())
    	{
    		bool readingNodes = false;
    		bool readingElems = false;
    		bool readingBdrys = false;
    		std::string markerTag;

    		while (getline(file, line))
    		{
    			// start by selecting reader
    			std::vector<std::string> vect;
    			this->String2VectorOfString(line, vect);
    			SelectReader(vect, readingNodes, readingElems, readingBdrys, markerTag);

    			if (vect[0] == "NMARK=")
    			{}
    			else if (vect[0] == "MARKER_TAG=")
    			{}
    			else if (vect[0] == "MARKER_ELEMS=")
    			{}
    			else if (vect[0] == "NPERIODIC=")
    			{
    				readingElems = false;
    				readingNodes = false;
    				readingBdrys = false;
    			}
    			// case reading elements
    			else if (readingElems && !readingNodes && !readingBdrys && vect[0] != "NELEM=")
    			{
    				std::vector<int> vect;
    				this->String2VectorOfNumerics(line, vect);
    				// check that the type of element is supported
    				bool isElementTypeSupported = vect[0] == 9 || vect[0] == 5;
    				if (!isElementTypeSupported)
    					std::runtime_error("Found not supported element type in the mesh.");
    				_meshElements.push_back(Element(vect));
    			}
    			// case reading nodes
    			else if (!readingElems && readingNodes && !readingBdrys && vect[0] != "NPOIN=")
    			{
    				std::vector<double> vect;
    				this->String2VectorOfNumerics(line, vect);
    				_meshNodes.push_back(Node(vect[0], vect[1], vect[2]));
    			}
    			// case reading boundaries
    			else if (!readingElems && !readingNodes && readingBdrys)
    			{
    				std::vector<int> vect;
    				this->String2VectorOfNumerics(line, vect);
    				_meshBoundaries.push_back(Edge(vect, markerTag));
    			}
    			else
    			{
    				// nothing to do here...
    			}
    		}
    	}
    	else
    	{
    		std::runtime_error("Unable to read mesh");
    	}
    }

    /// Get the mesh element list.
    std::vector<Element> &GetElementList()
    {
    	return this->_meshElements;
    }

    /// Get the mesh element list.
    const std::vector<Node> &GetNodeList()
    {
    	return this->_meshNodes;
    }

    /// Get the mesh element list.
    const std::vector<Edge> &GetBoundaryList()
    {
    	return this->_meshBoundaries;
    }

protected:
    /// Members.
    std::vector<Element> _meshElements;
    std::vector<Node> _meshNodes;
    std::vector<Edge> _meshBoundaries;

    /// Selects the mesh entries being read
    void SelectReader(const std::vector<std::string> vect,
    	bool &readingNodes, bool &readingElems, bool &readingBdrys, std::string &markerTag)
    {
    	// check whether dimensionality is correct
    	if (vect[0]=="NDIME=" && vect[1]!="2")
    		std::runtime_error("Mesh dimensionality is not 2D");

    	// identify whether reading elements
    	if (vect[0]=="NELEM=")
    	{
    		readingElems = true;
    		readingNodes = false;
    		readingBdrys = false;
    	}

    	// identify whether reading nodes
    	if (vect[0]=="NPOIN=")
    	{
    		readingElems = false;
    		readingNodes = true;
    		readingBdrys = false;
    	}

    	// identify whether reading markers
    	if (vect[0]=="MARKER_TAG=")
    	{
    		markerTag = vect[1];
    		readingElems = false;
    		readingNodes = false;
    		readingBdrys = true;
    	}
    }
};

/// Export grid nodal coordinates in unstructured tecplot format together with the solution.
/*
 * Requires a file name, the mesh and the associated state
 */
void ExportSolutionToTecplot(const char *filename, Mesh &mesh, State &solution, State &residual)
{
	// retrieve the mesh nodes
	const auto &nodes = mesh.GetNodeContainer();
	// retrieve the mesh cells
	const auto &elements = mesh.GetElementContainer();

	// define the number of fields to be output and their names
	std::vector<std::string> FieldNames {"ElementIndex", "Surface",
		"Density", "VelocityX", "VelocityY", "Pressure",
		"Residual.Density", "Residual.MomentumX", "Residual.MomentumY", "Residual.StagnationEnergy"};
	int NFields = FieldNames.size();

	// open file to plot
	std::fstream outfile;
	outfile.open(filename, std::ios::out);

	// line 1
	outfile << "VARIABLES = \"X\", \"Y\" ";
	for (int i = 0; i < NFields; i++)
		outfile << ", \"" << FieldNames[i] << "\" ";
	outfile << "\n";

	// line 2
	outfile << "ZONE t=\"t:0\",  STRANDID=0, SOLUTIONTIME=0";
	outfile << "\n";
	outfile << " Nodes=" << nodes.size() << ", Elements=" << elements.size() << ", ZONETYPE=FEQuadrilateral";
	outfile << "\n";
	outfile << " DATAPACKING=BLOCK";
	outfile << "\n";
	outfile << " VARLOCATION=([";
	for (int i = 2; i < 2 + NFields ; i++)
	{
		outfile << i + 1;
		// only write a comma if not the last element of our list
		if (i != 2 + NFields - 1)
			outfile << ",";
	}
	outfile << "]=CELLCENTERED)";
	outfile << "\n";
	outfile << "DT=(";
	for (int i = 0; i < 2 + NFields ; i++)
	{
		outfile << "SINGLE ";
	}
	outfile << ")";
	outfile << "\n";

	//std::cout.precision(16);

	// write nodal coordinates
	for (int i = 0; i < nodes.size(); i++)
	{
		outfile << std::scientific << nodes[i].GetXCoord() << " ";
		if (i % 10 == 1) // avoid tecplot saying line is too long...
			outfile << "\n";
	}
	outfile << "\n";

	for (int i = 0; i < nodes.size(); i++)
	{
		outfile << std::scientific << nodes[i].GetYCoord() << " ";
		if (i % 10 == 1) // avoid tecplot saying line is too long...
			outfile << "\n";
	}
	outfile << "\n";

	// write the cell centered solution: cell index
	for (int e = 0; e < mesh.GetCellContainer().size(); ++e)
	{
		outfile << std::scientific << mesh.GetCellContainer()[e].GetCenter().GetNodeID() << " ";
		if (e % 10 == 1) // avoid tecplot saying line is too long...
			outfile << "\n";
	}
	outfile << "\n";

	// write the cell centered solution: surface
	for (int e = 0; e < mesh.GetCellContainer().size(); ++e)
	{
		outfile << std::scientific << mesh.GetCellContainer()[e].GetSurface() << " ";
		if (e % 10 == 1) // avoid tecplot saying line is too long...
			outfile << "\n";
	}
	outfile << "\n";

	// write the cell centered solution: density
	int i = 0;
	for (FieldVector &fieldVector : solution)
	{
		outfile << std::scientific << fieldVector.ConservativeToPrimitive()[0] << " ";
		if (i % 10 == 1) // avoid tecplot saying line is too long...
			outfile << "\n";
		++i;
	}
	outfile << "\n";

	// write the cell centered solution: velocity x
	i = 0;
	for (FieldVector &fieldVector : solution)
	{
		outfile << std::scientific << fieldVector.ConservativeToPrimitive()[1] << " ";
		if (i % 10 == 1) // avoid tecplot saying line is too long...
			outfile << "\n";
		++i;
	}
	outfile << "\n";

	// write the cell centered solution: velocity y
	i = 0;
	for (FieldVector &fieldVector : solution)
	{
		outfile << std::scientific << fieldVector.ConservativeToPrimitive()[2] << " ";
		if (i % 10 == 1) // avoid tecplot saying line is too long...
			outfile << "\n";
		++i;
	}
	outfile << "\n";

	// write the cell centered solution: pressure
	i = 0;
	for (FieldVector &fieldVector : solution)
	{
		outfile << std::scientific << fieldVector.ConservativeToPrimitive()[3] << " ";
		if (i % 10 == 1) // avoid tecplot saying line is too long...
			outfile << "\n";
		++i;
	}
	outfile << "\n";

	// write the cell centered residual: density
	i = 0;
	for (FieldVector &residualVector : residual)
	{
		outfile << std::scientific << residualVector.GetDensity() << " ";
		if (i % 10 == 1) // avoid tecplot saying line is too long...
			outfile << "\n";
		++i;
	}
	outfile << "\n";

	// write the cell centered residual: velocity x
	i = 0;
	for (FieldVector &residualVector : residual)
	{
		outfile << std::scientific << residualVector.GetMomentumX() << " ";
		if (i % 10 == 1) // avoid tecplot saying line is too long...
			outfile << "\n";
		++i;
	}
	outfile << "\n";

	// write the cell centered residual: velocity y
	i = 0;
	for (FieldVector &residualVector : residual)
	{
		outfile << std::scientific << residualVector.GetMomentumY() << " ";
		if (i % 10 == 1) // avoid tecplot saying line is too long...
			outfile << "\n";
		++i;
	}
	outfile << "\n";

	// write the cell centered residual: stagnation energy
	i = 0;
	for (FieldVector &residualVector : residual)
	{
		outfile << std::scientific << residualVector.GetStagnationEnergy() << " ";
		if (i % 10 == 1) // avoid tecplot saying line is too long...
			outfile << "\n";
		++i;
	}
	outfile << "\n";


	// write the connectivity information
	for(int e = 0; e < elements.size(); e++)
	{
		// get the element node IDs
		const auto nodeIDs = elements[e].GetNodeIDs();

		// loop over the element nodes
		for (int j = 0; j < nodeIDs.size(); ++j)
			outfile << nodeIDs[j] + 1 << " ";

		outfile << "\n";
	}

	// end
	outfile.close();
}

void PrintLogo()
{
  std::cout<< std::endl;
  std::cout<< "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout<< "+                                                                  +" << std::endl;
  std::cout<< "+                                                                  +" << std::endl;
  std::cout<< "+                       ** M E R I D I A N **                      +" << std::endl;
  std::cout<< "+                                                                  +" << std::endl;
  std::cout<< "+                                                                  +" << std::endl;
  std::cout<< "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout<< std::endl;
}
