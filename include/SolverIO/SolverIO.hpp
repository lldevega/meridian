/*
* Solver Input/Output functionalities
*/

#include <string.h>

///
/*
 *
 */
class SU2Reader
{
public:
    /// Constructor.
    SU2Reader(std::string name) : _filename(name) {}

public:
    void ReadMesh()
    {
    	std::ifstream file;
    	std::string line;

    	file.open(this->_filename.c_str());
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
    			String2VectorOfString(line, vect);
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
    				String2VectorOfNumerics(line, vect);
    				_meshElements.push_back(Element(vect));
    			}
    			// case reading nodes
    			else if (!readingElems && readingNodes && !readingBdrys && vect[0] != "NPOIN=")
    			{
    				std::vector<double> vect;
    				String2VectorOfNumerics(line, vect);
    				_meshNodes.push_back(Node(vect[0], vect[1], vect[2]));
    			}
    			// case reading boundaries
    			else if (!readingElems && !readingNodes && readingBdrys)
    			{
    				std::vector<int> vect;
    				String2VectorOfNumerics(line, vect);
    				_meshBoundaries.push_back(Edge(vect, markerTag));
    			}
    			else {}
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
    std::string _filename;
    std::vector<Element> _meshElements;
    std::vector<Node> _meshNodes;
    std::vector<Edge> _meshBoundaries;

    /// Method reading and storing numeric entries from a file.
    template <class container>
    void String2VectorOfNumerics(const std::string& str, container& cont, char delimiter = ' ')
    {
    	std::stringstream ss(str);
    	std::string token;
    	while(std::getline(ss, token, delimiter))
    	{
    		cont.push_back(std::stod(token.c_str()));
    	}
    }

    /// Method reading and storing string entries from a file.
    template <class container>
    void String2VectorOfString(const std::string& str, container& cont, char delimiter = ' ')
    {
    	std::stringstream ss(str);
    	std::string token;
    	while(std::getline(ss, token, delimiter))
    	{
    		cont.push_back(token.c_str());
    	}
    }

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

    outfile << nodeIDs[0] + 1 << " " << nodeIDs[1] + 1 << " " << nodeIDs[2] + 1 << " " << nodeIDs[3] + 1;
    outfile << "\n";
  }

  // end
  outfile.close();
}

void PlotTecTimeIntegration(const char *filename, std::vector<std::vector<double>> &convergenceMonitorVector)
{
  // Open file to plot
  std::fstream outfile;
  outfile.open(filename, std::ios::out);

  // Line 1:
  outfile<<"VARIABLES="<< "\t"<<"Iteration"<< "\t"<<"DensityResiduals"<<"\t"<<"MomentumXResiduals"<< "\t"
      <<"MomentumYResiduals"<< "\t"<<"TotalEnergyResiduals"<< "\n";
  // Write lines
  for(int i=0; i<=convergenceMonitorVector[0].size()-1; i++)
  {
    outfile<<convergenceMonitorVector[0][i]<<"\t"<<convergenceMonitorVector[1][i]<<"\t"
        <<convergenceMonitorVector[2][i]<<"\t"<<convergenceMonitorVector[3][i]<<"\t"
        <<convergenceMonitorVector[4][i]<<"\n";
  }
  outfile<<"\n";

  // End
  outfile.close();
}

void PrintLogo()
{
  std::cout<< "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout<< "+                                                                  +" << std::endl;
  std::cout<< "+                                                                  +" << std::endl;
  std::cout<< "+                       ** M E R I D I A N **                      +" << std::endl;
  std::cout<< "+                                                                  +" << std::endl;
  std::cout<< "+                                                                  +" << std::endl;
  std::cout<< "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

}
