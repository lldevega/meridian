/*
* Mesh representation classes.
*/

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cassert>
#include <memory>

/// Node class.
/*
* A node represents a location in the 2D (x, y) space
*/
class Node
{
public:
  /// Constructor.
  /*
  * Requires the location (x, y) coordinates and the node identifier
  */
  Node(double x, double y, int ID): _x(x), _y(y), _ID(ID) {}

  /// Constructor from another node.
  Node(const Node &other)
  {
     this->_x = other.GetXCoord();
     this->_y = other.GetYCoord();
     this->_ID = other.GetNodeID();
  }

  /// Default constructor
  Node() {}

  /// x-coordinate.
  double GetXCoord() const
  {
      return this->_x;
  }

  /// y-coordinate.
  double GetYCoord() const
  {
      return this->_y;
  }

  /// ID.
  int GetNodeID() const
  {
      return this->_ID;
  }

  /// Operator == returns true if the coordinates x and y are the same (within some tolerance)
  bool operator==(const Node &other)
  {
      double tol = 1.0e-12;
      double dx = this->_x - other.GetXCoord();
      double dy = this->_y - other.GetYCoord();
      double dist = std::sqrt(dx * dx + dy * dy);
      return dist < tol;
  }


protected:
  /// Members.
  double _x;
  double _y;
  int _ID;
};

/// Cell class.
/*
* A cell is a quadrilateral in the 2D (x, y) space
* It is defined by 4 nodes, has 4 faces and a cell center
* which is a node itself.
*/
class Cell
{
public:
  /// Constructor.
  /*
  * Requires 4 nodes
  */
  Cell(const Node &node0, const Node &node1, const Node &node2, const Node &node3, int elemID):
      _nodes({node0, node1, node2, node3}), _center(ComputeCenter(elemID)), _surface(ComputeSurface(elemID))
  { }

  /// Computes the centroid of 4 nodes and stores it as a node.
  const Node ComputeCenter(int elemID)
  {
      double x = 0;
      double y = 0;
      for (int i = 0; i < _nodes.size(); ++i)
      {
          x += 0.25 * _nodes[i].GetXCoord();
          y += 0.25 * _nodes[i].GetYCoord();
      }

      Node center(x, y, elemID);

      return center;
  }

  /// Compute cell surface.
  double ComputeSurface(int elemID)
  {
    double v1x = _nodes[1].GetXCoord() - _nodes[0].GetXCoord();
    double v1y = _nodes[1].GetYCoord() - _nodes[0].GetYCoord();

    double v2x = _nodes[2].GetXCoord() - _nodes[0].GetXCoord();
    double v2y = _nodes[2].GetYCoord() - _nodes[0].GetYCoord();

    double v3x = _nodes[3].GetXCoord() - _nodes[0].GetXCoord();
    double v3y = _nodes[3].GetYCoord() - _nodes[0].GetYCoord();

    double S1 = 0.5 * (v1x * v2y - v2x * v1y);
    double S2 = 0.5 * (v2x * v3y - v3x * v2y);

    double surf = S1 + S2;

    assert (S1 >= 0.0 && S2 >= 0.0 && "Element has negative surface!"); /// TODO: give ID

    return surf;
  }

  /// Get center.
  const Node &GetCenter() const
  {
      return this->_center;
  }

  /// Get surface.
  double GetSurface() const
  {
      return this->_surface;
  }

protected:
    /// Members.
    const std::vector<Node> _nodes;
    const Node _center;
    double _surface;
};

/// Face class.
/*
* A face is a segment linking two nodes
* It is defined by 2 nodes, has at least 1 attached element
*/
class Face
{
public:
  /// Constructor.
  /*
  * Requires 2 nodes
  */
  Face(const Node &node0, const Node &node1, int faceID, const std::string tag): _nodes({node0, node1}),
      _center(ComputeCenter(faceID)), _length(ComputeLength(faceID)), _normal(ComputeNormal(faceID)),
      _tag(tag) { }

  /// Operator == return true if the nodes defining a face are the same
  bool operator ==(const Face &other)
  {
      return this->_nodes[0] == other.GetNodes()[0] && this->_nodes[1] == other.GetNodes()[1] ||
          this->_nodes[0] == other.GetNodes()[1] && this->_nodes[1] == other.GetNodes()[0];
  }

  /// Compute the face center.
  const Node ComputeCenter(int faceID)
  {
      double x;
      double y;
      for (int i = 0; i < this->_nodes.size(); ++i)
      {
          x += 0.5 * this->_nodes[i].GetXCoord();
          y += 0.5 * this->_nodes[i].GetYCoord();
      }

      Node center(x, y, faceID);

      return center;
  }

  /// Compute the face length.
  double ComputeLength(int faceID)
  {
      double dx = this->_nodes[1].GetXCoord() - this->_nodes[0].GetXCoord();
      double dy = this->_nodes[1].GetYCoord() - this->_nodes[0].GetYCoord();
      double length = sqrt(dx * dx + dy * dy);

      return length;
  }

  /// Compute the face normal.
  std::vector<double> ComputeNormal(int faceID)
  {
      double dx = this->_nodes[1].GetXCoord() - this->_nodes[0].GetXCoord();
      double dy = this->_nodes[1].GetYCoord() - this->_nodes[0].GetYCoord();
      double nx = dy / sqrt(dx * dx + dy * dy);
      double ny = -dx / sqrt(dx * dx + dy * dy);

      std::vector<double> normal{nx, ny};

      return normal;
  }

  /// Get center.
  const Node &GetCenter() const
  {
      return this->_center;
  }

  /// Get length.
  double GetLength() const
  {
      return this->_length;
  }

  /// Get normal.
  const std::vector<double> &GetNormal() const
  {
      return this->_normal;
  }

  /// Get tag.
  const std::string &GetTag() const
  {
      return this->_tag;
  }

  /// Get nodes.
  std::vector<Node> GetNodes() const
  {
      return this->_nodes;
  }

  /// Members
protected:
    std::vector<Node> _nodes;
    const Node _center;
    double _length;
    const std::vector<double> _normal;
    const std::string _tag;
};

/// Element class.
/*
* An element is a collection of 4 node IDs that make up a cell
*/
class Element
{
public:
    /// Constructor.
    Element(std::vector<int> data): _elemID(data[0]),
        _nodeIDs({data[1], data[2], data[3], data[4]}) {}

    /// Get element ID.
    int GetElementID()
    {
        return this->_elemID;
    }

    /// Get node ID list.
    std::vector<int> GetNodeIDs() const
    {
        return this->_nodeIDs;
    }

    /// Reorder nodes to ensure positive cell surface.
    void ReorderNodes(const std::vector<int> &newOrder)
    {
    	std::vector<int> oldNodeIDs = this->_nodeIDs;
    	this->_nodeIDs[0] = oldNodeIDs[newOrder[0]];
    	this->_nodeIDs[1] = oldNodeIDs[newOrder[1]];
    	this->_nodeIDs[2] = oldNodeIDs[newOrder[2]];
    	this->_nodeIDs[3] = oldNodeIDs[newOrder[3]];
    }

protected:
    /// Members
    int _elemID;
    std::vector<int> _nodeIDs;
};

/// Element and face stencil class.
/*
*
*/
class Stencil
{
public:

    /// Constructor.
    Stencil(int nEntries) : _stencil(InitStencil(nEntries))
    { }

    /// Default constructor.
    Stencil(): _stencil() { }

    /// Init stencil
    std::vector<std::vector<int>> InitStencil(const int nEntries)
	{
    	std::vector<std::vector<int>> tmpStcl;
    	tmpStcl.resize(nEntries);
        return tmpStcl;
	}

    /// Fills the `idx` element or face stencil.
    /*
    * For elements, the n-th entry contains the list of element
    *     indices the n-th element is attached to.
    * For faces, the the n-th entry contains the list of element
    *     indices the n-th face is attached to.
    */
    void FillAt(int idx, std::vector<int> otherIdx)
    {
        // resize the stencil
        const int stencilSize = otherIdx.size();
        this->_stencil[idx].resize(stencilSize);

        // the rest of the entries are the direct neighbors
        for (int i = 0; i < otherIdx.size(); ++i)
        {
            this->_stencil[idx][i] = otherIdx[i];
        }
    }

     /// Gets the stencil (non-const mode)
    std::vector<std::vector<int>> &GetStencil()
    {
        return this->_stencil;
    }

protected:
    /// Members.
    std::vector<std::vector<int>> _stencil;
};

/// Edge class.
/*
* An edge is a collection of 2 node IDs that make up a face
*/
class Edge
{
public:
    ///Constructor
    Edge(std::vector<int> data, std::string &tag): _nodeIDs({data[1], data[2]}),
        _tag(tag) {}

    /// Get element ID
    const std::string &GetTag()
    {
        return this->_tag;
    }

    /// Get node ID list
    const std::vector<int> &GetNodeIDs()
    {
        return this->_nodeIDs;
    }

protected:
    /// Members.
    const std::string _tag;
    const std::vector<int> _nodeIDs;
};

/// Mesh class.
/*
* The mesh is defined from:
*    1. A list of nodes.
*    2. A list of elements, defined from the previous nodes.
*    3. A list of boundary faces, delimiting the mesh extent.
*/
class Mesh
{
public:
    // define types
    using CellContainer = typename std::vector<Cell>;
    using FaceContainer = typename std::vector<Face>;

    /// Constructor
    /*
    * The mesh is created from a list of nodes, elements and boundaries
    */
    Mesh(const std::vector<Node> nodes, std::vector<Element> elements, const std::vector<Edge> boundaries):
      _nodes(nodes), _elements(elements), _boundaries(boundaries), _cellContainer(nullptr), _elementStencil(nullptr),
	  _internalFaceContainer(nullptr), _internalFaceStencil(nullptr)
    { }

    /// Compute mesh.
    /*
    * This method computes:
    *     1. A container of cells
    *     2. A container of faces
    *     3. The stencil of cells
    *     4. The stencil of faces
    */
    void Compute()
    {
        // get the number of elements in the mesh
        int nbElements = this->_elements.size();

        // initialize the internal face index
        int internalFaceIdx = 0;

        // initialize the containers
        CellContainer cellContainer;
        FaceContainer internalFaceContainer;

        // initialize the element stencil
        Stencil elementStencil(nbElements);
        // initialize the face stencil
        Stencil internalFaceStencil(10 * nbElements);

        // loop over all elements
        for (int i = 0; i < nbElements; i++)
        {
            // retrieve the list of node IDs of the current element
            const auto &nodeIDs = this->_elements[i].GetNodeIDs();

            // get the nodes
            const Node &node0 = this->_nodes[nodeIDs[0]];
            const Node &node1 = this->_nodes[nodeIDs[1]];
            const Node &node2 = this->_nodes[nodeIDs[2]];
            const Node &node3 = this->_nodes[nodeIDs[3]];

            // put them in an array
            std::vector<Node> tmpCellNodes{node0, node1, node2, node3};

            // get the order the indices in an anticlock-wise manner to get right oriented cells
            const auto &orderedIndices = OtherNodesAntiClockWise(node0, node1, node2, node3);

            // reorder the nodes in the element object
            this->_elements[i].ReorderNodes(orderedIndices);

            // compute the associated cell
            Cell cell(tmpCellNodes[orderedIndices[0]], tmpCellNodes[orderedIndices[1]],
                tmpCellNodes[orderedIndices[2]], tmpCellNodes[orderedIndices[3]], i);

            // and push it back to the container
            cellContainer.push_back(cell);

            /// The following handles the connectivity

            // initialize the list of potential neighbor elements
            std::vector<int> neighborElements;
            neighborElements.reserve(4);

            // loop for all elements that are different from the current element
            for (int j = 0; j < nbElements; j++)
            {
                if (i != j)
                {
                    // retrieve the list of node IDs of the (potential) neighbor element
                    const auto &neighborNodeIDs = this->_elements[j].GetNodeIDs();

                    // initialize a counter of common nodes between 2 elements
                    int nbCommonNodes = 0;

                    // initialize a vector containing the indices of the common nodes
                    std::vector<int> commonNodeIndices;

                    // The node IDs of elements i and j are compared to check whether 2 of them are the same
                    for (int k = 0; k < neighborNodeIDs.size(); k++)
                    {
                        //std::vector<int>::iterator it;
                        auto it = std::find(nodeIDs.begin(), nodeIDs.end(), neighborNodeIDs[k]);

                        if (it != nodeIDs.end())
                        {
                            nbCommonNodes++;
                            int index = std::distance(nodeIDs.begin(), it);
                            commonNodeIndices.push_back(nodeIDs[index]);
                        }
                    }

                    // test everything is ok
                    assert (nbCommonNodes < 3 && "Two different quad elements can only share up to 2 nodes");

                    // if we have 2 common nodes, the element j is a neighbor of element i
                    if (nbCommonNodes == 2)
                    {
                        // also, an internal face must be created **if it does not exist yet**
                    	std::string identifier = "internal";

                        Face face(this->_nodes[commonNodeIndices[0]], this->_nodes[commonNodeIndices[1]],
                            internalFaceIdx, identifier);

                        bool faceExists = false;
                        for (int l = 0; l < internalFaceContainer.size(); ++l)
                        {
                            faceExists = face == internalFaceContainer[l];
                            if (faceExists)
                                break;
                        }

                        if (!faceExists)
                        {
                            internalFaceContainer.push_back(face);

                            // and we know its stencil
                            internalFaceStencil.FillAt(internalFaceIdx /*the face index*/,
                                std::vector<int>{i, j} /*the left and right elements*/);

                            ++internalFaceIdx;
                        }

                        // element j is one of element i's neighbors
                        neighborElements.push_back(j);
                    }
                }
            }

            // fill in the element stencil
            elementStencil.FillAt(i, neighborElements);

            // define the right size of the face stencil
            int nInternalFaces = internalFaceIdx + 1;
            auto &internalFaceStencilContainer = internalFaceStencil.GetStencil();
            //internalFaceStencilContainer.resize(nInternalFaces);
            //assert(internalFaceStencilContainer.size() == nInternalFaces);
        }

        /// The following handles the boundaries

        // create the boundary face container
        FaceContainer boundaryFaceContainer;
        // initialize the boundary face stencil
        int nbBoundaryFaces = this->_boundaries.size();
        Stencil boundaryFaceStencil(nbBoundaryFaces);

        for (int bdryIdx = 0; bdryIdx < nbBoundaryFaces; ++bdryIdx)
        {
            // retrieve the boundary edge
            auto bdry = this->_boundaries[bdryIdx];

            // get the node indices
            const auto bdryNodeIDs = bdry.GetNodeIDs();

            // get the nodes of the current boundary
            const auto &node0 = this->_nodes[bdryNodeIDs[0]];
            const auto &node1 = this->_nodes[bdryNodeIDs[1]];

            // create a face
            Face bdryFace(node0, node1, bdryIdx, bdry.GetTag());

            // and push it back to the container
            boundaryFaceContainer.push_back(bdryFace);

            // now we do again a loop over all elements to check to which one this boundary is attached to
            for (int i = 0; i < nbElements; ++i)
            {
                const auto &elementNodeIDs = this->_elements[i].GetNodeIDs();
                int nbCommonNodes = 0;

                // The node IDs of elements i and j are compared to check whether 2 of them are the same
                for (int m = 0; m < elementNodeIDs.size(); m++)
                {
                    //std::vector<int>::iterator it;
                    auto it = std::find(bdryNodeIDs.begin(), bdryNodeIDs.end(), elementNodeIDs[m]);

                    if (it != bdryNodeIDs.end())
                        nbCommonNodes++;
                }

                // if the 2 nodes of the boundary also make up the element, then the boundary is attached to this element
                if (nbCommonNodes == 2)
                {
                    boundaryFaceStencil.FillAt(bdryIdx, std::vector<int>{i});
                    break; // no need to check other elements since a boundary face can only be attached to one element
                }
            }
        }

        /// Store all the compute information
        _cellContainer = std::make_shared<CellContainer> (cellContainer);
        _internalFaceContainer = std::make_shared<FaceContainer> (internalFaceContainer);
        _boundaryFaceContainer = std::make_shared<FaceContainer> (boundaryFaceContainer);
        _elementStencil = std::make_shared<Stencil> (elementStencil);
        _internalFaceStencil = std::make_shared<Stencil>(internalFaceStencil);
        _boundaryFaceStencil = std::make_shared<Stencil> (boundaryFaceStencil);
    }

    /// Get node container.
    const std::vector<Node> &GetNodeContainer() const
    {
        return this->_nodes;
    }

    /// Get element container.
    const std::vector<Element> &GetElementContainer() const
    {
        return this->_elements;
    }

    /// Get element container.
    const std::vector<Edge> &GetBoundaryContainer() const
    {
        return this->_boundaries;
    }

    /// Get element container.
    CellContainer &GetCellContainer() const
    {
        return *this->_cellContainer;
    }

    /// Get internal face container.
    FaceContainer &GetInternalFaceContainer()
    {
        return *this->_internalFaceContainer;
    }

    /// Get boundary face container.
    FaceContainer &GetBoundaryFaceContainer()
    {
        return *this->_boundaryFaceContainer;
    }

    /// Get element stencil.
    Stencil &GetElementStencil() const
    {
        return *this->_elementStencil;
    }

    /// Get element container.
    Stencil &GetInternalFaceStencil()
    {
        return *this->_internalFaceStencil;
    }

    /// Get element container.
    Stencil &GetBoundaryFaceStencil()
    {
        return *this->_boundaryFaceStencil;
    }

    void DisplayInfo()
    {
    	// get access to the containers
    	auto &cellCntr = GetCellContainer();
    	auto &nodeCntr = GetNodeContainer();
    	auto &internalFaceCntr = GetInternalFaceContainer();
    	auto &bcFaceCntr = GetBoundaryFaceContainer();
    	auto &elemStencil = GetElementStencil().GetStencil();
    	auto &innerFaceStencil = GetInternalFaceStencil().GetStencil();
    	auto &bcFaceStencil = GetBoundaryFaceStencil().GetStencil();

        std::cout << "+---------------------------------------+" << std::endl;
        std::cout << "+            Mesh global info           +" << std::endl;
        std::cout << "Nb cells: " << cellCntr.size() << std::endl;
        std::cout << "Nb nodes: " << nodeCntr.size() << std::endl;
        std::cout << "Nb internal faces: " << internalFaceCntr.size() << std::endl;
        std::cout << "Nb boundary faces: " << bcFaceCntr.size() << std::endl;
        std::cout << "+---------------------------------------+" << std::endl;

        // display cell info
        for (int i = 0; i < GetCellContainer().size(); i++)
        {
            std::cout << "Cell ID: " << cellCntr[i].GetCenter().GetNodeID() << std::endl;
            std::cout << "  Surface = " << cellCntr[i].GetSurface() << std::endl;
            std::cout << "  Center = (" << cellCntr[i].GetCenter().GetXCoord() << ", "
              <<  cellCntr[i].GetCenter().GetYCoord() << ")" << std::endl;
            std::cout << "  Stencil = ";
            for (int j = 0; j < elemStencil[i].size(); ++ j)
            {
                std::cout << elemStencil[i][j] << " ";
            }
            std::cout << std::endl;
        }

        // display internal face info
        for (int i = 0; i < internalFaceCntr.size(); ++i)
        {
            std::cout << "Internal face ID: " << internalFaceCntr[i].GetCenter().GetNodeID() << std::endl;
            std::cout << "  Tag = " << internalFaceCntr[i].GetTag()<< std::endl;
            std::cout << "  Length = " << internalFaceCntr[i].GetLength()<< std::endl;
            std::cout << "  Normal = " << internalFaceCntr[i].GetNormal()[0]<<
              " " << internalFaceCntr[i].GetNormal()[1] << std::endl;
            std::cout << "  Stencil = ";
            for (int j = 0; j < innerFaceStencil[i].size(); ++ j)
            {
                std::cout << innerFaceStencil[i][j] << " ";
            }
            std::cout << std::endl;
        }

        // display boundary face info
        for (int i = 0; i < bcFaceCntr.size(); ++i)
        {
            std::cout << "Boundary face ID: " << bcFaceCntr[i].GetCenter().GetNodeID() << std::endl;
            std::cout << "  Tag = " << bcFaceCntr[i].GetTag()<< std::endl;
            std::cout << "  Length = " << bcFaceCntr[i].GetLength()<< std::endl;
            std::cout << "  Normal = " << bcFaceCntr[i].GetNormal()[0]<<
              " " << bcFaceCntr[i].GetNormal()[1] << std::endl;
            std::cout << "  Stencil = ";
            for (int j = 0; j < bcFaceStencil[i].size(); ++ j)
            {
                std::cout << bcFaceStencil[i][j] << " ";
            }
            std::cout << std::endl;
        }
      std::cout << "+---------------------------------------+" << std::endl;
    }

    /// Export grid nodal coordinates in unstructured tecplot format
    /*
     * Requires a file name and the mesh
    */
    void ExportToTecplot(const char *filename)
    {
      // retrieve the mesh nodes
      const auto &nodes = GetNodeContainer();
      // retrieve the mesh cells
      const auto &elements = GetElementContainer();

      // open file to plot
      std::fstream outfile;
      outfile.open(filename, std::ios::out);

      // line 1
      outfile << "VARIABLES = \"X\", \"Y\" ";
      outfile << "\n";

      // line 2
      outfile << "ZONE t=\"t:0\", N=" << nodes.size() << ", E=" << elements.size() << ", F=FEPOINT, ET=";
      outfile << "QUADRILATERAL";
      outfile << "\n";

      //std::cout.precision(16);

      // write nodal coordinates
      for (int i = 0; i < nodes.size(); i++)
      {
        outfile << std::scientific << nodes[i].GetXCoord() << " ";
        outfile << std::scientific << nodes[i].GetYCoord() << " ";
        outfile << "\n";
      }

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

protected:
    /// Members.
    const std::vector<Node> _nodes;
    std::vector<Element> _elements;
    const std::vector<Edge> _boundaries;
    ///
    std::shared_ptr<CellContainer> _cellContainer;
    std::shared_ptr<FaceContainer> _internalFaceContainer;
    std::shared_ptr<FaceContainer> _boundaryFaceContainer;
    ///
    std::shared_ptr<Stencil> _elementStencil;
    std::shared_ptr<Stencil> _internalFaceStencil;
    std::shared_ptr<Stencil> _boundaryFaceStencil;

    /// Sets a list of nodes in anti-clock wise order
    /*
    *
    */
    std::vector<int> OtherNodesAntiClockWise(const Node &node0, const Node &node1,
      const Node &node2, const Node &node3)
    {
        // define the baricenter
        double xc = 0.25 * (node0.GetXCoord() + node1.GetXCoord() + node2.GetXCoord() + node3.GetXCoord());
        double yc = 0.25 * (node0.GetYCoord() + node1.GetYCoord() + node2.GetYCoord() + node3.GetYCoord());
        Node center(xc, yc, 0);

        // compute the angular position
        std::vector<Node> nodeVect{node0, node1, node2, node3};
        std::vector<double> angularPos;

        for (int i = 0; i < nodeVect.size(); i++)
        {
           double dx = nodeVect[i].GetXCoord() - center.GetXCoord();
           double dy = nodeVect[i].GetYCoord() - center.GetYCoord();
           double th = atan2(dy, dx);
           th = fmod((th + M_PI / 2.), (2 * M_PI));
           bool quadrant = std::signbit(th);

           if (quadrant)
               th += 2 * M_PI;

           angularPos.push_back(th * 180 / M_PI);
       }

       std::vector<int> sortedIndices = SortIndices<double>(angularPos);
       return sortedIndices;
    }

    /// Sorts indices.
    template <typename T>
    std::vector<int> SortIndices(const std::vector<T> &v)
    {

        // initialize original index locations
        std::vector<int> idx(v.size());
        std::iota(idx.begin(), idx.end(), 0);

        // sort indexes based on comparing values in v
        // using std::stable_sort instead of std::sort
        // to avoid unnecessary index re-orderings
        // when v contains elements of equal values
        std::stable_sort(idx.begin(), idx.end(),
            [&v](int i1, int i2) {return v[i1] < v[i2];});

        return idx;
    }
};
