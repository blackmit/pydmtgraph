#include <cstdlib>
#include <utility> // pair
#include <limits>
#include <vector>
#include <queue>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace p = boost::python;
namespace np = boost::python::numpy;

//------------------------/ TODO /------------------------//
// write documentation
// define own method for returning unstable manifold
// define methods for returning persistence pairs
// define own class for vertex hash
// look into more robust typing for img
// rename package and module in setup.py
// create union find class

//------------------------/ Enums /------------------------//
// enum to label an edge as belonging to a vertex-edge
// persistence pair or an edge-triangle persistent pair
enum EdgePairType
{
  VERTEX_EDGE_PAIR,
  EDGE_TRIANGLE_PAIR,
  UNKNOWN_PAIR_TYPE
};

//------------------------/ Classes /------------------------//

struct Vertex; struct Edge; class DMTGraph;

struct Vertex
{
  int index;
  float value;
  int x; int y;
  // `parent` is used for computing persistence with union find.
  Vertex * parent;
  // `morseParent` is used for finding the path to the
  // critical vertex  in the delta forest.
  Vertex * morseParent;
  // `neighbors` stores the list of neighbors
  // in the delta forest.
  std::vector<Vertex *> neighbors;

  Vertex(int iindex, int ix, int iy, float ivalue)
  {
    index = iindex;
    x = ix; y = iy;
    value = ivalue;
    parent = this;
    morseParent = NULL;
  }

  void addNeighbor(Vertex* n)
  {
    neighbors.push_back(n);
  }
};


namespace std
{
  template<> struct hash<Vertex*>
  {
    // hash function for a vertex pointer
    std::size_t operator()(Vertex* const& v) const
    {
      return v->index;
    }
  };
}

struct Edge
{
  int index;
  double value;
  double persistence;
  EdgePairType pairType;
  // primal vertices
  Vertex * v1, * v2;
  // dual vertices/triangles
  Vertex * dv1, * dv2;

  Edge(int iindex, Vertex * iv1, Vertex * iv2, Vertex * idv1, Vertex * idv2)
  {
    index = iindex;
    v1 = iv1; v2 = iv2;
    dv1 = idv1; dv2 = idv2;
    value = fmax(v1->value, v2->value);
    pairType = UNKNOWN_PAIR_TYPE;
    persistence = std::numeric_limits<double>::infinity();
  }
};

class DMTGraph
{
public:
  bool persistenceComputed;
  std::vector<Edge *> edges;
  std::vector<Vertex *> vertices;
  std::vector<Vertex *> dualVertices;
  std::vector< std::pair<double, double> > VEPairs;
  std::vector< std::pair<double, double> > ETPairs;
  std::unordered_set<Vertex *> unstableManifoldVertices;
  std::vector< std::pair<Vertex *, Vertex *> >unstableManifoldEdges;

  DMTGraph(np::ndarray const &);
  void computePersistence();
  p::tuple computeGraph(double, double);

private:
  void createSimplices(np::ndarray const &);
  void clear();
  void collectTree(double);
  void cancelMorsePairs();
  void collectPathToMin(Vertex *);
  void collectUnstableManifold(double, double);
  p::tuple returnUnstableManifold();
  np::ndarray persistentIntervalsInDimension(int);
};

//------------------------/ Comparators /------------------------//
// Edge comparators for sorting list of edges
bool edgeCompare(Edge * e1, Edge * e2)
{
  return e1->value < e2->value
        || ((e1->value == e2->value) && (e1->index < e2->index));
}

bool edgeCompareDual(Edge * e1, Edge * e2)
{
  return e1->value > e2->value
        || ((e1->value == e2->value) && (e1->index > e2->index));
}

// Vertex comparators for use in union find
inline bool vertexCompare(Vertex * v1, Vertex * v2)
{
  return (v1->value < v2->value)
       || ((v1->value == v2->value) && (v1->index < v2->index));
}

inline bool vertexCompareDual(Vertex * dv1, Vertex * dv2)
{
  return (dv1->value > dv2->value)
       || ((dv1->value == dv2->value) && (dv1->index > dv2->index));
}

inline Vertex * vertexMin(Vertex * v1, Vertex * v2)
{
  return vertexCompare(v1, v2) ? v1 : v2;
}

//------------------------/ Union Find /------------------------//

Vertex * find(Vertex* v)
{
  if (v->parent == v) {
    return v;
  } else {
    Vertex * parent = find(v->parent);
    v->parent = parent;
    return parent;
  }
}

double merge(Vertex * v1, Vertex * v2, bool (* vertexLT)(Vertex *, Vertex *))
{
  Vertex * parent1 = find(v1);
  Vertex * parent2 = find(v2);

  if (parent1 == parent2)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  else if (vertexLT(parent1, parent2))
  {
    parent2->parent = parent1;
    return parent2->value;
  }
  else
  {
    parent1->parent = parent2;
    return parent1->value;
  }
}

//------------------------/ DMTGraph Methods /------------------------//

DMTGraph::DMTGraph(np::ndarray const & img)
{
  //------------------------/ Preconditions /------------------------//
  if(img.get_nd() != 2)
  {
    PyErr_SetString(PyExc_ValueError, "Input array must be 2-dimensional.");
    p::throw_error_already_set();
  }
  if(img.get_dtype() != np::dtype::get_builtin<double>())
  {
    PyErr_SetString(PyExc_TypeError, "Input Array must be of dtype double");
    p::throw_error_already_set();
  }
  //------------------------/ initialize /------------------------//
  persistenceComputed = false;
  createSimplices(img);
}

void DMTGraph::createSimplices(np::ndarray const & img)
{
  /********************************************************************
   * Creates the simplices in the simplicial complex representing the image.
   * Triangles in the plane are represented as dual vertices.
   *
   * The plane is triangulated like so:
   *
   *      O---O---O .  .  .
   *      |  /|  /|
   *      | / | / |
   *      |/  |/  |
   *      O---O---O
   *      |  /|  /|
   *      | / | / |
   *      |/  |/  |
   *      O---O---O
   *      .         .
   *      .           .
   *      .             .
   *
   ********************************************************************/
  double * img_ptr = reinterpret_cast<double*>(img.get_data());
  //------------------------/ Create Vertices /------------------------//
  int nRows = img.shape(0);
  int nCols = img.shape(1);
  int nVertices = nRows*nCols;
  vertices.reserve(nVertices);
  for(int r=0; r<nRows; r++)
  {
    for(int c=0; c<nCols; c++)
    {
      int index = r*nCols + c;
      vertices.push_back(new Vertex(index, r, c, -1*img_ptr[index]));
    }
  }
  //------------------------/ Create Dual Vertices (Triangles) /------------------------//
  int nDualVertices = 2*(nRows-1)*(nCols-1)+1;
  dualVertices.reserve(nDualVertices);
  for(int r=0; r<nRows-1; r++)
  {
    for(int c=0; c<nCols-1; c++)
    {
      int index = r*nCols + c;

      double value = fmax(-1*img_ptr[index], fmax(-1*img_ptr[index+1], -1*img_ptr[index+nCols]));
      dualVertices.push_back(new Vertex(2*index, r, 2*c, value));

      value = fmax(-1*img_ptr[index+1], fmax(-1*img_ptr[index+nCols], -1*img_ptr[index+nCols+1]));
      dualVertices.push_back(new Vertex(2*index+1, r, 2*c+1, value));
    }
  }
  // create a dual vertex which is connected to every edge on the boundary of the image
  dualVertices.push_back(new Vertex(-1, -1, -1, std::numeric_limits<double>::infinity()));
  Vertex * dvInf = dualVertices.back();

  //------------------------/ Create Edges /------------------------//
  int nEdges = 2*nRows*nCols - nRows - nCols + (nRows-1)*(nCols-1);
  edges.reserve(nEdges);
  // create vertical edges
  for(int r=0; r<nRows-1; r++)
  {
    for(int c=0; c<nCols; c++)
    {
      // retrieve primal vertices
      int index = r*nCols + c;
      Vertex * v1 = vertices[index];
      Vertex * v2 = vertices[index+nCols];

      // retrieve dual vertices
      //
      // if the edge is the first or last in the row
      // one of its dual vertices is the infinite vertex
      int dualIndex = r*2*(nCols-1) + c*2;
      Vertex * dv1, * dv2;
      if (c == 0)
      {
        dv1 = dvInf;
        dv2 = dualVertices[dualIndex];
      }
      else if (c == nCols-1)
      {
        dv1 = dualVertices[dualIndex-1];
        dv2 = dvInf;
      } else
      {
        dv1 = dualVertices[dualIndex-1];
        dv2 = dualVertices[dualIndex];
      }

      // create edge
      int edgeIndex = edges.size();
      edges.push_back(new Edge(index, v1, v2, dv1, dv2));
    }
  }
  // create horizontal edges
  for(int r=0; r<nRows; r++)
  {
    for(int c=0; c<nCols-1; c++)
    {
      // retrieve primal vertices
      int index = r*nCols + c;
      Vertex * v1 = vertices[index];
      Vertex * v2 = vertices[index+1];

      // retrieve dual vertices
      int dualIndex = r*2*(nCols-1) + c*2;
      Vertex * dv1, * dv2;
      if (r == 0)
      {
        dv1 = dvInf;
        dv2 = dualVertices[dualIndex];
      }
      else if (r == nRows-1)
      {
        dv1 = dualVertices[dualIndex - 2*(nCols-1) + 1];
        dv2 = dvInf;
      }
      else
      {
        dv1 = dualVertices[dualIndex - 2*(nCols-1) + 1];
        dv2 = dualVertices[dualIndex];
      }

      // create edge
      int edgeIndex = edges.size();
      edges.push_back(new Edge(index, v1, v2, dv1, dv2));
    }
  }

  // create diagonal edges
  for(int r=0; r<nRows-1; r++)
  {
    for(int c=0; c<nCols-1; c++)
    {
        // retrieve primal vertices
        int index = r*nCols + c;
        Vertex * v1 = vertices[index + 1];
        Vertex * v2 = vertices[index + nCols];

        // retrieve dual vertices
        int dualIndex = r*2*(nCols-1) + c*2;
        Vertex * dv1 = dualVertices[dualIndex];
        Vertex * dv2 = dualVertices[dualIndex+1];

        // create edge
        int edgeIndex = edges.size();
        edges.push_back(new Edge(index, v1, v2, dv1, dv2));
    }
  }
} // end of create simplices

void DMTGraph::computePersistence()
{
  /*
   *  Compute the persistent homology of the image.
   *
   *  This method computes persistence using union find on the
   *  primal and dual graphs.
   */
  // if persistence already computed, skip this method
  if(persistenceComputed)
  {
    return;
  } else {
    persistenceComputed = true;
  }
  //------------------------/ Compute 0-dimensional persistence /------------------------//
  std::sort(edges.begin(), edges.end(), edgeCompare);
  for(Edge * e : edges)
  {
    double death = e->value;
    double birth = merge(e->v1, e->v2, vertexCompare);
    // `merge` return a NaN if the edge created a 1-dimensional cycle.
    // Thus, this if-statement handles the case the edge killed a 0-dimensional cycle.
    if(!std::isnan(birth))
    {
      e->persistence = death - birth;
      e->pairType = VERTEX_EDGE_PAIR;
      if(birth != death)
      {
        VEPairs.push_back(std::make_pair(birth, death));
      }
    }
  }
  //------------------------/ Compute 1-dimensional persistence /------------------------//
  std::sort(edges.begin(), edges.end(), edgeCompareDual);
  for(Edge * e : edges)
  {
    double birth = e->value;
    double death;
    if(e->pairType == UNKNOWN_PAIR_TYPE)
    {
      death = merge(e->dv1, e->dv2, vertexCompareDual);
    } else {
      death = std::numeric_limits<double>::quiet_NaN();
    } 
    // `merge` return a NaN if the edge killed a 1-dimensional dual cycle.
    // Thus, by duality, this if-statement handles the case the edge
    // created a 1-dimensional primal cycle.
    if(!std::isnan(death))
    {
      e->persistence = death - birth;
      e->pairType = EDGE_TRIANGLE_PAIR;
      if(birth != death)
      {
        ETPairs.push_back(std::make_pair(birth, death));
      }
    }
  }
}

void DMTGraph::clear()
{
  unstableManifoldEdges.clear();
  unstableManifoldVertices.clear();
  for(Vertex * v : vertices)
  {
    v->neighbors.clear();
    v->morseParent = NULL;
  }
}

void DMTGraph::collectTree(double delta1)
{
  /*
   * Build the tree induced by all edges that
   *  - are in vertex-edge pairs
   *  - have persistence less than delta
   * This is the algorithm 'PerSimpTree' in [Dey, Wang, Wang 2018]
   */
  for(Edge * e : edges)
  {
    if(e->pairType == VERTEX_EDGE_PAIR && e->persistence < delta1)
    {
      e->v1->addNeighbor(e->v2);
      e->v2->addNeighbor(e->v1);
    }
  }
}

void DMTGraph::cancelMorsePairs()
{
  /*
   * Compute the Morse vector field by performing Morse cancelation
   * on all vertex-edge persistence pairs with persistence < delta
   *
   * This uses the simplified algorithm presented in [Dey, Wang, Wang 2018]
   */
  for(Vertex * v : vertices)
  {
    if(!v->morseParent)
    {
      //------------/ find min in v's connected component /------------//
      std::unordered_set<Vertex *> explored;
      std::queue<Vertex *> queue;
      queue.push(v);
      Vertex * min = v;
      while(!queue.empty())
      {
        Vertex * curr = queue.front();
        queue.pop();
        explored.insert(curr);
        min = vertexMin(min, curr);
        for(Vertex * n : curr->neighbors)
        {
          if(!explored.count(n))
          {
            queue.push(n);
          }
        }
      }
      //------------/ set morseParent of each node in connected component  /------------//
      min->morseParent = min;
      queue.push(min);
      // perform a bfs from the root and set parents to the previous node
      // here, we test if a node has been explored by seeing if its parent has been set
      while(!queue.empty())
      {
        Vertex * curr = queue.front();
        queue.pop();
        for(Vertex * n : curr->neighbors)
        {
          if(!n->morseParent)
          {
            n->morseParent = curr;
            queue.push(n);
          }
        }
      }
    }
  }
}

void DMTGraph::collectUnstableManifold(double delta1, double delta2)
{
  for (Edge * e : edges)
  {
    if(e->persistence > delta1 && e->value < -1*delta2)
    {
      collectPathToMin(e->v1);
      collectPathToMin(e->v2);
      unstableManifoldEdges.push_back(std::make_pair(e->v1, e->v2));
    }
  }
}

void DMTGraph::collectPathToMin(Vertex * v)
{
  /*
   * A helper method for collectUnstableManifold
   *
   * Add the path between a vertex 'v' and its critical point
   * to the 1-unstable manifold of the image
   */
  Vertex * curr = v;
  while(!unstableManifoldVertices.count(curr) && curr->morseParent != curr)
  {
    unstableManifoldVertices.insert(curr);
    unstableManifoldEdges.push_back(std::make_pair(curr, curr->morseParent));
    curr = curr->morseParent;
  }
}

p::tuple DMTGraph::returnUnstableManifold()
{
  /*
    Return the vertices and edges of the 1-unstable manifold
    formatted as a tuple of numpy arrays

    Retval:
      tuple:
        - First element of tuple is numpy array of vertex positions of
          1-unstable manifold vertices using image coordinates.
        - Second element of tuple is numpy array of the endpoints
          of each edge in the 1-unstable manifold. Each endpoint is represented
          by its index in the array of vertex positions.
  */
  std::vector< std::pair<int,int> > vertexPositions;
  std::vector< std::pair<int,int> > edgeIndices;
  // indexMap stores the index in the return array of each vertex
  std::unordered_map<Vertex *, int> indexMap;
  // unstableManifoldEdges stores pointers to vertices
  // we now add the coordinates of these pointers to the return array
  // TODO: change order of coorindates
  for(auto e : unstableManifoldEdges)
  {
    // If this is the first time encoutering one of e's endpoints,
    // add it to the list of vertices and store its index.
    if(!indexMap.count(e.first))
    {
      indexMap[e.first] = vertexPositions.size();
      vertexPositions.push_back(std::make_pair(e.first->x, e.first->y));
    }
    if(!indexMap.count(e.second))
    {
      indexMap[e.second] = vertexPositions.size();
      vertexPositions.push_back(std::make_pair(e.second->x, e.second->y));
    }
    edgeIndices.push_back(std::make_pair(indexMap[e.first], indexMap[e.second]));
  }

  np::dtype dtype = np::dtype::get_builtin<int>();
  p::tuple stride = p::make_tuple(2*sizeof(int), sizeof(int));   // stride = (stride to next row, stride to next column)

  p::tuple shapeVertices = p::make_tuple(vertexPositions.size(), 2);
  np::ndarray npVertices = np::from_data((int *) vertexPositions.data(), dtype, shapeVertices, stride, p::object());

  p::tuple shapeEdges = p::make_tuple(edgeIndices.size(), 2);
  np::ndarray npEdges = np::from_data((int *) edgeIndices.data(), dtype, shapeEdges, stride, p::object());

  return p::make_tuple(npVertices.copy(), npEdges.copy());
}

p::tuple DMTGraph::computeGraph(double delta1, double delta2=0)
{
  clear();
  computePersistence();
  collectTree(delta1);
  cancelMorsePairs();
  collectUnstableManifold(delta1, delta2);
  return returnUnstableManifold();
}

p::tuple computeDMTGraph(np::ndarray const & img, double delta1, double delta2=0)
{
  /*
    Returns the DMTGraph for `img`
    with parameters `delta1` and `delta2`.
    Use this method if you want a DMTGraph
    without having to create a DMTGraph class.
  */
  DMTGraph G(img);
  return G.computeGraph(delta1, delta2);
}

np::ndarray DMTGraph::persistentIntervalsInDimension(int dim)
{
  // TODO: This isn't returning the right intervals.
  // I think the problem is that is that we are using reverse persistence.
  //------------------------/ Preconditions /------------------------//
  if(dim != 0 && dim != 1)
  {
    PyErr_SetString(PyExc_ValueError, "Dimension must be 0 or 1");
    p::throw_error_already_set();
  }
  //------------------------/ Compute homology /------------------------//
  computePersistence();
  //------------------------/ Prepare output /------------------------//
  np::dtype dtype = np::dtype::get_builtin<double>();
  p::tuple stride = p::make_tuple(2*sizeof(double), sizeof(double));
  if(dim == 0) {
    // As we calculated the persistence using the reverse filtration,
    // the 0-dimensional classes are actually the edge-triangle pairs.
    p::tuple shape = p::make_tuple(ETPairs.size(), 2);
    np::ndarray retArray = np::from_data((double *) ETPairs.data(), dtype, shape, stride, p::object());
    return retArray.copy();
  } else {
    // As we calculated the persistence using the reverse filtration,
    // the 1-dimensional classes are actually the vertex-edge pairs.
    p::tuple shape = p::make_tuple(VEPairs.size(), 2);
    np::ndarray retArray = np::from_data((double *) VEPairs.data(), dtype, shape, stride, p::object());
    return retArray.copy();
  }
}

//------------------------/ C++/Python Compatibility Boilerplater /------------------------//
BOOST_PYTHON_MODULE(dmtgraph)
{
  // initialize the boost numpy namespace
  Py_Initialize();
  np::initialize();
  // define the functions so we can call them in Python
  p::def("computeDMTGraph", computeDMTGraph);
  // define the DMT graph class
  p::class_<DMTGraph>("DMTGraph", p::init<np::ndarray>())
    .def("computeGraph", &DMTGraph::computeGraph);
  ;
}
