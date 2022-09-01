#ifndef IGL_DIJKSTRA
#define IGL_DIJKSTRA
#include "igl_inline.h"

#include <Eigen/Core>
#include <vector>
#include <set>

namespace igl {

    /*
       Dijkstra's algorithm for vertex-weighted shortest paths, with multiple targets.
       Adapted from http://rosettacode.org/wiki/Dijkstra%27s_algorithm .
  
       Inputs:
         source           index of source vertex
         targets          target vector set
         adjList               #vers list of lists of incident vertices (adjacency list), e.g.  as returned by igl::adjacency_list
         weights          #vers list of scalar vertex weights
  
       Output:
         min_distance     #vers by 1 list of the minimum distances from source to all vertices
         previous           list of the previous visited vertices (for each vertex) - used for backtracking
  
    */
  template <typename IndexType, typename DerivedD, typename DerivedP>
  IGL_INLINE int dijkstra(
    const IndexType &source,
    const std::set<IndexType> &targets,
    const std::vector<std::vector<IndexType> >& adjList,
    const std::vector<double>& weights,
    Eigen::PlainObjectBase<DerivedD> &min_distance,
    Eigen::PlainObjectBase<DerivedP> &previous);


  // Dijkstra's algorithm for shortest paths, with multiple targets.
  // Adapted from http://rosettacode.org/wiki/Dijkstra%27s_algorithm .
  //
  // Inputs:
  //   source           index of source vertex
  //   targets          target vector set
  //   adjList               #vers list of lists of incident vertices (adjacency list), e.g.
  //                    as returned by igl::adjacency_list
  //
  // Output:
  //   min_distance     #vers by 1 list of the minimum distances from source to all vertices
  //   previous         #vers by 1 list of the previous visited vertices (for each vertex) - used for backtracking
  //
  template <typename IndexType, typename DerivedD, typename DerivedP>
  IGL_INLINE int dijkstra(
    const IndexType &source,
    const std::set<IndexType> &targets,
    const std::vector<std::vector<IndexType> >& adjList,
    Eigen::PlainObjectBase<DerivedD> &min_distance,
    Eigen::PlainObjectBase<DerivedP> &previous);


  // Backtracking after Dijkstra's algorithm, to find shortest path.
  //
  // Inputs:
  //   vertex           vertex to which we want the shortest path (from same source as above)
  //   previous         #vers by 1 list of the previous visited vertices (for each vertex) - result of Dijkstra's algorithm
  //
  // Output:
  //   path             #P by 1 list of vertex indices in the shortest path from vertex to source
  //
  template <typename IndexType, typename DerivedP>
  IGL_INLINE void dijkstra(
    const IndexType &vertex,
    const Eigen::MatrixBase<DerivedP> &previous,
    std::vector<IndexType> &path);


  
  /*
       Dijkstra's algorithm for shortest paths on a mesh, with multiple targets, using edge length
  
       Inputs:
         vers                #vers by 3 list of vertex positions
         adjList               #vers list of lists of incident vertices (adjacency list), e.g.
                          as returned by igl::adjacency_list, will be generated if empty.
         source           index of source vertex
         targets          target vector set
  
       Output:
         min_distance     #vers by 1 list of the minimum distances from source to all vertices
         previous         #vers by 1 list of the previous visited vertices (for each vertex) - used for backtracking
  */
  template <typename IndexType, typename DerivedV,
  typename DerivedD, typename DerivedP>
  IGL_INLINE int dijkstra(
    const Eigen::MatrixBase<DerivedV> &vers,
    const std::vector<std::vector<IndexType> >& adjList,
    const IndexType &source,
    const std::set<IndexType> &targets,
    Eigen::PlainObjectBase<DerivedD> &min_distance,
    Eigen::PlainObjectBase<DerivedP> &previous);

}

#ifndef IGL_STATIC_LIBRARY
#include "dijkstra.cpp"
#endif


#endif /* defined(IGL_DIJKSTRA) */
