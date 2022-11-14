#include "march_cube.h"

// march_cube()——生成marching cubes算法中的单个方块；
/*
     Something bad is happening when I made this a function.
     Maybe something is not inlining? 
     It ends up 1.25× slower than if the code is pasted into the respective functions in igl::marching_cubes

     Even if I make it a lambda with no arguments (all capture by reference [&]) and call it immediately I get a 1.25× slow-down. 
 
     Maybe keeping it out of a function allows the compiler to optimize with the loop? 
     But then I guess that measn this function is not getting inlined? 
     Or that it's not getting optimized after inlining?
*/
template <
  typename DerivedGV, 
  typename Scalar, 
  typename Index, 
  typename DerivedV, 
  typename DerivedF>
IGL_INLINE void igl::march_cube(
  const DerivedGV & gridCenters, 
  const Eigen::Matrix<Scalar, 8, 1> & cornerSDF, 
  const Eigen::Matrix<Index, 8, 1> & cornerIdx, 
  const Scalar & isovalue, 
  Eigen::PlainObjectBase<DerivedV> &versResult, 
  Index & curVersCount, 
  Eigen::PlainObjectBase<DerivedF> &tris, 
  Index & curTrisCount, 
  std::unordered_map<int64_t, int> & E2V)
{

// These consts get stored reasonably
#include "marching_cubes_tables.h"

  // Seems this is also successfully inlined
  const auto ij2vertex = [&E2V, &versResult, &curVersCount, &gridCenters](const Index & i,  const Index & j,  const Scalar & t)->Index
  {
    // Seems this is successfully inlined.
    const auto ij2key = [](int32_t i, int32_t j)
    {
      if(i>j){ std::swap(i, j); }
      std::int64_t ret = 0;
      ret |= i;
      ret |= static_cast<std::int64_t>(j) << 32;
      return ret;
    };

    const auto key = ij2key(i, j);
    const auto it = E2V.find(key);
    int v = -1;

    if(it == E2V.end())
    {
      // 生成新的顶点：
      if(curVersCount==versResult.rows())
          versResult.conservativeResize(versResult.rows()*2+1, versResult.cols()); 
      versResult.row(curVersCount) = gridCenters.row(i) + t*(gridCenters.row(j) - gridCenters.row(i));
      v = curVersCount;
      E2V[key] = v;
      curVersCount++;
    }
    else
      v = it->second;
    
    return v;
  };

    int c_flags = 0;
    for(int c = 0; c < 8; c++)
      if(cornerSDF(c) > isovalue)
        c_flags |= 1<<c; 

    // Find which edges are intersected by the surface
    int e_flags = aiCubeEdgeFlags[c_flags];
    
    // If the cube is entirely inside or outside of the surface,  then there will be no intersections
    if(e_flags == 0)
        return; 

    // Find the point of intersection of the surface with each edge. Then find the normal to the surface at those points
    Eigen::Matrix<Index, 12, 1> edge_vertices;
    for(int e = 0; e < 12; e++)
    {
#ifndef NDEBUG
      edge_vertices[e] = -1;
#endif

      //if there is an intersection on this edge
      if(e_flags & (1<<e))
      {
        // find crossing point assuming linear interpolation along edges
        const Scalar & a = cornerSDF(a2eConnection[e][0]);
        const Scalar & b = cornerSDF(a2eConnection[e][1]);
        Scalar t;
        {
          const Scalar delta = b-a;
          if(delta == 0)
              t = 0.5; 

          t = (isovalue - a)/delta;
        };

        // record global index into local table
        edge_vertices[e] = ij2vertex(cornerIdx(a2eConnection[e][0]), cornerIdx(a2eConnection[e][1]), t);
        assert(edge_vertices[e] >= 0);
        assert(edge_vertices[e] < curVersCount);
      }
    }

    // Insert the triangles that were found.  There can be up to five per cube
    for(int f = 0; f < 5; f++)
    {
      if(a2fConnectionTable[c_flags][3*f] < 0) 
          break;
      
      if(curTrisCount==tris.rows())
          tris.conservativeResize(tris.rows()*2+1, tris.cols()); 

      assert(edge_vertices[a2fConnectionTable[c_flags][3*f+0]]>=0);
      assert(edge_vertices[a2fConnectionTable[c_flags][3*f+1]]>=0);
      assert(edge_vertices[a2fConnectionTable[c_flags][3*f+2]]>=0);

      tris.row(curTrisCount) <<
        edge_vertices[a2fConnectionTable[c_flags][3*f+0]], 
        edge_vertices[a2fConnectionTable[c_flags][3*f+1]], 
        edge_vertices[a2fConnectionTable[c_flags][3*f+2]];
      curTrisCount++;
    }
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::march_cube<Eigen::MatrixBase<Eigen::Matrix<float,  -1,  -1,  0,  -1,  -1> >,  float,  unsigned int,  Eigen::Matrix<float,  -1,  3,  1,  -1,  3>,  Eigen::Matrix<int,  -1,  3,  1,  -1,  3> >(Eigen::MatrixBase<Eigen::Matrix<float,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::Matrix<float,  8,  1,  0,  8,  1> const&,  Eigen::Matrix<unsigned int,  8,  1,  0,  8,  1> const&,  float const&,  Eigen::PlainObjectBase<Eigen::Matrix<float,  -1,  3,  1,  -1,  3> >&,  unsigned int&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  3,  1,  -1,  3> >&,  unsigned int&,  std::unordered_map<int64_t,  int,  std::hash<int64_t>,  std::equal_to<int64_t>,  std::allocator<std::pair<int64_t const,  int> > >&);
template void igl::march_cube<Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >,  double,  unsigned int,  Eigen::Matrix<double,  -1,  3,  1,  -1,  3>,  Eigen::Matrix<int,  -1,  3,  1,  -1,  3> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::Matrix<double,  8,  1,  0,  8,  1> const&,  Eigen::Matrix<unsigned int,  8,  1,  0,  8,  1> const&,  double const&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  3,  1,  -1,  3> >&,  unsigned int&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  3,  1,  -1,  3> >&,  unsigned int&,  std::unordered_map<int64_t,  int,  std::hash<int64_t>,  std::equal_to<int64_t>,  std::allocator<std::pair<int64_t const,  int> > >&);
template void igl::march_cube<Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >,  double,  long,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::Matrix<double,  8,  1,  0,  8,  1> const&,  Eigen::Matrix<long,  8,  1,  0,  8,  1> const&,  double const&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >&,  long&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  long&,  std::unordered_map<int64_t,  int,  std::hash<int64_t>,  std::equal_to<int64_t>,  std::allocator<std::pair<int64_t const,  int> > >&);
template void igl::march_cube<Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >,  double,  unsigned int,  Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1>,  Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >(Eigen::MatrixBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> > const&,  Eigen::Matrix<double,  8,  1,  0,  8,  1> const&,  Eigen::Matrix<unsigned int,  8,  1,  0,  8,  1> const&,  double const&,  Eigen::PlainObjectBase<Eigen::Matrix<double,  -1,  -1,  0,  -1,  -1> >&,  unsigned int&,  Eigen::PlainObjectBase<Eigen::Matrix<int,  -1,  -1,  0,  -1,  -1> >&,  unsigned int&,  std::unordered_map<int64_t,  int,  std::hash<int64_t>,  std::equal_to<int64_t>,  std::allocator<std::pair<int64_t const,  int> > >&);
#ifdef WIN32
template void __cdecl igl::march_cube<class Eigen::MatrixBase<class Eigen::Matrix<double, -1, -1, 0, -1, -1> >, double, __int64, class Eigen::Matrix<double, -1, -1, 0, -1, -1>, class Eigen::Matrix<int, -1, -1, 0, -1, -1> >(class Eigen::MatrixBase<class Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, class Eigen::Matrix<double, 8, 1, 0, 8, 1> const &, class Eigen::Matrix<__int64, 8, 1, 0, 8, 1> const &, double const &, class Eigen::PlainObjectBase<class Eigen::Matrix<double, -1, -1, 0, -1, -1> > &, __int64 &, class Eigen::PlainObjectBase<class Eigen::Matrix<int, -1, -1, 0, -1, -1> > &, __int64 &, class std::unordered_map<__int64, int, struct std::hash<__int64>, struct std::equal_to<__int64>, class std::allocator<struct std::pair<__int64 const , int> > > &);
#endif
#endif 
