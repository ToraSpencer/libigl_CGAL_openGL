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
  Eigen::PlainObjectBase<DerivedF> &trisResult, 
  Index & curTrisCount, 
  std::unordered_map<int64_t, int> & edgeIsctMap)
{
    /*
    const DerivedGV& gridCenters,															栅格数据
    const Eigen::Matrix<Scalar, 8, 1>& cornerSDF,									当前立方体八个顶点的SDF值
    const Eigen::Matrix<Index, 8, 1>& cornerIdx,									当前立方体八个顶点在栅格中的索引；
    const Scalar& isovalue,																		需要提取的等值面的SDF值
    Eigen::PlainObjectBase<DerivedV>& versResult,								输出网格的顶点
    Index& curVersCount,																			当前累计生成的输出网格顶点数
    Eigen::PlainObjectBase<DerivedF>& trisResult,									输出网格的三角片
    Index& curTrisCount,																			当前累计生成的输出网格三角片数
    std::unordered_map<int64_t, int>& edgeIsctMap								边编码-边交点索引的哈希表；
*/


// These consts get stored reasonably
#include "marching_cubes_tables.h"

  // Seems this is also successfully inlined
  const auto ij2vertex = [&edgeIsctMap, &versResult, &curVersCount, &gridCenters](const Index & i,  const Index & j,  const Scalar & t)->Index
  {
    // Seems this is successfully inlined.
    const auto genMCedgeCode = [](int32_t vaIdx, int32_t vbIdx)
    {
      if(vaIdx > vbIdx)
          std::swap(vaIdx, vbIdx);
      std::int64_t edgeCode = 0;
      edgeCode |= vaIdx;
      edgeCode |= static_cast<std::int64_t>(vbIdx) << 32;
      return edgeCode;
    };

    const auto edgeCode = genMCedgeCode(i, j);
    const auto iter = edgeIsctMap.find(edgeCode);
    int v = -1;

    if(iter == edgeIsctMap.end())
    {
      // 生成新的顶点：
      if(curVersCount==versResult.rows())
          versResult.conservativeResize(versResult.rows()*2+1, versResult.cols()); 
      versResult.row(curVersCount) = gridCenters.row(i) + t*(gridCenters.row(j) - gridCenters.row(i));
      v = curVersCount;
      edgeIsctMap[edgeCode] = v;
      curVersCount++;
    }
    else
      v = iter->second;
    
    return v;
  };

  // 1. 计算当前立方体的顶点状态编码，即8个顶点在等值面的内外状态1；
  Eigen::Matrix<Index, 12, 1> isctVerIdxes;           // 立方体边上的交点的绝对索引——是在最终输出网格中的索引；
  int cornerState = 0;											// 立方体顶点状态编码；256种情形；
  for(int c = 0; c < 8; c++)
      if(cornerSDF(c) > isovalue)
          cornerState |= 1<<c; 

    // 2. 确定当前立方体中和等值面相交的边；
    int edgeState = aiCubeEdgeFlags[cornerState];
    
    // If the cube is entirely inside or outside of the surface,  then there will be no intersections
    if(edgeState == 0)
        return; 

    // 3. 确定等值面和当前立方体的边的交点； Find the point of intersection of the surface with each edge. Then find the normal to the surface at those points
    for(int i = 0; i < 12; i++)
    {
#ifndef NDEBUG
      isctVerIdxes[i] = -1;
#endif

      //if there is an intersection on this edge
      if(edgeState & (1<<i))
      {
        // find crossing point assuming linear interpolation along edges
        const Scalar & SDFa = cornerSDF(a2eConnection[i][0]);
        const Scalar & SDFb = cornerSDF(a2eConnection[i][1]);
        Scalar t;
        const Scalar delta = SDFb-SDFa;
        t = (isovalue - SDFa)/delta;

        // record global index into local table
        isctVerIdxes[i] = ij2vertex(cornerIdx(a2eConnection[i][0]), cornerIdx(a2eConnection[i][1]), t);
        assert(isctVerIdxes[i] >= 0);
        assert(isctVerIdxes[i] < curVersCount);
      }
    }

    // 4. 生成当前立方体中的三角片，一个立方体中最多生成5个三角片；
    for(int i = 0; i < 5; i++)
    {
      if(a2fConnectionTable[cornerState][3*i] < 0) 
          break;
      
      if(curTrisCount==trisResult.rows())
          trisResult.conservativeResize(trisResult.rows()*2+1, trisResult.cols()); 

      assert(isctVerIdxes[a2fConnectionTable[cornerState][3*i+0]]>=0);
      assert(isctVerIdxes[a2fConnectionTable[cornerState][3*i+1]]>=0);
      assert(isctVerIdxes[a2fConnectionTable[cornerState][3*i+2]]>=0);

      trisResult.row(curTrisCount) <<
        isctVerIdxes[a2fConnectionTable[cornerState][3*i+0]], 
        isctVerIdxes[a2fConnectionTable[cornerState][3*i+1]], 
        isctVerIdxes[a2fConnectionTable[cornerState][3*i+2]];
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
