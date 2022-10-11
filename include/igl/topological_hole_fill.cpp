#include "topological_hole_fill.h"


  template <
  typename DerivedF,
  typename Derivedb,
  typename VectorIndex,
  typename DerivedF_filled>
IGL_INLINE void igl::topological_hole_fill(
  const Eigen::MatrixBase<DerivedF> & tris,
  const Eigen::MatrixBase<Derivedb> & b,
  const std::vector<VectorIndex> & holes,
  Eigen::PlainObjectBase<DerivedF_filled> &trisOut)
{
  int filledCount = 0;
  int holesCount = holes.size();
  int trisCount = tris.rows();
  const int versCount = tris.maxCoeff()+1;

  for (int i = 0; i < holesCount; i++)
    filledCount += holes[i].size();

  trisOut.resize(filledCount + trisCount, 3);
  trisOut.topRows(trisCount) = tris;

  int newVerIdx = versCount;
  int newTriIdx = trisCount;

  // 补每一个洞的循环
  for (int i = 0; i < holesCount; i++, newVerIdx++)
  {
    int holeVersCount = holes[i].size();                // 当前洞所含顶点数；
    int index = 0;

    // 洞中生成新的三角片；
    trisOut.row(newTriIdx++) << holes[i][index], holes[i][holeVersCount - 1], newVerIdx;
    while (index != holeVersCount - 1)
    {
          trisOut.row(newTriIdx++) << holes[i][(index + 1)], holes[i][(index)], newVerIdx;
          index++;
    }
  }
  assert(newTriIdx == trisOut.rows());
  assert(newVerIdx == versCount + holesCount);

}


#ifdef IGL_STATIC_LIBRARY
template void igl::topological_hole_fill<\
            Eigen::Matrix<int, -1, -1, 0, -1, -1>, \
            Eigen::Matrix<int, -1, 1, 0, -1, 1>, \
            std::vector<int, std::allocator<int> > >
            (
                Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, \
                Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, \
                std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > const&, \
                Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&
            );
 
#endif
