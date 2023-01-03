#ifndef IGL_CONNECTED_COMPONENTS_H
#define IGL_CONNECTED_COMPONENTS_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <Eigen/Sparse>


namespace igl
{
    /*
       Determine the connected components of a graph described by the input adjacency matrix (similar to MATLAB's graphconncomp).
  
       Inputs:
          adjSM  #adjSM by #adjSM adjacency matrix (treated as describing an undirected graph)
                                     �����������ڽӾ���

       Outputs:
          connectedLabels  #adjSM list of component indices into [0,#connectedCount-1]
                                        ͬһ��ͨ�����ڵĶ����ǩ������Ϊ0, 1, 2,...
                                        ����Ϊi�Ķ���ı�ǩΪconnectedLabels(i);

          connectedCount  #connectedCount list of sizes of each component
                                        ÿ����ǩ��Ӧ�ĵ���ͨ�����ڰ����Ķ�����
                                        ��ǩΪi�ĵ���ͨ��������Ķ�����ΪconnectedCount(i)
       
       Returns number of connected components
    */
  template < typename Atype, typename DerivedC, typename DerivedK>
  IGL_INLINE int connected_components(
    const Eigen::SparseMatrix<Atype> & adjSM,
    Eigen::PlainObjectBase<DerivedC> & connectedLabels,
    Eigen::PlainObjectBase<DerivedK> & connectedCount);
}

#ifndef IGL_STATIC_LIBRARY
#  include "connected_components.cpp"
#endif

#endif
