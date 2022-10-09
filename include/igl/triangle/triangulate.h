#pragma once
#ifndef IGL_TRIANGLE_TRIANGULATE_H
#define IGL_TRIANGLE_TRIANGULATE_H
 
#include "../igl_inline.h"
#include <string>
#include <Eigen/Core>
#include "triangle.h"

// for debug
#include <vector>


// 封装了triangle.h中接口的三角剖分接口；
namespace igl
{
  namespace triangle
  {
      /*
         Triangulate the interior of a polygon using the triangle library.

     Inputs:
       vers #vers by 2 list of 2D vertex positions
       E #E by 2 list of vertex ids forming unoriented edges of the boundary of the polygon
       H #H by 2 coordinates of points contained inside holes of the polygon
           M #vers list of markers for input vertices
       flags  string of options pass to triangle (see triangle documentation)
     Outputs:
       V2  #V2 by 2  coordinates of the vertives of the generated triangulation
       F2  #F2 by 3  list of indices forming the faces of the generated triangulation
           M2  #V2 list of markers for output vertices

     TODO: expose the option to prevent Steiner points on the boundary
*/
      template <typename DerivedV, typename DerivedE, typename DerivedH, typename DerivedVM, typename DerivedEM,
          typename DerivedV2, typename DerivedF2, typename DerivedVM2, typename DerivedEM2>
          IGL_INLINE void triangulate(
              const Eigen::MatrixBase<DerivedV>& vers,
              const Eigen::MatrixBase<DerivedE>& E,
              const Eigen::MatrixBase<DerivedH>& H,
              const Eigen::MatrixBase<DerivedVM>& VM,
              const Eigen::MatrixBase<DerivedEM>& EM,
              const std::string flags,
              Eigen::PlainObjectBase<DerivedV2>& V2,
              Eigen::PlainObjectBase<DerivedF2>& F2,
              Eigen::PlainObjectBase<DerivedVM2>& VM2,
              Eigen::PlainObjectBase<DerivedEM2>& EM2)
      {
          using namespace std;
          using namespace Eigen;

          assert((VM.size() == 0 || vers.rows() == VM.size()) &&
              "Vertex markers must be empty or same size as vers");
          assert((EM.size() == 0 || E.rows() == EM.size()) &&
              "Segment markers must be empty or same size as E");
          assert(vers.cols() == 2);
          assert(E.size() == 0 || E.cols() == 2);
          assert(H.size() == 0 || H.cols() == 2);

          // Prepare the flags
          string full_flags = flags + "pz" + (EM.size() || VM.size() ? "" : "B");

          // for debug
          full_flags = "pY";

          typedef Map< Matrix<double, Dynamic, Dynamic, RowMajor> > MapXdr;
          typedef Map< Matrix<int, Dynamic, Dynamic, RowMajor> > MapXir;

          // Prepare the input struct
          triangulateio in;
          in.numberofpoints = vers.rows();
          in.pointlist = (double*)calloc(vers.size(), sizeof(double));             // 二维顶点数组的首地址
          {
              MapXdr inpl(in.pointlist, vers.rows(), vers.cols());
              inpl = vers.template cast<double>();
          }

          // for debug:
          std::vector<double> numVec;
          for (unsigned i = 0; i < vers.cols() * vers.rows(); ++i)
              numVec.push_back(in.pointlist[i]);

          in.numberofpointattributes = 0;
          in.pointmarkerlist = (int*)calloc(vers.size(), sizeof(int));
          for (unsigned i = 0; i < vers.rows(); ++i) in.pointmarkerlist[i] = VM.size() ? VM(i) : 1;

          in.trianglelist = NULL;
          in.numberoftriangles = 0;
          in.numberofcorners = 0;
          in.numberoftriangleattributes = 0;
          in.triangleattributelist = NULL;

          in.numberofsegments = E.size() ? E.rows() : 0;            // 环路边界中的边数
          in.segmentlist = (int*)calloc(E.size(), sizeof(int));         // 环路边界中的边数据
          {
              MapXir insl(in.segmentlist, E.rows(), E.cols());
              insl = E.template cast<int>();
          }
          in.segmentmarkerlist = (int*)calloc(E.rows(), sizeof(int));
          for (unsigned i = 0; i < E.rows(); ++i) in.segmentmarkerlist[i] = EM.size() ? EM(i) : 1;

          in.numberofholes = H.size() ? H.rows() : 0;
          in.holelist = (double*)calloc(H.size(), sizeof(double));
          {
              MapXdr inhl(in.holelist, H.rows(), H.cols());
              inhl = H.template cast<double>();
          }
          in.numberofregions = 0;

          // Prepare the output struct
          triangulateio out;
          out.pointlist = NULL;
          out.trianglelist = NULL;
          out.segmentlist = NULL;
          out.segmentmarkerlist = NULL;
          out.pointmarkerlist = NULL;

          // Call triangle
          triangulate(const_cast<char*>(full_flags.c_str()), &in, &out, 0);

          // Return the mesh
          V2 = MapXdr(out.pointlist, out.numberofpoints, 2).cast<typename DerivedV2::Scalar>();
          F2 = MapXir(out.trianglelist, out.numberoftriangles, 3).cast<typename DerivedF2::Scalar>();
          if (VM.size())
          {
              VM2 = MapXir(out.pointmarkerlist, out.numberofpoints, 1).cast<typename DerivedVM2::Scalar>();
          }
          if (EM.size())
          {
              EM2 = MapXir(out.segmentmarkerlist, out.numberofsegments, 1).cast<typename DerivedEM2::Scalar>();
          }

          // Cleanup in
          free(in.pointlist);
          free(in.pointmarkerlist);
          free(in.segmentlist);
          free(in.segmentmarkerlist);
          free(in.holelist);
          // Cleanup out
          free(out.pointlist);
          free(out.trianglelist);
          free(out.segmentlist);
          free(out.segmentmarkerlist);
          free(out.pointmarkerlist);
      }


          /*
             Triangulate the interior of a polygon using the triangle library.
    
             Inputs:
               vers #vers by 2 list of 2D vertex positions
               E #E by 2 list of vertex ids forming unoriented edges of the boundary of the polygon
               H #H by 2 coordinates of points contained inside holes of the polygon
               flags  string of options pass to triangle (see triangle documentation)
             Outputs:
               V2  #V2 by 2  coordinates of the vertives of the generated triangulation
               F2  #F2 by 3  list of indices forming the faces of the generated triangulation
        */
        template <
          typename DerivedV, typename DerivedE,  typename DerivedH, typename DerivedV2,  typename DerivedF2>
        IGL_INLINE void triangulate(
          const Eigen::MatrixBase<DerivedV> & vers,
          const Eigen::MatrixBase<DerivedE> & E,
          const Eigen::MatrixBase<DerivedH> & H,
          const std::string flags,
          Eigen::PlainObjectBase<DerivedV2> & V2,
          Eigen::PlainObjectBase<DerivedF2> & F2)
        {
            Eigen::VectorXi VM, EM, VM2, EM2;
            return triangulate(vers, E, H, VM, EM, flags, V2, F2, VM2, EM2);
        }
  }
}
 
#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template void igl::triangle::triangulate<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::basic_string<char, std::char_traits<char>, std::allocator<char> >, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::triangle::triangulate<Eigen::Matrix<double, -1, -1, 1, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 1, -1, -1>, Eigen::Matrix<double, -1, -1, 1, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 1, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 1, -1, -1> > const&, std::basic_string<char, std::char_traits<char>, std::allocator<char> >, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 1, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::triangle::triangulate<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, std::basic_string<char, std::char_traits<char>, std::allocator<char> >, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
#endif
 

#endif
