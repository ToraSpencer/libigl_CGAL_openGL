// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "writeWRL.h"
#include <iostream>
#include <fstream>
template <typename DerivedV, typename DerivedF>
IGL_INLINE bool igl::writeWRL(
  const std::string & str,
  const Eigen::MatrixBase<DerivedV> & vers,
  const Eigen::MatrixBase<DerivedF> & tris)
{
  using namespace std;
  using namespace Eigen;
  assert(vers.cols() == 3 && "vers should have 3 columns");
  assert(tris.cols() == 3 && "tris should have 3 columns");
  ofstream s(str);
  if(!s.is_open())
  {
    cerr<<"IOError: writeWRL() could not open "<<str<<endl;
    return false;
  }
  // Append column of -1 to tris
  Matrix<typename DerivedF::Scalar,Dynamic,4> FF(tris.rows(),4);
  FF.leftCols(3) = tris;
  FF.col(3).setConstant(-1);

  s<<R"(#VRML V2.0 utf8
DEF default Transform {
translation 0 0 0
children [
Shape {
geometry DEF default-FACES IndexedFaceSet {
ccw TRUE
)"<<
    vers.format(
      IOFormat(
        FullPrecision,
        DontAlignCols,
        " ",",\n","","",
        "coord DEF default-COORD Coordinate { point [ \n","]\n}\n"))<<
    FF.format(
      IOFormat(
        FullPrecision,
        DontAlignCols,
        ",","\n","","",
        "coordIndex [ \n"," ]\n"))<<
    "}\n}\n]\n}\n";
  return true;
}

// write mesh and colors-by-vertex to an ascii off file
template <typename DerivedV, typename DerivedF, typename DerivedC>
IGL_INLINE bool igl::writeWRL(
  const std::string & str,
  const Eigen::MatrixBase<DerivedV> & vers,
  const Eigen::MatrixBase<DerivedF> & tris,
  const Eigen::MatrixBase<DerivedC> & C)
{
  using namespace std;
  using namespace Eigen;
  assert(vers.cols() == 3 && "vers should have 3 columns");
  assert(tris.cols() == 3 && "tris should have 3 columns");
  ofstream s(str);
  if(!s.is_open())
  {
    cerr<<"IOError: writeWRL() could not open "<<str<<endl;
    return false;
  }
  // Append column of -1 to tris
  Matrix<typename DerivedF::Scalar,Dynamic,4> FF(tris.rows(),4);
  FF.leftCols(3) = tris;
  FF.col(3).setConstant(-1);


  //Check if RGB values are in the range [0..1] or [0..255]
  double rgbScale = (C.maxCoeff() <= 1.0)?1.0:1.0/255.0;
  Eigen::MatrixXd RGB = rgbScale * C;

  s<<R"(#VRML V2.0 utf8
DEF default Transform {
translation 0 0 0
children [
Shape {
geometry DEF default-FACES IndexedFaceSet {
ccw TRUE
)"<<
    vers.format(
      IOFormat(
        FullPrecision,
        DontAlignCols,
        " ",",\n","","",
        "coord DEF default-COORD Coordinate { point [ \n","]\n}\n"))<<
    FF.format(
      IOFormat(
        FullPrecision,
        DontAlignCols,
        ",","\n","","",
        "coordIndex [ \n"," ]\n"))<<
    RGB.format(
      IOFormat(
        FullPrecision,
        DontAlignCols,
        ",","\n","","",
        "colorPerVertex TRUE\ncolor Color { color [ \n"," ] }\n"))<<
    "}\n}\n]\n}\n";
  return true;
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template bool igl::writeWRL<Eigen::Matrix<double, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 3, 1, -1, 3> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&);
// generated by autoexplicit.sh
template bool igl::writeWRL<Eigen::Matrix<float, -1, 3, 1, -1, 3>, Eigen::Matrix<int, -1, 3, 1, -1, 3> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&);
// generated by autoexplicit.sh
template bool igl::writeWRL<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&);
// generated by autoexplicit.sh
template bool igl::writeWRL<Eigen::Matrix<double, 8, 3, 0, 8, 3>, Eigen::Matrix<int, 12, 3, 0, 12, 3> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 8, 3, 0, 8, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, 12, 3, 0, 12, 3> > const&);
template bool igl::writeWRL<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
template bool igl::writeWRL<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(std::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&);
#endif
