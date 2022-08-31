#include "AABB.h"
#include "EPS.h"
#include "barycenter.h"
#include "colon.h"
#include "doublearea.h"
#include "point_simplex_squared_distance.h"
#include "project_to_line_segment.h"
#include "sort.h"
#include "volume.h"
#include "ray_box_intersect.h"
#include "parallel_for.h"
#include "ray_mesh_intersect.h"
#include <iostream>
#include <iomanip>
#include <limits>
#include <list>
#include <queue>
#include <stack>


// init()重载1
template <typename DerivedV, int DIM>
template <typename DerivedEle, typename Derivedbb_mins, typename Derivedbb_maxs, typename Derivedelements>
IGL_INLINE void igl::AABB<DerivedV,DIM>::init(
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedEle> & Ele,
    const Eigen::MatrixBase<Derivedbb_mins> & bb_mins,
    const Eigen::MatrixBase<Derivedbb_maxs> & bb_maxs,
    const Eigen::MatrixBase<Derivedelements> & elements,
    const int i)
{
  using namespace std;
  using namespace Eigen;

  deinit();
  if(bb_mins.size() > 0)
  {
    assert(bb_mins.rows() == bb_maxs.rows() && "Serial tree arrays must match");
    assert(bb_mins.cols() == vers.cols() && "Serial tree array dim must match vers");
    assert(bb_mins.cols() == bb_maxs.cols() && "Serial tree arrays must match");
    assert(bb_mins.rows() == elements.rows() &&
        "Serial tree arrays must match");
    // construct from serialization
    m_box.extend(bb_mins.row(i).transpose());
    m_box.extend(bb_maxs.row(i).transpose());
    m_primitive = elements(i);
    // Not leaf then recurse
    if(m_primitive == -1)
    {
      m_left = new AABB();
      m_left->init( vers,Ele,bb_mins,bb_maxs,elements,2*i+1);
      m_right = new AABB();
      m_right->init( vers,Ele,bb_mins,bb_maxs,elements,2*i+2);
      //m_depth = std::max( m_left->m_depth, m_right->m_depth)+1;
    }
  }else
  {
    VectorXi allI = colon<int>(0,Ele.rows()-1);
    MatrixXDIMS BC;
    if(Ele.cols() == 1)
    {
      // points
      BC = vers;
    }else
    {
      // Simplices
      barycenter(vers,Ele,BC);
    }
    MatrixXi SI(BC.rows(),BC.cols());
    {
      MatrixXDIMS _;
      MatrixXi IS;
      igl::sort(BC,1,true,_,IS);
      // Need SI(i) to tell which place i would be sorted into
      const int dim = IS.cols();
      for(int i = 0;i<IS.rows();i++)
      {
        for(int d = 0;d<dim;d++)
        {
          SI(IS(i,d),d) = i;
        }
      }
    }
    init(vers,Ele,SI,allI);
  }
}


// init()重载1.1——直接读取网格，构造BVH对象；
template <typename DerivedV, int DIM>
template <typename DerivedEle>
void igl::AABB<DerivedV,DIM>::init(
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedEle> & Ele)
{
  using namespace Eigen;
  // deinit will be immediately called...
  return init(vers,Ele,MatrixXDIMS(),MatrixXDIMS(),VectorXi(),0);
}


// init()重载2
  template <typename DerivedV, int DIM>
template <
  typename DerivedEle,
  typename DerivedSI,
  typename DerivedI>
IGL_INLINE void igl::AABB<DerivedV,DIM>::init(
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedEle> & Ele,
    const Eigen::MatrixBase<DerivedSI> & SI,
    const Eigen::MatrixBase<DerivedI> & I)
{
  using namespace Eigen;
  using namespace std;
  deinit();
  if(vers.size() == 0 || Ele.size() == 0 || I.size() == 0)
  {
    return;
  }
  assert(DIM == vers.cols() && "vers.cols() should matched declared dimension");
  //const Scalar inf = numeric_limits<Scalar>::infinity();
  m_box = AlignedBox<Scalar,DIM>();
  // Compute bounding box
  for(int i = 0;i<I.rows();i++)
  {
    for(int closestVer = 0;closestVer<Ele.cols();closestVer++)
    {
      m_box.extend(vers.row(Ele(I(i),closestVer)).transpose());
      m_box.extend(vers.row(Ele(I(i),closestVer)).transpose());
    }
  }
  switch(I.size())
  {
    case 0:
      {
        assert(false);
      }
    case 1:
      {
        m_primitive = I(0);
        break;
      }
    default:
      {
        // Compute longest direction
        int max_d = -1;
        m_box.diagonal().maxCoeff(&max_d);
        // Can't use median on BC directly because many may have same value,
        // but can use median on sorted BC indices
        VectorXi SIdI(I.rows());
        for(int i = 0;i<I.rows();i++)
        {
          SIdI(i) = SI(I(i),max_d);
        }
        // Pass by copy to avoid changing input
        const auto median = [](VectorXi A)->int
        {
          size_t n = (A.size()-1)/2;
          nth_element(A.data(),A.data()+n,A.data()+A.size());
          return A(n);
        };
        const int med = median(SIdI);
        VectorXi LI((I.rows()+1)/2),RI(I.rows()/2);
        assert(LI.rows()+RI.rows() == I.rows());
        // Distribute left and right
        {
          int li = 0;
          int ri = 0;
          for(int i = 0;i<I.rows();i++)
          {
            if(SIdI(i)<=med)
            {
              LI(li++) = I(i);
            }else
            {
              RI(ri++) = I(i);
            }
          }
        }
        //m_depth = 0;
        if(LI.rows()>0)
        {
          m_left = new AABB();
          m_left->init(vers,Ele,SI,LI);
          //m_depth = std::max(m_depth, m_left->m_depth+1);
        }
        if(RI.rows()>0)
        {
          m_right = new AABB();
          m_right->init(vers,Ele,SI,RI);
          //m_depth = std::max(m_depth, m_right->m_depth+1);
        }
      }
  }
}


template <typename DerivedV, int DIM>
IGL_INLINE bool igl::AABB<DerivedV,DIM>::is_leaf() const
{
  return m_primitive != -1;
}


template <typename DerivedV, int DIM>
template <typename DerivedEle, typename Derivedq>
IGL_INLINE std::vector<int> igl::AABB<DerivedV,DIM>::find(
    const Eigen::MatrixBase<DerivedV> & vers,
    const Eigen::MatrixBase<DerivedEle> & Ele,
    const Eigen::MatrixBase<Derivedq> & q,
    const bool first) const
{
  using namespace std;
  using namespace Eigen;
  assert(q.size() == DIM &&
      "Query dimension should match aabb dimension");
  assert(Ele.cols() == vers.cols()+1 &&
      "AABB::find only makes sense for (d+1)-simplices");
  const Scalar epsilon = igl::EPS<Scalar>();
  // Check if outside bounding box
  bool inside = m_box.contains(q.transpose());
  if(!inside)
  {
    return std::vector<int>();
  }
  assert(m_primitive==-1 || (m_left == NULL && m_right == NULL));
  if(is_leaf())
  {
    // Initialize to some value > -epsilon
    Scalar a1=0,a2=0,a3=0,a4=0;
    switch(DIM)
    {
      case 3:
        {
          // Barycentric coordinates
          typedef Eigen::Matrix<Scalar,1,3> RowVector3S;
          const RowVector3S V1 = vers.row(Ele(m_primitive,0));
          const RowVector3S V2 = vers.row(Ele(m_primitive,1));
          const RowVector3S V3 = vers.row(Ele(m_primitive,2));
          const RowVector3S V4 = vers.row(Ele(m_primitive,3));
          a1 = volume_single(V2,V4,V3,(RowVector3S)q);
          a2 = volume_single(V1,V3,V4,(RowVector3S)q);
          a3 = volume_single(V1,V4,V2,(RowVector3S)q);
          a4 = volume_single(V1,V2,V3,(RowVector3S)q);
          break;
        }
      case 2:
        {
          // Barycentric coordinates
          typedef Eigen::Matrix<Scalar,2,1> Vector2S;
          const Vector2S V1 = vers.row(Ele(m_primitive,0));
          const Vector2S V2 = vers.row(Ele(m_primitive,1));
          const Vector2S V3 = vers.row(Ele(m_primitive,2));
          // Hack for now to keep templates simple. If becomes bottleneck
          // consider using std::enable_if_t
          const Vector2S q2 = q.head(2);
          a1 = doublearea_single(V1,V2,q2);
          a2 = doublearea_single(V2,V3,q2);
          a3 = doublearea_single(V3,V1,q2);
          break;
        }
      default:assert(false);
    }
    // Normalization is important for correcting sign
    Scalar sum = a1+a2+a3+a4;
    a1 /= sum;
    a2 /= sum;
    a3 /= sum;
    a4 /= sum;
    if(
        a1>=-epsilon &&
        a2>=-epsilon &&
        a3>=-epsilon &&
        a4>=-epsilon)
    {
      return std::vector<int>(1,m_primitive);
    }else
    {
      return std::vector<int>();
    }
  }
  std::vector<int> left = m_left->find(vers,Ele,q,first);
  if(first && !left.empty())
  {
    return left;
  }
  std::vector<int> right = m_right->find(vers,Ele,q,first);
  if(first)
  {
    return right;
  }
  left.insert(left.end(),right.begin(),right.end());
  return left;
}

template <typename DerivedV, int DIM>
IGL_INLINE int igl::AABB<DerivedV,DIM>::subtree_size() const
{
  // 1 for self
  int n = 1;
  int n_left = 0,n_right = 0;
  if(m_left != NULL)
  {
    n_left = m_left->subtree_size();
  }
  if(m_right != NULL)
  {
    n_right = m_right->subtree_size();
  }
  n += 2*std::max(n_left,n_right);
  return n;
}


template <typename DerivedV, int DIM>
template <typename Derivedbb_mins, typename Derivedbb_maxs, typename Derivedelements>
IGL_INLINE void igl::AABB<DerivedV,DIM>::serialize(
    Eigen::PlainObjectBase<Derivedbb_mins> & bb_mins,
    Eigen::PlainObjectBase<Derivedbb_maxs> & bb_maxs,
    Eigen::PlainObjectBase<Derivedelements> & elements,
    const int i) const
{
  using namespace std;
  using namespace Eigen;
  // Calling for root then resize output
  if(i==0)
  {
    const int m = subtree_size();
    //cout<<"m: "<<m<<endl;
    bb_mins.resize(m,DIM);
    bb_maxs.resize(m,DIM);
    elements.resize(m,1);
  }
  //cout<<i<<" ";
  bb_mins.row(i) = m_box.min();
  bb_maxs.row(i) = m_box.max();
  elements(i) = m_primitive;
  if(m_left != NULL)
  {
    m_left->serialize(bb_mins,bb_maxs,elements,2*i+1);
  }
  if(m_right != NULL)
  {
    m_right->serialize(bb_mins,bb_maxs,elements,2*i+2);
  }
}


// squared_distance()重载1.1
template <typename DerivedV, int DIM>
template <typename DerivedEle>
IGL_INLINE typename igl::AABB<DerivedV,DIM>::Scalar
igl::AABB<DerivedV,DIM>::squared_distance(
  const Eigen::MatrixBase<DerivedV> & vers,
  const Eigen::MatrixBase<DerivedEle> & Ele,
  const RowVectorDIMS & ver0,
  int & i,
  Eigen::PlainObjectBase<RowVectorDIMS> & closestVer) const
{
  return squared_distance(vers,Ele,ver0, std::numeric_limits<Scalar>::infinity(), i, closestVer);
}


// squared_distance()重载1
template <typename DerivedV, int DIM>
template <typename DerivedEle>
IGL_INLINE typename igl::AABB<DerivedV,DIM>::Scalar
igl::AABB<DerivedV,DIM>::squared_distance(
  const Eigen::MatrixBase<DerivedV> & vers,
  const Eigen::MatrixBase<DerivedEle> & Ele,
  const RowVectorDIMS & ver0,
  Scalar low_sqr_d,             // 距离下界，一般取0.0
  Scalar up_sqr_d,              // 距离上界，一般取std::numeric_limits<Scalar>::infinity()
  int & i,                              // 和ver0点距离最近的网格面片索引；
  Eigen::PlainObjectBase<RowVectorDIMS> & closestVer) const
{
  using namespace Eigen;
  using namespace std;
 
  if(low_sqr_d > up_sqr_d)
    return low_sqr_d;

  Scalar sqr_d = up_sqr_d;

  assert((Ele.cols() == 3 || Ele.cols() == 2 || Ele.cols() == 1)
    && "Code has only been tested for simplex sizes 3,2,1");

  assert(m_primitive==-1 || (m_left == NULL && m_right == NULL));

  if(is_leaf())
    leaf_squared_distance(vers, Ele, ver0, low_sqr_d, sqr_d, i, closestVer);
  else
  {
    bool looked_left = false;
    bool looked_right = false;

    // lambda
    const auto & look_left = [&]()
    {
      int i_left;
      RowVectorDIMS c_left = closestVer;
      Scalar sqr_d_left = m_left->squared_distance(vers,Ele,ver0,low_sqr_d,sqr_d,i_left,c_left);
      this->set_min(ver0,sqr_d_left,i_left,c_left,sqr_d, i, closestVer);
      looked_left = true;
    };

    // lambda
    const auto & look_right = [&]()
    {
      int i_right;
      RowVectorDIMS c_right = closestVer;
      Scalar sqr_d_right = m_right->squared_distance(vers,Ele,ver0,low_sqr_d,sqr_d,i_right,c_right);
      this->set_min(ver0,sqr_d_right,i_right,c_right,sqr_d,i,closestVer);
      looked_right = true;
    };

    // must look left or right if in box
    if(m_left->m_box.contains(ver0.transpose()))
      look_left();
    if(m_right->m_box.contains(ver0.transpose()))
      look_right();

    // if haven't looked left and could be less than current min, then look
    Scalar left_up_sqr_d = m_left->m_box.squaredExteriorDistance(ver0.transpose());
    Scalar right_up_sqr_d = m_right->m_box.squaredExteriorDistance(ver0.transpose());

    if(left_up_sqr_d < right_up_sqr_d)
    {
      if(!looked_left && left_up_sqr_d<sqr_d)
        look_left();
      if( !looked_right && right_up_sqr_d<sqr_d)
        look_right();
    }
    else
    {
      if( !looked_right && right_up_sqr_d<sqr_d)
        look_right();
      if(!looked_left && left_up_sqr_d<sqr_d)
        look_left();
    }
  }

  return sqr_d;
}


// squared_distance()重载1.2
template <typename DerivedV, int DIM>
template <typename DerivedEle>
IGL_INLINE typename igl::AABB<DerivedV,DIM>::Scalar
igl::AABB<DerivedV,DIM>::squared_distance(
  const Eigen::MatrixBase<DerivedV> & vers,
  const Eigen::MatrixBase<DerivedEle> & Ele,
  const RowVectorDIMS & ver0,
  Scalar up_sqr_d,
  int & i,
  Eigen::PlainObjectBase<RowVectorDIMS> & closestVer) const
{
  return squared_distance(vers, Ele, ver0, 0.0, up_sqr_d, i, closestVer);
}


template <typename DerivedV, int DIM>
template <
  typename DerivedEle,
  typename DerivedP,
  typename DerivedsqrD,
  typename DerivedI,
  typename DerivedC>
IGL_INLINE void igl::AABB<DerivedV,DIM>::squared_distance(
  const Eigen::MatrixBase<DerivedV> & vers,
  const Eigen::MatrixBase<DerivedEle> & Ele,
  const Eigen::MatrixBase<DerivedP> & vers0,
  Eigen::PlainObjectBase<DerivedsqrD> & minSqrDis,
  Eigen::PlainObjectBase<DerivedI> & I,
  Eigen::PlainObjectBase<DerivedC> & C) const
{
  assert(vers0.cols() == vers.cols() && "cols in vers0 should match dim of cols in vers");
  minSqrDis.resize(vers0.rows(),1);
  I.resize(vers0.rows(),1);
  C.resizeLike(vers0);

  // O( #vers0 * log #Ele ), where log #Ele is really the depth of this AABB hierarchy

  igl::parallel_for(vers0.rows(),[&](int num)
        {
          RowVectorDIMS Pp = vers0.row(num), closestVer;
          int Ip;
          minSqrDis(num) = squared_distance(vers, Ele, Pp, Ip, closestVer);
          I(num) = Ip;
          C.row(num).head(DIM) = closestVer;
        },\
    10000);
}


// squared_distance()重载2
template <typename DerivedV, int DIM>
template <
  typename DerivedEle,
  typename Derivedother_V,
  typename Derivedother_Ele,
  typename DerivedsqrD,
  typename DerivedI,
  typename DerivedC>
IGL_INLINE void igl::AABB<DerivedV,DIM>::squared_distance(
  const Eigen::MatrixBase<DerivedV> & vers,
  const Eigen::MatrixBase<DerivedEle> & Ele,
  const AABB<Derivedother_V,DIM> & other,
  const Eigen::MatrixBase<Derivedother_V> & other_V,
  const Eigen::MatrixBase<Derivedother_Ele> & other_Ele,
  Eigen::PlainObjectBase<DerivedsqrD> & minSqrDis,
  Eigen::PlainObjectBase<DerivedI> & I,
  Eigen::PlainObjectBase<DerivedC> & C) const
{
  assert(other_Ele.cols() == 1 &&
    "Only implemented for other as list of points");
  assert(other_V.cols() == vers.cols() && "other must match this dimension");
  minSqrDis.setConstant(other_Ele.rows(),1,std::numeric_limits<double>::infinity());
  I.resize(other_Ele.rows(),1);
  C.resize(other_Ele.rows(),other_V.cols());

  // All points in other_V currently think they need to check against root of
  // this. The point of using another AABB is to quickly prune chunks of
  // other_V so that most points just check some subtree of this.

  // This holds a conservative estimate of max(sqr_D) where sqr_D is the current best minimum squared distance for all points in this subtree
  double up_sqr_d = std::numeric_limits<double>::infinity();
  squared_distance_helper(vers, Ele, &other, other_V, other_Ele, 0, up_sqr_d, minSqrDis, I, C);
}



template <typename DerivedV, int DIM>
template <
  typename DerivedEle,
  typename Derivedother_V,
  typename Derivedother_Ele,
  typename DerivedsqrD,
  typename DerivedI,
  typename DerivedC>
IGL_INLINE typename igl::AABB<DerivedV,DIM>::Scalar
  igl::AABB<DerivedV,DIM>::squared_distance_helper(
  const Eigen::MatrixBase<DerivedV> & vers,
  const Eigen::MatrixBase<DerivedEle> & Ele,
  const AABB<Derivedother_V,DIM> * other,
  const Eigen::MatrixBase<Derivedother_V> & other_V,
  const Eigen::MatrixBase<Derivedother_Ele> & other_Ele,
  const Scalar /*up_sqr_d*/,
  Eigen::PlainObjectBase<DerivedsqrD> & minSqrDis,
  Eigen::PlainObjectBase<DerivedI> & I,
  Eigen::PlainObjectBase<DerivedC> & C) const
{
  using namespace std;
  using namespace Eigen;

  // This implementation is a bit disappointing. There's no major speed up. Any
  // performance gains seem to come from accidental cache coherency and
  // diminish for larger "other" (the opposite of what was intended).

  // Base case
  if(other->is_leaf() && this->is_leaf())
  {
    Scalar sqr_d = minSqrDis(other->m_primitive);
    int i = I(other->m_primitive);
    RowVectorDIMS closestVer = C.row(      other->m_primitive);
    RowVectorDIMS ver0 = other_V.row(other->m_primitive);
    leaf_squared_distance(vers,Ele,ver0,sqr_d,i,closestVer);
    minSqrDis( other->m_primitive) = sqr_d;
    I(other->m_primitive) = i;
    C.row(other->m_primitive) = closestVer;
    return sqr_d;
  }

  if(other->is_leaf())
  {
    Scalar sqr_d = minSqrDis(other->m_primitive);
    int i = I(other->m_primitive);
    RowVectorDIMS closestVer = C.row(      other->m_primitive);
    RowVectorDIMS ver0 = other_V.row(other->m_primitive);
    sqr_d = squared_distance(vers,Ele,ver0,sqr_d,i,closestVer);
    minSqrDis( other->m_primitive) = sqr_d;
    I(other->m_primitive) = i;
    C.row(other->m_primitive) = closestVer;
    return sqr_d;
  }

  if(this->is_leaf())
  {
    if(true)
    {
      this->squared_distance_helper(vers,Ele,other->m_left,other_V,other_Ele,0,minSqrDis,I,C);
      this->squared_distance_helper(vers,Ele,other->m_right,other_V,other_Ele,0,minSqrDis,I,C);
    }
    return 0;
  }

  // FORCE DOWN TO OTHER LEAF EVAL
  if(true)
  {
    if(true)
    {
      this->squared_distance_helper(vers,Ele,other->m_left,other_V,other_Ele,0,minSqrDis,I,C);
      this->squared_distance_helper(vers,Ele,other->m_right,other_V,other_Ele,0,minSqrDis,I,C);
    }
    else // this direction never seems to be faster
    {
      this->m_left->squared_distance_helper(vers,Ele,other,other_V,other_Ele,0,minSqrDis,I,C);
      this->m_right->squared_distance_helper(vers,Ele,other,other_V,other_Ele,0,minSqrDis,I,C);
    }
  }
  return 0;

}


template <typename DerivedV, int DIM>
template <typename DerivedEle>
IGL_INLINE void igl::AABB<DerivedV,DIM>::leaf_squared_distance(
  const Eigen::MatrixBase<DerivedV> & vers,
  const Eigen::MatrixBase<DerivedEle> & Ele,
  const RowVectorDIMS & ver0,
  const Scalar low_sqr_d,
  Scalar & sqr_d,
  int & i,
  Eigen::PlainObjectBase<RowVectorDIMS> & closestVer) const
{
  using namespace Eigen;
  using namespace std;
  if(low_sqr_d > sqr_d)
  {
    sqr_d = low_sqr_d;
    return;
  }
  RowVectorDIMS c_candidate;
  Scalar sqr_d_candidate;
  igl::point_simplex_squared_distance<DIM>(
    ver0,vers,Ele,m_primitive,sqr_d_candidate,c_candidate);
  set_min(ver0,sqr_d_candidate,m_primitive,c_candidate,sqr_d,i,closestVer);
}

template <typename DerivedV, int DIM>
template <typename DerivedEle>
IGL_INLINE void igl::AABB<DerivedV,DIM>::leaf_squared_distance(
  const Eigen::MatrixBase<DerivedV> & vers,
  const Eigen::MatrixBase<DerivedEle> & Ele,
  const RowVectorDIMS & ver0,
  Scalar & sqr_d,
  int & i,
  Eigen::PlainObjectBase<RowVectorDIMS> & closestVer) const
{
  return leaf_squared_distance(vers,Ele,ver0,0,sqr_d,i,closestVer);
}


template <typename DerivedV, int DIM>
IGL_INLINE void igl::AABB<DerivedV,DIM>::set_min(
  const RowVectorDIMS &
#ifndef NDEBUG
  ver0
#endif
  ,
  const Scalar sqr_d_candidate,
  const int i_candidate,
  const RowVectorDIMS & c_candidate,
  Scalar & sqr_d,
  int & i,
  Eigen::PlainObjectBase<RowVectorDIMS> & closestVer) const
{
#ifndef NDEBUG
  //std::cout<<matlab_format(c_candidate,"c_candidate")<<std::endl;
  //// This doesn't quite make sense to check with bounds
  // const Scalar pc_norm = (ver0-c_candidate).squaredNorm();
  // const Scalar diff = fabs(sqr_d_candidate - pc_norm);
  // assert(diff<=1e-10 && "distance should match norm of difference");
#endif
  if(sqr_d_candidate < sqr_d)
  {
    i = i_candidate;
    closestVer = c_candidate;
    sqr_d = sqr_d_candidate;
  }
}


template <typename DerivedV, int DIM>
template <typename DerivedEle>
IGL_INLINE bool
igl::AABB<DerivedV,DIM>::intersect_ray(
  const Eigen::MatrixBase<DerivedV> & vers,
  const Eigen::MatrixBase<DerivedEle> & Ele,
  const RowVectorDIMS & origin,
  const RowVectorDIMS & dir,
  std::vector<igl::Hit> & hits) const
{
  hits.clear();
  const Scalar t0 = 0;
  const Scalar t1 = std::numeric_limits<Scalar>::infinity();
  {
    Scalar _1,_2;
    if(!ray_box_intersect(origin,dir,m_box,t0,t1,_1,_2))
    {
      return false;
    }
  }
  if(this->is_leaf())
  {
    // Actually process elements
    assert((Ele.size() == 0 || Ele.cols() == 3) && "Elements should be triangles");
    // Cheesecake way of hitting element
    bool ret = ray_mesh_intersect(origin,dir,vers,Ele.row(m_primitive),hits);
    // Since we only gave ray_mesh_intersect a single face, it will have set
    // any hits to id=0. Set these to this primitive's id
    for(auto & hit : hits)
    {
      hit.id = m_primitive;
    }
    return ret;
  }
  std::vector<igl::Hit> left_hits;
  std::vector<igl::Hit> right_hits;
  const bool left_ret = m_left->intersect_ray(vers,Ele,origin,dir,left_hits);
  const bool right_ret = m_right->intersect_ray(vers,Ele,origin,dir,right_hits);
  hits.insert(hits.end(),left_hits.begin(),left_hits.end());
  hits.insert(hits.end(),right_hits.begin(),right_hits.end());
  return left_ret || right_ret;
}

template <typename DerivedV, int DIM>
template <typename DerivedEle>
IGL_INLINE bool
igl::AABB<DerivedV,DIM>::intersect_ray(
  const Eigen::MatrixBase<DerivedV> & vers,
  const Eigen::MatrixBase<DerivedEle> & Ele,
  const RowVectorDIMS & origin,
  const RowVectorDIMS & dir,
  igl::Hit & hit) const
{
#if false
  // BFS
  std::queue<const AABB *> Q;
  // Or DFS
  //std::stack<const AABB *> Q;
  Q.push(this);
  bool any_hit = false;
  hit.t = std::numeric_limits<Scalar>::infinity();
  while(!Q.empty())
  {
    const AABB * tree = Q.front();
    //const AABB * tree = Q.top();
    Q.pop();
    {
      Scalar _1,_2;
      if(!ray_box_intersect(
        origin,dir,tree->m_box,Scalar(0),Scalar(hit.t),_1,_2))
      {
        continue;
      }
    }
    if(tree->is_leaf())
    {
      // Actually process elements
      assert((Ele.size() == 0 || Ele.cols() == 3) && "Elements should be triangles");
      igl::Hit leaf_hit;
      if(
        ray_mesh_intersect(origin,dir,vers,Ele.row(tree->m_primitive),leaf_hit)&&
        leaf_hit.t < hit.t)
      {
        // correct the id
        leaf_hit.id = tree->m_primitive;
        hit = leaf_hit;
      }
      continue;
    }
    // Add children to queue
    Q.push(tree->m_left);
    Q.push(tree->m_right);
  }
  return any_hit;
#else
  // DFS
  return intersect_ray(
    vers,Ele,origin,dir,std::numeric_limits<Scalar>::infinity(),hit);
#endif
}

template <typename DerivedV, int DIM>
template <typename DerivedEle>
IGL_INLINE bool
igl::AABB<DerivedV,DIM>::intersect_ray(
  const Eigen::MatrixBase<DerivedV> & vers,
  const Eigen::MatrixBase<DerivedEle> & Ele,
  const RowVectorDIMS & origin,
  const RowVectorDIMS & dir,
  const Scalar _min_t,
  igl::Hit & hit) const
{
  //// Naive, slow
  //std::vector<igl::Hit> hits;
  //intersect_ray(vers,Ele,origin,dir,hits);
  //if(hits.size() > 0)
  //{
  //  hit = hits.front();
  //  return true;
  //}else
  //{
  //  return false;
  //}
  Scalar min_t = _min_t;
  const Scalar t0 = 0;
  {
    Scalar _1,_2;
    if(!ray_box_intersect(origin,dir,m_box,t0,min_t,_1,_2))
    {
      return false;
    }
  }
  if(this->is_leaf())
  {
    // Actually process elements
    assert((Ele.size() == 0 || Ele.cols() == 3) && "Elements should be triangles");
    // Cheesecake way of hitting element
    bool ret = ray_mesh_intersect(origin,dir,vers,Ele.row(m_primitive),hit);
    hit.id = m_primitive;
    return ret;
  }

  // Doesn't seem like smartly choosing left before/after right makes a
  // differnce
  igl::Hit left_hit;
  igl::Hit right_hit;
  bool left_ret = m_left->intersect_ray(vers,Ele,origin,dir,min_t,left_hit);
  if(left_ret && left_hit.t<min_t)
  {
    // It's scary that this line doesn't seem to matter....
    min_t = left_hit.t;
    hit = left_hit;
    left_ret = true;
  }else
  {
    left_ret = false;
  }
  bool right_ret = m_right->intersect_ray(vers,Ele,origin,dir,min_t,right_hit);
  if(right_ret && right_hit.t<min_t)
  {
    min_t = right_hit.t;
    hit = right_hit;
    right_ret = true;
  }else
  {
    right_ret = false;
  }
  return left_ret || right_ret;
}

// This is a bullshit template because AABB annoyingly needs templates for bad
// combinations of 3D vers with DIM=2 AABB
//
// _Define_ as a no-op rather than monkeying around with the proper code above
//
// Meanwhile, GCC seems to have a bug. Let's see if GCC likes using explicit
// namespace block instead. https://stackoverflow.com/a/25594681/148668
namespace igl
{
  template<> template<> IGL_INLINE float AABB<Eigen::Matrix<float, -1, 3, 1, -1, 3>, 2>::squared_distance( Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&, Eigen::Matrix<float, 1, 2, 1, 1, 2> const&, int&, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> >&) const { assert(false);return -1;};
  template<> template<> IGL_INLINE float igl::AABB<Eigen::Matrix<float, -1, 3, 1, -1, 3>, 2>::squared_distance( Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&, Eigen::Matrix<float, 1, 2, 1, 1, 2> const&, float, float, int&, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> >&) const { assert(false);return -1;};
  template<> template<> IGL_INLINE void igl::AABB<Eigen::Matrix<float, -1, 3, 1, -1, 3>, 2>::init (Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&) { assert(false);};
  template<> template<> IGL_INLINE double AABB<Eigen::Matrix<double, -1, 3, 1, -1, 3>, 2>::squared_distance( Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&, Eigen::Matrix<double, 1, 2, 1, 1, 2> const&, int&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&) const { assert(false);return -1;};
  template<> template<> IGL_INLINE double igl::AABB<Eigen::Matrix<double, -1, 3, 1, -1, 3>, 2>::squared_distance( Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&, Eigen::Matrix<double, 1, 2, 1, 1, 2> const&, double, double, int&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&) const { assert(false);return -1;};
  template<> template<> IGL_INLINE void igl::AABB<Eigen::Matrix<double, -1, 3, 1, -1, 3>, 2>::init (Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&) { assert(false);};
  template<> template<> IGL_INLINE void igl::AABB<Eigen::Matrix<float, -1, 3, 0, -1, 3>, 2>::init(Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&) {assert(false);};
  template<> template<> IGL_INLINE float igl::AABB<Eigen::Matrix<float, -1, 3, 0, -1, 3>, 2>::squared_distance(Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::Matrix<float, 1, 2, 1, 1, 2> const&, float, float, int&, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 2, 1, 1, 2> >&) const { assert(false);return -1;};
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template void igl::AABB<Eigen::Matrix<double, -1, 3, 1, -1, 3>, 3>::squared_distance<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 3, 1, -1, 3>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 3, 1, -1, 3> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> >&) const;
// generated by autoexplicit.sh
template void igl::AABB<Eigen::Matrix<double, -1, 3, 1, -1, 3>, 3>::init<Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
template bool igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 3>::intersect_ray<Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::Matrix<double, 1, 3, 1, 1, 3> const&, Eigen::Matrix<double, 1, 3, 1, 1, 3> const&, igl::Hit&) const;
template double igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 2>::squared_distance<Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::Matrix<double, 1, 2, 1, 1, 2> const&, int&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&) const;
template double igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 3>::squared_distance<Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::Matrix<double, 1, 3, 1, 1, 3> const&, double, int&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&) const;
template double igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 3>::squared_distance<Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::Matrix<double, 1, 3, 1, 1, 3> const&, int&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&) const;
template double igl::AABB<Eigen::Matrix<double, -1, 3, 1, -1, 3>, 3>::squared_distance<Eigen::Matrix<int, -1, 3, 1, -1, 3> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&, Eigen::Matrix<double, 1, 3, 1, 1, 3> const&, double, double, int&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&) const;
template float igl::AABB<Eigen::Matrix<float, -1, 3, 0, -1, 3>, 3>::squared_distance<Eigen::Matrix<int, -1, 3, 0, -1, 3> >(Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::Matrix<float, 1, 3, 1, 1, 3> const&, float, float, int&, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 3, 1, 1, 3> >&) const;
template float igl::AABB<Eigen::Matrix<float, -1, 3, 1, -1, 3>, 3>::squared_distance<Eigen::Matrix<int, -1, 3, 1, -1, 3> >(Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&, Eigen::Matrix<float, 1, 3, 1, 1, 3> const&, int&, Eigen::PlainObjectBase<Eigen::Matrix<float, 1, 3, 1, 1, 3> >&) const;
template std::vector<int, std::allocator<int> > igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 2>::find<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, 1, -1, 1, 1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, bool) const;
template std::vector<int, std::allocator<int> > igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 3>::find<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Block<Eigen::Matrix<double, -1, -1, 0, -1, -1> const, 1, -1, false> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Block<Eigen::Matrix<double, -1, -1, 0, -1, -1> const, 1, -1, false> > const&, bool) const;
template std::vector<int, std::allocator<int> > igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 3>::find<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, 1, -1, 1, 1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 1, -1, 1, 1, -1> > const&, bool) const;
template void igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 2>::init<Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
template void igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 2>::init<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, int);
template void igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 2>::init<Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&);
template void igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 2>::serialize<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, int) const;
template void igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 2>::squared_distance<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&) const;
template void igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 2>::squared_distance<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&) const;
template void igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 2>::squared_distance<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&) const;
template void igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 2>::squared_distance<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 3, 0, -1, 3> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&) const;
template void igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 2>::squared_distance<Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&) const;
template void igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 3>::init<Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
template void igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 3>::init<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, int);
template void igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 3>::init<Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&);
template void igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 3>::serialize<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, int) const;
template void igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 3>::squared_distance<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&) const;
template void igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 3>::squared_distance<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&) const;
template void igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 3>::squared_distance<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&) const;
template void igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 3>::squared_distance<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&) const;
template void igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 3>::squared_distance<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, 3, 0, -1, 3> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> >&) const;
template void igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 3>::squared_distance<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, 2, 3, 0, 2, 3>, Eigen::Matrix<double, 2, 1, 0, 2, 1>, Eigen::Matrix<int, 2, 1, 0, 2, 1>, Eigen::Matrix<double, 2, 3, 0, 2, 3> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, 2, 3, 0, 2, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, 2, 1, 0, 2, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, 2, 1, 0, 2, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, 2, 3, 0, 2, 3> >&) const;
template void igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 3>::squared_distance<Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&) const;
template void igl::AABB<Eigen::Matrix<double, -1, 3, 1, -1, 3>, 3>::init<Eigen::Matrix<int, -1, 3, 1, -1, 3> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&);
template void igl::AABB<Eigen::Matrix<float, -1, 3, 0, -1, 3>, 3>::init<Eigen::Matrix<int, -1, 3, 0, -1, 3> >(Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 0, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&);
template void igl::AABB<Eigen::Matrix<float, -1, 3, 1, -1, 3>, 3>::init<Eigen::Matrix<int, -1, 3, 1, -1, 3> >(Eigen::MatrixBase<Eigen::Matrix<float, -1, 3, 1, -1, 3> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 3, 1, -1, 3> > const&);
template double igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 3>::squared_distance<Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::Matrix<double, 1, 3, 1, 1, 3> const&, int&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 3, 1, 1, 3> >&) const;
template double igl::AABB<Eigen::Matrix<double, -1, -1, 0, -1, -1>, 2>::squared_distance<Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, Eigen::Matrix<double, 1, 2, 1, 1, 2> const&, int&, Eigen::PlainObjectBase<Eigen::Matrix<double, 1, 2, 1, 1, 2> >&) const;
#ifdef WIN32
template void igl::AABB<class Eigen::Matrix<double,-1,-1,0,-1,-1>,2>::squared_distance<class Eigen::Matrix<int,-1,-1,0,-1,-1>,class Eigen::Matrix<double,-1,-1,0,-1,-1>,class Eigen::Matrix<double,-1,1,0,-1,1>,class Eigen::Matrix<__int64,-1,1,0,-1,1>,class Eigen::Matrix<double,-1,3,0,-1,3> >(class Eigen::MatrixBase<class Eigen::Matrix<double,-1,-1,0,-1,-1> > const &,class Eigen::MatrixBase<class Eigen::Matrix<int,-1,-1,0,-1,-1> > const &,class Eigen::MatrixBase<class Eigen::Matrix<double,-1,-1,0,-1,-1> > const &,class Eigen::PlainObjectBase<class Eigen::Matrix<double,-1,1,0,-1,1> > &,class Eigen::PlainObjectBase<class Eigen::Matrix<__int64,-1,1,0,-1,1> > &,class Eigen::PlainObjectBase<class Eigen::Matrix<double,-1,3,0,-1,3> > &)const;
template void igl::AABB<class Eigen::Matrix<double,-1,-1,0,-1,-1>,3>::squared_distance<class Eigen::Matrix<int,-1,-1,0,-1,-1>,class Eigen::Matrix<double,-1,-1,0,-1,-1>,class Eigen::Matrix<double,-1,1,0,-1,1>,class Eigen::Matrix<__int64,-1,1,0,-1,1>,class Eigen::Matrix<double,-1,3,0,-1,3> >(class Eigen::MatrixBase<class Eigen::Matrix<double,-1,-1,0,-1,-1> > const &,class Eigen::MatrixBase<class Eigen::Matrix<int,-1,-1,0,-1,-1> > const &,class Eigen::MatrixBase<class Eigen::Matrix<double,-1,-1,0,-1,-1> > const &,class Eigen::PlainObjectBase<class Eigen::Matrix<double,-1,1,0,-1,1> > &,class Eigen::PlainObjectBase<class Eigen::Matrix<__int64,-1,1,0,-1,1> > &,class Eigen::PlainObjectBase<class Eigen::Matrix<double,-1,3,0,-1,3> > &)const;
#endif
#endif
