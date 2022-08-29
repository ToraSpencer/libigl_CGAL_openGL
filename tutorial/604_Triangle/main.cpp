#include <igl/opengl/glfw/Viewer.h>
#include <igl/triangle/triangulate.h>
#include <igl/writeOBJ.h>
#include <igl/readOBJ.h>

// Input polygon
Eigen::MatrixXd vers;
Eigen::MatrixXi bdry;
Eigen::MatrixXd hole;

// Triangulated interior
Eigen::MatrixXd V2;
Eigen::MatrixXi F2;


// original main();
void test0(int argc, char* argv[])
{
    using namespace Eigen;
    using namespace std;

    // Create the boundary of a square
    vers.resize(8, 2);
    bdry.resize(8, 2);
    hole.resize(1, 2);

    // create two squares, one with edge length of 4, one with edge length of 2

    // both centered at origin
    vers << -1, -1, 1, -1, 1, 1, -1, 1,
        -2, -2, 2, -2, 2, 2, -2, 2;

    // add the bdry of the squares
    bdry << 0, 1, 1, 2, 2, 3, 3, 0,
        4, 5, 5, 6, 6, 7, 7, 4;

    // specify a point that is inside a closed shape , where we do not want triangulation to happen
    hole << 0, 0;

    // for debug
    Eigen::MatrixXd vers3D = Eigen::MatrixXd::Zero(8, 3);
    vers3D.leftCols(2) = vers;
    igl::writeOBJ("E:/vers3D.obj", vers3D, Eigen::MatrixXi{});

    igl::triangle::triangulate(vers, bdry, hole, "a0.005q", V2, F2);

    // Plot the generated mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V2, F2);
    viewer.launch();
}


int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;

  Eigen::MatrixXd vers3D;
  Eigen::MatrixXi tris;
  igl::readOBJ("E:/circleVers.obj", vers3D, tris);
  unsigned versCount = vers3D.rows();

  // Create the boundary of a square
  vers = vers3D.leftCols(2);

  // 2. ���ɶ�ά����ı�Ե��Ϣ��
 bdry.resize(versCount, 2);				// ���ɱպϻ�·��һϵ�еı�
  for (int i = 0; i < versCount; ++i)
  {
      bdry(i, 0) = i + 1;								// ������ע��triangle���ж���������1��ʼ��
      bdry(i, 1) = (i + 1) % versCount + 1;
  }
 
  igl::triangle::triangulate(vers, bdry, hole,"a0.005q",V2, F2);

  // Plot the generated mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V2,F2);
  viewer.launch();
}
