#include <iostream>
#include <string>
#include "Eigen/Dense"


// 导入libigl的各个静态库：
#pragma comment(lib,"igl.lib")		
#pragma comment(lib,"glad.lib")	
#pragma comment(lib,"glfw3.lib")		
#pragma comment(lib,"igl_opengl.lib")		
#pragma comment(lib,"igl_opengl_glfw.lib")		


#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/opengl/glfw/Viewer.h>


// 注：本项目禁用了vcpkg




////////////////////////////////////////////////////////////////////////////////////////////// DEBUG 接口
namespace MY_DEBUG
{
	static std::string g_debugPath = "E:/";


	static void debugDisp()			// 递归终止
	{						//		递归终止设为无参或者一个参数的情形都可以。
		std::cout << std::endl;
		return;
	}

	template <typename T, typename... Types>
	static void debugDisp(const T& firstArg, const Types&... args)
	{
		std::cout << firstArg << " ";
		debugDisp(args...);
	}


	template <typename T, int M, int N>
	static void dispData(const Eigen::Matrix<T, M, N>& m)
	{
		auto dataPtr = m.data();
		unsigned elemsCount = m.size();

		for (unsigned i = 0; i < elemsCount; ++i)
			std::cout << dataPtr[i] << ", ";

		std::cout << std::endl;
	}


	template <typename Derived>
	static void dispData(const Eigen::PlainObjectBase<Derived>& m)
	{
		int m0 = m.RowsAtCompileTime;
		int n0 = m.ColsAtCompileTime;

		auto dataPtr = m.data();
		unsigned elemsCount = m.size();

		for (unsigned i = 0; i < elemsCount; ++i)
			std::cout << dataPtr[i] << ", ";

		std::cout << std::endl;
	}


	template <typename Derived>
	static void dispElem(const Eigen::MatrixBase<Derived>& m)
	{
		const Derived& mm = m.derived();
		std::cout << mm(1, 1) << std::endl;
	}


	template<typename DerivedV>
	static void debugWriteVers(const char* name, const Eigen::PlainObjectBase<DerivedV>& vers)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteVerticesMat(path, vers);
	}

	template<typename DerivedV>
	static void debugWriteVers2D(const char* name, const Eigen::PlainObjectBase<DerivedV>& vers)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteVerticesMat2D(path, vers);
	}


	template<typename T>
	static void debugWriteMesh(const char* name, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteMeshMat(path, vers, tris);
	}


	template<typename DerivedV>
	static void debugWriteEdges(const char* name, const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& vers)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
		objWriteEdgesMat(path, edges, vers);
	}
}
using namespace MY_DEBUG;



int main(int argc, char** argv) 
{ 
	Eigen::MatrixXd vers;
	Eigen::MatrixXi tris;


	igl::readOBJ("E:/材料/tooth.obj", vers, tris);
	debugDisp("versCount == ", vers.rows());
	debugDisp("trisCount == ", tris.rows());

	igl::opengl::glfw::Viewer viewer;				// libigl中的基于glfw的显示窗口；

	// 1. 窗口中装载数据；
	viewer.data().set_mesh(vers, tris);

	// 2. 设定三轴旋转；默认是两轴旋转，set_rotation_type()方法可以指定其他旋转风格
	viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);

	// 3. show_lines指定是否画出网格线；
	viewer.data().show_lines = 0;

	// 4. 进入渲染循环：
	viewer.launch();


	debugDisp("main finished.");


	return 0;
}