#include <chrono>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "spline.hpp"

int main()
{
	using std::cout;
	using std::endl;

	const int Npieces = 8;      // Number of polynomials in the spline
	constexpr unsigned Deg = 2; // Polynomial degree

	const int Npoints = 10240;
	const float pi = atan(1.0f)*4.0f;
    
	//Discretize the sine function for [-pi/2..pi/2]
	Eigen::Matrix2Xf xy(2, Npoints);
	xy.row(0) = Eigen::VectorXf::LinSpaced(Npoints, -pi*0.5f, pi*0.5f);
	xy.row(1) = xy.row(0).array().sin();

	auto start = std::chrono::steady_clock::now();
	//initialize the spline
	Spline<Deg> spline = Spline<Deg>::fromSorted(xy, Npieces);
	auto diff = std::chrono::steady_clock::now() - start;
	double dtMs = std::chrono::duration<double, std::milli>(diff).count();
	//Takes 0.13 ms on i7-4930K CPU
	cout << "Spline<>::fromSorted() took: " << dtMs << " ms" << endl;

	start = std::chrono::steady_clock::now();
	Eigen::VectorXf approx = xy.row(0).unaryExpr([&](float f){return spline.value_unsafe(f);});
	float sum = approx.sum();
	diff = std::chrono::steady_clock::now() - start;
	dtMs = std::chrono::duration<double, std::milli>(diff).count();
	//Takes 0.041 ms on i7-4930K CPU
	cout << "spline.value_unsafe() took: " << dtMs << " ms, sum = " << sum
	     << endl;

	float maxAbsDiff = (approx.transpose() - xy.row(1)).cwiseAbs().maxCoeff();
	cout << "maximum absolute diference = " << std::setprecision(8)
	     << maxAbsDiff << endl;
	return 0;
}
