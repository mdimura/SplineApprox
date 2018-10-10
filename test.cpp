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

	// Discretize the sine function for [-pi/2..pi/2]
	const int Npoints = 1024;
	Eigen::Matrix2Xf xy(2, Npoints);
	const float pi = 3.14159f;
	for (int i = 0; i < Npoints; ++i) {
		xy(0, i) = -pi * 0.5f + pi * i / (Npoints - 1);
		xy(1, i) = std::sin(xy(0, i));
	}

	auto start = std::chrono::steady_clock::now();
	// initialize the spline
	Spline<Deg> spline = Spline<Deg>::fromSorted(xy, Npieces);
	auto diff = std::chrono::steady_clock::now() - start;
	double dtMs = std::chrono::duration<double, std::milli>(diff).count();
	// Takes 0.02 ms on i7-4930K CPU
	cout << "Spline<>::fromSorted() took: " << dtMs << " ms" << endl;

	start = std::chrono::steady_clock::now();
	float sum = 0.0;
	for (int j = 0; j < 1000; ++j) {
		for (int i = 0; i < Npoints; i += 1) {
			sum += spline.value_unsafe(xy(0, i));
		}
	}
	diff = std::chrono::steady_clock::now() - start;
	dtMs = std::chrono::duration<double, std::milli>(diff).count();
	// Takes 2.8 ms on i7-4930K CPU
	cout << "spline.value_unsafe() took: " << dtMs << " ms, sum = " << sum
	     << endl;


	// cout << "x\ty\tspline(x)\tdiff*100\n";
	cout << std::setprecision(4) << std::fixed;
	float maxAbsDiff = 0.0f;
	for (int i = 0; i < Npoints; i++) {
		float appr = spline.value_unsafe(xy(0, i));
		float diff = appr - xy(1, i);
		// cout<<xy(0,i)<<'\t'<<xy(1,i)<<'\t'<<appr<<'\t'<<diff*100.0f<<'\n';
		maxAbsDiff = std::max(maxAbsDiff, std::fabs(diff));
	}
	cout << "maximum absolute diference = " << std::setprecision(8)
	     << maxAbsDiff << endl;
	return 0;
}
