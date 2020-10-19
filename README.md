# SplineApprox
Trivial header-only implementation of polynomial spline approximation in C++.

## Usage:
```C++
#include "spline.hpp"

...

const int Npieces = 8;      // Number of polynomials in the spline
constexpr unsigned Deg = 2; // Polynomial degree

const int Npoints = 10240; //Number of data points
const float pi = atan(1.0f)*4.0f;

//Generate reference data
Eigen::Matrix2Xf xy(2, Npoints);
xy.row(0) = Eigen::VectorXf::LinSpaced(Npoints, -pi*0.5f, pi*0.5f);
xy.row(1) = xy.row(0).array().sin();

// initialize the spline
Spline<Deg> spline = Spline<Deg>::fromSorted(xy, Npieces);

//Calculate approximated values
float x = pi/2.0f+1.0f;
cout << spline.value_or(x, 999.9f) << '\n'; // with boundary checks

float x = pi/3.0f;
cout << spline.value_unsafe(x) << '\n'; // no boundary checks 
```
