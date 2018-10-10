#SplineApprox
Trivial header-only implementation of polynomial spline approximation in C++.

##Usage:
```C++
#include "spline.hpp"

...

const int Npoints = 1024;
const int Npieces = 8;
constexpr unsigned Deg = 2; #Polynomial degree

//Discretize the sine function for [-pi/2..pi/2]
const int Npoints = 1024;
Eigen::Matrix2Xf xy(2, Npoints);
const float pi = 3.14159f;
for (int i = 0; i < Npoints; ++i) {
    xy(0, i) = -pi * 0.5f + pi * i / (Npoints - 1);
    xy(1, i) = std::sin(xy(0, i));
}

// initialize the spline
Spline<Deg> spline = Spline<Deg>::fromSorted(xy, Npieces);

// with boundary checks 
float x = pi/2.0f+1.0f;
cout << spline.value_or(x, 999.9f) << '\n';

// no boundary checks 
float x = pi/3.0f;
cout << spline.value_unsafe(x) << '\n';
```