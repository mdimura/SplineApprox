#ifndef SPLINE_H
#define SPLINE_H

#include "polynomial.hpp"
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <iostream>

template <unsigned Deg> class Spline
{
private:
	std::vector<Polynomial<Deg>> polynomials;
	float xmin, xmax;
	float invdx;
	Spline() = default;
	void setFromSorted(const Eigen::Matrix2Xf &sortedXY,
			   const unsigned NPieces);

public:
	float value_unsafe(const float &x) const
	{
		int i = (x - xmin) * invdx;
		assert(i >= 0);
		assert(i < polynomials.size());
		return polynomials[i].value(x);
	}
	float value_or(const float &x, const float &def) const
	{
		if (x < xmin || x > xmax) {
			return def;
		}
		int i = (x - xmin) * invdx;
		return polynomials[i].value(x);
	}
	static Spline<Deg> fromSorted(const Eigen::Matrix2Xf &sortedXY,
				      const unsigned NPieces);

	template <unsigned Nd>
	friend std::ostream &operator<<(std::ostream &os,
					const Spline<Nd> &spline);
};
template <unsigned Deg>
std::ostream &operator<<(std::ostream &os, const Spline<Deg> &spline)
{
	for (int i = 0; i < spline.polynomials.size(); ++i) {
		float dx = 1.0 / spline.invdx;
		float from = spline.xmin + dx * i;
		float to = from + dx;
		os << '[' << from << ".." << to << "] " << spline.polynomials[i]
		   << '\n';
	}
	return os;
}

template <unsigned int Deg>
Spline<Deg> Spline<Deg>::fromSorted(const Eigen::Matrix2Xf &sortedXY,
				    const unsigned NPieces)
{
	Spline<Deg> sp;
	sp.setFromSorted(sortedXY, NPieces);
	return sp;
}
template <unsigned int Deg>
void Spline<Deg>::setFromSorted(const Eigen::Matrix2Xf &sortedXY,
				const unsigned NPieces)
{
	const int nPoints = sortedXY.cols();
	assert(nPoints > Deg);

	xmin = sortedXY(0, 0);
	xmax = sortedXY(0, nPoints - 1);
	constexpr float inf = std::numeric_limits<float>::infinity();
	const float dx = nextafterf((xmax - xmin) / NPieces, inf);
	invdx = 1.0f / dx;

	assert((xmax - xmin) * invdx < NPieces);

	polynomials.resize(NPieces);
	int iStart = 0;
	for (int piece = 0; piece < NPieces - 1; ++piece) {
		float xEnd = dx * (piece + 1) + xmin;
		int iEnd = iStart;
		while (sortedXY(0, iEnd) < xEnd) {
			++iEnd;
		}
		int len = std::max(iEnd - iStart, int(Deg + 1));
		Polynomial<Deg> p(sortedXY.middleCols(iStart, len));
		polynomials[piece] = std::move(p);
		if (iEnd > iStart + Deg) {
			iStart = std::min(int(nPoints - Deg - 1), iEnd);
		}
	}

	Polynomial<Deg> p(sortedXY.middleCols(iStart, nPoints - iStart));
	polynomials[NPieces - 1] = std::move(p);
}

#endif // SPLINE_H
