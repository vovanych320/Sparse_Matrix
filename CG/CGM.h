#pragma once


#include <iostream>
#include <cmath>


#include "Matrix.h"


template<typename T> void prec(const vector<T>& x, vector<T>& y)
{
	size_t size = y.size();

	y[0] = x[0];
	for (size_t i = 0; i < size; i++)
	{
		y[i] = 0.5 * x[i];
	}
};


template<typename T> double vector_sum(const vector<T>& a)
{
	size_t size = a.size();
	double sum = 0;

	for (size_t i = 0; i < size; i++)
	{
		sum += fabs(a[i]);
	}

	return sum;
};


template<typename T>double dot(const vector<T>& a, const vector<T>& b)
{
	double sum = 0;
	size_t size = a.size();

	for (size_t i = 0; i < size; i++)
	{
		sum += a[i] * b[i];
	}

	return sum;
};


template <typename Matrix, typename Vector>int conj_grad_method(const Matrix& A, Vector& x, const Vector& b, const double& eps = 0.00001)
{
	double rho = 0;
	double rho_1 = 0;
	double alpha = 0;

	size_t uni_size = x.size();

	Vector p(uni_size), q(uni_size), r(uni_size), z(uni_size);

	r = b - A * x;
	int iter = 0;

	while (vector_sum(r) >= eps)
	{
		prec(r, z);

		rho = dot(r, z);

		if (!iter)
			p = z;
		else
			p = z + (rho / rho_1) * p;
		q = A * p;
		alpha = rho / dot(p, q);
		x += alpha * p;
		r -= alpha * q;
		rho_1 = rho;
		++iter;
	}

	return iter;
};



template <typename Matrix, typename Vector>int CGM(const Matrix& A, Vector& x, const Vector& b)
{
	if (isSymetric(A))
	{
		return conj_grad_method(A, x, b);
	}
	else
	{
		return -1;
	}
}