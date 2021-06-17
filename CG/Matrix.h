#pragma once

#ifndef __MATRIX_H__

#define	__MATRIX_H__

#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

#include "Exceptions.h"

template<typename T> class Matrix
{
private:
	int m, n;

	vector<T>* vals;
	vector<int>* rows, * cols;


public:
	Matrix(const int& n);
	Matrix(const int& r, const int& c);

	Matrix(const Matrix<T>& m);
	Matrix<T>& operator = (const Matrix<T>& m);

	~Matrix(void);

	int getAmountOfRows() const;
	int getAmountOfColumns() const;


	T get(const int& row, const int& col) const;
	Matrix& set(const T& val, const int& row, const int& col);


	vector<T> operator* (const vector<T>& x) const;
	Matrix<T> operator* (const Matrix<T>& m) const;
	Matrix<T> operator+ (const Matrix<T>& m) const;
	Matrix<T> operator- (const Matrix<T>& m) const;
};

template<typename T> bool isSymetric(const Matrix<T>& a);
template<typename X> ostream& operator<< (ostream& os, const Matrix<X>& matrix);
template<typename T> ostream& operator<< (ostream& os, const vector<T>& a);
template<typename T>vector<T> operator- (const vector<T>& a, const vector<T>& b);
template<typename T>vector<T> operator* (const T& a, const vector<T>& b);
template<typename T>vector<T> operator+ (const vector<T>& a, const vector<T>& b);
template<typename T>vector<T> operator+= (vector<T>& a, const vector<T>& b);
template<typename T>vector<T> operator-= (vector<T>& a, const vector<T>& b);



//-----------------------Implementaation-----------------------//


template<typename T>Matrix<T>::Matrix(const int& n)
{
	int rows = n;
	int columns = n;


	if (rows < 1 || columns < 1)
	{
		throw InvalidDimensionsException("The dimensions of the Matrix must be more then zero.");
	}

	this->m = rows;
	this->n = columns;

	this->vals = NULL;
	this->cols = NULL;
	this->rows = new vector<int>(rows + 1, 1);

};


template<typename T>Matrix<T>::Matrix(const int& r, const int& c)
{
	int rows = r;
	int columns = c;


	if (rows < 1 || columns < 1)
	{
		throw InvalidDimensionsException("The dimensions of the Matrix must be more then zero.");
	}

	this->m = rows;
	this->n = columns;

	this->vals = NULL;
	this->cols = NULL;
	this->rows = new vector<int>(rows + 1, 1);
};


template<typename T>Matrix<T>::Matrix(const Matrix<T>& matrix)
{
	this->m = matrix.m;
	this->n = matrix.n;
	this->rows = new vector<int>(*(matrix.rows));

	if (matrix.vals != NULL)
	{
		this->cols = new vector<int>(*(matrix.cols));
		this->vals = new vector<T>(*(matrix.vals));
	}
};


template<typename T>Matrix<T>& Matrix<T>::operator = (const Matrix<T>& matrix)
{
	if (&matrix != this)
	{
		if (this->vals != NULL)
		{
			delete this->vals;
			delete this->cols;
		}
		delete this->rows;


		this->m = matrix.m;
		this->n = matrix.n;
		this->rows = new vector<int>(*(matrix.rows));
		if (matrix.vals != NULL)
		{
			this->cols = new vector<int>(*(matrix.cols));
			this->vals = new vector<T>(*(matrix.vals));
		}
	}

	return *this;
};


template<typename T>Matrix<T>::~Matrix()
{
	if (this->vals != NULL)
	{
		delete this->vals;
		delete this->cols;
	}

	delete this->rows;
};


template<typename T>int Matrix<T>::getAmountOfRows() const
{
	return this->m;
};


template<typename T>int Matrix<T>::getAmountOfColumns() const
{
	return this->n;
};


template<typename T>T Matrix<T>::get(const int& row, const int& col) const
{
	if (row < 1 || col < 1 || row > this->m || col > this->n)
	{
		throw InvalidCoordinatesException("Coordinates out of range.");
	}


	int currCol;

	for (int pos = (*(this->rows))[row - 1] - 1; pos < (*(this->rows))[row] - 1; ++pos)
	{
		currCol = (*(this->cols))[pos];

		if (currCol == col)
		{
			return (*(this->vals))[pos];
		}
		else if (currCol > col)
		{
			break;
		}
	}

	return T();
};


template<typename T>Matrix<T>& Matrix<T>::set(const T& val, const int& row, const int& col)
{
	if (row < 1 || col < 1 || row > this->m || col > this->n)
	{
		throw InvalidCoordinatesException("Coordinates out of range.");
	}

	int pos = (*(this->rows))[row - 1] - 1;
	int currCol = 0;

	for (; pos < (*(this->rows))[row] - 1; pos++)
	{
		currCol = (*(this->cols))[pos];

		if (currCol >= col)
		{
			break;
		}
	}

	if (currCol != col)
	{
		if (!(val == T()))
		{
			int index = pos;

			if (this->vals == NULL)
			{
				this->vals = new vector<T>(1, val);
				this->cols = new vector<int>(1, col);
			}
			else
			{
				this->vals->insert(this->vals->begin() + index, val);
				this->cols->insert(this->cols->begin() + index, col);
			}
			for (int i = row; i <= this->m; i++)
			{
				(*(this->rows))[i] += 1;
			}
		}
	}
	else if (val == T())
	{
		int index = pos;
		this->vals->erase(this->vals->begin() + index);
		this->cols->erase(this->cols->begin() + index);

		for (int i = row; i <= this->m; i++)
		{
			(*(this->rows))[i] -= 1;
		}
	}
	else
	{
		(*(this->vals))[pos] = val;
	}

	return *this;
};





//-----------------------Utylitues-----------------------//

template<typename T> bool isSymetric(const Matrix<T>& a)
{
	size_t cols = a.getAmountOfColumns();
	size_t rows = a.getAmountOfRows();

	for (size_t i = 1; i < cols; i++)
	{
		for (size_t j = 1; j < rows; j++)
		{
			if (i != j)
			{

				if (a.get(j, i) != a.get(i, j))
				{
					return false;
				}

			}
		}

	}

	return true;
};


template<typename T>vector<T> Matrix<T>::operator* (const vector<T>& x) const
{
	if (this->n != (int)x.size())
	{
		throw InvalidDimensionsException("Cannot multiply: Matrix column count and vector size don't match.");
	}

	vector<T> result(this->m, T());

	if (this->vals != NULL)
	{
		for (int i = 0; i < this->m; i++)
		{
			T sum = T();
			for (int j = (*(this->rows))[i]; j < (*(this->rows))[i + 1]; j++)
			{
				sum = sum + (*(this->vals))[j - 1] * x[(*(this->cols))[j - 1] - 1];
			}

			result[i] = sum;
		}
	}

	return result;
};


template<typename T>Matrix<T> Matrix<T>::operator* (const Matrix<T>& m) const
{
	if (this->n != m.m)
	{
		throw InvalidDimensionsException("Cannot multiply: Left matrix column count and right matrix row count don't match.");
	}

	Matrix<T> result(this->m, m.n);

	T a;

	for (int i = 1; i <= this->m; i++)
	{
		for (int j = 1; j <= m.n; j++)
		{
			a = T();

			for (int k = 1; k <= this->n; k++)
			{
				a = a + this->get(i, k) * m.get(k, j);
			}

			result.set(a, i, j);
		}
	}

	return result;
};


template<typename T>Matrix<T> Matrix<T>::operator+ (const Matrix<T>& m) const
{
	if (this->m != m.m || this->n != m.n)
	{
		throw InvalidDimensionsException("Cannot add: matrices dimensions don't match.");
	}

	Matrix<T> result(this->m, this->n);


	for (int i = 1; i <= this->m; i++)
	{
		for (int j = 1; j <= this->n; j++)
		{
			result.set(this->get(i, j) + m.get(i, j), i, j);
		}
	}

	return result;
};


template<typename T>Matrix<T> Matrix<T>::operator- (const Matrix<T>& m) const
{
	if (this->m != m.m || this->n != m.n)
	{
		throw InvalidDimensionsException("Cannot subtract: matrices dimensions don't match.");
	}

	Matrix<T> result(this->m, this->n);


	for (int i = 1; i <= this->m; i++)
	{
		for (int j = 1; j <= this->n; j++)
		{
			result.set(this->get(i, j) - m.get(i, j), i, j);
		}
	}

	return result;
};


template<typename T>ostream& operator << (ostream& os, const Matrix<T>& matrix)
{
	for (int i = 1; i <= matrix.getAmountOfRows(); i++)
	{
		os << "(";
		for (int j = 1; j <= matrix.getAmountOfColumns(); j++)
		{
			if (j != 1)
			{
				os << " ";
			}

			os << matrix.get(i, j);
		}

		if (i < matrix.getAmountOfRows())
		{
			os << ")" << endl;
		}
	}
	os << ")";

	return os;
};


template<typename T> ostream& operator<<(ostream& os, const vector<T>& a)
{
	size_t size = a.size();
	cout << "(";
	for (size_t i = 0; i < size; i++)
	{
		cout << a[i] << (i != (size - 1) ? "," : ")");
	}
	return os;
};


template<typename T>vector<T> operator-(const vector<T>& a, const vector<T>& b)
{
	size_t size = a.size();

	vector<T> res(size);

	for (size_t i = 0; i < size; i++)
	{
		res[i] = a[i] - b[i];
	}

	return res;
};


template<typename T>vector<T> operator*(const T& a, const vector<T>& b)
{
	size_t size = b.size();
	vector<T> res(size);


	for (size_t i = 0; i < size; i++)
	{
		res[i] = a * b[i];
	}

	return res;
};


template<typename T>vector<T> operator+ (const vector<T>& a, const vector<T>& b)
{
	size_t size = a.size();
	vector<T> res(size);

	for (size_t i = 0; i < size; i++)
	{
		res[i] = a[i] + b[i];
	}

	return res;
};


template<typename T>vector<T> operator+= (vector<T>& a, const vector<T>& b)
{
	size_t size = a.size();


	for (size_t i = 0; i < size; i++)
	{
		a[i] += b[i];
	}

	return a;
};


template<typename T>vector<T> operator-= (vector<T>& a, const vector<T>& b)
{
	size_t size = a.size();

	for (size_t i = 0; i < size; i++)
	{
		a[i] -= b[i];
	}

	return a;

};



#endif