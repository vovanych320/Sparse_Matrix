#pragma once
#ifndef __MATRIX_EXCEPTIONS_H__
#define	__MATRIX_EXCEPTIONS_H__

#include <string>
#include <exception>

using namespace std;


class Exception : public exception
{
protected:
	string message;

public:
	explicit Exception(const string& message) : exception(), message(message)
	{

	};


	virtual ~Exception() throw ()
	{

	};


	string getMessage(void) const
	{
		return this->message;
	};
};


class InvalidDimensionsException : public Exception
{
public:
	InvalidDimensionsException(string message) : Exception(message) {};
};


class InvalidCoordinatesException : public Exception
{
public:
	InvalidCoordinatesException(const string& message) : Exception(message) {};
};


#endif