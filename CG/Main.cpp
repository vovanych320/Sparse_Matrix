#include "CGM.h"

int main()
{

	vector<double> x(3);
	
	vector<double> b(3);
	b[0] = 1;
	b[1] = 5;
	b[2] = 9;

	vector<double> r(3);

	Matrix<double> A(3);
	A.set(3, 1, 1);
	A.set(4, 1, 2);
	A.set(0, 1, 3);

	A.set(4, 2, 1);
	A.set(-3, 2, 2);
	A.set(0, 2, 3);

	A.set(0, 3, 1);
	A.set(0, 3, 2);
	A.set(5, 3, 3);
	
	cout << "Vector b " << b << endl;

	cout << A << endl << endl;

	int con = CGM(A, x, b);

	


	for (size_t i = 0; i < 3; i++)
	{
		cout << x[i] << " ";
	}


	cout << "dfdfv " << con << endl;


	Matrix<double> two(2);
	two.set(4, 1, 1);
	two.set(1, 1, 2);
	two.set(1, 2, 1);
	two.set(3, 2, 2);

	cout << two << endl << endl;

	vector<double> x1(2);

	vector<double> b1(2);
	b1[0] = 1;
	b1[1] = 2;


	cout << endl;

	for (size_t i = 0; i < 2; i++)
	{
		cout << x1[i] << " ";
	}

	
	

	

	return 0;
}



