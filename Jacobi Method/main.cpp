#include <iostream>
#include <fstream> //read and write file
#include <iomanip> //format output to the screen
#include <cmath> //for absolute value functions
/*
The solution is obtained iteratively via
x_(k+1) = D^-1(b - (L+U)x_k)
Do this formula until convergence is reached

*/
using namespace std;
int main()
{
	//---------- Set up ----------//

	//Variables
	int n = 0;
	int iteration = 0;
	int setWidth = 10; //used for output formatting
	int setPrecision = 20; //used for output formatting
	double actualValue, checkValue = 0.0;
	double error = 0.0;
	double tolerableError = 0.001;
	double sumOfb = 0.0; //used to figure out the check value

	//Matrices
	double **A; //LHS
	double **L_U; //L + U
	double **D; //Diagonal

	//Vectors
	double *b; //RHS
	double *newB; //used to calculate actual value
	double *x; //solution
	double *workingVector;

	//File streams
	std::ifstream fin;
	std::ofstream fout;

	//Open file containing the matrix
	fin.open("q12.txt");

	//Find dimensions
	fin >> n;
    cout << "DEBUG: n = " << n << endl;

	//Finish constructing matrices
	A = new double *[n]; //LHS
	for (int i = 0; i < n; i++)
		A[i] = new double[n];

    D = new double *[n]; //LHS
	for (int i = 0; i < n; i++)
		D[i] = new double[n];

    L_U = new double *[n];
	for (int i = 0; i < n; i++)
		L_U[i] = new double[n];

	//Fill A, D, L_U
	for(int row = 0; row < n; row++)
	{
		for (int col = 0; col < n; col++)
		{
			fin >> A[row][col];
			if(row == col)
			{
			    if(A[row][col] != 0)
                    D[row][col] = 1/A[row][col]; //D gets initialized inverted
                else
                    D[row][col] = 0.0;

			    L_U[row][col] = 0.0;
			}
			else
			{
			    L_U[row][col] = A[row][col];
			    D[row][col] = 0.0;
			}
		}
	}

	//Finish constructing vectors
	b = new double[n]; //RHS
	newB = new double[n];
	x = new double[n]; //solution vector
	workingVector = new double[n];

	//Fill RHS,workingVector, and solution vector
	for(int row = 0; row < n; row++)
	{
		fin >> b[row];
		x[row] = 0.0; //make first guess here
		newB[row] = 0.0;
	}
	fin.close();

	//Initialize check value
	for(int i = 0; i < n; i++)
    {
        sumOfb += b[i] * b[i];
    }
    checkValue = sqrt(sumOfb);

	//---------- Jacobi Iteration ----------//
	do
    {
        //set working vector to x
        for(int i = 0; i < n; i++){workingVector[i] = x[i];}

        //(L+U)x_i-1
        for(int i = 0; i < n; i ++)
        {
            double sumOfThis = 0.0;
            for(int j = 0; j < n; j++)
            {
                sumOfThis = sumOfThis + (L_U[i][j] * workingVector[j]);
            }
            x[i] = sumOfThis;
        }

        //b-((L+U)x_i-1)
        for(int i = 0; i < n; i++)
        {
            x[i] = b[i] - x[i];
        }

        //set working vector to x
        for(int i = 0; i < n; i++){workingVector[i] = x[i];}

        //D^-1(b-((L+U)x_i-1))
        for(int i = 0; i < n; i++)
        {
            double sumOfThis = 0.0;
            for(int j = 0; j < n; j++)
            {
                sumOfThis += D[i][j] * workingVector[j];
            }
            x[i] = sumOfThis;
        }

        //calculate actualValue
        for(int i = 0; i < n; i++)
        {
            double sumOfThis = 0.0;
            for(int j = 0; j < n; j++)
            {
                sumOfThis += A[i][j] * x[j];
            }
            newB[i] = sumOfThis;
        }
        double sumOfNewB = 0.0;
        for(int i = 0; i < n; i++)
        {
            sumOfNewB = sumOfNewB + newB[i] * newB[i];
        }
        actualValue = sqrt(sumOfNewB);
        //calculate error
        error = fabs(checkValue - actualValue);

        //-------Print out-------//
        std::cout << "checkValue = " << checkValue << std::endl;
        std::cout << "actual = " << actualValue << std::endl;
        iteration++;
        std::cout << "error = " << error << std::endl;
        std::cout << "----------------- iteration = " << iteration << "-----------------" << std::endl;
        std::cout << std::endl;

        //fail stop
        if(iteration == 1000)
        {
            goto STOP;
        }
    }
    while(error > tolerableError);

    STOP:
    int rootN = sqrt(n);
    fout.open("q12_Solution.txt");
    fout << x[0] << " ";
    for(int i = 1; i < n; i++)
    {
        if(i%rootN == 0)
            fout << std::endl;
        fout << x[i] << " ";
    }

    return 0;
}
