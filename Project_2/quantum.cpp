#include <iostream>
#include "armadillo"
#include <cmath>
#include "catch.hpp"
#include "time.h"

using namespace std;
using namespace arma;

vector<int> readvalues(string file);
mat Toeplitz_Matrix_vec(int N, vec d, double a);
vector<int> offdiagmaximum(mat T, int N);
mat eigenvalues(mat T, double eps);

/*

TEST_CASE("Test for N = 5"){
    int N = 5;
    vec X = linspace(0,1,N+1);
    double h = (X.back() - X[0])/N;
    double eps = pow(10,-10);

    // Construct Toeplitz Matrix
    double d = 2/pow(h,2);
    double a = -1/pow(h,2);
    mat T = Toeplitz_Matrix(N,d,a);
    // mat T = mat("0 0 1;5 2 13;0 0 1");

    vec D = eigenvalues(T, eps);

    vec exact(N-1);

    //Analytical Eigenvalues
    for (int i = 1;i < N;i++){
        exact(i) = d + 2*a*cos((i*M_PI)/(N));
    }
    // Check Eigenvalues
    REQUIRE(D(0)==Approx(exact(0)).epsilon(0.0000000001));
    REQUIRE(D(1)==Approx(exact(1)).epsilon(0.0000000001));
    REQUIRE(D(2)==Approx(exact(2)).epsilon(0.0000000001));
    REQUIRE(D(3)==Approx(exact(3)).epsilon(0.0000000001));

    //Check Frobenius Norm:
    double Frobnorm1 = 0;
    // Froebenious Norm
    for (uword i = 0; i < N-1; i ++){
    for (uword j = 0; j < N-1; j ++){
          Frobnorm1 += pow(T(i,j),2);
            }

        }
        Frobnorm1 = sqrt(Frobnorm1);
    double Frobnorm2 = 0;
    for (uword i = 0;i<N-1;i++){
    Frobnorm2 += pow(D(i),2);
    }
    Frobnorm2 = sqrt(Frobnorm2);
    REQUIRE(Frobnorm1 == Approx(Frobnorm2).epsilon(0.0000000001));
}

*/


int main()
{
    // Values for N
    vector<int> N_values;
    N_values = readvalues("N_values.txt");
    vector<int> Xmax_values = readvalues("Xmax_values.txt");


    string save_transformations;
    string save_eigenpairs;

    cout << "do you want to save transformations? y or n" << endl;
    cin >> save_transformations;
    cout << "do you want to save eigenpairs? y or n" << endl;
    cin >> save_eigenpairs;
    cout << "-------------------------------" << endl;

    vec N_transformations(N_values.size());

    // Run segment for all values of N
    for (int i=0;i<N_values.size();i++){
        for (int j = 0; j<Xmax_values.size();j++){
        int N = N_values[i];
        int Xmax = Xmax_values[j];
        vec X = linspace(0,Xmax,N+1);
        double h = (X.back() - X[0])/N;
        double eps = pow(10,-8);

        // Construct Toeplitz Matrix
        vec d = vec(X.n_elem);
        for (uword i = 0; i<N+1;i++) {d(i) = 2/pow(h,2) + X(i)*X(i);}
        double a = -1/pow(h,2);
        mat T = Toeplitz_Matrix_vec(N,d,a);
        // mat T = mat("0 0 1;5 2 13;0 0 1");

        //cout << T << endl;

        time_t start, finish;
        start = clock();
        mat D = eigenvalues(T, eps);
        finish = clock();

        cout << (double)(finish-start)/CLOCKS_PER_SEC << endl;


        vec lambdas = reshape(D.row(N-1),N,1);
        N_transformations(i) = lambdas(N-1);
        lambdas.resize(N-1); vec sorted_lambdas = sort(lambdas);
        D.resize(N-1,N-1); mat sorted_D = mat(N-1,N-1,fill::zeros);

        for (uword i = 0; i < N-1;i++){
            for (uword j = 0;j < N-1;j++){
                if (sorted_lambdas(i) == lambdas(j)){
                    sorted_D.col(i) = D.col(j);
                }
            }
        }

        mat eigenvecs = mat(N+1,N-1,fill::zeros);
        for (uword i = 1;i<N;i++){for (uword j = 0;j<N-1;j++){
                eigenvecs(i,j) = sorted_D(i-1,j);
            }}

        // Printing stuff
        /*
        cout << "final prints" << endl;
        cout << sorted_D << endl;
        cout << sorted_lambdas << endl;
        cout << "---" << endl;
        cout << N_transformations(i) << endl;
        cout << "next N" << endl;
        */


        if (save_eigenpairs == "y"){
            string filenameeigen = "Quantum_Eigenpairs_N_" + to_string(N) + "_X_" + to_string(Xmax) + "_test" + ".txt";
            ofstream output;
            output.open(filenameeigen,ios::out);
            output << sorted_lambdas << endl;
            for (uword i = 0; i < N-1;i++){
                output << eigenvecs.col(i) << endl;
            }
            output.close();
            }
            else{}

        cout << "N = " + to_string(N) << " , Xmax = " + to_string(Xmax) << endl << "----" << endl;
    }} // End of N-loop



    if (save_transformations == "y"){
    string filenametrans = "Quantum_Number of transformations.txt";
    ofstream output;
    output.open(filenametrans,ios::out);
    for (int i = 0;i<N_values.size();i++){
        output << N_values[i] << endl;
    }
    output << endl;
    output << N_transformations << endl;
    output.close();
    }
    else{}



    return 0;
}
