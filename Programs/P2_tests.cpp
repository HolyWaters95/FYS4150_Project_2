#include <iostream>
#include "armadillo"
#include <cmath>
#include "catch.hpp"
#include "time.h"

using namespace std;
using namespace arma;

vector<int> readvalues(string file);
mat Toeplitz_Matrix(int N, double d, double a);
mat Toeplitz_Matrix_vec(int N, vec d, double a);
vector<int> offdiagmaximum(mat T, int N);
mat eigenvalues(mat T, double eps);

// Tests the eigenvalue solver for a matrix of size 5x5
// The test is performed by comparing with analytical eigenvalues
// And comparing Frobenious norms of original and final matrix
void eigtest_N_5(){
    cout << "Running tests" << endl;
    int N = 5;
    vec X = linspace(0,1,N+1);
    double h = (X.back() - X[0])/N;
    double eps = pow(10,-10);

    // Construct Toeplitz Matrix
    double d = 2/pow(h,2);
    double a = -1/pow(h,2);
    mat T = Toeplitz_Matrix(N,d,a);
    // mat T = mat("0 0 1;5 2 13;0 0 1");

    mat D = eigenvalues(T, eps);

    vec lambdas = reshape(D.row(N-1),N,1);
    lambdas.resize(N-1); vec sorted_lambdas = sort(lambdas);
    D.resize(N-1,N-1); mat sorted_D = mat(N-1,N-1,fill::zeros);

    for (uword i = 0; i < N-1;i++){
        for (uword j = 0;j < N-1;j++){
            if (sorted_lambdas(i) == lambdas(j)){
                sorted_D.col(i) = D.col(j);
            }
        }
    }

    vec exact(N-1);
    double EPS = pow(10,-8);

    //Analytical Eigenvalues
    for (int i = 1;i < N;i++){
        exact(i-1) = d + 2*a*cos((i*M_PI)/(N));
    }
    // Check Eigenvalues
    for (uword i = 0;i<4;i++){
    if (sorted_lambdas(i) > exact(i) - EPS and sorted_lambdas(i) < exact(i)+EPS){
        cout << "Lambda " << i+1 << " is correct" << endl;
    }

    else{cout << "LAMBDA " << i+1 << " NOT CORRECT" << endl << "Exact value: " << exact(i) << " | " << "Calculated value: " << sorted_lambdas(i) << endl;}
    }


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
    Frobnorm2 += pow(sorted_lambdas(i),2);
    }
    Frobnorm2 = sqrt(Frobnorm2);
    if (Frobnorm1 > Frobnorm2 - EPS and Frobnorm1 < Frobnorm2 + EPS){
        cout << "Frobenious norm is preserved" << endl;
    }
    else{
        cout << "FROBENIOUS NORM IS NOT PRESERVED" << endl;
        cout << "ORIGINAL NORM: " << Frobnorm1 << " | " << "FINAL NORM: " << Frobnorm2 << endl;
    }

}


// Tests eigenvalue solver for a matrix of size 5x5
// The test is performed by comparing with the Armadillo solver eigenvalues
void armadillotest_quantum_N_5(){
    cout << "Running tests" << endl;
    int N = 5;
    vec X = linspace(0,1,N+1);
    double h = (X.back() - X[0])/N;
    double eps = pow(10,-10);
    double w = 1;

    // Construct Toeplitz Matrix
    double a = -1/pow(h,2);
    vec d = vec(X.n_elem);
    for (uword i = 0; i<N+1;i++) {d(i) = 2/pow(h,2) + w*w*X(i)*X(i) + 1/X(i);}
    mat T = Toeplitz_Matrix_vec(N,d,a);
    // mat T = mat("0 0 1;5 2 13;0 0 1");

    mat D = eigenvalues(T, eps);

    vec lambdas = reshape(D.row(N-1),N,1);
    lambdas.resize(N-1); vec sorted_lambdas = sort(lambdas);
    D.resize(N-1,N-1); mat sorted_D = mat(N-1,N-1,fill::zeros);

    for (uword i = 0; i < N-1;i++){
        for (uword j = 0;j < N-1;j++){
            if (sorted_lambdas(i) == lambdas(j)){
                sorted_D.col(i) = D.col(j);
            }
        }
    }

    vec armaeigval;
    mat armaeigvec;
    eig_sym(armaeigval,armaeigvec,T);


    // Check Eigenvalues
    double EPS = pow(10,-8);
    for (uword i = 0;i<4;i++){
    if (sorted_lambdas(i) > armaeigval(i) - EPS and sorted_lambdas(i) < armaeigval(i)+EPS){
        cout << "Lambda " << i+1 << " is correct according to arma::eigsym" << endl;
    }

    else{cout << "LAMBDA " << i+1 << " NOT CORRECT ACCORDING TO ARMA::EIGSYM" << endl << "Armadillo value: " << armaeigval(i) << " | " << "Calculated value: " << sorted_lambdas(i) << endl;}
    }
    //cout << D << endl << endl << armaeigvec << endl;

    mat diffmat = abs(sorted_D) - abs(armaeigvec);
    if (norm(diffmat,"fro") < EPS){
        cout << "Eigenvectors are correct according to arma::eigsym" << endl;
    }
    else{cout << "EIGENVECTORS NOT CORRECT ACCORDING TO ARMA::EIGSYM" << endl << "CALCULATED MATRIX:" << endl << sorted_D << endl << "ARMA MATRIX:" << endl << armaeigvec << endl;}
}
