#include <iostream>
#include "armadillo"
#include <cmath>
#include "catch.hpp"
#include "time.h"

using namespace std;
using namespace arma;

vector<int> readvalues(string file);
mat Toeplitz_Matrix(int N, double d, double a);
vector<int> offdiagmaximum(mat T, int N);
mat eigenvalues(mat T, double eps);
void eigtest_N_5();




int main(){

    // Run test of eigenvalue solver
    eigtest_N_5();


    // Values for N
    vector<int> N_values;
    N_values = readvalues("N_values.txt");

    // To save or not to save?
    string save_transformations;
    string save_eigenpairs;
    string save_eigexact;


    cout << "do you want to save transformations? y or n" << endl;
    cin >> save_transformations;
    cout << "do you want to save eigenpairs? y or n" << endl;
    cin >> save_eigenpairs;
    cout << "do you want to save analytical eigenvalues? y or n" << endl;
    cin >> save_eigexact;
    cout << "-------------------------------" << endl;

    vec N_transformations(N_values.size());
    vec Runtimes(N_values.size());

    // Run segment for all values of N
    for (int i=0;i<N_values.size();i++){
        // Setting up boundary conditions and step length
        int N = N_values[i];
        vec X = linspace(0,1,N+1);
        double h = (X.back() - X[0])/N;
        double eps = pow(10,-10);

        // Construct Toeplitz Matrix
        double d = 2/pow(h,2);
        double a = -1/pow(h,2);
        mat T = Toeplitz_Matrix(N,d,a);
        // mat T = mat("0 0 1;5 2 13;0 0 1");

        // Run eigenvalue solver and print runtime
        time_t start, finish;
        start = clock();
        mat D = eigenvalues(T, eps);
        finish = clock();

        cout << (double)(finish-start)/CLOCKS_PER_SEC << endl;
        Runtimes(i) = (double)(finish-start)/CLOCKS_PER_SEC;

        // Reorganizing: Get sorted eigenvalues in a separate vector, and sort eigenvectors accordingly
        // Also, separate out N_transformations
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


        // Calculate analytical eigenvalues
        vec analytical_eigenvalues = vec(N-1);
        for (int i = 1;i < N;i++){
            analytical_eigenvalues(i-1) = d + 2*a*cos((i*M_PI)/(N));
        }

        // Save eigenpairs (both computed and exact) to file
        if (save_eigenpairs == "y"){
            string filenameeigen = "Eigenpairs_" + to_string(N) + ".txt";
            ofstream output;
            output.open(filenameeigen,ios::out);
            output << sorted_lambdas << endl;
            for (uword i = 0; i < N-1;i++){
                output << sorted_D.col(i) << endl;
            }
            output.close();
            }
            else{}
        if (save_eigexact == "y"){
            string filenameeigexact = "Eigenvalues_exact_" + to_string(N) + ".txt";
            ofstream output;
            output.open(filenameeigexact,ios::out);
            output << analytical_eigenvalues << endl;
            output.close();
            }
            else{}
    } // End of N-loop


    // Save transformation number and runtimes to file
    if (save_transformations == "y"){
    string filenametrans = "Transformations and Runtimes.txt";
    ofstream output;
    output.open(filenametrans,ios::out);
    for (int i = 0;i<N_values.size();i++){
        output << N_values[i] << endl;
    }
    output << endl;
    output << N_transformations << endl;
    output << Runtimes << endl;
    output.close();
    }
    else{}



    return 0;
}
