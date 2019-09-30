#include<iostream>
#include "armadillo"
#include <vector>
#include "time.h"

using namespace std;
using namespace arma;

// Read a Coloumn Vector from a file
vector<int> readvalues(string file){
    string line;
    vector<int> N_values;
    ifstream values;
    values.open(file,ios::in);
    while (getline(values,line)){
        if (atoi(line.c_str()) != 0){
        N_values.push_back(atoi(line.c_str()));
        }
       else{}
   }
    values.close();
    return N_values;
    }

// Construct a Toeplitz Matrix with equal diagonal elements
mat Toeplitz_Matrix(int N, double d, double a){
    mat T = mat(N-1,N-1,fill::zeros);
    for (int i = 0; i < N-1;i++){
        T.at(i,i) = d;
        if (i != 0 and i != N-2){
            T.at(i,i+1) = a;
            T.at(i,i-1) = a;
        }
        else if(i == 0){
            T.at(i,i+1) = a;
        }
        else if(i == N-2){
            T.at(i,i-1) = a;
        }
    }
    return T;
}

// Construct a Toeplitz Matrix with varying diagonal elements
mat Toeplitz_Matrix_vec(int N, vec d, double a){
    mat T = mat(N-1,N-1,fill::zeros);
    for (int i = 0; i < N-1;i++){
        T.at(i,i) = d(i+1);
        if (i != 0 and i != N-2){
            T.at(i,i+1) = a;
            T.at(i,i-1) = a;
        }
        else if(i == 0){
            T.at(i,i+1) = a;
        }
        else if(i == N-2){
            T.at(i,i-1) = a;
        }
    }
    return T;
}

// Find the offdiagonal maximum of a matrix
vector<int> offdiagmaximum(mat T, int N){

    vector<int> A;
    double offdiagmax = 0;
    uword k = 0;
    uword l = 0;

    for (uword i = 0; i < N-1; i++){
        T(i,i) = 0;
        for (uword j = 0; j < N-1; j++){
            if (T(i,j)<0){T(i,j) = -T(i,j);}
            //T(i,j) = abs(T(i,j));
            if (T(i,j) > offdiagmax){
                offdiagmax = T(i,j);
                k = i;
                l = j;
            }
        }
    }

    if (k!=l){
        A.push_back(k);
        A.push_back(l);
    }
    else {
        A.push_back(0);
        A.push_back(1);}
    return A;
}


// Jacobi's solver for finding eigenvalues of a matrix
mat eigenvalues(mat T, double eps){
    int N = T.col(0).n_elem + 1;
    double kk = 0;
    double ll = 0;
    double ik = 0;
    double il = 0;

    // Orthogonal basis for eigenvectors
    mat R = mat(N-1,N-1,fill::eye);

    // Find offdiagonal maximum value indices
    vector<int> indices = offdiagmaximum(T,N);

    uword k = indices[0];
    uword l = indices[1];
    double offdiagmax = T(k,l);
    int transtimer = 0;

    while (abs(offdiagmax) > eps){

        // Calculating trigonometrics
        double tau = (T(l,l)-T(k,k))/(2*T(k,l));
        double t = 0;
        if(tau>=0){
            t=1/(tau+sqrt(1+tau*tau));}
        else{
            t=-1/(-tau+sqrt(1+tau*tau));}

        double c = 1/sqrt(1 + t*t);
        double s = t*c;

        // Transforming matrix elements
        kk = T(k,k)*c*c - 2*T(k,l)*c*s + T(l,l)*s*s;
        ll = T(l,l)*c*c + 2*T(k,l)*c*s + T(k,k)*s*s;

        T(k,k) = kk;
        T(l,l) = ll;
        T(k,l) = 0;
        T(l,k) = 0;

        double rik = 0;
        double ril = 0;

        for (uword i = 0; i < N-1; i ++){
            if (i != k and i != l){
                ik = c*T(i,k)-s*T(i,l);
                il = c*T(i,l)+s*T(i,k);
                T(i,k) = ik;
                T(k,i) = ik;
                T(i,l) = il;
                T(l,i) = il;
            }
            rik = R(i,k);
            ril = R(i,l);
            R(i,k) = c * rik - s * ril;
            R(i,l) = c * ril + s * rik;
        }

        // Find maximum value of transformed matrix
        indices = offdiagmaximum(T,N);

        k = indices[0];
        l = indices[1];
        offdiagmax = T(k,l);
        transtimer++;

    }

    // Adjust R to accomodate eigenvalues and transformation timer
    R.resize(N,N);

    cout << transtimer << endl << "---" << endl << "---" << endl;

    for (uword i = 0; i < N-1;i++){
        R(N-1,i) = T(i,i);
    }

    R(N-1,N-1) = transtimer;
    return R;

}



