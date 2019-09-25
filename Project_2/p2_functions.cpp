#include<iostream>
#include "armadillo"
#include <vector>
#include "time.h"

using namespace std;
using namespace arma;

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

//_________________________________________________________________
vector<int> offdiagmaximum(mat T, int N){

    vector<int> A;
    //mat T_zerodiag;
    double offdiagmax = 0;
    uword k = 0;
    uword l = 0;

    //time_t startzero,finishzero;
    //time_t startindex,finishindex;

    //startzero = clock();
    //mat T_zerodiag = abs(T);
    //for (uword i = 0;i<N-1;i++){
        //T(i,i) = 0;
    //}
    //finishzero = clock();

    //startindex = clock();
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
    //finishindex = clock();
    //cout << "-----" << endl;
    //cout << finishzero - startzero << endl;
    //cout << finishindex - startindex << endl;
    //cout << "-----" << endl;

    if (k!=l){
        A.push_back(k);
        A.push_back(l);
    }
    else {
        A.push_back(0);
        A.push_back(1);}
    return A;
}
//----------------------------------------------------------------------


mat eigenvalues(mat T, double eps){
    int N = T.col(0).n_elem + 1;
    double kk = 0;
    double ll = 0;
    double ik = 0;
    double il = 0;
    //double Frobnorm = 0;

    time_t startMax, finishMax;
    double maxtimer = 0;

    //time_t starttrans, finishtrans;
    //double transftimer = 0;

    mat R = mat(N-1,N-1,fill::eye);

    // Find offdiagonal maximum value indices
    startMax = clock();
    vector<int> indices = offdiagmaximum(T,N);
    finishMax = clock();
    maxtimer += (double) (finishMax - startMax)/CLOCKS_PER_SEC;
    //cout << maxtimer << endl << "---" << "---" << endl;

    uword k = indices[0];
    uword l = indices[1];
    double offdiagmax = T(k,l);
    int transtimer = 0;

    //mat test = mat("1 2 3;4 5 13;7 8 9");
    //cout << test(0,1) << endl;
    //cout << test(1,0) << endl;

    while (abs(offdiagmax) > eps){

        // Froebenious Norm
        /*
        for (uword i = 0; i < N-1; i ++){
            for (uword j = 0; j < N-1; j ++){
                Frobnorm += pow(T(i,j),2);
            }

        }
        Frobnorm = sqrt(Frobnorm);
        */

        /* printing
        if (true){
        cout << "****" << endl;
        cout << T << endl;
        cout << "Frobnorm = " << Frobnorm << endl;

        cout << offdiagmax << endl;
        cout << k << "   " << l << endl;
        cout << "****" << endl;

        // Calculate tau, t, c, s
        cout << "k,l = " << T(k,l) << endl;
        cout << "k,k = " << T(k,k) << endl;
        cout << "l,l = " << T(l,l) << endl;
        }
        */

        // Calculating trigonometrics
        double tau = (T(l,l)-T(k,k))/(2*T(k,l));
        double t = 0;
        if(tau>=0){
            t=1/(tau+sqrt(1+tau*tau));}
        else{
            t=-1/(-tau+sqrt(1+tau*tau));}

        double c = 1/sqrt(1 + t*t);
        double s = t*c;

        /* printing
        if (true){
        cout << "^^^^^^^^" << endl;
        cout << "tau = " << tau << endl;
        cout << "t = " << t << endl;
        cout << "c = " << c << endl;
        cout << "s = " << s << endl;
        cout << "^^^^^^^^" << endl;
        }
        */

        //starttrans = clock();
        // Matrix transformation
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
        //finishtrans = clock();
        //transftimer += (double) (finishtrans - starttrans)/CLOCKS_PER_SEC;


        startMax = clock();
        indices = offdiagmaximum(T,N);
        finishMax = clock();
        maxtimer += (double) (finishMax - startMax)/CLOCKS_PER_SEC;

        k = indices[0];
        l = indices[1];
        offdiagmax = T(k,l);
        transtimer++;

    }


    time_t startresize, finishresize;
    double resizetimer = 0;

    startresize = clock();
    R.resize(N,N);
    finishresize = clock();

    //cout << endl << "loop end" << endl;
    //cout << transftimer << endl << "---" << endl;
    cout << maxtimer << endl << "---" << endl;
    //cout << "Resize timer: " << (double)(finishresize - startresize)/CLOCKS_PER_SEC << endl << "---" << endl;

    cout << transtimer << endl << "---" << endl << "---" << endl;

    for (uword i = 0; i < N-1;i++){
        R(N-1,i) = T(i,i);
    }

    R(N-1,N-1) = transtimer;
    return R;

}



