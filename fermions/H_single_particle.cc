#include "H_single_particle.h"
#include <itensor/all.h>
#include <iostream>

// Build single particle Hamiltonian H

// \hat{H} = \dag{c} H c
// with c = (c_1,\sigma c_1,\sigma' , c_2,\sigma , c_2,\sigma' , ...)

// ITensor
// H_number_conserving_fermions(const int N , const vector<double> J, const vector<double> h, const bool spinful = false)
// {
//     int L = 2 * N;
//     vector<double> Jp = J;
//     vector<double> hp = h; 

//     if(spinful)
//     {
//         L = 4 * N;
//         if(length(J) != 2*(N-1) && length(J) !=1)
//         {
//             cerr <<  "Hopping not specified for both species.\n";
//             exit(0);
//         }
//         if(length(J) == 1)
//         {
//             Jp = {};
//             for(int i = 1 ; i <= 2 * N ; i++) Jp.push_back(J[0]);
//         }


    
//         if(length(h) != 2 *N && length(h) !=1)
//         {
//             cerr <<  "Local fields not specified for both species.\n";
//             exit(0);
//         }

//         if(length(h) == 1)
//         {
//             hp = {};
//             for(int i = 1 ; i <= 2 * N ; i++) hp.push_back(h[0]);
//         }


//     }

//     Index r = Index(L);
//     Index l = Index(L);

//     ITensor H = ITensor(r,l);

//     for(int i=1 ; i <= dim(l) ; i++)
//     {
//         for(int j=1 ; j <= dim(r) ; j++)
//         {
//             if(i==j) H.set(l(i),r(j),hp[(i-1)%int(L/2)]);
//             if(abs(i-j)==1)  H.set(l(i),r(j), Jp[(min(i,j)-1)%int(L/2-1)]);
//         }
//     } 

//     PrintData(H);

//     return H;

// }



// ITensor
// H_number_conserving_fermions(const int N , const vector<double> J, const vector<double> h, const bool spinful = false)
// {
//     int L = 2 * N;
//     vector<double> Jp = J;
//     vector<double> hp = h; 

//     if(spinful) L = 2 * N;

//     Index r = Index(L);
//     Index l = prime(r);

//     ITensor H = ITensor(r,l);

//     for(int i=1 ; i <= dim(l) ; i++)
//     {
//         for(int j=1 ; j <= dim(r) ; j++)
//         {   
//             if(i==j) H.set(l(i),r(j),hp[0]);
//             if(abs(i-j)==1)  H.set(l(i),r(j), Jp[0]);
//         }
//     } 

//     return H;

// }

ITensor
H_number_conserving_fermions(const int N , const vector<double> J, const vector<double> h, const bool spinful = false)
{
    int L = N;
    vector<double> Jp = J;
    vector<double> hp = h; 

    if(spinful) L = 2 * N;

    Index r = Index(L);
    Index l = prime(r);

    ITensor H = ITensor(r,l);

    for(int i=1 ; i <= dim(l) ; i++)
    {
        for(int j=1 ; j <= dim(r) ; j++)
        {   
            if(i==j) H.set(l(i),r(j),hp[0]);
            if(abs(i-j)==1)  H.set(l(i),r(j), Jp[0]);
        }
    } 

    return H;

}

ITensor
H_number_conserving_fermions_impurity(const int N , const vector<double> J, const vector<double> h, const bool spinful = false)
{
    int L = N;
    vector<double> Jp = J;
    vector<double> hp = h; 

    if(spinful) L = 2 * N;

    Index r = Index(L);
    Index l = prime(r);

    ITensor H = ITensor(r,l);

    for(int i=1 ; i <= dim(l) ; i++)
    {
        for(int j=1 ; j <= dim(r) ; j++)
        {   
            if(i==j) H.set(l(i),r(j),hp[0]);
            if(abs(i-j)==1)  
            {
                if(j==1 || i==1) H.set(l(i),r(j), Jp[0]);
                else H.set(l(i),r(j), Jp[1]);
            }

        }
    } 

    return H;

}
