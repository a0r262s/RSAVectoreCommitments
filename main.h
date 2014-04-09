#ifndef MAIN_H
#define MAIN_H

#include <vector>
#include <string>
#include <map>
#include <utility>
#include "gmp.h"
#include <iostream>
#include<stdlib.h>
#include<string.h>
#include<sstream>
#include <malloc.h>
using namespace std;

template <class public_parameter>
class rsa_vector_commitment
{
    public:

    virtual public_parameter rsa_gen(int para,mpz_t p,mpz_t q,mpz_t w[],mpz_t N,mpz_t a,mpz_t S[])=0;
    virtual string rsa_com(int para, vector<string> messagesino,public_parameter p)=0;

    virtual string rsa_open(int para,vector<string> messagesino,int index,public_parameter pp)=0;

    virtual bool rsa_verify(string Commitment,string message,string Proof,int index,public_parameter pp)=0;

    virtual string rsa_update(string Commitment,string message,string newmessage,int index,public_parameter pp)=0;

    virtual string rsa_proof_update(int para,string Proof,string newmsg,
                          int i,int j,string msg,public_parameter pp)=0;
};

/**
* Vector Commitments based on the RSA assumption.
*/
class RSA_public_parameter

{
    public:


    RSA_public_parameter( int param,string a_param,string n_param, vector<string>primes,vector<string> S)
    {

        rsa_p=param;
        nins = n_param;//N
        ains=a_param;//a
        for (int i=0 ; i <param;i++)
            Ss.push_back(S.at(i));
        for (int i=0 ; i <param;i++)
                    primess.push_back(primes.at(i));

    };
    unsigned int rsa_p;
    string nins;
    string ains;
    vector <string> Ss;
    vector <string> primess;
};


class RSA_vector_commitment : public rsa_vector_commitment <RSA_public_parameter>
{
    public:
    RSA_public_parameter rsa_gen(int para,mpz_t p,mpz_t q,mpz_t w[],mpz_t N,mpz_t a,mpz_t S[]);
    string rsa_com(int para, vector<string> messagesino,  RSA_public_parameter pp);
    string rsa_open(int para,vector<string> messagesino,int index,RSA_public_parameter pp);
    bool rsa_verify(string Commitment,string message,string Proof,int index,RSA_public_parameter pp);
    string rsa_update(string Commitment,string message,string newmessage,int index,RSA_public_parameter pp);
    string rsa_proof_update(int para,string Proof,string newmsg,
                          int i,int j,string msg,RSA_public_parameter pp);
};
#endif
