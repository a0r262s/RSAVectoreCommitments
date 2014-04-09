#include "main.h"
#include <stdio.h>
#include "gmp.h"
#include <iostream>
#include <utility>
#include <map>
#include<stdlib.h>
#include<string.h>
#include <sstream>
#include <malloc.h>
#include <unistd.h>
#include <math.h>
#include <gmp.h>
#include <gmpxx.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <vector>
#include <algorithm>
using namespace std;
#define MODULUS_SIZE 1024                    /* This is the number of bits we want in the modulus */
#define BUFFER_SIZE_P ((MODULUS_SIZE/8)/16) /* This is the number of *BYTES* in primes p and q */
#define BUFFER_SIZE_E ((MODULUS_SIZE/8)/8) /* This is the number of *BYTES* in e_i to e_j */

RSA_public_parameter RSA_vector_commitment:: rsa_gen(int para,mpz_t p,mpz_t q,mpz_t w[],mpz_t N,mpz_t a,mpz_t S[]){
    char buf_p[BUFFER_SIZE_P];
    char buf_e[BUFFER_SIZE_E];
    int i;
    mpz_t tmp1; mpz_init(tmp1);
    mpz_t tmp2; mpz_init(tmp2);
    mpz_t tmp_array[para];
    for(int i= 0; i< para; i++)
        mpz_init( tmp_array[i]);
    srand(time(NULL));
    /* Select p and q */
    /* Start with p */
    // Set the bits of tmp randomly
    for(int i = 0; i < BUFFER_SIZE_P; i++)
        buf_p[i] = rand() % 0xFF;
    // Set the top two bits to 1 to ensure int(tmp) is relatively large
    buf_p[0] |= 0xC0;
    // Set the bottom bit to 1 to ensure int(tmp) is odd (better for finding primes)
    buf_p[BUFFER_SIZE_P - 1] |= 0x01;
    // Interpret this char buffer as an int
    mpz_import(tmp1, BUFFER_SIZE_P, 1, sizeof(buf_p[0]), 0, 0, buf_p);
    // Pick the next prime starting from that random number
    mpz_nextprime(p, tmp1);
    /*Select prime q*/
    for(int i = 0; i < BUFFER_SIZE_P; i++)
        buf_p[i] = rand() % 0xFF;
    // Set the top two bits to 1 to ensure int(tmp) is relatively large
    buf_p[0] |= 0xC0;
    // Set the bottom bit to 1 to ensure int(tmp) is odd
    buf_p[BUFFER_SIZE_P - 1] |= 0x01;
    // Interpret this char buffer as an int
    mpz_import(tmp1, (BUFFER_SIZE_P), 1, sizeof(buf_p[0]), 0, 0, buf_p);
    mpz_nextprime(q, tmp1);
    /*Selecting primes e_1 to e_q*/
    for(int j= 0; j< para; j++){
        for(i = 0; i < BUFFER_SIZE_E; i++){
            buf_e[i] = rand() % 0xFF;
            // Set the top two bits to 1 to ensure int(tmp) is relatively large
            buf_e[0] |= 0xC0;
            // Set the bottom bit to 1 to ensure int(tmp) is odd
            buf_e[BUFFER_SIZE_E - 1] |= 0x01;
            // Interpret this char buffer as an int
            mpz_import(tmp_array[j], (BUFFER_SIZE_E), 1, sizeof(buf_e[0]), 0, 0, buf_e);
        }
    }
    mpz_mul(N,q,p);//N=p*q
    mpz_sub_ui(tmp1, p, 1);//p-1
    mpz_sub_ui(tmp2, q, 1);//q-1
    mpz_mul(tmp1,tmp1,tmp2);//Phi(N)=(p-1)(q-1)
    mpz_t(res1);mpz_init(res1);
    mpz_t(res2);mpz_init(res2);
    MP_INT ONE;
    mpz_init_set_str (&ONE, "1", 10);
    for(int i= 0; i< para; i++){
        mpz_nextprime(w[i],tmp_array[i]);
        if(mpz_cmp(res1,&ONE)!=0){  //if res1!=&ONE
            mpz_nextprime(w[i],tmp_array[i]);
            mpz_gcd(res1,tmp1,w[i]);
        }
        if(mpz_cmp(res2,&ONE)!=0){//if res1!=&ONE
            mpz_nextprime(w[i],tmp_array[i]);
            mpz_gcd(res2,w[i],N);
        }
    }
    /*generate random a which is blongs to Z_N */
    gmp_randstate_t rs;
    //gmp_randinit_mt(rs);
    gmp_randinit_default(rs);
    gmp_randseed_ui(rs, 1233367890);
    mpz_urandomm(a,rs,N);
    gmp_randclear(rs);
    /*Computing S_i*/
    mpz_t ex[para];
    for(i= 0; i< para; i++)
        //mpz_init_set_ui( ex[i],usi);
        mpz_init_set_str (ex[i], "1", 10);
    for (int i =0; i<para;i++){
        for (int j =0;j<para;j++){
            if (j != i){
                mpz_mul(ex[i],ex[i],w[j]);
            }
        }

        mpz_powm(S[i],a,ex[i],N);//pow mod N. cf. sec. 3.2
    }
    mpz_class a_class;
    a_class.get_mpz_t();
    mpz_set(a_class.get_mpz_t(),a);
    string ains;
    ains =a_class.get_str();
    cout << ains << "<-- a\n";

    mpz_class N_class;
    N_class.get_mpz_t();
    mpz_set(N_class.get_mpz_t(),N);
    string Nins;
    Nins =N_class.get_str();
    cout << Nins << " <--N\n";
    //
    vector<string> Ss;
    //mpz_class S_class[para];
    mpz_class *S_class = new mpz_class[para];
    for(int i= 0; i< para; i++){
        S_class[i].get_mpz_t();
        mpz_set(S_class[i].get_mpz_t(),S[i]);
        Ss.push_back(S_class[i].get_str());
        cout << Ss[i] << " <-- 'S'\n";
    }

    vector <string> primess;
    //mpz_class w_class[para];
    mpz_class *w_class = new mpz_class[para];

    for(int i= 0; i< para; i++){
        w_class[i].get_mpz_t();
        mpz_set(w_class[i].get_mpz_t(),w[i]);
        primess.push_back( w_class[i].get_str());
        cout << w[i] << " <--Primes\n";
    }
    cout <<  "\n";
    //clears
    mpz_clear(tmp1);
    mpz_clear(tmp2);

    mpz_clear(res1);
    mpz_clear(res2);
    for(int i= 0; i< para; i++){
        mpz_clear(tmp_array[i]);
        mpz_clear(ex[i]);
    }
    RSA_public_parameter pp(para,ains,Nins,primess,Ss);
    return pp;
}

string  RSA_vector_commitment:: rsa_com(int para, vector<string> messagesin,RSA_public_parameter pp){
    //initializations
    mpz_t N;
    mpz_set_str(N,pp.nins.c_str(),10);
    mpz_class C;
    mpz_init_set_str (C.get_mpz_t(), "1", 10);
    mpz_t S[para];
    for(int i= 0; i< para; i++)
        mpz_init( S[i]);
    for(int i= 0; i< para; i++)
        mpz_set_str(S[i],pp.Ss.at(i).c_str(),10);
    mpz_t tmp[para];
    for(int i= 0; i< para; i++)
        mpz_init(tmp[i]);
    mpz_t(com_int);mpz_init(com_int);
    //computations
    for(int i = 0; i <para; i++){
        mpz_set_str(com_int,messagesin.at(i).c_str(),2);
        //mpz_set_str(S[i],pp.Ss.at(i).c_str(),10);
        mpz_powm(tmp[i],S[i],com_int,N);//pow mod N. cf. sec. 3.2
        mpz_mul(C.get_mpz_t(),tmp[i],C.get_mpz_t());
    }
    mpz_mod(C.get_mpz_t(),C.get_mpz_t(),N);//ATT
    cout << C.get_mpz_t() << " <--Commitment from rsa_com\n";
    //clears
    for(int i= 0; i< para; i++)
        mpz_clear(tmp[i]);
    mpz_clear(com_int);
    return C.get_str();
}

string RSA_vector_commitment:: rsa_open(int para,vector<string> messages,int index,RSA_public_parameter pp){
    mpz_t w[para];
    for(int i= 0; i< para; i++)
        mpz_init( w[i]);
    mpz_t a;mpz_init(a);
    mpz_t N;mpz_init(N);
    mpz_class proof;
    mpz_t com_int;mpz_init(com_int);
    mpz_t tmpz; mpz_init(tmpz);
    mpz_t tmp[para];
    for(int i= 0; i< para; i++)
        mpz_init_set_str (tmp[i], "1", 10);

    for (int i=0;i<para;i++){
        if(i!=index){
            for (int j=0;j<para;j++){
                if(j!=i){
                    if(j!=index){
                        mpz_set_str(w[j],pp.primess.at(j).c_str(),10);
                        mpz_mul(tmp[i],w[j],tmp[i]);

                    }
                }
            }
            mpz_set_str(com_int,messages.at(i).c_str(),2);
            mpz_mul(tmp[i],tmp[i],com_int);
            mpz_add(tmpz,tmpz,tmp[i]);
        }

    }
    mpz_set_str(a,pp.ains.c_str(),10);
    mpz_set_str(N,pp.nins.c_str(),10);
    mpz_powm(tmpz,a,tmpz,N);
    mpz_set(proof.get_mpz_t(),tmpz);
    cout << proof.get_str()<< " <--Proof from rsa_open\n";
    return proof.get_str();
}

bool RSA_vector_commitment::rsa_verify(string Commitment,string message,string Proof,int index,RSA_public_parameter pp){
    //initializations
    bool b=true;
    mpz_t stcommitment; mpz_init(stcommitment);
    mpz_set_str(stcommitment,Commitment.c_str(),10);
    mpz_t stproof;
    mpz_set_str(stproof,Proof.c_str(),10);
    cout << stproof<< " <--Proof in rsa_verify\n";
    mpz_t com_int;mpz_init(com_int);
    mpz_t(tmp1);mpz_init(tmp1);
    mpz_t(tmp2);mpz_init(tmp2);
    mpz_set_str(com_int,message.c_str(),2);
    //computations
    mpz_t S; mpz_init(S);
    mpz_set_str(S,pp.Ss.at(index).c_str(),10);
    mpz_t w; mpz_init(w);
    mpz_set_str(w,pp.primess.at(index).c_str(),10);
    mpz_t N; mpz_init(N);
    mpz_set_str(N,pp.nins.c_str(),10);
    mpz_powm(tmp1,S,com_int,N);
    mpz_powm(tmp2,stproof,w,N);
    mpz_mul(tmp1,tmp2,tmp1);
    mpz_mod(tmp1,tmp1,N);//ATT
    cout <<tmp1<< " <--Right side from rsa_verify\n";
    cout <<stcommitment<< " <--Left side from rsa_verify\n";
    //return
    if(mpz_cmp(stcommitment,tmp1)==0)//if C==tmp1
        b=true;
    else
        b=false;

    printf("\nBoolaen value of rsa_verify function:%s\n",(b)?"true":"false");
    cout<<"\n"<<"";
    //clears
    mpz_clear(tmp2);
    mpz_clear(tmp1);
    mpz_clear(com_int);
    return true;
}

string RSA_vector_commitment::rsa_update(string Commitment,string message,string newmessage,int index,RSA_public_parameter pp){
    mpz_class outupdate;
    mpz_t stcommitment;mpz_init(stcommitment);
    mpz_set_str(stcommitment,Commitment.c_str(),10);
    //initializations
    mpz_t com_int;mpz_init(com_int);
    mpz_t msg_int;mpz_init(msg_int);
    mpz_t tmp1;mpz_init(tmp1);
    //computations
    mpz_set_str(com_int,newmessage.c_str(),2);
    mpz_set_str(msg_int,message.c_str(),2);
    mpz_sub(msg_int,msg_int,com_int);
    mpz_t S; mpz_init(S);
    mpz_set_str(S,pp.Ss.at(index).c_str(),10);
    mpz_t N; mpz_init(N);
    mpz_set_str(N,pp.nins.c_str(),10);
    mpz_powm(tmp1,S,msg_int,N);
    mpz_mul(outupdate.get_mpz_t(),stcommitment,tmp1);
    mpz_mod(outupdate.get_mpz_t(),outupdate.get_mpz_t(),N);//ATT
    cout << outupdate.get_mpz_t()<< " <--Commitment from rsa_update\n";
    //clears
    mpz_clear(tmp1);
    mpz_clear(msg_int);
    mpz_clear(com_int);
    return outupdate.get_str();
}

string RSA_vector_commitment:: rsa_proof_update(int para,string Proof,string newmsg,
                                                int i,int j,string msg,RSA_public_parameter pp){

    mpz_class proofout;

    mpz_t N; mpz_init(N);
    mpz_set_str(N,pp.nins.c_str(),10);
    mpz_t a; mpz_init(a);
    mpz_set_str(a,pp.ains.c_str(),10);

    mpz_t w[para];
    for(int m= 0; m< para; m++)
        mpz_init( w[m]);
    for(int m= 0; m< para; m++)
        mpz_set_str(w[m],pp.primess.at(m).c_str(),10);
    mpz_t stproof;mpz_init(stproof);
    mpz_set_str(stproof,Proof.c_str(),10);
    //initializations
    mpz_t com_int;mpz_init(com_int);
    mpz_t msg_int;mpz_init(msg_int);
    mpz_t tmp1;mpz_init(tmp1);
    //computations: updated commitment
    mpz_set_str(com_int,msg.c_str(),2);
    mpz_set_str(msg_int,newmsg.c_str(),2);
    mpz_sub(msg_int,msg_int,com_int);
    mpz_t(tmp2);
    mpz_init_set_str (tmp2, "1", 10);
    if(i!=j){

        for(int k=0;k<para;k++){
            if(k!=i){
                if(k!=j){
                    mpz_mul(tmp2,w[k],tmp2);
                }
            }
        }

        mpz_mul(tmp2,tmp2,msg_int);

        mpz_powm(tmp2,a,tmp2,N);

        mpz_mul(tmp2,stproof,tmp2);
        mpz_mod(tmp2,tmp2,N);
        mpz_set(proofout.get_mpz_t(),tmp2);



    }
    if (i==j){
        mpz_set(proofout.get_mpz_t(),stproof);
    }
    //clears
    mpz_clear(tmp2);
    mpz_clear(tmp1);
    mpz_clear(com_int);
    mpz_clear(msg_int);
    cout<<proofout.get_str()<<"<-- proof from rsa_proof_update\n";
    return proofout.get_str();
}


int main(){
    //inputs
    vector<string> messagesin;//Not Alphabetical

    messagesin.push_back    ("1111010010001010110111111111110");
    messagesin.push_back    ("1111011011001010110111111110000");
    messagesin.push_back    ("1111010111101010110111111010000");

    messagesin.push_back    ("1111011110001011111111111010010");

    messagesin.push_back    ("1000110010001010110111111111111");
    messagesin.push_back    ("1111010010001010110110100010100");

    messagesin.push_back    ("1111111111111010110110011011111");
    messagesin.push_back    ("1101011111111111110110001110000");
    messagesin.push_back    ("1111111111111010110110001010000");

    messagesin.push_back    ("1100010111111111111110011010011");

    string message_verify=   "1111010010001010110111111111110"; // equal to messagein[0], argument of rsa_verify
    string message_update=   "1000110010001010110111000010011"; // is used for rsa_updat
    //string message_verify= "1111110010001010110110000010000";//another message not equal to [i] in vector

    int parameter=10;
    int indexin=1;
    int indexinin=2;
    //init
    int count=1;
    int i;
    mpz_t *qin;
    mpz_t *pin;
    mpz_t ain, Nin;
    mpz_init(ain);
    mpz_init(Nin);
    mpz_t Cprimein;
    mpz_init_set_str (Cprimein, "1", 10);
    mpz_t proofin; mpz_init(proofin);
    mpz_t proofpuin; mpz_init(proofpuin);
    mpz_t Cpuin; mpz_init(Cpuin);
    pin = (mpz_t *)malloc(count*sizeof(mpz_t));
    qin = (mpz_t *)malloc(count*sizeof(mpz_t));
    if(pin == NULL) printf("Failed to allocate memory\n");
    if(qin == NULL) printf("Failed to allocate memory\n");
    mpz_t win[parameter];//q prime numbers q is here 10
    for(i= 0; i< parameter; i++)
        mpz_init( win[i]);
    mpz_t Sin[parameter];
    for(i= 0; i< parameter; i++)
        mpz_init( Sin[i]);
    //functions
    RSA_vector_commitment vc;
    RSA_public_parameter pp = vc.rsa_gen(parameter,*pin,*qin,win,Nin,ain,Sin);

    string commitment =vc.rsa_com(parameter,messagesin,pp);
    string proof=vc.rsa_open(parameter,messagesin,indexin,pp);
    vc.rsa_verify(commitment,messagesin.at(indexin),proof,indexin,pp);
    string update_comm =vc.rsa_update(commitment,message_verify,message_update,indexinin,pp);
    string newproof = vc.rsa_proof_update(parameter,proof,message_verify,indexinin,indexin,message_update,pp);
    vc.rsa_verify(update_comm, messagesin.at(indexin), newproof,indexin,  pp);//proved

    return 0;
}
