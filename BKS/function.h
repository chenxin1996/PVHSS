#include "HSS.h"
#include "GenData.h"
#include <algorithm>
#include <cstdlib>
using namespace std;

void Eval_Poly(int deg, int num_data) {
    //para
    PKE_Para pkePara;
    vec_ZZ_pX pkePk, pkeSk, hssEk_1, hssEk_2;
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    hssEk_1.SetLength(2);
    hssEk_2.SetLength(2);
    pkePara.msg_bit = 32;
    pkePara.d = deg;
    pkePara.num_data = num_data;
    PKE_Gen(pkePara, pkePk, pkeSk,1);
    ZZ_pXModulus modulus(pkePara.xN);
    HSS_Gen(hssEk_1, hssEk_2, pkePara, pkeSk);
    //data
    Data data;
    GenData(data, pkePara, pkePk);
    ZZ_pX y, y1, y2;

    cout << "******************** Server 1 Evaluating ********************" << "\n";
    double time = GetTime();
    HSS_Eval(y1, 1, pkePara, modulus, hssEk_1, data.C_X, data.PRF);
    time = GetTime() - time;
    cout << "running time Server 1: " << time * 1000 << "ms\n";

    cout << "******************** Server 2 Evaluating ********************" << "\n";
    time = GetTime();
    HSS_Eval(y2, 2, pkePara, modulus, hssEk_2, data.C_X, data.PRF);
    time = GetTime() - time;
    cout << "running time Server 2: " << time * 1000 << "ms\n";

    cout << "******************** Client Decrypting ********************" << "\n";
    time = GetTime();
    y = y1 + y2;
    cout << "The value of y= " << y << "\n";
    time = GetTime() - time;
    cout << "running time Client: " << time * 1000 << "ms\n";

    ZZ y_ZZ;
    NativeEval(y_ZZ, pkePara.d, pkePara.num_data, data.X);
    cout << "native comnputing y=" << y_ZZ << "\n";
}

void Run_Gen(int msg_bit, int num_data) {
    PKE_Para pkePara;
    vec_ZZ_pX pkePk, pkeSk, hssEk_1, hssEk_2;
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    hssEk_1.SetLength(2);
    hssEk_2.SetLength(2);
    pkePara.msg_bit = msg_bit;
    pkePara.num_data = 2;
    PKE_Gen(pkePara, pkePk, pkeSk,0);
    ZZ_pXModulus modulus(pkePara.xN);
    HSS_Gen(hssEk_1, hssEk_2, pkePara, pkeSk);
}

void Time_Gen(int msg_bit, int cyctimes) {
    auto *Time = new double[cyctimes];
    double time,mean,stdev;
    for (int i = 0; i < cyctimes; i++) {
        time = GetTime();
        Run_Gen(msg_bit, 2);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean,stdev,Time,cyctimes);
    cout << "gen algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";
}


void Time_Enc(int msg_bit,  int cyctimes) {
    //para
    PKE_Para pkePara;
    vec_ZZ_pX pkePk, pkeSk, hssEk_1, hssEk_2;
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    hssEk_1.SetLength(2);
    hssEk_2.SetLength(2);
    pkePara.msg_bit = msg_bit;
    pkePara.num_data = 2;
    PKE_Gen(pkePara, pkePk, pkeSk,0);
    ZZ_pXModulus modulus(pkePara.xN);
    HSS_Gen(hssEk_1, hssEk_2, pkePara, pkeSk);

    vec_ZZ_pX C;
    ZZ x;
    C.SetLength(4);
    auto *Time = new double[cyctimes];
    double time,mean,stdev;
    for (int i = 0; i < cyctimes; i++) {
        x = RandomBits_ZZ(pkePara.msg_bit);
        time = GetTime();
        HSS_Enc(C, pkePara, modulus, pkePk, x);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean,stdev,Time,cyctimes);
    cout << "Enc algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";
}

void Time_Dec(int msg_bit, int cyctimes)
{
    //para
    PKE_Para pkePara;
    vec_ZZ_pX pkePk, pkeSk, hssEk_1, hssEk_2;
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    hssEk_1.SetLength(2);
    hssEk_2.SetLength(2);
    pkePara.msg_bit = msg_bit;
    pkePara.num_data = 2;
    PKE_Gen(pkePara, pkePk, pkeSk,0);
    ZZ_pXModulus modulus(pkePara.xN);
    HSS_Gen(hssEk_1, hssEk_2, pkePara, pkeSk);

    ZZ_pX y_ZZ_pZ,y1,y2;
    ZZ_p y_ZZ_p;
    ZZ y_ZZ;
    auto *Time = new double[cyctimes];
    double time,mean,stdev;
    for (int i = 0; i < cyctimes; i++) {
        RandomBits(y_ZZ,msg_bit);
        conv(y_ZZ_p,y_ZZ);
        SetCoeff(y_ZZ_pZ,0,y_ZZ_p);
        Random_ZZ_pX(y1,pkePara.N,pkePara.q_bit);
        y2=y_ZZ_pZ-y1;
        time = GetTime();
        y1=y1+y2;
        //cout << "y: " << y1[0] << "\n";
        Time[i] = GetTime() - time;
    }
    DataProcess(mean,stdev,Time,cyctimes);
    cout << "Dec algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";
}

void Time_Eval_Subalgo(int msg_bit,  int cyctimes) {
    //para
    PKE_Para pkePara;
    vec_ZZ_pX pkePk, pkeSk, hssEk_1, hssEk_2;
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    hssEk_1.SetLength(2);
    hssEk_2.SetLength(2);
    pkePara.msg_bit = msg_bit;
    pkePara.num_data = 2;
    PKE_Gen(pkePara, pkePk, pkeSk,0);
    ZZ_pXModulus modulus(pkePara.xN);
    HSS_Gen(hssEk_1, hssEk_2, pkePara, pkeSk);
    //data
    Data data;
    GenData(data, pkePara, pkePk);

    //load
    vec_ZZ_pX db1,db2;
    db1.SetLength(2);
    db2.SetLength(2);
    auto *Time = new double[cyctimes];
    double time,mean,stdev;
    for (int i = 0; i < cyctimes; i++) {
        int index=rand()%10;
        time = GetTime();
        HSS_Mult(db1,pkePara,modulus,hssEk_1,data.C_X[rand()%pkePara.num_data]);
        db1[0]=db1[0]+data.PRF[index][0];
        db1[1]=db1[1]+data.PRF[index][1];
        Time[i] = GetTime() - time;
    }
    DataProcess(mean,stdev,Time,cyctimes);
    cout << "Load algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";


    //add1
    for (int i = 0; i < cyctimes; i++) {
        Random_ZZ_pX(db1[0],pkePara.N,pkePara.q_bit);
        Random_ZZ_pX(db1[1],pkePara.N,pkePara.q_bit);
        Random_ZZ_pX(db2[0],pkePara.N,pkePara.q_bit);
        Random_ZZ_pX(db2[1],pkePara.N,pkePara.q_bit);
        int index=rand()%10;
        time = GetTime();
        db1[0]=db1[0]+db2[0]+data.PRF[index][0];
        db1[1]=db1[1]+db2[1]+data.PRF[index][1];
        Time[i] = GetTime() - time;
    }
    DataProcess(mean,stdev,Time,cyctimes);
    cout << "Add1 algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";

    //add2
    vec_ZZ_pX C;
    C.SetLength(4);
    for (int i = 0; i < cyctimes; i++) {
        time = GetTime();
        C[0]=data.C_X[0][0]+data.C_X[1][0];
        C[1]=data.C_X[0][1]+data.C_X[1][1];
        C[2]=data.C_X[0][2]+data.C_X[1][2];
        C[3]=data.C_X[0][3]+data.C_X[1][3];
        Time[i] = GetTime() - time;
    }
    DataProcess(mean,stdev,Time,cyctimes);
    cout << "Add2 algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";

    //mult
    for (int i = 0; i < cyctimes; i++) {
        HSS_Mult(db1,pkePara,modulus,hssEk_1,data.C_X[rand()%pkePara.num_data]);
        int index=rand()%10;
        time = GetTime();
        HSS_Mult(db1,pkePara,modulus,db1,data.C_X[rand()%pkePara.num_data]);
        db1[0]=db1[0]+data.PRF[index][0];
        db1[1]=db1[1]+data.PRF[index][1];
        Time[i] = GetTime() - time;
    }
    DataProcess(mean,stdev,Time,cyctimes);
    cout << "Mult algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";


    //out
    ZZ_pX yb_ZZ_pX;
    ZZX yb_ZZX;
    ZZ coeff;
    ZZ msg_size= power_ZZ(2,msg_bit);
    for (int i = 0; i < cyctimes; i++) {
        Random_ZZ_pX(yb_ZZ_pX,pkePara.N,pkePara.q_bit);
        time = GetTime();
        conv(yb_ZZX,yb_ZZ_pX);
        for(int j=0;j<pkePara.N;j++)
        {
            GetCoeff(coeff,yb_ZZX,i);
            coeff=coeff%msg_size;
            SetCoeff(yb_ZZX,i,coeff);
        }
        Time[i] = GetTime() - time;
    }
    DataProcess(mean,stdev,Time,cyctimes);
    cout << "Output algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";

}

