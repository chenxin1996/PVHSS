#include "VHSS.h"
#include "GenData.h"
#include <algorithm>
#include <cstdlib>
using namespace std;

void Eval_Poly(int deg, int num_data) {
    //para
    PKE_Para pkePara;
    VHSS_Para vhssPara;
    vec_ZZ_pX pkePk, pkeSk;
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    pkePara.msg_bit = 32;
    pkePara.d = deg;
    pkePara.num_data = num_data;
    PKE_Gen(pkePara, pkePk, pkeSk,1);
    ZZ_pXModulus modulus(pkePara.xN);
    VHSS_Gen(vhssPara,pkePara, modulus,pkeSk);
    //data
    Data data;
    GenData(data, pkePara, pkePk);
    ZZ_pX y, tau, y1, y2,tau1,tau2,y_alpha;


    cout << "******************** Server 1 Evaluating ********************" << "\n";
    double time = GetTime();
    VHSS_Eval(y1, 1, pkePara, modulus, vhssPara.vhssEk_1, data.C_X, data.PRF);
    VHSS_Eval(tau1, 1, pkePara, modulus, vhssPara.vhssEk_3, data.C_X, data.PRF);
    time = GetTime() - time;
    cout << "running time Server 1: " << time * 1000 << "ms\n";

    cout << "******************** Server 2 Evaluating ********************" << "\n";
    time = GetTime();
    VHSS_Eval(y2, 2, pkePara, modulus, vhssPara.vhssEk_2, data.C_X, data.PRF);
    VHSS_Eval(tau2, 2, pkePara, modulus, vhssPara.vhssEk_4, data.C_X, data.PRF);
    time = GetTime() - time;
    cout << "running time Server 2: " << time * 1000 << "ms\n";

    cout << "******************** Client Verifying ********************" << "\n";
    time = GetTime();
    y = y1 + y2;
    tau=tau1+tau2;
    y_alpha=y*vhssPara.alpha;
    if(y_alpha==tau)
    {
        cout << "Verification passed!\nThe value of y= " << y[0] << "\n";
    } else {
        printf("ERROR\n");
    }
    time = GetTime() - time;
    cout << "running time Client: " << time * 1000 << "ms\n";

    ZZ y_ZZ;
    NativeEval(y_ZZ,pkePara.d,pkePara.num_data,data.X);
    cout << "native comnputing y=" << y_ZZ << "\n";
}

void Run_Gen(int msg_bit) {
    //para
    PKE_Para pkePara;
    VHSS_Para vhssPara;
    vec_ZZ_pX pkePk, pkeSk;
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    pkePara.msg_bit = msg_bit;
    pkePara.num_data = 2;
    PKE_Gen(pkePara, pkePk, pkeSk,0);
    ZZ_pXModulus modulus(pkePara.xN);
    VHSS_Gen(vhssPara,pkePara, modulus,pkeSk);
}

void Time_Gen(int msg_bit, int cyctimes) {
    auto *Time = new double[cyctimes];
    double time,mean,stdev;
    for (int i = 0; i < cyctimes; i++) {
        time = GetTime();
        Run_Gen(msg_bit);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean,stdev,Time,cyctimes);
    cout << "Gen algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";
}


void Time_Enc(int msg_bit, int cyctimes) {
    //para
    PKE_Para pkePara;
    VHSS_Para vhssPara;
    vec_ZZ_pX pkePk, pkeSk;
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    pkePara.msg_bit = msg_bit;
    pkePara.num_data = 2;
    PKE_Gen(pkePara, pkePk, pkeSk,0);
    ZZ_pXModulus modulus(pkePara.xN);
    VHSS_Gen(vhssPara,pkePara, modulus,pkeSk);

    vec_ZZ_pX C;
    ZZ x;
    C.SetLength(4);
    auto *Time = new double[cyctimes];
    double time,mean,stdev;
    for (int i = 0; i < cyctimes; i++) {
        x = RandomBits_ZZ(pkePara.msg_bit);
        time = GetTime();
        VHSS_Enc(C, pkePara, modulus, pkePk, x);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean,stdev,Time,cyctimes);
    cout << "Enc algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";
}

void Time_Ver(int msg_bit, int cyctimes)
{
    //para
    PKE_Para pkePara;
    VHSS_Para vhssPara;
    vec_ZZ_pX pkePk, pkeSk;
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    pkePara.msg_bit = msg_bit;
    pkePara.num_data = 2;
    PKE_Gen(pkePara, pkePk, pkeSk,0);
    ZZ_pXModulus modulus(pkePara.xN);
    VHSS_Gen(vhssPara,pkePara, modulus,pkeSk);

    ZZ_pX y1,y2,tau1,tau2,alpha_y,y_ZZ_pX,tau;
    ZZ y;
    ZZ_p y_ZZ_p;

    auto *Time = new double[cyctimes];
    double time,mean,stdev;
    for (int i = 0; i < cyctimes; i++) {

        RandomBits(y,640);
        conv(y_ZZ_p,y);
        SetCoeff(y_ZZ_pX,0,y_ZZ_p);
        alpha_y=y_ZZ_pX*vhssPara.alpha;
        Random_ZZ_pX(y1,pkePara.N,pkePara.q_bit);
        Random_ZZ_pX(tau1,pkePara.N,pkePara.q_bit);
        y2=y_ZZ_pX-y1;
        tau2=alpha_y-tau1;

        time = GetTime();
        y_ZZ_pX=y1+y2;
        tau=tau1+tau2;
        alpha_y=vhssPara.alpha*y_ZZ_pX;
        if(alpha_y==tau)
        {
            //cout << "Verification passed!\nThe value of function: " << y1[0] << "\n";
        } else {
            //cout << "ERROR << "\n";
        }
        Time[i] = GetTime() - time;
    }
    DataProcess(mean,stdev,Time,cyctimes);
    cout << "Ver algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";
}

void Time_Eval_Subalgo(int msg_bit, int cyctimes) {
    //para
    PKE_Para pkePara;
    VHSS_Para vhssPara;
    vec_ZZ_pX pkePk, pkeSk;
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    pkePara.msg_bit = msg_bit;
    pkePara.num_data = 2;
    PKE_Gen(pkePara, pkePk, pkeSk,0);
    ZZ_pXModulus modulus(pkePara.xN);
    VHSS_Gen(vhssPara,pkePara, modulus,pkeSk);
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
        VHSS_Mult(db1,pkePara,modulus,vhssPara.vhssEk_1,data.C_X[rand()%pkePara.num_data]);
        VHSS_Mult(db2,pkePara,modulus,vhssPara.vhssEk_3,data.C_X[rand()%pkePara.num_data]);
        db1[0]=db1[0]+data.PRF[index][0];
        db1[1]=db1[1]+data.PRF[index][1];
        db2[0]=db2[0]+data.PRF[index][0];
        db2[1]=db2[1]+data.PRF[index][1];
        Time[i] = GetTime() - time;
    }
    DataProcess(mean,stdev,Time,cyctimes);
    cout << "Load algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";

    //add1
    vec_ZZ_pX db3,db4;
    db3.SetLength(2);
    db4.SetLength(2);
    for (int i = 0; i < cyctimes; i++) {
        Random_ZZ_pX(db1[0],pkePara.N,pkePara.q_bit);
        Random_ZZ_pX(db1[1],pkePara.N,pkePara.q_bit);
        Random_ZZ_pX(db2[0],pkePara.N,pkePara.q_bit);
        Random_ZZ_pX(db2[1],pkePara.N,pkePara.q_bit);
        Random_ZZ_pX(db3[0],pkePara.N,pkePara.q_bit);
        Random_ZZ_pX(db3[1],pkePara.N,pkePara.q_bit);
        Random_ZZ_pX(db4[0],pkePara.N,pkePara.q_bit);
        Random_ZZ_pX(db4[1],pkePara.N,pkePara.q_bit);
        int index=rand()%10;
        time = GetTime();
        db1[0]=db1[0]+db2[0]+data.PRF[index][0];
        db1[1]=db1[1]+db2[1]+data.PRF[index][1];
        db3[0]=db3[0]+db4[0]+data.PRF[index][0];
        db3[1]=db3[1]+db4[1]+data.PRF[index][1];
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

    //cmult
    ZZ constant= RandomBits_ZZ(msg_bit);
    for (int i = 0; i < cyctimes; i++) {

        VHSS_Mult(db1,pkePara,modulus,vhssPara.vhssEk_1,data.C_X[rand()%pkePara.num_data]);
        VHSS_Mult(db2,pkePara,modulus,vhssPara.vhssEk_3,data.C_X[rand()%pkePara.num_data]);
        int index=rand()%10;
        time = GetTime();
        ZZ_pX_ScaleMul_ZZ(db1[0],db1[0],constant);
        ZZ_pX_ScaleMul_ZZ(db1[1],db1[1],constant);
        VHSS_Mult(db1,pkePara,modulus,db1,data.C_X[rand()%pkePara.num_data]);

        ZZ_pX_ScaleMul_ZZ(db2[0],db2[0],constant);
        ZZ_pX_ScaleMul_ZZ(db2[1],db2[1],constant);
        VHSS_Mult(db2,pkePara,modulus,db2,data.C_X[rand()%pkePara.num_data]);
        db1[0]=db1[0]+data.PRF[index][0];
        db1[1]=db1[1]+data.PRF[index][1];
        db2[0]=db2[0]+data.PRF[index][0];
        db2[1]=db2[1]+data.PRF[index][1];
        Time[i] = GetTime() - time;
    }
    DataProcess(mean,stdev,Time,cyctimes);
    cout << "cMult algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";


    //mult
    for (int i = 0; i < cyctimes; i++) {
        VHSS_Mult(db1,pkePara,modulus,vhssPara.vhssEk_1,data.C_X[rand()%pkePara.num_data]);
        VHSS_Mult(db2,pkePara,modulus,vhssPara.vhssEk_3,data.C_X[rand()%pkePara.num_data]);
        int index=rand()%10;
        time = GetTime();
        VHSS_Mult(db1,pkePara,modulus,db1,data.C_X[rand()%pkePara.num_data]);
        VHSS_Mult(db2,pkePara,modulus,db2,data.C_X[rand()%pkePara.num_data]);
        db1[0]=db1[0]+data.PRF[index][0];
        db1[1]=db1[1]+data.PRF[index][1];
        db2[0]=db2[0]+data.PRF[index][0];
        db2[1]=db2[1]+data.PRF[index][1];
        Time[i] = GetTime() - time;
    }
    DataProcess(mean,stdev,Time,cyctimes);
    cout << "Mult algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";


    //out
    ZZ_pX yb_ZZ_pX,taub_ZZ_pX;
    ZZX yb_ZZX,taub_ZZX;
    ZZ coeff;
    ZZ msg_size= power_ZZ(2,msg_bit);
    for (int i = 0; i < cyctimes; i++) {
        Random_ZZ_pX(yb_ZZ_pX,pkePara.N,pkePara.q_bit);
        Random_ZZ_pX(taub_ZZ_pX,pkePara.N,pkePara.q_bit);
        time = GetTime();
        conv(yb_ZZX,yb_ZZ_pX);
        conv(taub_ZZX,taub_ZZ_pX);
        for(int j=0;j<pkePara.N;j++)
        {
            GetCoeff(coeff,yb_ZZX,i);
            coeff=coeff%msg_size;
            SetCoeff(yb_ZZX,i,coeff);

            GetCoeff(coeff,taub_ZZX,i);
            coeff=coeff%msg_size;
            SetCoeff(taub_ZZX,i,coeff);
        }
        Time[i] = GetTime() - time;
    }
    DataProcess(mean,stdev,Time,cyctimes);
    cout << "Output algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";

}