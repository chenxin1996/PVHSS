#include <iostream>
#include <helib/helib.h>

using namespace std;
using namespace helib;
using namespace NTL;

void Depth2Tree(Ctxt &res, Ctxt p, Ctxt lc, Ctxt rc, Ctxt cipher_1);

void Depth2Tree(Ctxt &res, Ctxt p, Ctxt lc, Ctxt rc, Ctxt cipher_1) {
    cipher_1 -= p;
    cipher_1 *= lc;
    p *= rc;
    cipher_1 += p;
    res = cipher_1;
}

void DataProcess(double &mean, double &stdev, double *Time, int cyctimes) {
    double temp;
    double sum = 0;
    for (int i = 0; i < cyctimes; i++) {
        sum = sum + Time[i];
    }
    mean = sum / cyctimes;

    double temp_sum = 0;
    for (int i = 0; i < cyctimes; i++) {
        temp = mean - Time[i];
        temp = temp * temp;
        temp_sum = temp_sum + temp;
    }

    stdev = sqrt(temp_sum / cyctimes);
}

void DT(Ctxt &res,int model, int d, vector<Ctxt> b, vector<Ctxt> c,Ctxt cipher_1) {
    if(model==1&&d==4)
    {
        Ctxt temp1_cipher=cipher_1;
        temp1_cipher-=b[3];//1-b4
        temp1_cipher*=c[3];//(1-b4)*c4
        Ctxt temp2_cipher=b[3];
        temp2_cipher*=c[4];//b4*c5
        temp1_cipher+=temp2_cipher;//b4*c5+(1-b4)*c4
        temp1_cipher*=b[2];//b3*(b4*c5+(1-b4)*c4)
        Ctxt temp3_cipher=cipher_1;
        temp3_cipher-=b[2];//1-b3
        temp3_cipher*=c[2];//(1-b3)*c3
        temp1_cipher+=temp3_cipher;//b3*(b4*c5+(1-b4)*c4)+(1-b3)*c3
        temp1_cipher*=b[0];//b1*(b3*(b4*c5+(1-b4)*c4)+(1-b3)*c3)
        temp2_cipher=cipher_1;
        temp2_cipher-=b[1];//1-b2
        temp2_cipher*=c[0];//(1-b2)*c1
        temp3_cipher=b[1];
        temp3_cipher*=c[1];//b2*c2
        temp2_cipher*= temp3_cipher;//b2*c2+(1-b2)*c1
        temp3_cipher=cipher_1;
        temp3_cipher-=b[0];//1-b1
        temp2_cipher*=temp3_cipher;//(1-b1)*(b2*c2+(1-b2)*c1)
        temp1_cipher+=temp2_cipher;//res
        res=temp1_cipher;
    }

    if(model==2&&d==4) {
        Ctxt temp1 = cipher_1;
        temp1 -= b[3];
        temp1 *= c[2];
        Ctxt temp2 = b[3];
        temp2 *= c[3];
        temp1 += temp2;//(1-b4)c3+b4c4
        temp2 = cipher_1;
        temp2 -= b[1];
        temp1 *= temp2;//(1-b2)((1-b4)c3+b4c4)
        temp2 = cipher_1;
        temp2 -= b[4];
        temp2 *= c[4];
        Ctxt temp3 = b[4];
        temp3 *= c[5];
        temp2 += temp3;
        temp2 *= b[1];//b2[(1-b5)c5+c5c6]
        temp1 += temp2;
        temp1 *= b[0];
        temp2 = cipher_1;
        temp2 -= b[2];
        temp2 *= c[0];
        temp3 = b[2];
        temp3 *= c[1];
        temp2 += temp3;
        temp3 = cipher_1;
        temp3 -= b[0];
        temp2 *= temp3;
        temp1 += temp2;
        res = temp1;
    }

    if(model==1&&d==6) {
        Ctxt temp1 = cipher_1;
        Depth2Tree(temp1, b[7], c[6], c[7], cipher_1);
        Ctxt temp2 = cipher_1;
        Depth2Tree(temp2, b[5], c[4], temp1, cipher_1);
        Ctxt temp3 = cipher_1;
        Depth2Tree(temp3, b[8], c[8], c[9], cipher_1);
        Ctxt temp4 = cipher_1;
        Depth2Tree(temp4, b[6], temp3, c[5], cipher_1);
        Ctxt temp5 = cipher_1;
        Depth2Tree(temp5, b[3], temp2, temp4, cipher_1);
        temp1 = cipher_1;
        Depth2Tree(temp1, b[4], c[2], c[3], cipher_1);
        temp2 = cipher_1;
        Depth2Tree(temp2, b[2], temp5, temp1, cipher_1);
        temp3 = cipher_1;
        Depth2Tree(temp3, b[1], c[0], c[1], cipher_1);
        Depth2Tree(res, b[0], temp3, temp2, cipher_1);
    }

    if(model==2&&d==6) {
        Ctxt temp1 = cipher_1;
        Depth2Tree(temp1, b[5], c[5], c[6], cipher_1);
        Ctxt temp2 = cipher_1;
        Depth2Tree(temp2, b[4], c[4], temp1, cipher_1);
        temp1 = cipher_1;
        Depth2Tree(temp1, b[3], c[3], temp2, cipher_1);
        temp2 = cipher_1;
        Depth2Tree(temp2, b[2], c[1], c[2], cipher_1);
        Ctxt temp3 = cipher_1;
        Depth2Tree(temp3, b[1], temp1, temp2, cipher_1);
        Depth2Tree(res, b[0], temp3, c[0], cipher_1);
    }

    if(model==1&&d==8) {
        Ctxt temp1 = cipher_1;
        Depth2Tree(temp1, b[12], c[11], c[12], cipher_1);
        Depth2Tree(temp1, b[10], temp1, c[9], cipher_1);
        Depth2Tree(temp1, b[8], temp1, c[7], cipher_1);
        Ctxt temp2 = cipher_1;
        Depth2Tree(temp2, b[13], c[13], c[14], cipher_1);
        Depth2Tree(temp2, b[11], temp2, c[10], cipher_1);
        Depth2Tree(temp2, b[9], temp2, c[8], cipher_1);
        Depth2Tree(temp1, b[6], temp1, temp2, cipher_1);
        Depth2Tree(temp2, b[7], c[5], c[6], cipher_1);
        Depth2Tree(temp1, b[5], temp1, temp2, cipher_1);
        Depth2Tree(temp2, b[4], c[3], c[4], cipher_1);
        Depth2Tree(temp1, b[2], temp2, temp1, cipher_1);
        Depth2Tree(temp2, b[3], c[1], c[2], cipher_1);
        Depth2Tree(temp2, b[1], temp2, c[0], cipher_1);
        Depth2Tree(res, b[0], temp2, temp1, cipher_1);
    }

    if(model==2&&d==8) {
        Ctxt temp1 = cipher_1;
        Depth2Tree(temp1, b[8], c[8], c[9], cipher_1);
        Ctxt temp2 = cipher_1;
        Depth2Tree(temp2, b[7], c[7], temp1, cipher_1);
        Depth2Tree(temp1, b[6], c[6], temp2, cipher_1);
        Depth2Tree(temp2, b[5], c[5], temp1, cipher_1);
        Depth2Tree(temp1, b[3], temp2, c[2], cipher_1);
        Depth2Tree(temp2, b[4], c[3], c[4], cipher_1);
        Ctxt temp3 = cipher_1;
        Depth2Tree(temp3, b[2], temp2, c[1], cipher_1);
        Depth2Tree(temp2, b[1], temp3, temp1, cipher_1);
        Depth2Tree(res, b[0], temp2, c[0], cipher_1);
    }
}

void Eval_DT(int FHE, int model, int d, int cyctimes) {

    if (FHE == 1) {
        long m, bits;
        if (d == 4) {
            m = 10897;
            bits = 160;
        }
        if (d == 6) {
            m = 17711;
            bits = 240;
        }
        if (d == 8) {
            m = 18721;
            bits = 320;
        }

        helib::Context context = helib::ContextBuilder<helib::BGV>().m(m).bits(bits).r(1).c(2).build();
        SecKey secretKey(context);
        secretKey.GenSecKey();
        const PubKey &publicKey = secretKey;

        Ctxt single_cipher(publicKey);
        Ctxt cipher_1(publicKey);
        Ctxt res(publicKey);
        PtxtArray p(context, 1);
        p.encrypt(cipher_1);
        long num_not_leaf_node = 20;
        long num_leaf_node = 21;
        vector<long> b(num_not_leaf_node);
        vector<long> c(num_leaf_node);
        for (int i = 0; i < num_not_leaf_node; i++) {
            b[i] = random() % 2;
        }
        for (int i = 0; i < num_leaf_node; i++) {
            c[i] = random() % 2;
        }
        vector<Ctxt> cipher_b(num_not_leaf_node, single_cipher);
        vector<Ctxt> cipher_c(num_leaf_node, single_cipher);
        for (int i = 0; i < num_not_leaf_node; i++) {
            PtxtArray p(context, b[i]);
            p.encrypt(cipher_b[i]);
        }
        for (int i = 0; i < num_leaf_node; i++) {
            PtxtArray p(context, c[i]);
            p.encrypt(cipher_c[i]);
        }

        //eval
        cout << "******************** Server Evaluating ********************" << "\n";
        double time;
        auto Time = new double[cyctimes];
        for (int i = 0; i < cyctimes; i++) {
            time = GetTime();
            DT(res,model,d, cipher_b, cipher_c, cipher_1);
            time = GetTime() - time;
            Time[i] = time;
        }
        double eval_time = 0;
        double stdev = 0;
        DataProcess(eval_time, stdev, Time, cyctimes);
        cout << "running time Server: " << eval_time*1000<< " ms  RSD: "<<stdev*100<<"%\n";

        //dec
        cout << "******************** Client Evaluating ********************" << "\n";
        for (int i = 0; i < cyctimes; i++) {
            time = GetTime();
            PtxtArray temp_res(context);
            temp_res.decrypt(res, secretKey);
            vector<long> classes;
            temp_res.store(classes);
            time = GetTime() - time;
            Time[i] = time;
        }
        PtxtArray temp_res(context);
        temp_res.decrypt(res, secretKey);
        vector<long> classes;
        temp_res.store(classes);
        double dec_time = 0;
        stdev = 0;
        DataProcess(dec_time, stdev, Time, cyctimes);
        cout << "running time Client: " << dec_time*1000 << " ms  RSD: "<<stdev*100<<"%\n";
    }

    if (FHE == 2) {
        long m, bits;
        if (d == 4) {
            m = 32768;
            bits = 239;
        }
        if (d == 6) {
            m = 32768;
            bits = 239;
        }
        if (d == 8) {
            m = 32768;
            bits = 358;
        }

        helib::Context context = helib::ContextBuilder<helib::CKKS>().m(m).bits(bits).precision(20).c(6)
                .build();
        SecKey secretKey(context);
        secretKey.GenSecKey();
        const PubKey &publicKey = secretKey;

        Ctxt single_cipher(publicKey);
        Ctxt cipher_1(publicKey);
        Ctxt res(publicKey);
        PtxtArray p(context, 1);
        p.encrypt(cipher_1);
        long num_not_leaf_node = 20;
        long num_leaf_node = 21;
        vector<long> b(num_not_leaf_node);
        vector<long> c(num_leaf_node);
        for (int i = 0; i < num_not_leaf_node; i++) {
            b[i] = random() % 2;
        }
        for (int i = 0; i < num_leaf_node; i++) {
            c[i] = random() % 2;
        }
        vector<Ctxt> cipher_b(num_not_leaf_node, single_cipher);
        vector<Ctxt> cipher_c(num_leaf_node, single_cipher);
        for (int i = 0; i < num_not_leaf_node; i++) {
            PtxtArray p(context, b[i]);
            p.encrypt(cipher_b[i]);
        }
        for (int i = 0; i < num_leaf_node; i++) {
            PtxtArray p(context, c[i]);
            p.encrypt(cipher_c[i]);
        }

        //eval
        cout << "******************** Server Evaluating ********************" << "\n";
        double time;
        auto Time = new double[cyctimes];
        for (int i = 0; i < cyctimes; i++) {
            time = GetTime();
            DT(res,model,d, cipher_b, cipher_c, cipher_1);
            time = GetTime() - time;
            Time[i] = time;
        }
        double eval_time = 0;
        double stdev = 0;
        DataProcess(eval_time, stdev, Time, cyctimes);
        cout << "running time Server: " << eval_time*1000<< " ms  RSD: "<<stdev*100<<"%\n";

        //dec
        cout << "******************** Client Evaluating ********************" << "\n";
        for (int i = 0; i < cyctimes; i++) {
            time = GetTime();
            PtxtArray temp_res(context);
            temp_res.decrypt(res, secretKey);
            vector<long> classes;
            temp_res.store(classes);
            time = GetTime() - time;
            Time[i] = time;
        }
        PtxtArray temp_res(context);
        temp_res.decrypt(res, secretKey);
        vector<long> classes;
        temp_res.store(classes);
        double dec_time = 0;
        stdev = 0;
        DataProcess(dec_time, stdev, Time, cyctimes);
        cout << "running time Client: " << dec_time*1000 << " ms  RSD: "<<stdev*100<<"%\n";
    }
}