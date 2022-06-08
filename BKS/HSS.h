#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include "PKE.h"

using namespace std;
using namespace NTL;


void HSS_Gen(vec_ZZ_pX &hssEk_1, vec_ZZ_pX &hssEk_2,
             PKE_Para pkePara, vec_ZZ_pX pkeSk) {

    Random_ZZ_pX(hssEk_1[0], pkePara.N, pkePara.q_bit);
    Random_ZZ_pX(hssEk_1[1], pkePara.N, pkePara.q_bit);

    hssEk_2[0] = pkeSk[0] - hssEk_1[0];
    hssEk_2[1] = pkeSk[1] - hssEk_1[1];


}

void HSS_Enc(vec_ZZ_pX &C, const PKE_Para &pkePara, ZZ_pXModulus &modulus, vec_ZZ_pX &pkePk, const ZZ &x) {
    PKE_OKDM(C, pkePara, modulus, pkePk, x);
}

void HSS_Mult(vec_ZZ_pX &db, const PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX pkeSk, vec_ZZ_pX C) {
    PKE_DDec(db, pkePara, modulus, pkeSk, C);
}



int prfkey = 1;

void f(vec_ZZ_pX &tb, int b, int d, int num_data, int loop, int beg_ind, int *ind_var,
       PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek, Vec<vec_ZZ_pX> C_X, Vec<vec_ZZ_pX> PRF) {
    if (loop == d) {
        vec_ZZ_pX tb_temp;
        tb_temp.SetLength(2);
        if (d == 1) {
            prfkey = (prfkey + 1) % 10;
            HSS_Mult(tb_temp, pkePara, modulus, ek, C_X[ind_var[0]]);

            if (b == 1) {
                tb_temp[0] = tb_temp[0] + PRF[prfkey][0];
                tb_temp[1] = tb_temp[1] + PRF[prfkey][1];
            } else {
                tb_temp[0] = tb_temp[0] - PRF[prfkey][0];
                tb_temp[1] = tb_temp[1] - PRF[prfkey][1];
            }

            tb[0] = tb[0] + tb_temp[0];
            tb[1] = tb[1] + tb_temp[1];

        } else {
            HSS_Mult(tb_temp, pkePara, modulus, ek, C_X[ind_var[0]]);
            prfkey = (prfkey + 1) % 10;

            if (b == 1) {
                tb_temp[0] = tb_temp[0] + PRF[prfkey][0];
                tb_temp[1] = tb_temp[1] + PRF[prfkey][1];
            } else {
                tb_temp[0] = tb_temp[0] - PRF[prfkey][0];
                tb_temp[1] = tb_temp[1] - PRF[prfkey][1];
            }

            for (int i = 1; i < d; i++) {
                HSS_Mult(tb_temp, pkePara, modulus, tb_temp, C_X[ind_var[i]]);
                prfkey = (prfkey + 1) % 10;
                if (b == 1) {
                    tb_temp[0] = tb_temp[0] + PRF[prfkey][0];
                    tb_temp[1] = tb_temp[1] + PRF[prfkey][1];
                } else {
                    tb_temp[0] = tb_temp[0] - PRF[prfkey][0];
                    tb_temp[1] = tb_temp[1] - PRF[prfkey][1];
                }
            }
            tb[0] = tb[0] + tb_temp[0];
            tb[1] = tb[1] + tb_temp[1];
        }
    } else {
        loop = loop + 1;
        for (int i = beg_ind; i < num_data; i++) {
            ind_var[loop - 1] = i;
            f(tb, b, d, num_data, loop, i, ind_var,
              pkePara, modulus, ek, C_X, PRF);
        }
    }
}


void HSS_Eval(ZZ_pX &tb_y, int b, PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek, Vec<vec_ZZ_pX> C_X,
              Vec<vec_ZZ_pX> PRF) {
    int loop = 0;
    int beg_ind = 0;
    int *ind_var = (int *) malloc(sizeof(int) * pkePara.d);
    vec_ZZ_pX tb, tb_temp;
    tb.SetLength(2);
    tb_temp.SetLength(2);
    for (int i = 1; i < pkePara.d + 1; i++) {
        f(tb_temp, b, i, pkePara.num_data, loop, beg_ind, ind_var, pkePara, modulus, ek, C_X, PRF);
        prfkey = 1;
        tb[0] = tb[0] + tb_temp[0];
        tb[1] = tb[1] + tb_temp[1];
        tb_temp[0] = 0;
        tb_temp[1] = 0;
    }
    tb_y = tb[0];
}