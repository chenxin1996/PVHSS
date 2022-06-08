#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include "PKE.h"

using namespace std;
using namespace NTL;

struct VHSS_Para {
    vec_ZZ_pX vhssEk_1;
    vec_ZZ_pX vhssEk_2;
    vec_ZZ_pX vhssEk_3;
    vec_ZZ_pX vhssEk_4;
    ZZ_pX alpha;
};

void VHSS_Gen(VHSS_Para &vhssPara,PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX pkeSk) {
    Random_ZZ_pX(vhssPara.alpha,pkePara.N,1);
    vec_ZZ_pX alpha_pkeSk;
    alpha_pkeSk.SetLength(2);
    vhssPara.vhssEk_1.SetLength(2);
    vhssPara.vhssEk_2.SetLength(2);
    vhssPara.vhssEk_3.SetLength(2);
    vhssPara.vhssEk_4.SetLength(2);

    alpha_pkeSk[0]=vhssPara.alpha;
    MulMod(alpha_pkeSk[1],vhssPara.alpha,pkeSk[1],modulus);

    Random_ZZ_pX(vhssPara.vhssEk_1[0], pkePara.N, pkePara.q_bit);
    Random_ZZ_pX(vhssPara.vhssEk_1[1], pkePara.N, pkePara.q_bit);
    Random_ZZ_pX(vhssPara.vhssEk_3[0], pkePara.N, pkePara.q_bit);
    Random_ZZ_pX(vhssPara.vhssEk_3[1], pkePara.N, pkePara.q_bit);

    vhssPara.vhssEk_2[0] = pkeSk[0] - vhssPara.vhssEk_1[0];
    vhssPara.vhssEk_2[1] = pkeSk[1] - vhssPara.vhssEk_1[1];
    vhssPara.vhssEk_4[0] =  alpha_pkeSk[0] - vhssPara.vhssEk_3[0];
    vhssPara.vhssEk_4[1] =  alpha_pkeSk[1] - vhssPara.vhssEk_3[1];

}

void VHSS_Enc(vec_ZZ_pX &C, const PKE_Para &pkePara, ZZ_pXModulus &modulus, vec_ZZ_pX &pkePk, const ZZ &x) {
    PKE_OKDM(C, pkePara, modulus, pkePk, x);
}

void VHSS_Mult(vec_ZZ_pX &db, const PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX pkeSk, vec_ZZ_pX C) {
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
            VHSS_Mult(tb_temp, pkePara, modulus, ek, C_X[ind_var[0]]);

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
            VHSS_Mult(tb_temp, pkePara, modulus, ek, C_X[ind_var[0]]);
            prfkey = (prfkey + 1) % 10;

            if (b == 1) {
                tb_temp[0] = tb_temp[0] + PRF[prfkey][0];
                tb_temp[1] = tb_temp[1] + PRF[prfkey][1];
            } else {
                tb_temp[0] = tb_temp[0] - PRF[prfkey][0];
                tb_temp[1] = tb_temp[1] - PRF[prfkey][1];
            }

            for (int i = 1; i < d; i++) {
                VHSS_Mult(tb_temp, pkePara, modulus, tb_temp, C_X[ind_var[i]]);
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


void VHSS_Eval(ZZ_pX &tb_y, int b, PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek, Vec<vec_ZZ_pX> C_X,
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