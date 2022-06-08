#include <gmp.h>
extern "C" {
#include <relic/relic.h>
}
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include "PKE.h"

using namespace std;
using namespace NTL;
typedef struct {
    ZZ g1_order_ZZ;
    ZZ g2_order_ZZ;
    bn_t g1_order;
    bn_t g2_order;
    ep_t g1_gen;
    ep2_t g2_gen;
    fp12_t gT_gen;
    fp12_t A_gT;
    fp12_t const_gT;

} PVHSS_Para;

void PVHSS_Enc(vec_ZZ_pX &C, const PKE_Para &pkePara, ZZ_pXModulus &modulus, vec_ZZ_pX &pkePk, const ZZ_pX &x);

void PVHSS_Gen(vec_ZZ_pX &hssEk_1, vec_ZZ_pX &hssEk_2, vec_ZZ_pX &C_alpha, PVHSS_Para &pvhssPara,
               PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX pkePk, vec_ZZ_pX pkeSk) {
    // ek
    Random_ZZ_pX(hssEk_1[0], pkePara.N, pkePara.q_bit);
    Random_ZZ_pX(hssEk_1[1], pkePara.N, pkePara.q_bit);
    hssEk_2[0] = pkeSk[0] - hssEk_1[0];
    hssEk_2[1] = pkeSk[1] - hssEk_1[1];
    //alpha
    ZZ A;
    ZZ_pX alpha;
    A = RandomBits_ZZ(128);
    Decimal2Bin(alpha, A, 128);
    PVHSS_Enc(C_alpha, pkePara, modulus, pkePk, alpha);
    //pairing para
    bn_new(pvhssPara.g1_order);
    bn_new(pvhssPara.g2_order);
    ep_new(pvhssPara.g1_gen);
    ep2_new(pvhssPara.g2_gen);
    fp12_new(pvhssPara.gT_gen);
    fp12_new(pvhssPara.A_gT);
    fp12_new(pvhssPara.const_gt);
    core_init();
    ep_curve_init();
    ep_param_set(B12_P381);
    ep2_curve_init();
    ep2_curve_set_twist(EP_MTYPE);
    ep_curve_get_gen(pvhssPara.g1_gen);
    ep2_curve_get_gen(pvhssPara.g2_gen);
    pp_map_oatep_k12(pvhssPara.gT_gen, pvhssPara.g1_gen, pvhssPara.g2_gen);
    bn_t A_bn;
    ep_t temp;
    bn_new(A_bn);
    ep_new(temp);
    ZZ2bn(A_bn, A);
    ep_mul_gen(temp, A_bn);
    pp_map_oatep_k12(pvhssPara.A_gT, temp, pvhssPara.g2_gen);
    ep_curve_get_ord(pvhssPara.g1_order);
    ep2_curve_get_ord(pvhssPara.g2_order);
    int size = bn_size_str(pvhssPara.g1_order, 10);
    char *g1_order_str = new char[size];
    bn_write_str(g1_order_str, size, pvhssPara.g1_order, 10);
    ZZ g1_order_ZZ = conv<ZZ>(g1_order_str);
    pvhssPara.g1_order_ZZ = g1_order_ZZ;
    size = bn_size_str(pvhssPara.g2_order, 10);
    char *g2_order_str = new char[size];
    bn_write_str(g2_order_str, size, pvhssPara.g2_order, 10);
    ZZ g2_order_ZZ = conv<ZZ>(g2_order_str);
    pvhssPara.g2_order_ZZ = g2_order_ZZ;
    //compute e(g_1,g_2)^{q(2^N-1)}
    bn_t const_bn;
    bn_new(const_bn);
    ZZ temp1;
    temp1 = pkePara.q * (power_ZZ(2, pkePara.N) - 1);
    temp1 = temp1 % pvhssPara.g1_order_ZZ;
    ZZ2bn(const_bn, temp1);
    ep_mul_gen(temp, const_bn);
    pp_map_oatep_k12(pvhssPara.const_gT, temp, pvhssPara.g2_gen);
}

void PVHSS_Enc(vec_ZZ_pX &C, const PKE_Para &pkePara, ZZ_pXModulus &modulus, vec_ZZ_pX &pkePk, const ZZ_pX &x) {
    PKE_OKDM(C, pkePara, modulus, pkePk, x);
}

void PVHSS_Mult(vec_ZZ_pX &db, const PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX pkeSk, vec_ZZ_pX C) {
    PKE_DDec(db, pkePara, modulus, pkeSk, C);
}

void PVHSS_Add(vec_ZZ_pX &out, vec_ZZ_pX in1,vec_ZZ_pX in2) {
    out[0]=in1[0]+in2[0];
    out[1]=in1[1]+in2[1];
}

void PVHSS_Sub(vec_ZZ_pX &out, vec_ZZ_pX in1,vec_ZZ_pX in2) {
    out[0]=in1[0]-in2[0];
    out[1]=in1[1]-in2[1];
}

int prfkey = 1;

void f(vec_ZZ_pX &tb, int b, int d, int num_data, int loop, int beg_ind, int *ind_var,
       PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek, Vec<vec_ZZ_pX> C_X, Vec<vec_ZZ_pX> PRF) {
    if (loop == d) {
        vec_ZZ_pX tb_temp;
        tb_temp.SetLength(2);
        if (d == 1) {
            prfkey = (prfkey + 1) % 10;
            PVHSS_Mult(tb_temp, pkePara, modulus, ek, C_X[ind_var[0]]);

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
            PVHSS_Mult(tb_temp, pkePara, modulus, ek, C_X[ind_var[0]]);
            prfkey = (prfkey + 1) % 10;

            if (b == 1) {
                tb_temp[0] = tb_temp[0] + PRF[prfkey][0];
                tb_temp[1] = tb_temp[1] + PRF[prfkey][1];
            } else {
                tb_temp[0] = tb_temp[0] - PRF[prfkey][0];
                tb_temp[1] = tb_temp[1] - PRF[prfkey][1];
            }

            for (int i = 1; i < d; i++) {
                PVHSS_Mult(tb_temp, pkePara, modulus, tb_temp, C_X[ind_var[i]]);
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

void PVHSS_Eval(ZZ &Tb,ZZ_pX &tb0_test, ZZ_p &yb, ep_t g1T1, ep2_t g2T2, int b, PVHSS_Para pvhssPara, PKE_Para pkePara,
                ZZ_pXModulus modulus, vec_ZZ_pX ek, vec_ZZ_pX C_alpha,
                Vec<vec_ZZ_pX> C_X, Vec<vec_ZZ_pX> PRF) {
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
    tb0_test = tb[0];


    //output
    ZZX tb0;
    ZZ two_ZZ;
    ZZ_p two_ZZ_p;
    bn_t Tb_bn;
    two_ZZ_p = 2;
    bn_new(Tb_bn);
    eval(yb, tb[0], two_ZZ_p);
    PVHSS_Mult(tb, pkePara, modulus, tb, C_alpha);
    conv(tb0, tb[0]);
    eval_ZZX(Tb, tb0);
    if (b == 1) {
        Tb = Tb % pvhssPara.g1_order_ZZ;
        ZZ2bn(Tb_bn, Tb);
        ep_mul_gen(g1T1, Tb_bn);
    } else {
        Tb = Tb % pvhssPara.g2_order_ZZ;
        ZZ2bn(Tb_bn, Tb);
        ep2_mul_gen(g2T2, Tb_bn);
    }
}

void PVHSS_Ver(ZZ_p &y, ZZ_p y1, ZZ_p y2, ep_t g1T1, ep2_t g2T2, PVHSS_Para pvhssPara) {
    fp12_t gtT1;
    fp12_t gtT2;
    fp12_t right_side;
    fp12_t left_side;
    fp12_new(gtT1);fp12_new(gtT2);fp12_new(right_side);fp12_new(left_side);
    bn_t y_bn;
    bn_new(y_bn);
    ZZ y_ZZ;
    y = y1 + y2;
    conv(y_ZZ, y);
    y_ZZ = y_ZZ % pvhssPara.g1_order_ZZ;
    ZZ2bn(y_bn, y_ZZ);

    fp12_exp(right_side, pvhssPara.A_gT, y_bn);
    fp12_mul(right_side, right_side, pvhssPara.const_gT);

    pp_map_oatep_k12(gtT1, g1T1,pvhssPara.g2_gen);
    pp_map_oatep_k12(gtT2, pvhssPara.g1_gen, g2T2);
    fp12_mul(left_side, gtT1, gtT2);

    if (fp12_cmp(left_side, right_side) == RLC_EQ) {
    } else {
        printf("ERROR\n");
        cout <<  y << "\n";

        fp12_print(left_side);
        printf("\n\n");
        fp12_print(right_side);

    }
}