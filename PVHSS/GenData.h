#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/matrix.h>
#include <NTL/vec_vec_ZZ.h>



using namespace std;
using namespace NTL;
struct Data {
    vec_ZZ X_decimal;
    Vec<ZZ_pX> X;
    Vec<vec_ZZ_pX> C_X;
    Vec<vec_ZZ_pX> PRF;
};

void GenData(Data &data, PKE_Para pkePara,ZZ_pXModulus &modulus, vec_ZZ_pX pkePk) {
    data.X.SetLength(pkePara.num_data);
    vec_ZZ_pX C_x, prf;
    C_x.SetLength(4);
    prf.SetLength(2);
    data.X_decimal.SetLength(pkePara.num_data);
    double time;
    for (int i = 0; i < pkePara.num_data; i++) {
        data.X_decimal[i]= RandomBits_ZZ(pkePara.msg_bit);
        time = GetTime();
        Decimal2Bin(data.X[i],data.X_decimal[i],pkePara.msg_bit);
        time = GetTime() - time;
        //cout << "decimal2Bin time: " << time * 1000 <<"ms\n";
        PVHSS_Enc(C_x, pkePara, modulus, pkePk, data.X[i]);
        data.C_X.append(C_x);
    }
    for (int i = 0; i < 10; i++) {
        Random_ZZ_pX(prf[0], pkePara.N, pkePara.q_bit);
        Random_ZZ_pX(prf[1], pkePara.N, pkePara.q_bit);
        data.PRF.append(prf);
    }
}