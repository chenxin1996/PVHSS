#include <gmp.h>
extern "C"{
#include <relic/relic.h>
}
#include <sstream>

using namespace std;
using namespace NTL;
#define PI 3.141592654

void Random_ZZ_pX(ZZ_pX &a, int N, int q_bit) {
    ZZ_p coeff;
    for (int i = 0; i < N; i++) {
        conv(coeff, RandomBits_ZZ(q_bit));
        SetCoeff(a, i, coeff);
    }
}

void SecretKey(ZZ_pX &sk, int N, int hsk) {
    int interval = 0;
    interval = N / hsk;
    int index = rand() % interval;
    for (int i = 0; i < hsk; i++) {
        SetCoeff(sk, index, 1);
        index = index + rand() % interval;
    }

}

void GaussRand(ZZ_pX &e, int N) {
    double res_standard;
    int deviation = 8;
    int res;
    for (int i = 0; i < N; i++) {
        res_standard = sqrt(-2.0 * log(rand() / (RAND_MAX + 1.0))) * sin(2.0 * PI * rand() / (RAND_MAX + 1.0));
        res = res_standard * deviation;
        SetCoeff(e, i, res);
    }
}

void NativeEval_f(ZZ &y, int d, int num_data, int loop, int beg_ind, int *ind_var, vec_ZZ X) {
    if (loop == d) {
        ZZ temp;
        temp=1;
        for (int i = 0; i < d; i++) {
            temp = temp * X[ind_var[i]];
        }
        y = y + temp;
    } else {
        loop = loop + 1;
        for (int i = beg_ind; i < num_data; i++) {
            ind_var[loop - 1] = i;
            NativeEval_f(y, d, num_data, loop, i, ind_var, X);
        }
    }

}

void NativeEval(ZZ &y, int d, int num_data, vec_ZZ X) {
    ZZ temp;
    for (int i = 1; i < d + 1; i++) {
        int *ind_var = (int *) malloc(sizeof(int) * i);
        NativeEval_f(temp, i, num_data, 0, 0, ind_var, X);
        y = y + temp;
        temp = 0;
    }
}

void Decimal2Bin(ZZ_pX &out, ZZ in, int bits)
{
    ZZ rem;
    ZZX out_ZZ_pX;
    for(int i=0;i<bits;i++)
    {
        rem=in%2;
        SetCoeff(out_ZZ_pX,i, rem);
        in=(in-rem)/2;
    }
    conv(out,out_ZZ_pX);
}

void ZZ_pX_ScaleMul_ZZ(ZZ_pX &out, ZZ_pX in1,ZZ in2)
{
    ZZ_pX a;
    conv(a,in2);
    out=in1*a;
}


void ZZ2bn(bn_t out, ZZ in)
{
    if(in==0)
    {
        bn_read_str(out,"0",1,10);
    }else
    {
        std::stringstream buffer;
        buffer << in;
        const char* str=strdup(buffer.str().c_str());
        bn_read_str(out,str,floor(log(in)/log(10))+1,10);
    }

}

void eval_ZZX(ZZ &y,ZZX f)
{
    ZZ coeff,pow2;
    pow2=1;
    for(int i=0;i<deg(f)+1;i++)
    {
        GetCoeff(coeff,f,i);
        y=y+coeff*pow2;
        pow2= pow2*2;
    }

}

void DataProcess(double &mean, double &stdev,  double *Time,int cyctimes)
{
    double temp;
    double sum=0;
    for(int i=0;i<cyctimes;i++)
    {
        sum=sum+Time[i];
    }
    mean=sum/cyctimes;
    double temp_sum=0;
    for(int i=0;i<cyctimes;i++)
    {
        temp=mean-Time[i];
        temp=temp*temp;
        temp_sum=temp_sum+temp;
    }
    stdev=sqrt(temp_sum/cyctimes);
    stdev=stdev/mean;
}

