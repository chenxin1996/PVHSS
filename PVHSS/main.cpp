#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
extern "C" {
#include <relic/relic.h>
}
#include <cstring>
#include "function.h"
using namespace std;

int main() {
    cout << "******************************************************************************************" <<"\n";
    cout << "It is important to note that the generation algorithm, encryption algorithm" <<"\n";
    cout << "and verification algorithm are executed by resource-constrained clients in our model." <<"\n";
    cout << "If you want to reproduce the experimental results of the above three algorithms in" <<"\n";
    cout << "the paper, you need to limit the frequency of the CPU to about 800MHZ." <<"\n";
    cout << "There is no need to limit the CPU for the evaluation algorithm run by the servers." <<"\n";
    cout << "******************************************************************************************" <<"\n\n";

    cout << "************************Benchmarks (Section V-A in the paper)************************" <<"\n";
    // the size of the input (i.e. the value of Log(B_msg) in the paper),
    // and it set to 1,16,32,64,128,256,512,1024 in the experiments in our paper.
    int msg_bit=16;
    // the number of times the program is repeated
    int cyctimes=2;
    // test the running time of generation algorithm
    Time_Gen(msg_bit,  cyctimes);
    // test the running time of encryption algorithm
    Time_Enc(msg_bit,  cyctimes);
    // test the running time of verification algorithm
    Time_Ver(msg_bit,  cyctimes);
    // test the running time of each subroutine in evaluation algorithm
    Time_Eval_Subalgo(msg_bit,  cyctimes);

    cout << "\n\n**************Multivariate polynomial evaluation (Section V-B in the paper)**************" <<"\n";
    // the degree of the polynomial, and it set to 2,4,...,20 in the experiments in our paper.
    int deg=2;
    // the number of the variables, and it set to 4 in the experiments in our paper.
    int num_data=2;
    // evaluating a polynomial with the input size 32bit
    Eval_Poly(deg, num_data);


}
