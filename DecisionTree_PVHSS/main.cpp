#include <iostream>
#include "PVHSS.h"
#include "Eval_DT.h"
int main() {
    cout << "******************************************************************************************" <<"\n";
    cout << "It is important to note that the verification algorithm are executed by" <<"\n";
    cout << "resource-constrained clients in our model." <<"\n";
    cout << "If you want to reproduce the experimental results of the above three algorithms in" <<"\n";
    cout << "the paper, you need to limit the frequency of the CPU to about 800MHZ." <<"\n";
    cout << "There is no need to limit the CPU for the evaluation algorithm run by the servers." <<"\n";
    cout << "******************************************************************************************" <<"\n\n";

    // model=1: EEG Eye State Data Set
    // model=2: Bank Marketing Data Set
    int model=1;
    // the depth of decision tree
    int d=4;
    // the number of times the program is repeated
    int cyctimes=100;

    EvalDT(model,d,cyctimes);

    return 0;
}
