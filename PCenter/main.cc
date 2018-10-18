#include "pcenter_solver.h"
#include <iostream>

using namespace std;
using namespace pcenter_solver;

constexpr char help_info[] = {
"Input switch is as follows:\n\
 -i: instance name;\n\
 -p: facility number;\n\
 -k: consider range;\n\
 -t: time out seconds;\n\
 -s: random seed value;\n\
 -m: max iteration steps;\n\
Switch -i is necessary and other switchs is optional." };

int main(int argc, char *argv[]) {
    int facility_num = 10;
    char * instance_name = nullptr;
    PCenterSolver::Configure cfg;
    for (int i = 1; i < argc; ++i) {
        if (string(argv[i]) == "-i")
            instance_name = argv[++i];
        else if (string(argv[i]) == "-p")
            facility_num = atoi(argv[++i]);
        else if (string(argv[i]) == "-k")
            cfg.consider_range = atoi(argv[++i]);
        else if (string(argv[i]) == "-t")
            cfg.timeout_seconds = atoi(argv[++i]);
        else if (string(argv[i]) == "-s")
            cfg.random_seed = atol(argv[++i]);
        else if (string(argv[i]) == "-m")
            cfg.max_tabu_steps = atoi(argv[++i]);
    }
    if (instance_name != nullptr) {
        PCenterSolver solver(instance_name, facility_num);
        solver.Solve();
        solver.Record();
    } else {
        cout << help_info << endl;
    }
    return 0;
}
