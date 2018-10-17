#include <iostream>
#include "solver.h"

using namespace std;
using namespace pcenter_solver;

int main(int argc, char *argv[]) {
    Solver solver;
    solver.LoadGraph("instance/pmed2.txt");
    solver.Solve();
    system("pause");
}
