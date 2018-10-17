#include <iostream>
#include "solver.h"

using namespace std;
using namespace pcenter_solver;

int main(int argc, char *argv[]) {
    Solver solver(10);
    solver.LoadGraph("instance/pr226.tsp");
    solver.Solve();
    system("pause");
}
