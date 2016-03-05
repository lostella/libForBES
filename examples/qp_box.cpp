#include "ForBES.h"


int main(int argc, char** argv) {

    const size_t n = 4;
    double data_Q[] = {
        7, 2, -2, -1,
        2, 3, 0, -1,
        -2, 0, 3, -1,
        -1, -1, -1, 1
    };
    double data_q[] = {1, 2, 3, 4};
    Matrix Q(n, n, data_Q);
    Matrix q(n, 1, data_q);
    
    Function * f = new Quadratic(Q, q);
    
    double lb = -1;
    double ub = +1;
    Function * g = new IndBox(lb, ub);
    
    // Define the forward-backward problem specifications
    FBProblem prob(*f, *g);
    
    // Introduce a stopping criterion:
    const double rel_tolerance = 1e-3;    
    FBStoppingRelative sc(rel_tolerance);
    
    // Choose an initial guess
    Matrix x0(n, 1);
    
    // Construct an instance of a FB solver
    double gamma = 0.1;
    const int maxit = 100;
    FBSplitting * solver = new FBSplitting(prob, x0, gamma, sc, maxit);
    
    // Run the solver and get the minimizer
    int solver_status = solver->run();
    Matrix xstar = solver->getSolution();

    std::cout << "status : " << solver_status << std::endl
              << std::endl
              << "optimizer : " << xstar;

    delete solver;
    delete f;
    delete g;

    return EXIT_SUCCESS;
}
