#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 20;
double dt = 0.1;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

// Both the reference cross track and orientation errors are 0.
double v_ref = 30;

// The solver takes all the state variables and actuator
// variables in a singular vector. Thus, we should to establish
// when one variable starts and another ends to make our lifes easier.
size_t x_start = 0;
size_t y_start = x_start + N;
size_t p_start = y_start + N;
size_t v_start = p_start + N;
size_t cte_start = v_start + N;
size_t ep_start = cte_start + N;
size_t d_start = ep_start + N;
size_t a_start = d_start + N - 1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable
    // values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
    
    //---------------------
    // ERROR SETUP:
    //---------------------
    fg[0] = 0; // Init the error term
    
    double time_weight = 1;
    for (int t = 0; t < N; t++) {
      // TODO: maybe can be normal double
      AD<double> ev_t = vars[v_start + t] - v_ref;
      
      // Add in error terms: cte, angular err, velocity err
      fg[0] += CppAD::pow(vars[cte_start + t], 2) * 100 * time_weight;
      fg[0] += CppAD::pow(vars[ep_start + t], 2) * time_weight;
      fg[0] += CppAD::pow(ev_t, 2) * 100 * time_weight;
      
      if (t < N - 1) {
        // Discourage high actuation values (there are only N-1 actuations)
        fg[0] += CppAD::pow(vars[d_start + t], 2) * 80000 * time_weight;
        fg[0] += CppAD::pow(vars[a_start + t], 2) * time_weight;
        
        if (t < N - 2) {
          // Discourage high change in actuation values for smoother turning/acceleration.
          AD<double> d_diff = vars[d_start + t + 1] - vars[d_start + t];
          AD<double> a_diff = vars[a_start + t + 1] - vars[a_start + t];
          
          fg[0] += CppAD::pow(d_diff, 2) * 8000 * time_weight;
          fg[0] += CppAD::pow(a_diff, 2) * time_weight;
        }
      }
      time_weight *= 0.9;
    }
    
    //---------------------
    // CONSTRAINT SETUP:
    //---------------------
    
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`. This bumps up the position of all the other values.
    
    // Constraints at first timestep are just that all vars = initial value
    fg[1 + x_start]   = vars[x_start];
    fg[1 + y_start]   = vars[y_start];
    fg[1 + p_start]   = vars[p_start];
    fg[1 + v_start]   = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + ep_start]  = vars[ep_start];
    
    // Constraints at later timesteps are var_new = h(var_old) =>
    // var_new - h(var_old) = 0. Where h() is the transition function
    // for each variable.
    for (int t = 1; t < N - 1; t++) {
      // The state at time t - 1.
      AD<double> x_old   = vars[x_start + t - 1];
      AD<double> y_old   = vars[y_start + t - 1];
      AD<double> p_old   = vars[p_start + t - 1];
      AD<double> v_old   = vars[v_start + t - 1];
      AD<double> cte_old = vars[cte_start + t - 1];
      AD<double> ep_old  = vars[ep_start + t - 1];
      
      // Actuation at time t - 1
      AD<double> d_old = vars[d_start + t - 1];
      AD<double> a_old = vars[a_start + t - 1];
      
      // The desired trajectory at time t - 1
      AD<double> f_old = coeffs[0] +
                         coeffs[1] * x_old +
                         coeffs[2] * pow(x_old, 2) +
                         coeffs[3] * pow(x_old, 3);
      AD<double> p_des_old = CppAD::atan(coeffs[1] +
                                         2 * coeffs[2] * x_old +
                                         3 * coeffs[3] * pow(x_old, 2));
      
      // The state at time t
      AD<double> x_new   = vars[x_start + t];
      AD<double> y_new   = vars[y_start + t];
      AD<double> p_new   = vars[p_start + t];
      AD<double> v_new   = vars[v_start + t];
      AD<double> cte_new = vars[cte_start + t];
      AD<double> ep_new  = vars[ep_start + t];
      
      // Transition constraints:
      fg[2 + x_start + t] = x_new - (x_old + v_old * CppAD::cos(p_old) * dt);
      fg[2 + y_start + t] = y_new - (y_old + v_old * CppAD::sin(p_old) * dt);
      fg[2 + p_start + t] = p_new - (p_old + v_old * d_old / Lf * dt);
      fg[2 + v_start + t] = v_new - (v_old + a_old * dt);
      fg[2 + cte_start + t] =
        cte_new - ((f_old - y_old) + (v_old * CppAD::sin(ep_old) * dt));
      fg[2 + ep_start + t] =
        ep_new - ((p_old - p_des_old) + v_old * d_old / Lf * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;
  
  double x   = state[0];
  double y   = state[1];
  double p   = state[2];
  double v   = state[3];
  double cte = state[4];
  double ep  = state[5];

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = N * 6 + (N - 1) * 2;
  // TODO: Set the number of constraints
  size_t n_constraints = N * 6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
  vars[x_start] = x;
  vars[y_start] = y;
  vars[p_start] = p;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[ep_start] = ep;

  // Set lower and upper limits for variables.
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  
  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (int i = 0; i < d_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }
  
  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  for (int i = d_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }
  
  // Acceleration/decceleration upper and lower limits.
  for (int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[p_start] = p;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[ep_start] = ep;
  
  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[p_start] = p;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[ep_start] = ep;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  
  vector<double> results;
  
  // Add actuation to return vector
  results.push_back(solution.x[d_start]);
  results.push_back(solution.x[a_start]);
  
  // Add predicted locations to return vector
  for (int t = 0; t < N - 1; t++) {
    results.push_back(solution.x[x_start + t]);
  }
  for (int t = 0; t < N - 1; t++) {
    results.push_back(solution.x[y_start + t]);
  }

  return results;
}
