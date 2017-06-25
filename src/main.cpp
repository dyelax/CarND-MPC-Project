#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double x = j[1]["x"];
          double y = j[1]["y"];
          double p = j[1]["psi"];
          double v = j[1]["speed"];
          double d = j[1]["steering_angle"];
          double a = j[1]["throttle"];
          
          // Convert d to rad and v to m/s
          d *= -1;
          v *= 0.44704;
          
          // Simulate latency in the car state and transform into car space
          const double Lf = 2.67;
          double latency = 0.1;
          float x_car = 0;
          float y_car = 0;
          float p_car = 0;
          float v_car = v;
          
          float x_car_latency = x_car + v * cos(p_car) * latency;
          float y_car_latency = y_car + v * sin(p_car) * latency;
          float p_car_latency = p_car + (v / Lf) * d * latency;
          float v_car_latency = v_car + a * latency;

          // Convert world coords to car coords
          Eigen::VectorXd ptsx_car(ptsx.size());
          Eigen::VectorXd ptsy_car(ptsy.size());
          
          for (int i = 0; i < ptsx.size(); i++) {
            ptsx_car[i] = (ptsx[i] - x) * cos(p) + (ptsy[i] - y) * sin(p);
            ptsy_car[i] = -(ptsx[i] - x) * sin(p) + (ptsy[i] - y) * cos(p);
          }
          
          // Compute desired trajectory
          Eigen::VectorXd coeffs = polyfit(ptsx_car, ptsy_car, 3);
          
          // Compute errors in car space
          double cte_car = polyeval(coeffs, 0);
          double ep_car = -atan(coeffs[1]);
          
          // Simulate latency in the errors
          double cte_car_latency = cte_car - y_car + (v_car * sin(ep_car) * latency);
          double ep_car_latency = ep_car + (v / Lf) * d * latency;
          
          // Build the state
          Eigen::VectorXd state(6);
          state << x_car_latency,
                   y_car_latency,
                   p_car_latency,
                   v_car_latency,
                   cte_car_latency,
                   ep_car_latency;
//          state << 0, 0, 0, v, cte_car, ep_car;
          
          vector<double> solution = mpc.Solve(state, coeffs);
          
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          double steer_value = -solution[0] / deg2rad(25);
          double throttle_value = solution[1];
          
          json msgJson;
          
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          //Display the MPC predicted trajectory
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;
          
          int N = (solution.size() - 2) / 2;
          
          for (int i = 2; i < N + 2; i++) {
            mpc_x_vals.push_back(solution[i]);
            mpc_y_vals.push_back(solution[i+N]);
          }
          
          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

//          Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;
          
          for (int i = 0; i < ptsx.size(); i++) {
            next_x_vals.push_back(ptsx_car[i]);
            next_y_vals.push_back(ptsy_car[i]);
          }

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
