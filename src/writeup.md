# MPC Project

### Model
- The state of my model contains:
  - The car's x position
  - The car's y position
  - The car's orientation (psi)
  - The car's velocity
  - The cross track error
  - The orientation error
- My actuators are:
  - Acceleration / throttle
  - Steering angle
- The state prediction is updated in simulation via the standard kinematic motion models described in class, (e.g. x_t+1 = x_t + v_t * cos(psi_t) * delta_t).

### Timestep Length and Duration
I chose N and dt empirically to provide a balance between enough lookahead time and spatio-temporal resolution (for N and dt, respectively) vs wasting computational resources needed to run the controller in real-time.

### Preprocessing
I converted velocity input from MPH to m/s to match with units from other variables. I also converted world-coordinate measurements to car coordinates to fit the polynomial. This led to cleaner code and easier error computation.

### Latency
I handled latency by setting the initial state of the car in the optimizer to be the predicted state 100ms ahead. This enabled the optimizer to provide the correct actuations at the time when they were executed.

### Other notes
To get the model to work, I had to heavily weight the angle magnitude term relative to the cte term in the total error. Without this, the car oscilates wildly off the track. With more time, I'd like to experiment with tuning parameters further so that the car can stay under control while going faster.
