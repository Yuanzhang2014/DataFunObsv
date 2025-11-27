Target estimation of nonlinear neuron networks using limited sensors
===
We apply the data-driven observer framework to target estimation in neuron networks modeled by the classic Hindmarsh-Rose system.

usage
---
* `neuron_nonlinear_Koopman.m` : Main script for estimation performance for the nonlinear neuron networks using limited sensors by Koopman operator. It generates Figures 6 (b), (c).
* `neuron_nonlinear_Koopman_mc.m`: Monte Carlo evaluation of Koopman-based nonlinear observer performance by computing time-varying RMSE (Root Mean Square Error) across Nr experiments. It generates Figure 6 (d).
* `neuron_nonlinear.m`: Generate nonlinear neuron network system trajectories.

* `lift_koopman`: Create thin plate spline radial basis functions and lift original state data to higher-dimensional space.

* `Darouch_observer_with_trajectory_pure_data.m`: Design Darouch functional observer using pure trajectory data without system model.

* `check_consistent.m`: Check the consistency condition for observer design.

* `Simulation_agumentation_observer_with_trajectory_pure_data.m`: Simulate and evaluate the performance of a data-driven augmented state observer. It generates Figures. 6 (b).
