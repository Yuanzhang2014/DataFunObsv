function [A,B,B_affine] = discretized_lineared_neuron(x0,dt)

% =========================================================================
% Filename: discretized_lineared_neuron.m
% Purpose: Compute discretized linearized matrices for a single neuron model around an operating point
% Inputs:
%   x0 - Operating point for linearization
%   dt - Discretization time step
% Outputs:
%   A - Discrete-time state transition matrix (3×3)
%   B - Discrete-time input matrix (3×1)  
%   B_affine - Discrete-time affine term (constant vector)
% Main functions:
%   Discretize continuous-time linearized model using forward Euler method
% =========================================================================

a = 1;
b = 3;
c = 1;
d = 5;
r = 0.001;
s = 1;
x1 = -1.6;
A = [1+dt*(2*b*x0-3*a*x0*x0), dt, -1*dt; -1*dt*2*d*x0, 1-dt, 0; dt*r*s, 0, 1-r*dt];
B = [dt;0;0];
B_affine = [dt*(2*a*x0^3-b*x0^2); dt*(c+d*x0*x0); -1*dt*r*s*x1];
end

