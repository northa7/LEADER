% -------------------------------------------------------------------------
% Generates a coordinate transformation matrix (output: T)
% that transforms a vector from the inertial coordinate system
% to the asteroid's own coordinate system
%
% Recommended input parameters: phi0 = 0, omega = 2*pi/period, t0 = 0
% -------------------------------------------------------------------------

function T = transform_mat(phi0, omega, t, t0, beta, lambda)

arg1 = phi0 + omega*(t-t0);
RZ1 = [cos(arg1), sin(arg1), 0; -sin(arg1), cos(arg1), 0; 0, 0, 1];
RY1 = [cos(beta), 0, -sin(beta); 0, 1, 0; sin(beta), 0, cos(beta)];
RZ2 = [cos(lambda), sin(lambda), 0; -sin(lambda), cos(lambda), 0; 0, 0, 1];

T = RZ1*RY1*RZ2;