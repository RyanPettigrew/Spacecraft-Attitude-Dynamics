function [dq,dE] = kinematics(omega, q, E)
%#codegen

I = eye(3); % Identity

% unpack quaternion
eps1 = q(1);
eps2 = q(2);
eps3 = q(3);
eps = [eps1; eps2; eps3];

eda = q(4);

% Compute the derivative of the quaternion components here
             
% Quaternion rates Eq.1.56
eps_cross = [0,  -eps3, eps2;
             eps3,  0, -eps1;
             -eps2, eps1, 0];

deps = 0.5*(eda*I + eps_cross)*omega;
deda = -0.5*eps'*omega;

dq = [deps; deda];

% Compute the derivative of the Euler Angle components here
            
% Euler rates
phi = E(1);
theta = E(2);
R = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);
     0     cos(phi)           -sin(phi);
     0 sin(phi)*sec(theta) cos(phi)*sec(theta)];
    
dE = R*omega;

end
