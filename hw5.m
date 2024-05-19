% clear; clc;

addpath("Functions\");

% Compute Mass properties for normal operations phase

m = 640; %kg, mass of s/c
rvec_cm = [0; 0; 0.2344]; %m, position vector of s/c
% J = [812.04 0 0; 0 545.37 0; 0 0 627.71]; %kg*m^2, inertia matrix of s/c
J = diag([1200, 2000, 2800]);

% Part 2 - Torque Free Attitudue Simulation
Areas_1= 4*ones(6,1);
normals_1 = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
cps_1 = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];

% Append geometric properties for Solar Panel 1 (2 sides)
Areas_2 = 6*ones(2,1);
normals_2 = [0 0 1; 0 0 -1];
cps_2=[0 2.5 -0.025; 0 2.5 0.025];

% Append geometric properties for Solar Panel 2 (2 sides)
Areas_3 = 6*ones(2,1);
normals_3 = [0 0 1; 0 0 -1];
cps_3 = [0 -2.5 -0.025; 0 -2.5 0.025];

% Append geometric properties for Sensor (5 sides)
Areas_4 = [.25; .25; .25; .25; .0625];
normals_4 = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1];
cps_4 = [0.25 0 1.5;-0.25 0 1.5;0 0.25 1.5;0 -0.25 1.5;0 0 2];

% Now build the matrix
surfaceProperties = [Areas_1 cps_1 normals_1; Areas_2 cps_2 normals_2; Areas_3 cps_3 normals_3; Areas_4 cps_4 normals_4];

% Spacecraft Orbit Properties (given)
global mu
mu = 398600; % km^3/s^2
h = 53335.2; % km^2/s
ecc = 0; % none
raan = 0*pi/180; % radians
inc = 98.43*pi/180; % radians
aop = 0*pi/180; % radians
theta = 0*pi/180; % radians
r = h^2/mu/(1 - ecc^2);
orbital_period = 2*pi*sqrt(r^3/mu);

% Torque free scenario (Given)
T = [0;0;0];
% Set/Compute initial conditions
% intial orbital position and velocity
[rPeri,vPeri,rECI,vECI] = coes2rv(ecc,r,inc,raan,aop,theta,mu);
% Compute inital F_LVLH basis vectors in F_ECI components based on F_LVLH
% definition

z_LVLH = -rECI/norm(rECI);
y_LVLH = -cross(rECI,vECI)/norm(cross(rECI,vECI));
x_LVLH = cross(y_LVLH,z_LVLH);
% Initial Euler angles relating F_body and F_LVLH (given)
phi_0 = 0;
theta_0 = 0;
psi_0 = 0;
E_b_LVLH_0 = [phi_0; theta_0; psi_0];

% Initial Quaternion relating F_body and F_LVLH (given)
q_b_LVLH_0 = [0; 0; 0; 1];

% Compute initial C_LVLH_ECI_0, C_b_LHVL_0, and C_b_ECI_0 rotaiton matrices
C_b_ECI_0 = [x_LVLH y_LVLH z_LVLH]';
C_b_LVLH_0 = Euler2C(phi_0, theta_0, psi_0);

% Initial Euler angles relating body to ECI
E_b_ECI_0 = C2EulerAngles(C_b_ECI_0);

% Initial quaternion relating body to Eci
% q_b_ECI_0 = C2quat(C_b_ECI_0);

% % Initial body rates of spacecraft (given)
% w_b_ECI_0 = [0.01; -0.05; 0.05];

%% Finding Gain
close all;
Mp = 0.02;    % Max overshoot
ts = 100; % Settle time
zeta= 0.65;

wd = 4.4/(zeta*ts);
wn = wd/sqrt(1-zeta^2);
% wn = 4 / (zeta * ts);
% beta = atan(sqrt((1 - zeta^2)/zeta));

Kp = 2*J*wn^2.*[1 1 1];
Kd = J*2*zeta*wn.*[1 1 1];

epsilon_C = [0.2; -0.5; 0.3];
% q_C = [epsilon_C; sqrt(1-norm(epsilon_C)^2)];
% n0 = sqrt(1- epsilon_C' *epsilon_C)
% p = -3*wn + wn*sqrt(zeta^2 -1) *1i
epsilon_b_ECI_0 = [-.2; .4; .2];
q_b_ECI_0 = [epsilon_b_ECI_0; sqrt(1-norm(epsilon_b_ECI_0)^2)];
% q_C = sqrt(1 - epsilon_C'.*epsilon_C);
w_b_ECI_0 = [.1; -.05; .05];

% epsilon_C = [.1; -.3; .4];
% epsilon_C = [0; 0; 0];
q_C = [epsilon_C; sqrt(1-norm(epsilon_C)^2)]; 

% A = [zeros(3,3), eye(3);
%     -0.5* J\ Kp, -J\Kd];

% State-space representation
A = zeros(3,3);
B = 1\J;
C = eye(3);
D = zeros(3,3);
Q = eye(3);
R = eye(3);
K = lqr(A, B, Q, R);

T_control = [0; 0; 0]*2;
tspan = 150;

% out = sim('ProjectSim');
out = sim('ProjectSim.slx');


%% Plot Results
subplot(2,1,1); hold on;
plot(out.q_b_ECI.time, out.q_b_ECI.signals.values)
% plot(q_b_ECI.time, squeeze(q_b_ECI.signals.values))
% plot(out.tout, out.q_b_ECI(:,2), out.tout, out.q_b_ECI(:,3), out.tout, out.q_b_ECI(:,4), out.tout, out.q_b_ECI(:,5));
title('Quaternions');
ylabel('Quat Angles');
legend("q0", "q1", "q2", "q3");
hold off; grid on;

subplot(2,1,2); hold on;
plot(out.omega.time, out.omega.signals.values);
% plot(omega.time, omega.signals.values)
% plot(out.tout, out.omega(:,2), out.tout, out.w_b_ECI(:,3), out.tout, out.w_b_ECI(:,4));
title('Angular Velocities');
ylabel('Angular Velocity (rad/sec)');
legend("wx", "wy", "wz");
hold off; grid on;

% subplot(3,1,3); hold on;
% plot(out.tout, out.E_b_ECI(1:140,2),out.tout, out.E_b_ECI(1,3));
% title('Euler Angles');
% xlabel('Time (sec)'); ylabel('Angle (deg)');
% legend("phi", "theta", "psi");
% hold off; grid on;

figure; hold on; grid on;
plot(out.T_Control.time, out.T_Control.signals.values(1,:));
plot(out.T_Control.time, out.T_Control.signals.values(2,:));
plot(out.T_Control.time, out.T_Control.signals.values(3,:));
title('Torque Commands');
xlabel('Time (sec)'); ylabel('Torque');
legend("Tx", "Ty", "Tz");
hold off;
%% Functions Used

function [C] = principleRotations(inp,angle)
% Principle Rotation function returns Cx, Cy, Cz
if inp == 'x'
 C = [1 0 0;
 0 cos(angle) sin(angle);
 0 -sin(angle) cos(angle)];
elseif inp == 'y'
 C = [cos(angle) 0 -sin(angle);
 0 1 0;
 sin(angle) 0 cos(angle)];
elseif inp == 'z'
 C = [cos(angle) sin(angle) 0;
 -sin(angle) cos(angle) 0;
 0 0 1];
else
 error('Error: Input character must be x, y, or z')
end
end

function [rPeri,vPeri,rECI,vECI] = coes2rv(ecc,r,inc,raan,aop,theta,mu)
% ecc = eccentricity
% r = radius, km
% inc = inclination, degrees
% raan = RAAN, degrees
% aop = argument of perigee, degrees
% theta = true anomaly, degrees
% Put angles in radians for calculations:
aop = deg2rad(aop);
inc = deg2rad(inc);
raan = deg2rad(raan);
theta = deg2rad(theta);
h = sqrt(mu*(1+ecc*cos(theta))*r);
rPeri = h^2/(mu*(1+ecc*cos(theta)))*[cos(theta);sin(theta);0]; %km, radius
vPeri = (mu/h)*[-sin(theta);ecc+cos(theta);0]; %km/s, velocity
R1 = principleRotations('z',aop);
R2 = principleRotations('x',inc);
R3 = principleRotations('z',raan);
C_ECI_peri = R1*R2*R3;
rECI = C_ECI_peri'*rPeri; %km, radius in ECI
vECI = C_ECI_peri'*vPeri; %km/s, velocity in ECI
end

function euler = C2EulerAngles(C)
phi = atan(C(2,3)/C(3,3));
theta = -asin(C(1,3));
psi = atan(C(1,2)/C(1,1));
euler = [phi; theta; psi];
end

function quat = C2quat(C)
 % Rotation matrix to quaternion.
 % Outputs:
 % q = [qw, qx, qy, qz]
 if ~all(size(C) == [3, 3])
 error('Input a 3x3 matrix. This aint it');
 end
 q = zeros(1, 4);
 q(1) = 0.5*sqrt(1 + C(1,1) + C(2,2) + C(3,3));
 q(2) = (C(2,3)-C(3,2)) / (4*q(1));
 q(3) = (C(3,1)- C(1,3)) / (4*q(1));
 q(4) = (C(1,2)-C(2,1)) / (4*q(1));
 quat = q / norm(q);
end

function C = Euler2C(yaw, pitch, roll)
 % Rotation matrix C from Euler angles (yaw, pitch, roll)
 % yaw: Rotation about the Z-axis in degrees
 % pitch: Rotation about the Y-axis in degrees
 % roll : Rotation about the X-axis in degrees
 % Output:
 % C:rotation matrix
 yaw_rad = deg2rad(yaw);
 pitch_rad = deg2rad(pitch);
 roll_rad = deg2rad(roll);
 cy = cos(yaw_rad);
 sy = sin(yaw_rad);
 cp = cos(pitch_rad);
 sp = sin(pitch_rad);
 cr = cos(roll_rad);
 sr = sin(roll_rad);
 % Rotation Matrix
 C = [cy*cp, cy*sp*sr - sy*cr, cy*sp*cr + sy*sr;
 sy*cp, sy*sp*sr + cy*cr, sy*sp*cr - cy*sr;
 -sp, cp*sr, cp*cr];
end
