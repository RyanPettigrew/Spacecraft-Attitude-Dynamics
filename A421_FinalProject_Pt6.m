%% Final Simulation Project
% Part 5 - Detumble Simulation
% Aero 421

clc; clear; close all;

% Mass properties for normal operations phase
m = 640; %kg, mass of s/c
rvec_cm = [0; 0; 0.2344]; %m, position vector of s/c
J = [812.04 0 0; 0 545.37 0; 0 0 627.71]; %kg*m^2, inertia matrix of s/c

% Calculate the total mass, inertia and Center of Mass of the MehielSat
total_mass = m;
cm = rvec_cm;
J_sc = J;

fprintf('The spacecraft mass for the normal operations mode is:\n');
disp(total_mass);
fprintf('The spacecraft center of mass for the normal operations mode is:\n');
disp(cm);
fprintf('The Inertia Matrix for the normal operations mode is:\n')
disp(J_sc);

%% Geometric Properties of MehielSat during Normal Operations

% Define the sun in F_ECI and residual dipole moment in F_b
sun_ECI = [1;0;0];

% I constructed a matrix where the rows represent each surface of the
% MehielSat. The first column stores the Area of the surface, the next
% three columns define the normal vector of that surface in F_b, and the
% final three columns store the center of pressure of the surface (geometric
% center of the surface) in F_b.

% First get rho vectors with respect to the center of the spacecraft bus
% the MehielSat BUS is a box
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

%% Initialize Simulation States

% Current JD - has to be on the solar equinox, why? - we'll use 3/20/2024
% from https://aa.usno.navy.mil/data/JulianDate
% Need this so we can convert from F_ECEF to F_ECI and to F_b for the
% magnetic field model
JD_0 = 2460390;
m_b = [0;0;-0.5];

% Spacecraft Orbit Properties
mu = 398600;        % km^3/s^2
h = 53335.2;        % km^2/s
ecc = 0;            % none
raan = 0*pi/180;    % radians
inc = 98.43*pi/180; % radians
aop = 0*pi/180;     % radians
theta = 0*pi/180;   % radians

a = h^2/mu/(1 - ecc^2);
orbital_period = 2*pi*sqrt(a^3/mu);

% Set/Compute initial conditions
% intial orbital position and velocity
[rPeri,vPeri,rECI,vECI] = coes2rv(ecc,a,inc,raan,aop,theta,mu);
r_ECI_0 = rECI;
v_ECI_0 = vECI;

% No external command Torque and no disturbance torque
T = [0; 0; 0];     % Nm
T_c = [0; 0; 0]; % Nm
T_d = [0; 0; 0]; % Nm

% Compute inital F_LVLH basis vectors in F_ECI components based on F_LVLH
% definition
z_LVLH = -r_ECI_0/norm(r_ECI_0);
y_LVLH = -cross(r_ECI_0, v_ECI_0)/norm(cross(r_ECI_0, v_ECI_0));
x_LVLH = cross(y_LVLH, z_LVLH);

% Initial Euler angles relating F_body and F_LVLH (given)
phi_0 = 0;
theta_0 = 0;
psi_0 = 0;
E_b_LVLH_0 = [phi_0; theta_0; psi_0];

% Initial Quaternion relating F_body and F_LVLH (given)
q_b_LVLH_0 = [0; 0; 0; 1];

% Compute initial C_LVLH_ECI_0, C_b_LHVL_0, and C_b_ECI_0 rotaiton matrices
C_LVLH_ECI_0 = [x_LVLH'; y_LVLH'; z_LVLH'];
Cx = principleRotations('x',phi_0);
Cy = principleRotations('y',theta_0);
Cz = principleRotations('z',psi_0);
C_b_LVLH_0 = Cx*Cy*Cz;
C_b_ECI_0 = C_b_LVLH_0*C_LVLH_ECI_0;

% Initial Euler angles relating body to ECI
E_b_ECI_0 = C2EulerAngles(C_b_ECI_0);

% Initial quaternion relating body to E
q_b_ECI_0 = C2quat(C_b_ECI_0);

% Initial body rates of spacecraft (given)
w_b_ECI_0 = [0.001; -0.001; 0.002];

%% Finding Gain
Mp = 0.02;    % Max overshoot
ts = 100;     % Settle time

%% Single Axis Analysis
Mp_reg = .1;
% tr < 12;

zeta = 0.65;

wn = log(0.02*sqrt(1-zeta^2))/-zeta/ts;

beta = atan(sqrt(1-zeta^2)/zeta);
tr = (pi-beta)/wn/sqrt(1-zeta^2);

% Extend to each Axis
Kp = 2*J*wn^2.;
Kd = J*2*zeta*wn;

epsilon_b_ECI_0 = [-.2; .4; .2];
q_b_ECI_0 = [epsilon_b_ECI_0; sqrt(1-norm(epsilon_b_ECI_0)^2)];
w_0 = [.1; -.1; .2];

%epsilon_C = [.1; -.3; .4];
epsilon_C = [0; 0; 0];
q_C = [epsilon_C; sqrt(1-norm(epsilon_C)^2)]; 

%% Reaction Wheel Properties
m_w = 1;   % kg
I_s = 1.2; % kg/m2
I_t = 0.6; % kg/m2

I_tot = J_sc + (I_s + 2*I_t + 2*m_w)*eye(3);

% Is_inv = inv([1.2 0 0; 0 1.2 0; 0 0 1.2]);
Is = [1.2 0 0; 0 0.6 0; 0 0 1.2];

%% Simulate Results

% n_revs = 5; % revs
% tspan = n_revs * orbital_period;
tspan = 30000;
out = sim('Final_Project_Pt6_New.slx');  %CHANGE THIS TO MATCH YOUR SIM

%% Plot Results

subplot(3,1,1)
plot(out.tout, out.w_b_ECI.signals.values)
ylabel('angular velocity (rad/sec)')
legend('\omega_x','\omega_y','\omega_z')
grid on

subplot(3,1,2)
plot(out.tout, out.q_b_ECI.signals.values)
ylabel('Quaternion Parameter')
legend('\eta','\epsilon_1','\epsilon_2','\epsilon_3')
grid on

E_b_ECI = squeeze(out.E_b_ECI.signals.values);
subplot(3,1,3)
plot(out.tout, out.E_b_ECI.signals.values)
xlabel('time (seconds)')
ylabel('Angle (deg)')
legend('\phi','\theta','\psi')
grid on
sgtitle('Body to ECI')

figure
subplot(2,2,1)
plot(out.tout, out.T_a.signals.values)
ylabel('Atmospheric Drag Torque (N-m)')
legend('T_{ax}','T_{ay}','T_{az}')
grid on

subplot(2,2,2)
plot(out.tout, out.T_srp.signals.values)
ylabel('SRP Torque (N-m)')
legend('T_{sx}','T_{sy}','T_{sz}')
grid on

subplot(2,2,3)
plot(out.tout, out.T_gg.signals.values)
ylabel('Gravity Gradient Torque (N-m)')
legend('T_{ggx}','T_{ggy}','T_{ggz}')
grid on

subplot(2,2,4)
plot(out.tout, squeeze(out.T_b.signals.values))
ylabel('Magnetic Field Torque (N-m)')
xlabel('Time (sec)')
legend('T_{bx}','T_{by}','T_{bz}')
grid on
sgtitle('Disturbance Torques')

figure
subplot(3,1,1)
plot(out.tout, squeeze(out.T_Control.signals.values))
ylabel('CMD')
legend('X','Y','Z')
grid on

subplot(3,1,2)
plot(out.tout, squeeze(out.Moment_Wheel.signals.values))
ylabel('Moment');
legend('Mx','My','Mz')
grid on

subplot(3,1,3)
plot(out.tout, squeeze(out.reaction_wheel_vel.signals.values))
ylabel('CMD Speeds')
legend('Omega_w1','Omega_w2','Omega_w3')
grid on
sgtitle('Controller Responses')

% figure
% subplot(3,1,1)
% w_b_LVLH = squeeze(out.w_b_LVLH.signals.values);
% plot(out.tout, w_b_LVLH)
% title('Angular Velocities')
% ylabel('angular velocity (rad/sec)')
% legend('\omega_x','\omega_y','\omega_z')
% grid on
% 
% q_b_LVLH = squeeze(out.q_b_LVLH.signals.values);
% subplot(3,1,2)
% plot(out.tout, q_b_LVLH)
% title('Quaternions')
% ylabel('Quaternion Parameter')
% legend('\eta','\epsilon_1','\epsilon_2','\epsilon_3')
% grid on
% 
% E_b_LVLH = squeeze(out.E_b_LVLH.signals.values);
% subplot(3,1,3)
% plot(out.tout, 180/pi*E_b_LVLH)
% title('Euler Angles')
% xlabel('time (seconds)')
% ylabel('Angle (deg)')
% legend('\phi','\theta','\psi')
% grid on
% sgtitle('Body to LVLH')

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

rPeri = h^2/(mu*(1+ecc*cos(theta)))*[cos(theta);sin(theta);0]; %km, radius vector
vPeri = (mu/h)*[-sin(theta);ecc+cos(theta);0];                 %km/s, velocity vector

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
%   yaw: Rotation about the Z-axis in degrees
%   pitch: Rotation about the Y-axis in degrees
%   roll : Rotation about the X-axis in degrees
% Output:
%   C:rotation matrix

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
