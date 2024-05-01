%% Final Simulation Project
% Template for part 3 - Disturbance Modeling
% Aero 421
% Cal Poly, SLO

clear
close all
clc

%% Part 1 - Mass Properties

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

%% Part 2 - Geometric Properties of MehielSat during Normal Operations

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

%% Part 3 - Initialize Simulation States

% Current JD - has to be on the solar equinox, why? - we'll use 3/20/2024
% from https://aa.usno.navy.mil/data/JulianDate
% Need this so we can convert from F_ECEF to F_ECI and to F_b for the
% magnetic field model
JD_0 = 2460390;

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

% No external command Torque
T_c = [0; 0; 0]; % Nm

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

% Initial quaternion relating body to E
q_b_ECI_0 = C2quat(C_b_ECI_0);

% Initial body rates of spacecraft (given)
w_b_ECI_0 = [0.001; -0.001; 0.002];

%% Part 4 - Simulate Results

n_revs = 1; % revs
tspan = n_revs * orbital_period;
out = sim('Final_Project_Part_3_disturbances');

%% Part 5 - Plot Results

% Get variables
time = out.tout;

w = out.w_b_ECI.signals.values;
w1 = w(:,1);
w2 = w(:,2);
w3 = w(:,3);

EA = out.E_b_ECI.signals.values;
phi = (180/pi).*EA(:,1);
theta = (180/pi).*EA(:,2);
psi = (180/pi).*EA(:,3);

q = out.q_b_ECI.signals.values;
q1 = q(:,1);
q2 = q(:,2);
q3 = q(:,3);
eta = q(:,4);

Tg = out.T_gg.signals.values;
Tgx = Tg(:,1);
Tgy = Tg(:,2);
Tgz = Tg(:,3);

Ts = out.T_srp.signals.values;
Tsx = Ts(:,1);
Tsy = Ts(:,2);
Tsz = Ts(:,3);
Ta = out.T_a.signals.values;
Tax = Ta(:,1);
Tay = Ta(:,2);
Taz = Ta(:,3);

Tb = squeeze(out.T_b.signals.values)';
Tbx = Tb(:,1);
Tby = Tb(:,2);
Tbz = Tb(:,3);

% Plot Angular Velocities, Euler Angles and Quaternions
figure(2);
subplot(3,1,1);
plot(time, w1);
hold on;
plot(time, w2);
plot(time, w3);
ylabel('Angular Velocity [rad/s]', 'FontSize',8);
title('Angular Velocities');
legend('\omega_{x}', '\omega_{y}', '\omega_{z}',"Location","northeastoutside");
grid on;

subplot(3,1,2);
plot(time, q1);
hold on;
plot(time, q2);
plot(time, q3);
plot(time, eta);
ylabel('Quaternion Parameter', 'FontSize',8);
title('Quaternions');
legend('q_{1}', 'q_{2}', 'q_{3}', '\eta',"Location","northeastoutside");
grid on;

subplot(3,1,3);
plot(time, phi);
hold on;
plot(time, theta);
plot(time, psi);
xlabel('Time [s]');
ylabel('Euler Angles [deg]', 'FontSize',8);
title('Euler Angles');
legend('\phi', '\theta', '\psi',"Location","northeastoutside");
grid on;

% Plot Disturbance torques in F_b
figure(3);
subplot(2,1,1);
hold on;
grid on;
title('Atmospheric Drag Torque');

plot(time,Tax,time,Tay,time,Taz);
ylabel('Atmospheric Drag Torque (N-m)');
legend('T_{ax}','T_{ay}','T_{az}',"Location","northeastoutside");
subplot(2,1,2);
hold on;
grid on;
title('SRP Torque');

plot(time,Tsx,time,Tsy,time,Tsz);
ylabel('SRP Torque (N-m)');
legend('T_{sx}','T_{sy}','T_{sz}',"Location","northeastoutside");
figure(4);
subplot(2,1,1);
hold on;
grid on;
title('Gravity Gradient Torque');

plot(time,Tgx,time,Tgy,time,Tgz);
ylabel('Gravity Gradient Torque (N-m)');
legend('T_{ggx}','T_{ggy}','T_{ggz}',"Location","northeastoutside");
subplot(2,1,2);
hold on;
grid on;
title('Magnetic Field Torque');

plot(time,Tbx,time,Tby,time,Tbz);
ylabel('Magnetic Field Torque (N-m)');
xlabel('Time [s]');
legend('T_{bx}','T_{by}','T_{bz}',"Location","northeastoutside");


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


