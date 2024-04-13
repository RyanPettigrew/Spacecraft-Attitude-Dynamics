clear; clc; close all;

addpath('Functions\')

%Detumble
m=640; %kg
r_cm = [0;0;0];
J= [426.67 0 0; 0 426.67 0; 0 0 426.67];

%Normal Op
r_cm_normal_op = [0;0;0.2344];
J_normal_op= [812.04 0 0; 0 545.37 0; 0 0 627.71];

%% Part #2 - Torque Free Motion

wb_eci = [0.001 -0.001 0.002]; %rad/s

% Spacecraft Orbit Properties (given)
global mu
mu = 398600; % km^3/s^2
h = 53335.2; % km^2/s
e = 0; % none
Omega = 0*pi/180; % radians
inclination = 98.43*pi/180; % radians
omega = 0*pi/180; % radians
nu = 0*pi/180; % radians

a = h^2/mu/(1 - e^2);
orbital_period = 2*pi*sqrt(a^3/mu);

% Torque free scenario (Given)
T = [0;0;0];

% Set/Compute initial conditions
[r_ECI_0, v_ECI_0] = coe2rv(h, e, Omega, inclination, omega, nu, mu);

% Compute initial C_LVLH_ECI_0, C_b_LHVL_0, and C_b_ECI_0 rotaiton matrices
z_LVLH = -r_ECI_0 / norm(r_ECI_0);
y_LVLH = cross(r_ECI_0, v_ECI_0) / norm(cross(r_ECI_0, v_ECI_0));
x_LVLH = cross(y_LVLH, z_LVLH);
C_LVLH_ECI_0 = [x_LVLH, y_LVLH, z_LVLH]';

% Initial Euler angles relating F_body and F_LVLH (given)
phi_0 = 0;
theta_0 = 0;
psi_0 = 0;
E_b_LVLH_0 = [phi_0; theta_0; psi_0];

% Initial Quaternion relating F_body and F_LVLH (given)
q_b_LVLH_0 = [0; 0; 0; 1];

C_b_ECI_0 = Euler2C(phi_0, theta_0, psi_0);

% Initial Euler angles relating body to ECI
E_b_ECI_0 = C2EulerAngles(C_b_ECI_0);

% Initial quaternion relating body to E
q_b_ECI_0 = C2quat(C_b_ECI_0);

% Initial body rates of spacecraft (given)
w_b_ECI_0 = [0.001; -0.001; 0.002];

tspan = orbital_period;

% out = sim('Project_Solutions_Part2');
tspan = [0 orbital_period]; % From 0 to one orbital period
disp('Cb-eci')
disp(C_LVLH_ECI_0)


%% Plot Results


%% Functions
function [r_ECI, v_ECI] = coe2rv(h, e, Omega, inclination, omega, nu, mu)
    % Outputs:
    %   r_ECI - position vector in the ECI frame (km)
    %   v_ECI - velocity vector in the ECI frame (km/s)

    r = (h.^2 / mu) .* (1 ./ (1 + e * cos(nu)));
    r_pf = [r * cos(nu); r * sin(nu); 0];

    % Velocity in perifocal frame
    v_perif = (mu/h) * [-sin(nu); e + cos(nu); 0];

    % Rotation matrix perifocal frame to ECI 
    R3 = [cos(Omega), sin(Omega), 0; -sin(Omega), cos(Omega), 0; 0, 0, 1];
    R1 = [1, 0, 0; 0, cos(inclination), sin(inclination); 0, -sin(inclination), cos(inclination)];
    R2 = [cos(omega), sin(omega), 0; -sin(omega), cos(omega), 0; 0, 0, 1];

    Q_pX = R3 * (R1*R2);% Rotations

    % Change pos & vel vectors to the ECI frame
    r_ECI = Q_pX*r_pf;
    v_ECI = Q_pX*v_perif;
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

function q = C2quat(C)
    % Rotation matrix to quaternion.
    % Outputs:
    % q = [qw, qx, qy, qz]

    if ~all(size(C) == [3, 3])
        error('Input a 3x3 matrix. This aint it');
    end
    
    q = zeros(1, 4);
    q(1) = 0.5*sqrt(1 + C(1,1) + C(2,2) + C(3,3));
    q(2) = (C(3,2)-C(2,3)) / (4*q(1));
    q(3) = (C(1,3)- C(3,1)) / (4*q(1));
    q(4) = (C(2,1)-C(1,2)) / (4*q(1));
    q = q / norm(q);
end

function [yaw, pitch, roll] = C2EulerAngles(C)
    if ~all(size(C) == [3,3])
        error('Not a 3x3 matrix. This aint it');
    end
    pitch = atan2(-C(3,1), sqrt(C(1,1)^2 + C(2,1)^2));
    roll = atan2(C(3,2), C(3,3));
    yaw = atan2(C(2,1), C(1,1));

    roll = rad2deg(roll);
    pitch = rad2deg(pitch);
    yaw = rad2deg(yaw);
end

function dstate = EulerEOM(t, state, J, M)
 omega = state(1:3);
 euler = state(4:6);

 % Angular velocity rates
 omega_dot = J \ (M - cross(-omega, J*omega));

 % Euler rates
 phi = euler(1);
 theta = euler(2);
 R = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
 0, cos(phi), -sin(phi);
 0, sin(phi)/cos(theta), cos(phi)/cos(theta)];
 euler_dot = R*omea;
 dstate = [omega_dot; euler_dot];
end

function dstate = QuatEOM(t, state, J, M)
 omega = state(1:3);
 quat = state(4:7);
 % Angular velocity rates
 omega_dot = J \ (M - cross(omega, J*omega));
 % Quaternion rates
 q = 0.5*[0, -omega(1), -omega(2), -omega(3);
 omega(1), 0, omega(3), -omega(2);
 omega(2), -omega(3), 0, omega(1);
 omega(3), omega(2), -omega(1), 0];
 quat_dot = q * quat;
 dstate = [omega_dot; quat_dot];
end
