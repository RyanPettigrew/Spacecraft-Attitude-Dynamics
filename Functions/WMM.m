function Tb_b = WMM(r_ECI, q_b_ECI, JD_0, t)
% Function of World Magn. Model

% DGRF 2020 (given)
a = 6371.2;      %km
g_11 = -1450.9;  %nT
h_11 = 4652.5;   %nT
g_10 = -29404.8; %nT
m_ECEF = a^3*[g_11; h_11; g_10];
m_b = [0;0;-0.5];

I = eye(3);    % Identity

eps1 = q_b_ECI(1);
eps2 = q_b_ECI(2);
eps3 = q_b_ECI(3);
eps = [eps1; eps2; eps3];

eda = q_b_ECI(4);

eps_cross = [0,  -eps3, eps2;
             eps3,  0, -eps1;
             -eps2, eps1, 0];

C_b_ECI = (2*eda^2 - 1)*I + 2*(eps*eps') - 2*eda*eps_cross; % rotation w/ quat

r_ECEF = eci2ecef_manual(JD_0, r_ECI);
C_ECI_ECEF = r_ECI*r_ECEF';

% Solve
B_ECEF = (3*(m_ECEF'*r_ECEF)*r_ECEF - norm(r_ECEF)^2*m_ECEF)/(norm(r_ECEF))^5; % Mag. field vect in ECEF

B_ECI = C_ECI_ECEF*B_ECEF; 

B_b = C_b_ECI*B_ECI;

Tb_b = cross(m_b,B_b);

end

function r_ECEF = eci2ecef_manual(JD, r_ECI)
%     omega_earth = 7.2921159e-5; 
    T_UT1 = (JD - 2451545.0)/36525; % Time in Julian centuries from J2000
    GMST = mod(67310.54841 + (876600*3600 + 8640184.812866)*T_UT1 + 0.093104*T_UT1^2 - 6.2e-6*T_UT1^3, 86400) / 240; % in degrees
    GMST = deg2rad(GMST); 

    % Rotation matrix -ECI to ECEF
    R_ECI2ECEF = [cos(GMST) sin(GMST) 0; -sin(GMST) cos(GMST) 0; 0 0 1];
    r_ECEF = R_ECI2ECEF * r_ECI;
end
