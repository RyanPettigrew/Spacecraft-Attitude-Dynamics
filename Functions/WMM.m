function Tb_b = WMM(m_b, r_ECI, q_b_ECI, JD_0, t)
% Function of World Magn. Model

% DGRF 2020 (given)
a = 6371.2;      %km
g_11 = -1450.9;  %nT
h_11 = 4652.5;   %nT
g_10 = -29404.8; %nT
m_ECEF = a^3*[g_11; h_11; g_10];

I = eye(3);    % Identity

% Unpack Quat
eps1 = q_b_ECI(1);
eps2 = q_b_ECI(2);
eps3 = q_b_ECI(3);
eps = [eps1; eps2; eps3];

eda = q_b_ECI(4);

% Solve
eps_cross = [0,  -eps3, eps2;
             eps3,  0, -eps1;
             -eps2, eps1, 0];

C_b_ECI = (2*eda^2 - 1)*I + 2*(eps*eps') - 2*eda*eps_cross; % rotation w/ quat

[r_ECEF] = eci2ecef(JD_0, r_ECI);

% Solve
B_ECEF = (3*(m_ECEF'*r_ECEF)*r_ECEF - norm(r_ECEF)^2*m_ECEF)/(norm(r_ECEF))^5; % Mag. field vect in ECEF

B_ECI = B_ECEF; % NEED TO FIND THIS TRANSFORM

B_b = C_b_ECI*B_ECI;

Tb_b = m_b*B_b;

end
