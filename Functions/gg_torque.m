function T_gg = gg_torque(r_ECI, q_b_ECI, J)
% Function that calculates gravity gradient torque

% Constant(s)
muE = 398600;  % km^3/s^2, Mu earth
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

r_b = C_b_ECI*r_ECI;

r_b_cross = [0 -r_b(3,1) r_b(2,1);
           r_b(3,1) 0 -r_b(1,1);
           -r_b(2,1) r_b(1,1) 0];

T_gg = (3*muE)/(norm(r_ECI)^5)*r_b_cross*J*r_ECI;

end
