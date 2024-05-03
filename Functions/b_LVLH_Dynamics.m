function w_b_LVLH = b_LVLH_Dynamics(r_ECI, v_ECI, w_b_ECI, q_b_ECI)
%#codegen

% Compute the angular velocity of LVLH wrt ECI in Body components
I = eye(3);

% Unpack Quat
eps1 = q_b_ECI(1);
eps2 = q_b_ECI(2);
eps3 = q_b_ECI(3);
eps = [eps1; eps2; eps3];

eda = q_b_ECI(4);

eps_cross = [0,  -eps3, eps2;
             eps3,  0, -eps1;
             -eps2, eps1, 0];

C_b_ECI = (2*eda^2 - 1)*I + 2*(eps*eps') - 2*eda*eps_cross; % rotation w/ quat

h = cross(r_ECI,v_ECI);
w_LVLH_ECI = norm(v_ECI)/norm(r_ECI)*h/norm(h);
w_b_LVLH = w_b_ECI - C_b_ECI*w_LVLH_ECI;


end
