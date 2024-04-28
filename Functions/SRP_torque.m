function T_srp  = SRP_torque(sun_ECI, q_b_ECI, surfaceProperties)

p_sun = 4.5e-6; %N/m (given)


% Unpack Quat
eps1 = q_b_ECI(1);
eps2 = q_b_ECI(2);
eps3 = q_b_ECI(3);
eps = [eps1; eps2; eps3];

eda = q_b_ECI(4);

% Solve
A = surfaceProperties.Areas; % m^2

eps_cross = [0,  -eps3, eps2;
             eps3,  0, -eps1;
             -eps2, eps1, 0];

C_b_ECI = (2*eda^2 - 1)*I + 2*(eps*eps') - 2*eda*eps_cross; % rotation w/ quat

sun_b = C_b_ECI*sun_ECI;
sun_direction = sun_ECI / norm(sun_ECI);

T_srp = p_sun * A * cross(sun_b, sun_direction);

end
