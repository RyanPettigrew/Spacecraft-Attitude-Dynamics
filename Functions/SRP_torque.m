function T_srp  = SRP_torque(sun_ECI, q_b_ECI, surfaceProperties)

p_sun = 4.5e-6; %N/m (given)

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

s_b = C_b_ECI*sun_ECI;
s_d = sun_ECI / norm(sun_ECI);

% Solve

% Initialize torque
T_srp = [0; 0; 0];

% Iterate over surface properties
for i = 1:size(surfaceProperties, 1)
    % Extract surface properties
    S = surfaceProperties(i, 1);
    cp = surfaceProperties(i, 2:4)';
    normal = surfaceProperties(i, 5:7);  % Normal is now a row vector

    % Calculate relative velocity component normal to the surface
    A_wet = S*dot(s_b, normal);

    if A_wet > 0
        % Add torque contribution
        T_srp = p_sun * A_wet * cross(s_b, s_d);
    end
end

end
