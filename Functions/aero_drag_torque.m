function T_a = aero_drag_torque(r_ECI, v_ECI, q_b_ECI, surfaceProperties) 
% Constants
C_d = 2.5;
r_earth = 6478 * 1000;
rho0 = 1.225;
I = eye(3);    % Identity

% Extract quaternion components
eps = q_b_ECI(1:3);
eda = q_b_ECI(4);

% Quaternion to rotation matrix
eps_cross = [0, -eps(3), eps(2);
             eps(3), 0, -eps(1);
            -eps(2), eps(1), 0];
C_b_ECI = (2 * eda^2 - 1) * I + 2 * (eps * eps') - 2 * eda * eps_cross;

s_ECI = [1;0;0];

% Convert velocity from ECI to body frame
v_b = C_b_ECI * v_ECI;
s_b = C_b_ECI * s_ECI;

% Calculate atmospheric density
rho = rho0 * exp(-(norm(r_ECI*1000) - r_earth) / 8500); %NEED TO CHANGE VALUE IS TOO SMALL

% Initialize torque
T_a = [0; 0; 0];

% Iterate over surface properties
for i = 1:size(surfaceProperties, 1)
    % Extract surface properties
    S = surfaceProperties(i, 1);
    cp = surfaceProperties(i, 2:4)';
    normal = surfaceProperties(i, 5:7);

    % Calculate relative velocity component normal to the surface
    A_wet = S*dot(s_b, normal);

    if A_wet > 0
        % Calculate drag force
        F_drag = -0.5 * rho * norm(v_b*1000)^2 * A_wet * C_d * v_b/norm(v_b);

        % Add torque contribution
        T_a = T_a + cross(cp, F_drag);
    end
end

end


