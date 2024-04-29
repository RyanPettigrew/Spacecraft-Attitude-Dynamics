function T_a = aero_drag_torque(r_ECI, v_ECI, q_b_ECI, surfaceProperties) 
C_d = 2.5;
r_earth = 6478*1000;
rho0=1.225;

% Convert v_ECI to v_b
v_b = quatrotate(quatconj(q_b_ECI), v_ECI);

% Calculate density using standard exponetial decay model (or model of your choice)
rho = rho0*exp(-(norm(r_ECI)- r_earth) / 8500);


for i = 1:length(surfaceProperties)
    % Find wetted area
    A_wet = surfaceProperties(i).Areas;
    cp = surfaceProperties(i).cps;
    normal = surfaceProperties(i).normal;
    
    v_normal = dot(v_b, normal)*normal;

    if A_wet > 0
        % Add to total torque from atmospheric drag
        F_drag = -0.5*rho*norm(v_normal)^2 *A_wet*C_d *v_normal;
        T_a = cross(cp, F_drag);
    end
end
