function [r_ECI, v_ECI] = coe2rv(h, e, Omega, inclination, omega, nu, mu)
    % Outputs:
    %   r_ECI - position vector in the ECI frame (km)
    %   v_ECI - velocity vector in the ECI frame (km/s)

    r = (h.^2/mu) .* (1 ./ (1 + e * cos(nu)));
    r_perif = [r * cos(nu); r*sin(nu); 0];

    % Velocity in perifocal frame
    v_perif = (mu / h) * [-sin(nu); e + cos(nu); 0];

    % Rotation matrix: perifocal frame to ECI 
    R3 = [cos(Omega), sin(Omega), 0; -sin(Omega), cos(Omega), 0; 0, 0, 1];
    R1 = [1, 0, 0; 0, cos(inclination), sin(inclination); 0, -sin(inclination), cos(inclination)];
    R2 = [cos(omega), sin(omega), 0; -sin(omega), cos(omega), 0; 0, 0, 1];

    % Rotations
    Q_pX = R3 * (R1*R2);

    % Change pos & vel vectors to the ECI frame
    r_ECI = Q_pX * r_perif;
    v_ECI = Q_pX * v_perif;
end
