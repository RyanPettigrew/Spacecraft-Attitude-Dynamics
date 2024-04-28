function T_a = aero_drag_torque(r_ECI, v_ECI, q_b_ECI, surfaceProperties) 

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

% Convert v_ECI to v_b
v_b = C_b_ECI*v_ECI;

% Calculate density using standard exponetial decay model (or model of your
% choice)
%%% NEEDS MORE WORK

% Loop through all the surfaces of the MehielSat
for i = 1:length(surfaceProperties)
    % Find wetted area

    if A_wet > 0
        % Add to total torque from atmospheric drag
    end
end

T_a = c_pa * F_a;

end
