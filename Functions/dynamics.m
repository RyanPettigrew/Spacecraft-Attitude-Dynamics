function domega = dynamics(T, J, omega)
%#codegen

% Compute the derivative of the angular velocity $\omega$ here
domega = J \ (T - cross(-omega, J*omega));

end
