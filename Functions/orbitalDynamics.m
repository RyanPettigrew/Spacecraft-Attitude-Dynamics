function [dr,dv] = orbitalDynamics(r, v)
muE = 398600; % km^3/s^2, Mu Earth
rx = r(1);
ry = r(2);
rz = r(3);
vx = v(1); 
vy = v(2); 
vz = v(3); 

r_norm = norm([rx ry rz]);

dvx = -muE*rx / r_norm^3;
dvy = -muE*ry / r_norm^3;
dvz = -muE*rz / r_norm^3;

dr = [vx; vy; vz];
dv = [dvx; dvy; dvz];
end
