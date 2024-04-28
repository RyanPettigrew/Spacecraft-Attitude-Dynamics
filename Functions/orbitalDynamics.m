function [dr,dv] = orbitalDynamics(r, v)
%#codegen
muE = 398600;      %km^3/s^2, Mu Earth

rx = r(1);
ry = r(2);
rz = r(3);

vx = v(4);
vy = v(5);
vz = v(6);

r = norm([rx ry rz]);

dvx = -muE*rx/r^3;
dvy = -muE*ry/r^3;
dvz = -muE*rz/r^3;

dr = [vx; vy; vz];
dv = [dvx; dvy; dvz];

end
