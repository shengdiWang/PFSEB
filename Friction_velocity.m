%
function Ustar = Friction_velocity(wsi, wsHeight, RZ, Lstar)

zeta1 = wsHeight / Lstar;
zeta2 = RZ / Lstar;
VONK  = 0.4;

Ustar = (wsi * VONK) / (log(wsHeight / RZ) - psi_M(zeta1, zeta2));


end




