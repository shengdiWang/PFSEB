%
function QH = sensible_heat(roAir, CAir, TA, TS, wsi, RZ, Lstar)

VONK = 0.4;
wsHeight = 2;

QH = - roAir * CAir * VONK^2 * wsi / (log(wsHeight / RZ)- psi_M(wsHeight / Lstar, RZ / Lstar)) * (TA - TS) / (log(wsHeight / RZ)- psi_H(wsHeight / Lstar, RZ / Lstar));

end








