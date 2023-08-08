%
function QE = latent_heat(roAir, SH, ES0, PA, wsi, RZ, Lstar, root)

VONK = 0.4;
wsHeight = 2;
L_w = 1000 * (2500.8 - 2.36 * (root - 273.15));   % [J/kg] latent heat of evaporation of water
L_i = 1000 * 2835;                              % [J/kg] latent heat of sublimation

if root > 273.15
    QE = - roAir * L_w * VONK^2 * wsi / (log(wsHeight / RZ)- psi_M(wsHeight / Lstar, RZ / Lstar)) * (SH - 0.622 * ES0 / PA) ...
       / (log(wsHeight / RZ)- psi_H(wsHeight / Lstar, RZ / Lstar)) + 50 * VONK^2 * wsi / (log(wsHeight / RZ)- psi_M(wsHeight / Lstar, RZ / Lstar));
else
    QE = - roAir * L_i * VONK^2 * wsi / (log(wsHeight / RZ)- psi_M(wsHeight / Lstar, RZ / Lstar)) * (SH - 0.622 * ES0 / PA) ...
       / (log(wsHeight / RZ)- psi_H(wsHeight / Lstar, RZ / Lstar));
end










