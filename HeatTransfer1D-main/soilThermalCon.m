function soilK = soilThermalCon(kSolid, kWater,...
                                thetaIce, thetaWater, thetaSoild, thetaAir)

kIce = 2.29;   % ice thermal conductivity
kAir = 0.025;  % air thermal conductivity

soilK = kIce ^ 0.5 * thetaIce + kSolid ^ 0.5 * thetaSoild + kWater ^ 0.5 * thetaWater + kAir ^ 0.5 * thetaAir;
soilK = soilK ^ 2;
                        
end