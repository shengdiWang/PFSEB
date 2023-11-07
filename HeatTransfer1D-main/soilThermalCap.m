%
function soilCap = soilThermalCap(thetaU, thetaI, thetaS, thetaA, cU, cS)
% soil thermal capacity as Eq. 17

cI = 2.05;
cA = 0.001297;


soilCap = thetaU * cU + thetaI * cI + thetaS * cS + thetaA * cA; 
                        
end