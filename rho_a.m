function res = rho_a(Tz, p)

    R_a = 287.058;                                   % specific gas constant of air [ J/(kg K) ]
    
    res = p / (R_a*Tz);                               % air density [kg m^(-3)]

end