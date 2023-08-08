clc;
clear;

%% ========================================================================
%% ====  READ ME  ====
% A numerical model for surface energy balance and thermal regime of the
%   active layer and permafrost containing unfrozen water
%   paper: Https://linkinghub.elsevier.com/retrieve/pii/S0165232X03000570
% changes
%   [1] snow thermal conductivity and capacity are changed as HTESSEL
%   [2] constant values are revised as WIKI

%% ===== SETTING-UP =======================================================
% ---- directory ----
dir_mod = '/Users/shengdiwang/Desktop/JRA/SIM';
dir_foc = '/Users/shengdiwang/Desktop/JRA/FOC';
dir_var = '/Users/shengdiwang/Desktop/JRA/VAR';

% ---- metafile ----
metaf   = 'meta0.csv';
metaf   = fullfile(dir_var, metaf);
meta    = readtable(metaf, VariableNamingRule="preserve");
scSites = meta.site;

% ---- soil ----
% 5-layer snow, 16 layers for 0-1 m, 190 layers for 1-20 m [0.1 m interval]
% soil discretization for 0-1 m
soil1m  = [0, 3, 7, 11, 15, 22, 29, 35, 44, 50, 58, 66, 75, 84, 93, 100];
soil1m  = soil1m / 100; % in m
soilTck = 0.1; % thickness [m] for deep soil (> 1 m)

snowN     = 5;   % maximum snow layer
SLTck1    = 1; % soil layer 1 thickness (Table 1)
SLTck2    = 9; % soil layer 2 thickness (Table 1)
soilDepth = 10;  % total soil depth
NODE      = snowN + length(soil1m) + (soilDepth - 1)/soilTck;

RKKK(7,1) = 2.92 * 24.0 * 3600.0; % 0.564 is snow volumetric heat capacity

% soil node depth
for I = snowN + 1:NODE
    
    if (I <= snowN + length(soil1m)) % 5-layer snow & 0-1 m soil
        
        XYN(I,1) = soil1m(I - snowN);
        
    else % soil profile > 1 m
        
        XYN(I,1) = 1.0 + (I - snowN - length(soil1m)) * soilTck;
        
    end
end

% soil thickness
for I = snowN + 1:NODE-1
    
    DX(I,1) = XYN(I + 1) - XYN(I);
    
end


% ---- model run ----

% spin up
spinTime = 20;  % spin-up times

%% Call constants to calculate thermal parameters

Tzero   = 0; % frozen temperature
QQ      = 0; % the lower boundary condition

%% Assign all parameters needed in SURFACE ENERGY BALANCE METHOD
% AIR
CAir    = 1005.7;  % air thermal capacity [J m-3 C-1]

% WATER
roWater = 1000.0; % water density

% ICE
roIce = 920;   % ice density
cIce  = 2.05;  % ice volumetric heat capacity [MJ m-3 K-1]
kIce  = 2.29;  % thermal conductivity of ice  [W m-1 K-1]

% CONSTANT
GRAVIT     = 9.81;   % gravitational acceleration
VONK       = 0.4;     % Von Karman's constant
TF         = 273.15;  % unit C to K
latSub_i   = 2.838E6; % latent heat of ice sublimation [J kg-1]
latSub_w   = 2.5E6;   % latent heat of water [J kg-1]
wsHeight   = 2.0;     % REFERENCE HEIGHT for wind speed measure
albedoG    = 0.2;     % snow-free albedo
SIGMA      = 5.67e-8;

%% ---- Energy balance ----
kSL1 = 2.92;
kSL2 = 2.92;      % silt thermal conductivity

for si = 1:length(scSites)
% for si = 1

    clay = meta.clay_t(si);
    sand = meta.sand_t(si);
    bvsilt   = vb(clay);    % clay
    thaosilt = vthao(sand); % sand

    siteName = scSites(si);
    sitei = find(meta.site == siteName);
    
    % simulation period
    begYr = meta.beg(sitei);
    begDate1   = datetime(begYr, 1, 1);    % start date for run
    endDate1   = datetime(begYr + 5, 12, 31); % end date for run
    dateRange1 = begDate1:endDate1;
    
    endYr = meta.end(sitei);
    begDate2   = datetime(begYr, 1, 1); % start date for run
    endDate2   = datetime(endYr, 12, 31); % end date for run
    dateRange2 = begDate2:endDate2;

    %% ===== ENSEMBLE FORECASTING =============================================
    
    % theta water air & solid
    sm = meta.sm_t;
    dt = 1.0; % phase change

    thetaU1s     = sm(sitei);
    thetaU2s     = 0.2;
    thetaUStd1s  = 0;
    
    satSL1s      = thetaU1s + 0.05; % saturated water content soil layer 1
    satSL2s      = thetaU2s + 0.05; % saturated water content soil layer 2
    
    thetaS1s     = 1 - satSL1s; % solid particle soil layer 1
    thetaS2s     = 1 - satSL2s; % solid particle soil layer 2    
    
    satSL1       = satSL1s;
    satSL2       = satSL2s;

    thetaS1      = thetaS1s;
    thetaS2      = thetaS2s;

    thetaU1      = thetaU1s;
    thetaU2      = thetaU2s;

    thetaUStd    = thetaUStd1s;
    
    thetaU1Range = thetaU1-thetaUStd:0.01:thetaU1+thetaUStd;
    
    %% ===== IMPORT =======================================================
    
    % ---- import forcing ----
    focf        = sprintf('%d.csv', siteName); % forcing file
    outBaseName = 'CMA_'; % modeled outputs
    
    % ---- make file ----
    focf   = fullfile(dir_foc, focf);
    foc    = readtable(focf); % forcing
    
    DATES  = datestr(foc.date); % date
    AIRT   = foc.satFinal; % air temperature
    SND    = foc.snd;    % snow depth
    dewP   = foc.tempD;  % dew-point temperature
    SRin   = foc.sw;     % solar radiation
    TRin   = foc.lwin;   % theraml radiation
    WS     = foc.ws;     % wind speed
    PRE    = foc.Pa;     % pressure
    snowRo = foc.sdn;    % snow density
    emi    = foc.emi;    % emissivity
    rz     = foc.z;      % roughness length
    albe   = foc.albedo; % albedo
    obsGST = foc.tempS;  % observed soil temp at 0.01 m

    snowRo(SND > 0) = 200;
    emi(SND > 0)  = 0.99;
    emi(SND == 0) = 0.97;
    rz(SND > 0)   = 0.0005; % 0.005
    rz(SND == 0)  = 0.001; % 0.015
    
    EES     = 0.0001;
    TDAYS   = 3600.0 * 24.0; % 1 day in second
    RDELTAT = 1.0;
    
    % ASSIGN INITIAL SOIL TEMPERATURE
    intiST  = AIRT(1);  % initial soil temperature [C]
    for I = snowN + 1:NODE %
        RTT(I, 1) = intiST;
    end
    
    % ======= SPIN-UP =====================================================
    
    NTB1  = strmatch(datestr(begDate1), DATES); % forcing index of begDate

    for sdn = 200 % snow density
        
        for thetaU1 = thetaU1 % soil mositure layer 1

            root  = TF + AIRT(1); % initial root
            Lstar = -100000;      % initial Lstar
            
            for yr = 1:spinTime

                fprintf('%6.2f\n', yr);

                for NDAYI = 1:length(dateRange1)
%                 for NDAYI = 1:1
                    
                    dayi = NTB1 + NDAYI - 1;

                    % SNOW LAYER NUMBER AND THICKNESS
                    SNOWH  = SND(dayi);             % snow depth
                    ROSNOW = snowRo(dayi);          % snow density
                    KSNOW  = snowThermalCon(ROSNOW); % snow conductivity
                    CSNOW  = snowThermalCap(ROSNOW); % snow capacity
                    KSNOW  = KSNOW * 3600.0 * 24.0;
                    
                    if (SNOWH > 0)
                        % Albedo Eq. 31
                        albei = snowAlbedo(ROSNOW);
                    else
                        albei = albedoG;
                    end
                    
                    for I = 1:5 % snow layer [1-5]
                        RKKK(I, 1) = KSNOW;
                        RCCC(I, 1) = CSNOW;
                    end

                    %%
                    % snow free, no snow parameters
                    if (SNOWH < 0.003)
                        LMN = 6; % the first node, either for snow or soil
                        
                    elseif (SNOWH <= 0.07)
                        LMN       = 5; % 1-layer snow, node = 5
                        XYN(5)    = -SNOWH;

                        RTT(5)    = AIRT(dayi);
                        DX(LMN)   = SNOWH;
                        
                    elseif (SNOWH <= 0.15)
                        LMN       = 4;
                        XYN(4)    = -SNOWH;
                        XYN(5)    = -SNOWH / 2.0;
                        
                        RTT(4)    = AIRT(dayi);
                        RTT(5)    = AIRT(dayi);
                        DX(LMN)   = SNOWH / 2.0;
                        DX(LMN+1) = SNOWH / 2.0;
                        
                    elseif (SNOWH <= 0.24)
                        LMN = 3;
                        XYN(3)    = -SNOWH;
                        XYN(4)    = -SNOWH * 2.0 / 3.0;
                        XYN(5)    = -SNOWH / 3.0;

                        RTT(3)    = AIRT(dayi);
                        RTT(4)    = AIRT(dayi);
                        RTT(5)    = AIRT(dayi);
                        DX(LMN)   = SNOWH / 3.0;
                        DX(LMN+1) = SNOWH / 3.0;
                        DX(LMN+2) = SNOWH / 3.0;
                        
                    elseif (SNOWH <= 0.35)
                        LMN = 2;
                        XYN(2)    = -SNOWH;
                        XYN(3)    = -SNOWH * 3.0 / 4.0;
                        XYN(4)    = -SNOWH * 2.0 / 4.0;
                        XYN(5)    = -SNOWH / 4.0;

                        RTT(2)    = AIRT(dayi);
                        RTT(3)    = AIRT(dayi);
                        RTT(4)    = AIRT(dayi);
                        RTT(5)    = AIRT(dayi);
                        DX(LMN)   = SNOWH / 4.0;
                        DX(LMN+1) = SNOWH / 4.0;
                        DX(LMN+2) = SNOWH / 4.0;
                        DX(LMN+3) = SNOWH / 4.0;
                        
                    else
                        LMN = 1;
                        XYN(1)    = -SNOWH;
                        XYN(2)    = -SNOWH * 4.0 / 5.0;
                        XYN(3)    = -SNOWH * 3.0 / 5.0;
                        XYN(4)    = -SNOWH * 2.0 / 5.0;
                        XYN(5)    = -SNOWH / 5.0;

                        RTT(1)    = AIRT(dayi);
                        RTT(2)    = AIRT(dayi);
                        RTT(3)    = AIRT(dayi);
                        RTT(4)    = AIRT(dayi);
                        RTT(5)    = AIRT(dayi);
                        DX(LMN)   = SNOWH / 5.0;
                        DX(LMN+1) = SNOWH / 5.0;
                        DX(LMN+2) = SNOWH / 5.0;
                        DX(LMN+3) = SNOWH / 5.0;
                        DX(LMN+4) = SNOWH / 5.0;
                    end
                    
                    %% VARIABLES FOR ESTIMATING THE HEAT CONDUCTED
                    % THROUGH THE SOIL, QC. in SURFACE ENERGY BALANCE MODEL
                    XRTT7   = RTT(7) + TF; %
                    XXYN7   = XYN(7);   % depth of node 7
                    RKKPEAT = RKKK(7);  % thermal conductivity of node 7
                    %%
                    TA     = TF + AIRT(dayi);
                    tDew   = TF + dewP(dayi);
                    QSI    = SRin(dayi);
                    AQSI   = (1 - albei) * QSI;
                    wsi    = WS(dayi);
                    PA     = PRE(dayi);
                    roAir  = rho_a(TA, PA);
                    EA     = atmosphericVaporPressure(tDew);
                    SH     = 0.622 * EA / PA;
                    QLI    = TRin(dayi);
                    EMI    = emi(dayi);
                    RZ     = rz(dayi);
                    %%
                    NCONTRALER = 0;
                    REEE       = 1;
                    %%
                    syms f TS;
                    % CALL SUBQLE
                    Qle   = longWaveOut(emi(dayi), TS);
                    % CALL SUBES0
                    ES0   = atmosphericVaporPressure(TS);
                    % CALL SUBQH
                    QH    = sensible_heat(roAir, CAir, TA, TS, wsi, RZ, Lstar);
                    % CALL SUBQE
                    QE    = latent_heat(roAir, SH, ES0, PA, wsi, RZ, Lstar, root);
                    % CALL SUBQC
                    if (SNOWH < 0.003) % snow free
                        
                        RERKK =  XXYN7 / (RKKPEAT / TDAYS);
                        
                    else
                        
                        RERKK = SNOWH / ...
                            (KSNOW / TDAYS) + XXYN7 / (RKKPEAT / TDAYS);
                        
                    end
                    QC = - (TS - XRTT7) / RERKK;
                    %%
                    f(TS)   = AQSI + QLI - Qle - QH - QE + QC;
                    df(TS)  = diff(f);

                    %%
                    while EES < REEE
                        
                        NCONTRALER = NCONTRALER + 1;
                        RFTS0   = f(root);
                        DXRFTS0 = df(root);

                        TEMRTS0 = root - RFTS0 / DXRFTS0;
                        REEE    = abs(TEMRTS0 - root);
                        root    = round(TEMRTS0, 6);
%                         fprintf('%6.8f\n', REEE);

                    end
                    %%
                    Ustar   = Friction_velocity(wsi, wsHeight, RZ, Lstar);
                    Lstar   = Monin_obokhov(roAir, Ustar, TA, subs(QH, TS, root), subs(QE, TS, root));
                    %%
                    if (SNOWH > 0.003 && root > TF)
                        
                        root = TF;

                        QH   = sensible_heat(roAir, CAir, TA, root, wsi, RZ, Lstar);
                        ES0  = atmosphericVaporPressure(root);
                        QE   = latent_heat(roAir, SH, ES0, PA, wsi, RZ, Lstar, root);
                        Qle  = longWaveOut(emi(dayi), root);
                        QC   = -(root - XRTT7) / RERKK;

                        Ustar   = Friction_velocity(wsi, wsHeight, RZ, Lstar);
                        Lstar   = Monin_obokhov(roAir, Ustar, TA, QH, QE);
                                             
                    end

%                     fprintf('%6.2f\n', root);
%                     fprintf('%6.2f\n', Lstar);
%                     fprintf('%6.2f\n', subs(QH, TS, root));
%                     fprintf('%6.2f\n', subs(QE, TS, root));
%                     
                    %% SOIL THERMAL CONDUCTION
                    %% CALL SUBTHETA  WATER AND UNFROZEN WATER
                    
                    RTT(LMN) = root - TF; % Upper boundary
                    
                    for I = 6:NODE
                        if (XYN(I) <= SLTck1)    % peat layer
                            thetaS(I,1) = thetaS1; % solid
                            thetaW(I,1) = thetaU1;
                            satSL(I,1)  = satSL1;
                            thetaA(I,1) = satSL1 - thetaU1;
                        else
                            thetaS(I,1) = thetaS2;
                            thetaW(I,1) = thetaU2;
                            satSL(I,1)  = satSL2;
                            thetaA(I,1) = satSL2 - thetaU2;
                            
                        end
                        
                        if (RTT(I) < Tzero) % frozen
                            thetaU(I,1)   = min(thetaW(I), ...
                                NiuYang(satSL(I), (RTT(I) + TF), bvsilt, thaosilt));
                            thetaI(I,1)   = satSL(I) - thetaU(I) - thetaA(I);
                            RLL(I,1)      = roWater * ... % Eq. 19
                                (333.2 + 4.955 * RTT(I) + ...
                                0.02987 * RTT(I)^2) / 1000.0;
                            
                            if (RTT(I) < Tzero - dt)
                                DthetaU(I,1) = 0;
                            else
                                DthetaU(I,1)  = thetaI(I);
                            end
                            
                        else % unfrozen
                            thetaU(I,1)  = thetaW(I);
                            thetaI(I,1)  = 0.0;
                            RLL(I,1)     = 0.0;
                            DthetaU(I,1) = 0;
                        end
                    end
                    
                    %%  CALL SUBRCCC
                    for I = 6:NODE
                        
                        % Eq. 22
                        CCWATER = 4.20843 + 1.11362E-1 * RTT(I) + ...
                            5.12142E-3 * (RTT(I))^2 + 9.3482E-5 * (RTT(I))^3;
                        % Eq. 27
                        cSoil = 0.4091 + 5.433E-3 * (TF + RTT(I));
                        
                        RCV = soilThermalCap(thetaU(I), thetaI(I), ...
                            thetaS1, thetaA(I), CCWATER, cSoil);
                        % Eq. 14
                        RCCC(I,1) = (RCV + 334 * DthetaU(I)/dt) * 1.0E6 ;
                        
                    end
                    %% soil thermal conductivity
                    for I = 6:NODE
                        % Eq. 20
                        kWater = 0.11455 + 1.6318E-3 * (TF + RTT(I));
                        if (XYN(I) <= SLTck1)
                            kSolid = kSL1;
                        else
                            kSolid = kSL2;
                        end
                        RKKK(I,1) = soilThermalCon(kSolid, kWater, ...
                            thetaI(I), thetaU(I), thetaS(I), thetaA(I));
                        
                        % Eq. 16
                        RKKK(I,1) = RKKK(I,1) * 24. * 3600.;
                    end
                    %% CALL ABCD
                    
                    RKW = 2.0 * RKKK(LMN) * RKKK(LMN + 1) / (RKKK(LMN) + RKKK(LMN + 1));
                    RKE = 2.0 * RKKK(LMN + 1) * RKKK(LMN + 2) /...
                        (RKKK(LMN+1) + RKKK(LMN + 2));
                    A(LMN + 1) = RKW / DX(LMN);
                    C(LMN + 1) = RKE / DX(LMN + 1);
                    
                    RCW   = RCCC(LMN) * DX(LMN);
                    RCE   = RCCC(LMN + 2) * DX(LMN + 1);
                    RCCCP = (RCW + RCE) / (DX(LMN) + DX(LMN + 1));
                    
                    RDX      = (DX(LMN) + DX(LMN + 1)) / 2;
                    APZERO   = RCCCP * RDX / RDELTAT;
                    B(LMN+1) = -(A(LMN + 1) + C(LMN + 1) + APZERO);
                    D(LMN+1) = -APZERO * RTT(LMN + 1) - A(LMN + 1) * RTT(LMN);
                    %
                    for IK = LMN+2:NODE-1
                        RKW = 2 * RKKK(IK) * RKKK(IK - 1) / (RKKK(IK) + RKKK(IK - 1));
                        RKE = 2 * RKKK(IK) * RKKK(IK + 1) / (RKKK(IK) + RKKK(IK + 1));
                        A(IK) = RKW / DX(IK - 1);
                        C(IK) = RKE / DX(IK);
                        
                        RCW   = RCCC(IK - 1) * DX(IK - 1);
                        RCE   = RCCC(IK + 1) * DX(IK);
                        RCCCP = (RCW + RCE) / (DX(IK - 1) + DX(IK));
                        
                        RDX    = (DX(IK - 1) + DX(IK)) / 2.0;
                        APZERO = RCCCP * RDX / RDELTAT;
                        B(IK)  = -(A(IK) + C(IK) + APZERO);
                        D(IK)  = -APZERO * RTT(IK);
                    end
                    %
                    RKW     = RKKK(NODE - 1);
                    A(NODE) = RKW / DX(NODE - 1);
                    C(NODE) = 0.0;
                    B(NODE) = -A(NODE);
                    D(NODE) = -QQ;
                    
                    %% CALL TRAM
                    
                    RP(LMN+1) = -C(LMN + 1) / B(LMN + 1);
                    RQ(LMN+1) = D(LMN + 1) / B(LMN + 1);
                    
                    for I = LMN + 2:NODE
                        PP    = A(I) * RP(I - 1)+B(I);
                        RP(I) = -C(I)/PP;
                        RQQ   = D(I) - A(I) * RQ(I - 1);
                        RQ(I) = RQQ / PP;
                    end
                    
                    DTT(NODE,1) = RQ(NODE);
                    for I = NODE-1:-1:LMN+1
                        DTT(I) = RP(I) * DTT(I + 1) + RQ(I);
                    end
                    
                    for I = LMN + 1:NODE
                        RTT(I) = DTT(I);
                    end
                    
                end
            end
            
            % ======= SIMULATIONS =============================================
            
            % forcing index of begDate
            NTB2 = strmatch(datestr(begDate2), DATES);

            for NDAYI = 1:length(dateRange2)
                
                dayi = NTB1 + NDAYI - 1;

                % SNOW LAYER NUMBER AND THICKNESS
                SNOWH  = SND(dayi);             % snow depth
                ROSNOW = snowRo(dayi);          % snow density
                KSNOW  = snowThermalCon(ROSNOW); % snow conductivity
                CSNOW  = snowThermalCap(ROSNOW); % snow capacity
                KSNOW  = KSNOW * 3600.0 * 24.0;
                
                if (SNOWH > 0)
                    % Albedo Eq. 31
                    albei = snowAlbedo(ROSNOW);
                else
                    albei = albedoG;
                end
                
                for I = 1:5 % snow layer [1-5]
                    RKKK(I, 1) = KSNOW;
                    RCCC(I, 1) = CSNOW;
                end

                %%
                % snow free, no snow parameters
                if (SNOWH < 0.003)
                    LMN = 6; % the first node, either for snow or soil
                    
                elseif (SNOWH <= 0.07)
                    LMN       = 5; % 1-layer snow, node = 5
                    XYN(5)    = -SNOWH;

                    RTT(5)    = AIRT(dayi);
                    DX(LMN)   = SNOWH;
                    
                elseif (SNOWH <= 0.15)
                    LMN       = 4;
                    XYN(4)    = -SNOWH;
                    XYN(5)    = -SNOWH / 2.0;
                    
                    RTT(4)    = AIRT(dayi);
                    RTT(5)    = AIRT(dayi);
                    DX(LMN)   = SNOWH / 2.0;
                    DX(LMN+1) = SNOWH / 2.0;
                    
                elseif (SNOWH <= 0.24)
                    LMN = 3;
                    XYN(3)    = -SNOWH;
                    XYN(4)    = -SNOWH * 2.0 / 3.0;
                    XYN(5)    = -SNOWH / 3.0;

                    RTT(3)    = AIRT(dayi);
                    RTT(4)    = AIRT(dayi);
                    RTT(5)    = AIRT(dayi);
                    DX(LMN)   = SNOWH / 3.0;
                    DX(LMN+1) = SNOWH / 3.0;
                    DX(LMN+2) = SNOWH / 3.0;
                    
                elseif (SNOWH <= 0.35)
                    LMN = 2;
                    XYN(2)    = -SNOWH;
                    XYN(3)    = -SNOWH * 3.0 / 4.0;
                    XYN(4)    = -SNOWH * 2.0 / 4.0;
                    XYN(5)    = -SNOWH / 4.0;

                    RTT(2)    = AIRT(dayi);
                    RTT(3)    = AIRT(dayi);
                    RTT(4)    = AIRT(dayi);
                    RTT(5)    = AIRT(dayi);
                    DX(LMN)   = SNOWH / 4.0;
                    DX(LMN+1) = SNOWH / 4.0;
                    DX(LMN+2) = SNOWH / 4.0;
                    DX(LMN+3) = SNOWH / 4.0;
                    
                else
                    LMN = 1;
                    XYN(1)    = -SNOWH;
                    XYN(2)    = -SNOWH * 4.0 / 5.0;
                    XYN(3)    = -SNOWH * 3.0 / 5.0;
                    XYN(4)    = -SNOWH * 2.0 / 5.0;
                    XYN(5)    = -SNOWH / 5.0;

                    RTT(1)    = AIRT(dayi);
                    RTT(2)    = AIRT(dayi);
                    RTT(3)    = AIRT(dayi);
                    RTT(4)    = AIRT(dayi);
                    RTT(5)    = AIRT(dayi);
                    DX(LMN)   = SNOWH / 5.0;
                    DX(LMN+1) = SNOWH / 5.0;
                    DX(LMN+2) = SNOWH / 5.0;
                    DX(LMN+3) = SNOWH / 5.0;
                    DX(LMN+4) = SNOWH / 5.0;
                end
                
                %% VARIABLES FOR ESTIMATING THE HEAT CONDUCTED
                % THROUGH THE SOIL, QC. in SURFACE ENERGY BALANCE MODEL
                XRTT7   = RTT(7) + TF; %
                XXYN7   = XYN(7);   % depth of node 7
                RKKPEAT = RKKK(7);  % thermal conductivity of node 7
                %%
                TA     = TF + AIRT(dayi);
                tDew   = TF + dewP(dayi);
                QSI    = SRin(dayi);
                AQSI   = (1 - albei) * QSI;
                wsi    = WS(dayi);
                PA     = PRE(dayi);
                roAir  = rho_a(TA, PA);
                EA     = atmosphericVaporPressure(tDew);
                SH     = 0.622 * EA / PA;
                QLI    = TRin(dayi);
                EMI    = emi(dayi);
                RZ     = rz(dayi);
                %%
                NCONTRALER = 0;
                REEE       = 1;
                %%
                syms f TS;
                % CALL SUBQLE
                Qle   = longWaveOut(emi(dayi), TS);
                % CALL SUBES0
                ES0   = atmosphericVaporPressure(TS);
                % CALL SUBQH
                QH    = sensible_heat(roAir, CAir, TA, TS, wsi, RZ, Lstar);
                % CALL SUBQE
                QE    = latent_heat(roAir, SH, ES0, PA, wsi, RZ, Lstar, root); 
                % CALL SUBQC
                if (SNOWH < 0.003) % snow free
                    
                    RERKK =  XXYN7 / (RKKPEAT / TDAYS);
                    
                else
                    
                    RERKK = SNOWH / ...
                        (KSNOW / TDAYS) + XXYN7 / (RKKPEAT / TDAYS);
                    
                end
                QC = - (TS - XRTT7) / RERKK;
                %%
                f(TS)   = AQSI + QLI - Qle - QH - QE + QC;
                df(TS)  = diff(f);

                %%
                while EES < REEE
                    
                    NCONTRALER = NCONTRALER + 1;
                    RFTS0   = f(root);
                    DXRFTS0 = df(root);

                    TEMRTS0 = root - RFTS0 / DXRFTS0;
                    REEE    = abs(TEMRTS0 - root);
                    root    = round(TEMRTS0, 5);

                end
                %%
                Ustar   = Friction_velocity(wsi, wsHeight, RZ, Lstar);
                Lstar   = Monin_obokhov(roAir, Ustar, TA, subs(QH, TS, root), subs(QE, TS, root));
                %%
                if (SNOWH > 0.003 && root > TF)
                    
                    root = TF;

                    QH   = sensible_heat(roAir, CAir, TA, root, wsi, RZ, Lstar);
                    ES0  = atmosphericVaporPressure(root);
                    QE   = latent_heat(roAir, SH, ES0, PA, wsi, RZ, Lstar, root);
                    Qle  = longWaveOut(emi(dayi), root);
                    QC   = -(root - XRTT7) / RERKK;

                    Ustar   = Friction_velocity(wsi, wsHeight, RZ, Lstar);
                    Lstar   = Monin_obokhov(roAir, Ustar, TA, QH, QE);
                                         
                end

%                 fprintf('%6.2f\n', root);
                
                %% SOIL THERMAL CONDUCTION
                %% CALL SUBTHETA  WATER AND UNFROZEN WATER
                
                RTT(LMN) = root - TF; % Upper boundary

                for I = 6:NODE
                    
                    % soil contents
                    if (XYN(I) <= SLTck1)
                        thetaS(I,1) = thetaS1; % solid
                        thetaW(I,1) = thetaU1;
                        satSL(I,1)  = satSL1;
                        thetaA(I,1) = satSL1 - thetaU1;
                    else
                        thetaS(I,1) = thetaS2;
                        thetaW(I,1) = thetaU2;
                        satSL(I,1)  = satSL2;
                        thetaA(I,1) = satSL2 - thetaU2;
                        
                    end
                    
                    if (RTT(I) < Tzero) % frozen
                        thetaU(I,1) = min(thetaW(I), ...
                            NiuYang(satSL1, (RTT(I) + TF), bvsilt, thaosilt));
                        
                        thetaI(I,1) = satSL(I) - thetaU(I) - thetaA(I);
                        RLL(I,1)    = roWater * ... % Eq. 19
                            (333.2 + 4.955 * RTT(I) + ...
                            0.02987 * RTT(I)^2) / 1000.0;
                        
                        if (RTT(I) < Tzero - dt)
                            DthetaU(I,1) = 0;
                        else
                            DthetaU(I,1) = thetaI(I);
                        end
                        
                    else % unfrozen
                        thetaU(I,1)   = thetaW(I);
                        thetaI(I,1)   = 0.0;
                        RLL(I,1)      = 0.0;
                        DthetaU(I,1)  = 0;
                    end
                    
                end
                %% CALL SUBRCCC
                for I = 6:NODE
                    
                    % Eq. 22
                    CCWATER = 4.20843 + 1.11362E-1 * RTT(I) + ...
                        5.12142E-3 * (RTT(I))^2 + 9.3482E-5 * (RTT(I))^3;
                    
                    % Eq. 27
                    cSoil = 0.4091 + 5.433E-3 * (TF + RTT(I));
                    
                    RCV(I, 1) = soilThermalCap(thetaU(I), thetaI(I), ...
                        thetaS(I), thetaA(I), CCWATER, cSoil);
                    % Eq. 14
                    RCCC(I,1) = (RCV(I) + 334 * DthetaU(I)/dt) * 1.0E6;
                    
                end
                %% CALL SUBRKKK
                for I = 6:NODE
                    % Eq. 20
                    kWater = 0.11455 + 1.6318E-3 * (TF + RTT(I));
                    if (XYN(I) <= SLTck1)
                        kSolid = kSL1;
                    else
                        kSolid = kSL2;
                    end
                    RKKK(I,1) = soilThermalCon(kSolid, kWater,...
                        thetaI(I), thetaU(I), thetaS(I), thetaA(I));
                    % Eq. 16
                    RKKK(I,1) = RKKK(I,1) * 24. * 3600.;
                end
                %% CALL ABCD
                
                RKW = 2.0 * RKKK(LMN) * RKKK(LMN + 1) / (RKKK(LMN) + RKKK(LMN + 1));
                RKE = 2.0 * RKKK(LMN + 1) * RKKK(LMN + 2) /...
                    (RKKK(LMN+1) + RKKK(LMN + 2));
                A(LMN + 1) = RKW / DX(LMN);
                C(LMN + 1) = RKE / DX(LMN + 1);
                
                RCW   = RCCC(LMN) * DX(LMN);
                RCE   = RCCC(LMN + 2) * DX(LMN + 1);
                RCCCP = (RCW + RCE) / (DX(LMN) + DX(LMN + 1));
                
                RDX      = (DX(LMN) + DX(LMN + 1)) / 2;
                APZERO   = RCCCP * RDX / RDELTAT;
                B(LMN+1) = -(A(LMN + 1) + C(LMN + 1) + APZERO);
                D(LMN+1) = -APZERO * RTT(LMN + 1) - A(LMN + 1) * RTT(LMN);
                %
                for IK = LMN+2:NODE-1
                    RKW = 2 * RKKK(IK) * RKKK(IK - 1) / (RKKK(IK) + RKKK(IK - 1));
                    RKE = 2 * RKKK(IK) * RKKK(IK + 1) / (RKKK(IK) + RKKK(IK + 1));
                    A(IK) = RKW / DX(IK - 1);
                    C(IK) = RKE / DX(IK);
                    
                    RCW   = RCCC(IK - 1) * DX(IK - 1);
                    RCE   = RCCC(IK + 1) * DX(IK);
                    RCCCP = (RCW + RCE) / (DX(IK - 1) + DX(IK));
                    
                    RDX    = (DX(IK - 1) + DX(IK)) / 2.0;
                    APZERO = RCCCP * RDX / RDELTAT;
                    B(IK)  = -(A(IK) + C(IK) + APZERO);
                    D(IK)  = -APZERO * RTT(IK);
                end
                %
                RKW     = RKKK(NODE - 1);
                A(NODE) = RKW / DX(NODE - 1);
                C(NODE) = 0.0;
                B(NODE) = -A(NODE);
                D(NODE) = -QQ;
                
                %%  CALL TRAM
                
                RP(LMN+1) = -C(LMN + 1) / B(LMN + 1);
                RQ(LMN+1) = D(LMN + 1) / B(LMN + 1);
                
                for I = LMN + 2:NODE
                    PP    = A(I) * RP(I - 1)+B(I);
                    RP(I) = -C(I)/PP;
                    RQQ   = D(I) - A(I) * RQ(I - 1);
                    RQ(I) = RQQ / PP;
                end
                
                DTT(NODE,1) = RQ(NODE);
                for I = NODE-1:-1:LMN+1
                    DTT(I) = RP(I) * DTT(I + 1) + RQ(I);
                end
                
                for I = LMN + 1:NODE
                    RTT(I) = DTT(I);
                end
                
                resGST(:, NDAYI)    = RTT;
                resthetaU(:, NDAYI) = thetaU;
                resAlbedo(NDAYI) = albei;
                resEmi(NDAYI) = emi(dayi);
                resQli(NDAYI) = QLI;
                resQE(NDAYI) = subs(QE, TS, root);
                resQH(NDAYI) = subs(QH, TS, root);
                
            end

            % modeled outputs
            outf = strcat(outBaseName, string(siteName), '.txt');
            
            fprintf('%s\n', outf);

            outST  = fullfile(dir_mod, outf);
            [fileID, fmt] = fopen(outST, 'wt');
            fmt = '%s %s %s %s %s %s\n';
            fprintf(fileID, fmt, 'date', 'gst', 'gst2', 'wu', 'QE', 'QH');
            fmt = '%s %6.2f %6.2f %6.2f %6.2f %6.2f\n';
            
            begDate3   = datetime(begYr, 1, 1);
            endDate3   = datetime(endYr, 12, 31);
            dateRange3 = begDate3:days(1):endDate3;

            
            for k = 1:length(dateRange3)
                
                fprintf(fileID, fmt, char(dateRange3(k)), ...
                    resGST(7, k), resGST(15, k), resthetaU(7, k), resQE(k), resQH(k));
                
            end
            fclose(fileID);

        end % soil moisture
    end % snow density
end % site


