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
dir_mod = '/Users/shengdiwang/Library/CloudStorage/OneDrive-个人/TIB_test/scaled/SIM';
dir_foc = '/Users/shengdiwang/Library/CloudStorage/OneDrive-个人/TIB_test/scaled/FOC';


thetaU1s = [0.20, 0.27, 0.15, 0.13, 0.29, 0.18, 0.27, ...
    0.26, 0.34, 0.20, 0.16, 0.19, 0.20, 0.21, 0.11, 0.18, ... 
    0.18, 0.30, 0.22, 0.23, 0.28, 0.23, 0.20, 0.15];
thetaU2s = [0.30, 0.14, 0.24, 0.16, 0.16, 0.16, 0.17, ...
    0.23, 0.21, 0.18, 0.23, 0.14, 0.21, 0.20, 0.19, 0.12, ...
    0.17, 0.20, 0.20, 0.20, 0.22, 0.20, 0.19, 0.21];
thetaUStd1s = [0.02, 0.03, 0.02, 0.02, 0.03, 0.02, 0.02, ...
    0.03, 0.03, 0.02, 0.02, 0.02, 0.03, 0.01, 0.01, 0.02, ...
    0.01, 0.02, 0.03, 0.02, 0.03, 0.02, 0.03, 0.02];
 

% ---- soil ----

% 5-layer snow, 16 layers for 0-1 m, 190 layers for 1-20 m [0.1 m interval]
% soil discretization for 0-1 m
% discretization
soil1m = [0, 3, 7, 11, 15, 22, 29, 35, 44, 50, 58, 66, 75, 84, 93, 100];
soil1m = soil1m/100; % in m
soilDepth = 10;  % total soil depth
soilTck   = 0.1; % thickness [m] for deep soil (> 1 m)

snowN    = 5;   % maximum snow layer
NODE     = snowN + length(soil1m) + (soilDepth - 1)/soilTck;

SLTck1   = 2; % soil layer 1 thickness (Table 1)

clay = 3;
sand = 80;
bvsilt   = vb(clay);    % clay
thaosilt = vthao(sand); % sand
kSL2     = 2.92;      % silt thermal conductivity

intiST    = -8;  % initial soil temperature [C]
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

roWater = 1000.0; % water density

Tzero = 0; % frozen temperature
QQ    = 0; % the lower boundary condition

%% Assign all parameters needed in SURFACE ENERGY BALANCE METHOD
% AIR
roAir = 1.225;  % air density [kg m-3]
CAir  = 1003.5;  % air thermal capacity [J m-3 C-1]

% ICE
roIce = 920;   % ice density
cIce  = 2.05;  % ice volumetric heat capacity [MJ m-3 K-1]
kIce  = 2.29;  % thermal conductivity of ice  [W m-1 K-1]

% CONSTANT
GRAVIT   = 9.807;   % gravitational acceleration
VONK     = 0.4;     % Von Karman's constant
TF       = 273.15;  % unit C to K
latSub   = 2.838E6; % latent heat of sublimation [J kg-1]
wsHeight = 2.0;     % REFERENCE HEIGHT for wind speed measure
albedoG  = 0.2;     % snow-free albedo

%% ===== ENSEMBLE FORECASTING =============================================

dt = 1.0; % phase change
% thetaU2 = 0.10;

satSL1s = thetaU1s + 0.05; % saturated water content soil layer 1
satSL2s = thetaU2s + 0.05; % saturated water content soil layer 2

thetaS1s = 1 - satSL1s; % solid particle soil layer 1
thetaS2s = 1 - satSL2s; % solid particle soil layer 2

%% ---- Energy balance ----

% for kSL1 = 1.5:0.1:1.9  % thermal conductivity

kSL1 = 2.92;

    
% simulation period
begYr = 1900;
begDate1   = datetime(begYr, 1, 1);    % start date for run
endDate1   = datetime(begYr+5, 12, 30); % end date for run
dateRange1 = begDate1:endDate1;


begDate2   = datetime(begYr, 1, 1); % start date for run
endDate2   = datetime(2033, 11, 30); % end date for run
dateRange2 = begDate2:endDate2;


% soil mositure
satSL1 = 0.2;
satSL2 = 0.2;

thetaS1 = 0.2;
thetaS2 = 0.2;

thetaU1 = 0.2;
thetaU2 = 0.2;

thetaUStd = 0.2;

thetaU1Range = thetaU1-thetaUStd:0.01:thetaU1+thetaUStd;


% fprintf(sitei)

%% ===== IMPORT =======================================================

% ---- import forcing ----

Path = '/Users/shengdiwang/Library/CloudStorage/OneDrive-个人/TIB_test/scaled/FOC';
File = dir(fullfile(Path,'*.csv'));
FileNames = {File.name}';

for focf_name = 1:length(FileNames)
    focf_name = FileNames(focf_name);
    focf_name = char(focf_name);
    % focf_name        = 'Snow profile1_foc.csv'; % forcing file
    outBaseName = 'TIB_'; % modeled outputs
    
    % ---- make file ----
    focf       = fullfile(dir_foc, focf_name);
    
    foc = readtable(focf); % forcing
    
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
    RZ     = foc.z;      % roughness length
    albe   = foc.albedo; % albedo
    obsGST = foc.tempS;  % observed soil temp at 0.01 m
    
    EES     = 0.0001;
    TDAYS   = 3600.0; % 1 day in second
    RDELTAT = 1.0;
    
    % ASSIGN INITIAL SOIL TEMPERATURE
    for I = snowN + 1:NODE %
        
        RTT(I, 1) = intiST;
        
    end
    
    % ======= SPIN-UP =====================================================
    
    NTB1 = strmatch(datestr(begDate1), DATES); % forcing index of begDate
    
    for sdn = 300 % snow density
        
        snowRo(snowRo > 0) = sdn;
        
        %0.16:0.02:0.20 % soil mositure layer 1
        for thetaU1 = 0.2
            
            for yr = 1:spinTime
                
                for NDAYI = 1:length(dateRange1)
                    
                    dayi = NTB1 + NDAYI - 1;
                    RTS0 = TF + AIRT(dayi);
                    
                    % SNOW LAYER NUMBER AND THICKNESS
                    
                    SNOWH  = SND(dayi);             % snow depth
                    ROSNOW = snowRo(dayi);          % snow density
                    KSNOW = snowThermalCon(ROSNOW); % snow conductivity
                    CSNOW = snowThermalCap(ROSNOW); % snow capacity
                    KSNOW = KSNOW * 3600.0 * 24.0;
                    
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
                    if (SNOWH < 0.02)
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
                    %
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
                    EA     = atmosphericVaporPressure(tDew);
                    QLI    = TRin(dayi); %longWaveIn(EA, TA); %
                    NCONTRALER = 1;
                    
                    %%
                    REEE        = 1;
                    while EES < REEE
                        
                        NCONTRALER = NCONTRALER + 1;
                        
                        % CALL SUBQLE
                        Qle   = longWaveOut(emi(dayi), RTS0);
                        % CALL SUBDXQLE
                        DXQle = derivativeLongWaveOut(emi(dayi), RTS0);
                        % CALL SUBDHDE
                        % Eq.7
                        DH    = VONK^2 * wsi / (log(wsHeight / RZ(dayi)))^2;
                        DE    = DH; % Eq.7
                        % CALL SUBRIXI
                        % Eq.9
                        RI    = GRAVIT * wsHeight * (TA - RTS0) / (TA * wsi^2);
                        
                        if (RTS0 > TA)
                            
                            XI = 1.0;
                            
                        else
                            
                            XI = 1 / (1 + 10 * RI); % Eq.8
                            
                        end
                        
                        % CALL SUBES0
                        ES0   = atmosphericVaporPressure(tDew);
                        % CALL SUBQH
                        QH    = roAir * CAir * DH * XI * (TA - RTS0); % Eq. 5
                        % CALL SUBDXQH
                        
                        if(RTS0 <= TA)
                            YDXQH1 = roAir * CAir * ...
                                VONK^2 * wsi / (log(wsHeight / RZ(dayi)))^2;
                            YDXQH2 = GRAVIT * wsHeight*(TA-RTS0)/(TA * wsi^2);
                            DXQH1  = YDXQH1 * 10 * YDXQH2/(1 + 10 * YDXQH2)^2;
                            DXQH2  = YDXQH1 / (1 + 10 * YDXQH2);
                            DXQH   = DXQH1 - DXQH2;
                            
                        else
                            
                            DXQH1 = roAir * CAir * ...
                                VONK^2 * wsi / (log(wsHeight / RZ(dayi)))^2;
                            DXQH = -DXQH1;
                        end
                        % CALL SUBQE
                        % Eq.6 
                        QE = roAir * latSub * DE * XI * (0.622 * (EA - ES0)/PA);
                        
                        % CALL SUBDXQE
                        if(RTS0 <= TA)
                            
                            YDXDH1 = roAir * latSub * VONK^2 * wsi ...
                                / (log(wsHeight / RZ(dayi)))^2;
                            YDXQH2 = GRAVIT * wsHeight * (TA - RTS0)/(TA * wsi^2);
                            YDXQH3 = 0.622 * (EA - ES0) / PA;
                            DXQE0  = YDXDH1 * YDXQH3 * ...
                                10 * GRAVIT * wsHeight / (TA * wsi^2);
                            DXQE1  = DXQE0/(1 + 10 * YDXQH2)^2 ;
                            DXQE2  = YDXDH1 * 0.622 * ES0 / PA;
                            DXQE3  = 2353 * log(10) / RTS0^2 / (1 + 10 * YDXQH2);
                            DXQE   = DXQE1 - DXQE2 * DXQE3;
                            
                        else
                            
                            YDXDH1 = roAir * latSub * ...
                                VONK^2 * wsi / (log(wsHeight/RZ(dayi)))^2;
                            DXQE2 = 0.622 * ES0 / PA;
                            DXQE3 = 2353 * log(10)/ RTS0^2;
                            DXQE  = -YDXDH1 * DXQE2 * DXQE3;
                            
                        end
                        
                        % CALL SUBQC
                        if (SNOWH < 0.003) % snow free
                            
                            RERKK =  XXYN7 / (RKKPEAT / TDAYS);
                            
                        else
                            
                            RERKK = SNOWH / ...
                                (KSNOW / TDAYS) + XXYN7 / (RKKPEAT / TDAYS);
                            
                        end
                        
                        QC = -(RTS0 - XRTT7) / RERKK ;
                        %
                        % CALL SUBDXQC
                        if (SNOWH < 0.003)
                            RERKK = XXYN7 / (RKKPEAT / TDAYS);
                        else
                            RERKK = SNOWH / ...
                                (KSNOW / TDAYS) + XXYN7 / (RKKPEAT / TDAYS);
                        end
                        
                        DXQC = -1.0 / RERKK;
                        %%
                        RFTS0   = AQSI + QLI + Qle + QH + QE + QC;
                        DXRFTS0 = DXQle + DXQH + DXQE + DXQC;
                        
                        %%
                        TEMRTS0 = RTS0 - RFTS0 / DXRFTS0;
                        REEE    = abs(TEMRTS0 - RTS0);
                        RTS0    = TEMRTS0;
                        
                    end
                    
                    if (SNOWH > 0.003 && RTS0 > TF)
                        
                        RTS0 = TF;
                        QH  = roAir * CAir * DH * XI * (TA - RTS0);
                        QE  = roAir * latSub * DE * XI * (0.622 * (EA-ES0)/ PA);
                        Qle = longWaveOut(emi(dayi), RTS0);
                        QC = -(RTS0 - XRTT7) / RERKK;
                        
                    else
                        
                        RTS0 = TEMRTS0;
                        
                    end
                    
                    %% SOIL THERMAL CONDUCTION
                    %% CALL SUBTHETA  WATER AND UNFROZEN WATER
                    
                    RTT(LMN) = RTS0 - TF; % Upper boundary
                    
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
                
                dayi = NTB2 + NDAYI - 1;
                RTS0 = TF + AIRT(dayi);
                
                % INPUT iiTH DAY'S SNOW THICKNESS.
                % ESTABLISH THE LAYER NUMBER AND THICKNESS OF SNOW.
                
                SNOWH  = SND(dayi);    % snow depth
                ROSNOW = snowRo(dayi); % snow density
                KSNOW  = snowThermalCon(ROSNOW);
                CSNOW  = snowThermalCap(ROSNOW);
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
                    LMN       = 5;
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
                %
                %% VARIABLES FOR ESTIMATING THE HEAT CONDUCTED
                % THROUGH THE SOIL, QC. in SURFACE ENERGY BALANCE MODEL
                XRTT7   = RTT(7) + TF; %
                XXYN7   = XYN(7);  % depth of node 7
                RKKPEAT = RKKK(7); % thermal conductivity of node 7
                %%
                TA     = TF + AIRT(dayi);
                tDew   = TF + dewP(dayi);
                QSI    = SRin(dayi);
                AQSI   = (1 - albei) * QSI;
                wsi    = WS(dayi);
                PA     = PRE(dayi);
                EA     = atmosphericVaporPressure(tDew);
                QLI    = TRin(dayi); %longWaveIn(EA, TA); %
                NCONTRALER = 1;
                
                %%
                REEE        = 1;
                while EES < REEE
                    
                    NCONTRALER = NCONTRALER + 1;
                    
                    % CALL SUBQLE
                    Qle   = longWaveOut(emi(dayi), RTS0);
                    % CALL SUBDXQLE
                    DXQle = derivativeLongWaveOut(emi(dayi), RTS0);
                    % CALL SUBDHDE
                    % Eq.7
                    DH    = VONK^2 * wsi / (log(wsHeight / RZ(dayi)))^2;
                    DE    = DH;
                    % CALL SUBRIXI
                    % Eq.9
                    RI    = GRAVIT * wsHeight * (TA - RTS0) / (TA * wsi^2);
                    
                    if (RTS0 > TA)
                        
                        XI = 1.0;
                        
                    else
                        % Eq.8
                        XI = 1 / (1 + 10 * RI);
                        
                    end
                    
                    % CALL SUBES0
                    ES0   = atmosphericVaporPressure(tDew);
                    % CALL SUBQH
                    % Eq. 5
                    QH    = roAir * CAir * DH * XI * (TA - RTS0);
                    % CALL SUBDXQH
                    
                    if(RTS0 <= TA)
                        YDXQH1 = roAir * CAir * ...
                            VONK^2 * wsi / (log(wsHeight / RZ(dayi)))^2;
                        YDXQH2 = GRAVIT * wsHeight * (TA-RTS0) / (TA * wsi^2);
                        DXQH1  = YDXQH1 * 10 * YDXQH2 / (1 + 10 * YDXQH2)^2;
                        DXQH2  = YDXQH1 / (1 + 10 * YDXQH2);
                        DXQH   = DXQH1 - DXQH2;
                        
                    else
                        DXQH1 = roAir * CAir * ...
                            VONK^2 * wsi / (log(wsHeight / RZ(dayi)))^2;
                        DXQH = -DXQH1;
                    end
                    % CALL SUBQE [Eq.6]
                    
                    QE    = roAir * latSub * DE * XI * (0.622 * (EA - ES0) / PA );
                    
                    % CALL SUBDXQE
                    if(RTS0 <= TA)
                        
                        YDXDH1 = roAir * latSub * ...
                            VONK^2 * wsi / (log(wsHeight / RZ(dayi)))^2;
                        YDXQH2 = GRAVIT * wsHeight * (TA - RTS0) / (TA * wsi^2);
                        YDXQH3 = 0.622 * (EA - ES0) / PA;
                        DXQE0  = YDXDH1 * YDXQH3 * ...
                            10 * GRAVIT * wsHeight / (TA * wsi^2);
                        DXQE1  = DXQE0/(1 + 10 * YDXQH2)^2 ;
                        DXQE2  = YDXDH1 * 0.622 * ES0 / PA;
                        DXQE3  = 2353 * log(10) / RTS0^2 / (1 + 10 * YDXQH2);
                        DXQE   = DXQE1 - DXQE2 * DXQE3;
                        
                    else
                        
                        YDXDH1 = roAir * latSub * ...
                            VONK^2 * wsi / (log(wsHeight / RZ(dayi)))^2;
                        DXQE2 = 0.622 * ES0 / PA;
                        DXQE3 = 2353 * log(10)/ RTS0^2;
                        DXQE = -YDXDH1 * DXQE2 * DXQE3;
                        
                    end
                    
                    % CALL SUBQC
                    if (SNOWH < 0.003) % snow free
                        
                        RERKK =  XXYN7 / (RKKPEAT / TDAYS);
                        
                    else
                        
                        RERKK = SNOWH / ...
                            (KSNOW / TDAYS) + XXYN7 / (RKKPEAT / TDAYS);
                        
                    end
                    
                    QC = -(RTS0 - XRTT7) / RERKK ;
                    %
                    % CALL SUBDXQC
                    if (SNOWH < 0.003)
                        RERKK = XXYN7 / (RKKPEAT / TDAYS);
                    else
                        RERKK = SNOWH / ...
                            (KSNOW / TDAYS) + XXYN7 / (RKKPEAT / TDAYS);
                    end
                    
                    DXQC = -1.0 / RERKK;
                    %%
                    RFTS0   = AQSI + QLI + Qle + QH + QE + QC;
                    DXRFTS0 = DXQle + DXQH + DXQE + DXQC;
                    
                    %%
                    TEMRTS0 = RTS0 - RFTS0 / DXRFTS0;
                    REEE    = abs(TEMRTS0 - RTS0);
                    RTS0    = TEMRTS0;
                    
                end
                
                if (SNOWH > 0.003 && RTS0 > TF)
                    
                    RTS0 = TF;
                    QH   = roAir * CAir * DH * XI * (TA - RTS0);
                    QE   = roAir * latSub * DE * XI * (0.622 * (EA-ES0) / PA);
                    Qle  = longWaveOut(emi(dayi), RTS0);
                    QC   = -(RTS0 - XRTT7) / RERKK;
                    
                else
                    
                    RTS0 = TEMRTS0;
                    
                end
                
                %% SOIL THERMAL CONDUCTION
                %% CALL SUBTHETA  WATER AND UNFROZEN WATER
                
                RTT(LMN) = RTS0 - TF; % Upper boundary
                
                
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
                
            end
            
            % modeled outputs
            out_site = strsplit(focf_name,'.');
            out_site = out_site(1);
            outf = strcat(out_site, ...
                '_sdn', string(sdn), ...
                '_sm', string(thetaU1 * 100), ...
                '_sand', string(sand), ...
                '_clay', string(clay), ...
                '_tc', string(kSL1 * 100), ...   
                '.txt');
            
            fprintf('%s\n', outf);
    
            hold on
            plot(1:48912, resGST(7, 1:48912), 'r')
            ylim([-45, 45])
            xlim([1, 48912])
            plot(1:48912, SND*10, 'b')
            plot(1:48912, AIRT, 'k')
            plot(xlim, [0, 0])
            hold off
    
            outST  = fullfile(dir_mod, outf);
            [fileID, fmt] = fopen(outST, 'wt');
            fmt = '%s %s %s %s\n';
            fprintf(fileID, fmt, 'date', 'gst', 'wu', 'lwin');
            fmt = '%s %6.2f %6.2f %6.2f\n';
            
            begDate3   = datetime(2017, 1, 1, 0, 0, 0);
            endDate3   = datetime(2022, 7, 31, 23, 0, 0);
            dateRange3 = begDate3:hours(1):endDate3;

            
            for k = 1:length(dateRange2)
                
                fprintf(fileID, fmt, char(dateRange3(k)), ...
                    resGST(7, k), resthetaU(7, k), resQli(k));
                
            end
            fclose(fileID);
            
        end % soil moisture
    end % snow density
    
end
