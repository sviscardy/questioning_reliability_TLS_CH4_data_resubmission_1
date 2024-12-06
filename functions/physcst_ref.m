% Universal constants, standards, and factors for unit conversions
% 
% Author: sebastien.viscardy@aeronomie.be
%
%% Physical constants
rgas             = 8.314;            % Ideal gas constant         [J K-1 mol-1]
NA               = 6.022137e23;      % Avogadro constant          [particles mol-1]
kboltz           = 1.380649e-23;     % Boltzmann constant         [J K-1]
amu              = 1.660538921e-27;  % (unified) atomic mass unit [kg] (amu or Da)
hplanck          = 6.62607015e-34;   % Planck's constant          [kg m2 s-1]
clight           = 299792458;        % light speed                [m s-1]
Ggrav            = 6.67408e-11;      % Gravitational constant     [m3 kg-1 s-2]

%% Standards
vsmow = 155.76e-6;  % Vienna Standard Mean Ocean Water (D/H)     [/]
vpdb  = 1.12372e-2; % Vienna Pee Dee Belemnite (13C/12C)         [/]

%% Units conversion
prmum        = 3.34e18*1e4;
prmum_cm2    = 3.34e18;
prmum_m2     = 3.34e18*1e4;
DU           = 2.69e16*1e4;
DU_cm2       = 2.69e16;
DU_m2        = 2.69e16*1e4;
mumatm       = DU/10;
mumatm_cm2   = DU_cm2/10;
mumatm_m2    = DU_m2/10;

deg2rad      = pi/180;

%% Molar mass of key molecules [kg mol-1]
% ('M_xxx', where 'xxx' is the name of the molecules in capital letters)
M_CO2   = 44e-3;   % CO2 molar mass   [kg mol-1]
M_CO    = 18e-3;   % CO molar mass    [kg mol-1]
M_O2    = 32e-3;   % O2 molar mass    [kg mol-1]
M_H2O   = 18e-3;   % H2O molar mass   [kg mol-1]
M_HDO   = 19e-3;   % HDO molar mass   [kg mol-1]
M_H2SO4 = 98e-3;   % H2SO4 molar mass [kg mol-1]
M_CH4   = 16e-3;   % CH4 molar mass   [kg mol-1]
