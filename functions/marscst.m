% Key parameters for Mars
% 
% Author: sebastien.viscardy@aeronomie.be
%
%% Universal constants and standards
physcst_ref

%% Planet Mars
radius_mars  = 3389.5;           % Planet's radius     [km]
M_mars       = 6.4171e23;        % Planet's mass       [kg]
gmars        = 3.72;             % gravity at surface  [m/s2]

%% Time
% MY1, Ls = 0° : April 11, 1955
sol_h        = 24.6597;          % solar day           [h]
sol_s        = sol_h*3600;       % solar day           [s]
year_d       = 686.973;          % tropical year       [day]
year_s       = year_d*24*3600;   % tropical year       [s]
year_sol     = year_s/sol_s;     % # Sol per year      [Sol]
solhr        = 3698.969;

%% Atmosphere
mairmars     = 43.412;           % Mean molar air mass [g mol-1]
Mairmars     = 43.412e-3;        % Mean molar air mass [kg mol-1]
