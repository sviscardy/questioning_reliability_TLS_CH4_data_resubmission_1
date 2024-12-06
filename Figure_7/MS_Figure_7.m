% ========
% Figure 7
% ========
% 
% Illustration of the two-step experiment proposed in Section 7 of the
% manuscript
% 
% Author: sebastien.viscardy@aeronomie.be
%
%% Link to functions
addpath('../functions/');

%%
clearvars

tic

%% Physical constants
physcst_ref
marscst

%% Key parameters
psi_CH4        = 1.01e-14; % transport coeff. (CH4)                   [mol Pa-1 s-1]
select_CO2_CH4 = 10;       % CO2/CH4 selectivity (solution diffusion) [/]

%% Gas transport mechanisms
transport_mechanism = 'solution_diffusion';
% transport_mechanism = 'Knudsen_diffusion';

% select CO2/CH4 selectivity coefficient [/]
switch transport_mechanism
    case 'solution_diffusion'
        omega = select_CO2_CH4;
    case 'Knudsen_diffusion'
        omega = sqrt(M_CH4/M_CO2);
end
disp(' ')
disp(['type of transport   = ',transport_mechanism])
disp(['psi(CH4)            = ',num2str(psi_CH4),' mol Pa-1 s-1'])
disp(['CO2/CH4 selectivity = ',num2str(omega)])
disp(' ')

%% Mean cell pressure (empty-cell runs) / (E) Sol 2627
ref_prs_sol     = 2627;
ref_E_p_H       = 2.6;  % [Pa]
ref_E_delta_p_H = 12.6; % [Pa]

%% Figure: main infos
figid      = 611;
sfigtype   = 'png';
ftsz       = 12;
ftnm       = 'times';
savefig    = 1;
savedata   = 0;
collist_D  = [34 178 34]/220;
collist_EN = [178 34 34]/178;
collist_ED = [30 144 255]/255;

collistF   = [0.8 0 0.8];
collistE   = [0.2 0.8 0.6];

lab_panel  = {'(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)'};

%% Type of experiment
t_exp_list = {'D' 'D' 'D' 'D' 'D' 'D' 'D' 'D' 'D' 'D' 'E' 'E' 'D' 'E' 'E' 'E' 'E' 'E'};
sol_list   = [79 81 106 292 306 313 466 474 504 526 573 684 684 2442 2446 2615 2627 2644];
nsol       = length(sol_list);

%% Select experiments
isol        = 17; % E / Sol 2627
t_exp       = t_exp_list{isol};
sol_index   = sol_list(isol);

disp(' ')
disp(['Experiment: ',t_exp,' ',num2str(sol_index)])
disp(' ')

%% Load full data
if ( sol_index < 2442 )
    SS_MSL_full_data_Webster_2015
else
    SS_MSL_full_data_Webster_2021
end
if t_exp=='D', n_E = 1; elseif t_exp=='E', n_E = 25; end

%% duration of air sampling 'dt_samp' [s]
if t_exp=='D', dt_samp = 20*60; elseif t_exp=='E', dt_samp = 120*60; end

%% Parameters
l_HC    = 16.8;                      % HC optical length            [m]
l_FO    = 0.09;                      % FO optical length            [m]
V_HC    = 405e-6;                    % HC volume                    [m3]
V_FO    = 988.2e-6;                  % FO volume                    [m3]
nFpts   = length(F_Wefg_CH4);        % number of full-cell runs
nEpts   = length(E_Wefg_CH4);        % number of empty-cell runs
p_HC    = mean(F_HC_prs);            % TLS: mean HC pressure        [Pa]
p_FO    = mean(F_FO_prs) ...
    - min(0,F_FO_prs(1));            % TLS: mean FO pressure        [Pa]
T       = mean(F_FO_temp);           % TLS: mean temperature        [K]
eta_E   = mean(E_Wefg_CH4)*1e-9;     % TLS: mean empty-cell CH4 vmr [/]
eta_F   = mean(F_Wefg_CH4)*1e-9;     % TLS: mean empty-cell CH4 vmr [/]
eta_HC  = eta_F - eta_E;             % TLS: mean CH4 vmr in HC      [/]
eta_FO = p_HC/p_FO*l_HC/l_FO*eta_E;  % TLS: mean CH4 vmr in FO      [/]
N_FO    = p_HC/(rgas*T) ...
    *l_HC/l_FO*eta_E*V_FO;           % TLS: CH4 amount in FO        [mol]
N_HC    = eta_HC*p_HC*V_HC/(rgas*T); % TLS: CH4 amount in HC        [mol]

%% Errors
delta_p_FO    = std(F_FO_prs);                 % (full-cell) FO pressure [Pa]
delta_p_HC    = std(F_HC_prs);                 % (full-cell) HC pressure [Pa]
delta_T       = std(F_FO_temp);                % (full-cell) temperature [K]
delta_eta_E   = std(E_Wefg_CH4)*1e-9;          % (empty-cell) CH4 vmr    [/]
delta_eta_F   = std(F_Wefg_CH4)*1e-9;          % (full-cell) CH4 vmr     [/]
delta_eta_HC  = sqrt(delta_eta_F^2/nFpts + ...
    delta_eta_E^2/nEpts);                      % CH4 vmr in HC           [/]
delta_N_HC    = eta_HC*V_HC/(rgas*T)*delta_p_HC + ...
    p_HC*V_HC/(rgas*T)*delta_eta_HC + ...
    p_HC*eta_HC*V_HC/(rgas*T^2)*delta_T;       % CH4 amount in HC        [mol]

%% Key times [s]

% duration of each step
dtF      = F_elapsed_t(end) - F_elapsed_t(1); % duration full-cell runs          [s]
dtE      = E_elapsed_t(end) - E_elapsed_t(1); % duration empty-cell runs         [s]
dt_pump  = E_elapsed_t(1)-F_elapsed_t(end);   % duration of evacuation procedure [s]

% 1.1) first set of empty-cell runs
tE1_1    = 0;
tE1_2    = dtE;

% 1.2) first set of full-cell runs
tF1_1    = tE1_2 + dt_samp;
tF1_2    = tF1_1 + dtF;

% 2.1) second set of full-cell runs (1 sol after start of the two-step
% experiment)
tF2_1    = round(sol_s);
tF2_2    = tF2_1 + dtF;

% 2.2) second set of empty-cell runs
tE2_1    = tF2_2 + dt_pump;
tE2_2    = tE2_1 + dtE;

%% time grid
dtime    = 1;                        % time step            [s]
time     = 0:dtime:tE2_2;            % time grid            [s]
time_min = time/60;                  % time grid            [min.]
time_h   = time/3600;                % time grid            [h]
ntime    = length(time);

%% Variables
p_HC_t         = zeros(1,ntime); % HC pressure (t)                    [Pa]
N_CO2_HC_t     = zeros(1,ntime); % CO2 amount in HC (t)               [mol]
N_CH4_HC_t     = zeros(1,ntime); % CH4 amount in HC (t)               [mol]
dN_CO2_dt      = zeros(1,ntime); % CO2 time variation (t)             [mol s-1]
dN_CH4_dt      = zeros(1,ntime); % CH4 time variation (t)             [mol s-1]

dN_CO2_dt_diff = zeros(1,ntime); % CO2 time variation (t) (diffusion) [mol s-1]
dN_CH4_dt_diff = zeros(1,ntime); % CH4 time variation (t) (diffusion) [mol s-1]

%% Initial conditions
t_0_s         = 0;                                         % start of air sampling   [s]
N_CH4_HC_0    = 1e-40;                                     % initial CH4 amount      [mol]
N_CO2_HC_0    = 1e-40;                                     % Initial CO2 amount      [mol]
psi_CO2       = omega*psi_CH4;                             % transport coeff. (CO2)  [mol Pa-1 s-1]
p_HC_0        = (N_CO2_HC_0+N_CH4_HC_0)*rgas*T/V_HC;       % Initial pressure        [Pa]

disp(['N_CO2(0)  = ',num2str(N_CO2_HC_0*1e9,'%2.2e'),' nmol'])
disp(['N_CH4(0)  = ',num2str(N_CH4_HC_0*1e9,'%2.2e'),' nmol'])
disp(['p_H(0)    = ',num2str(p_HC_0,'%2.2f'),' Pa'])
disp(['psi_CH4   = ',num2str(psi_CH4,'%2.2e'),' mol Pa-1 s-1'])
disp(['psi_CO2   = ',num2str(psi_CO2,'%2.2e'),' mol Pa-1 s-1'])

p_HC_t(1)       = p_HC_0;
N_CO2_HC_t(1)   = N_CO2_HC_0;
N_CH4_HC_t(1)   = N_CH4_HC_0;

%% Time simulation
for it = 2:ntime
    tt = time(it);
    
    %% Cell pumped out: Parameters
    kk           = time==tF2_2;
    delta_p_pump = p_HC_0-p_HC_t(kk);          % pressure difference   [Pa]
    delta_N_pump = V_HC/(rgas*T)*delta_p_pump; % air amount difference [mol]
    
    %% CO2 amount in the Herriott cell N_CO2_t
    switch t_exp
        case 'D'
            tau_samp = 5*60;
            Gamma    = V_HC/(rgas*T*tau_samp); %
        case 'E'
            tau_samp = 5*60*(120/60)^2;
            Gamma    = V_HC/(rgas*T*tau_samp); %
    end
    
    dN_CO2_dt(it)      = psi_CO2*(p_FO-p_HC_t(it-1)) + ...
        Gamma*(p_HC - p_HC_t(it-1))*(SS_heaviside(tt-tE1_2)-SS_heaviside(tt-tF1_1)) + ...
        delta_N_pump/dt_pump*(SS_heaviside(tt-tF2_2)-SS_heaviside(tt-tE2_1));
    
    dN_CO2_dt_diff(it) = psi_CO2*(p_FO-p_HC_t(it-1));
    
    N_CO2_HC_t(it)     = N_CO2_HC_t(it-1)+dN_CO2_dt(it)*dtime;
    
    %% CH4 amount in the Herriott cell N_CH4_t
    eta_CH4        = N_CH4_HC_t(it-1)/(N_CO2_HC_t(it-1)+N_CH4_HC_t(it-1)); % CH4 vmr in cell [/]
    if ( N_CO2_HC_t(it-1)<1e-30 ), eta_CH4 = 1e-40; end
    
    dN_CH4_dt(it)      = psi_CH4*(eta_FO*p_FO - eta_CH4*p_HC_t(it-1)) + ...
        eta_CH4*delta_N_pump/dt_pump*(SS_heaviside(tt-tF2_2)-SS_heaviside(tt-tE2_1));
    
    dN_CH4_dt_diff(it) = psi_CH4*(eta_FO*p_FO - eta_CH4*p_HC_t(it-1));
    
    N_CH4_HC_t(it)     = N_CH4_HC_t(it-1)+dN_CH4_dt(it)*dtime;
    
    %% Cell pressure p_HC_t
    p_HC_t(it)     = (N_CO2_HC_t(it)+N_CH4_HC_t(it))*rgas*T/V_HC; % [Pa]
end

%% CH4 vmr in the Herriott cell eta_HC_t
eta_HC_t = N_CH4_HC_t./(N_CO2_HC_t+N_CH4_HC_t)*1e9; % [ppbv]
kk = N_CO2_HC_t < 1e-30;
eta_HC_t(kk) = NaN;

%% Simulation of data points
SS_sim_data

%% Make figure
SS_main_figure

toc