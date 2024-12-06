% ======================
% Figure 5 and Figure S3
% ======================
%
% Model for CO2 and CH4 diffusion through the O-ring during typical TLS
% experiments.
% 
% Options:
%   - transport_mechanism:
%       - Solution diffusion
%       - Knudsen diffusion
%       (main diff. only if cell not initially empty)
% 
%   - initial amounts of gas in the cell
%       - empty cell     (N_CH4(t0) and N_CO2(t0) = 0) ==> Psi calculated
%       - non-empty cell (N_CH4(t0) and N_CO2(t0) > 0) ==> Psi predefined
% 
% Author: sebastien.viscardy@aeronomie.be
%
%% Link to functions
addpath('../functions/');

%%
clearvars

%% Physical constants
physcst_ref

%% Gas transport mechanisms
transport_mechanism = 'solution_diffusion';
% transport_mechanism = 'Knudsen_diffusion';

% select CO2/CH4 selectivity coefficient [/]
switch transport_mechanism
    case 'solution_diffusion'      
        omega = 10;% CO2/CH4 selectivity (solution diffusion) [/]
    case 'Knudsen_diffusion'
        omega = sqrt(M_CH4/M_CO2);
end
disp(' ')
disp(['type of transport   = ',transport_mechanism])
disp(['CO2/CH4 selectivity = ',num2str(omega)])
disp(' ')

%% Initial conditions
init_cell = 'non_empty_cell';
% init_cell = 'empty_cell';

switch init_cell
    case 'non_empty_cell', disp('cell initially empty')
    case 'empty_cell'    , disp('cell initially non-empty')
end
disp(' ')

%% Figure: main infos
figid      = 621;
sfigtype   = 'png';
ftsz       = 11;
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
isol        = 12; % E / Sol 684
% isol        = 17; % E / Sol 2627
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
delta_p_FO    = std(F_FO_prs);               % (full-cell) FO pressure [Pa]
delta_p_HC    = std(F_HC_prs);              % (full-cell) HC pressure [Pa]
delta_T       = std(F_FO_temp);              % (full-cell) temperature [K]
delta_eta_E   = std(E_Wefg_CH4)*1e-9;          % (empty-cell) CH4 vmr    [/]
delta_eta_F   = std(F_Wefg_CH4)*1e-9;          % (full-cell) CH4 vmr     [/]
delta_eta_HC  = sqrt(delta_eta_F^2/nFpts + ...
    delta_eta_E^2/nEpts);                      % CH4 vmr in HC           [/]
delta_N_HC    = eta_HC*V_HC/(rgas*T)*delta_p_HC + ...
    p_HC*V_HC/(rgas*T)*delta_eta_HC + ...
    p_HC*eta_HC*V_HC/(rgas*T^2)*delta_T;       % CH4 amount in HC        [mol]

%% Key times
tF_1     = F_elapsed_t(1);           % first full-cell run  [s]
tF_2     = F_elapsed_t(end);         % last full-cell run   [s]
tE_1     = E_elapsed_t(1);           % first empty-cell run [s]
tE_2     = E_elapsed_t(end);         % last empty-cell run  [s]

%% Time grid [s]
dtime    = 1;                        % time step            [s]
time     = 0:dtime:tE_2;             % time grid            [s]
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
dtime_ini     = dt_samp + (tF_2 - tF_1)/2;                    % ini/mid-full time diff. [s]
t_0_s         = tF_1-dt_samp;                                 % start of air sampling   [s]

switch init_cell
    case 'non_empty_cell'
        psi_CH4    = 1.01e-14;                                % transport coeff. (CH4)  [mol Pa-1 s-1]
        N_CH4_HC_0 = N_HC-psi_CH4*N_FO*rgas*T/V_FO*dtime_ini; % initial CH4 amount      [mol]
        N_CO2_HC_0 = omega*N_CH4_HC_0/eta_FO;                 % Initial CO2 amount      [mol]
        eta_HC_0   = N_CH4_HC_0/N_CO2_HC_0*1e9;               % Initial CH4 vmr         [ppbv]
    case 'empty_cell'
        psi_CH4    = V_FO*N_HC/(N_FO*rgas*T*dtime_ini);       % transport coeff. (CH4)  [mol Pa-1 s-1]
        N_CH4_HC_0 = 1e-40;                                   % initial CH4 amount      [mol]
        N_CO2_HC_0 = 1e-40;                                   % Initial CO2 amount      [mol]
        eta_HC_0   = 1e-40;                                   % Initial CH4 vmr         [ppbv]
end
if (N_CH4_HC_0 < 0), error('initial CH4 amount < 0 nmol'), end
psi_CO2  = omega*psi_CH4;                                     % transport coeff. (CO2)  [mol Pa-1 s-1]
p_HC_0   = (N_CO2_HC_0+N_CH4_HC_0)*rgas*T/V_HC;               % Initial pressure        [Pa]
t_HC     = N_HC/(psi_CH4*eta_FO*p_FO);                        % time to get N_H in cell [s] 

disp(['psi_CH4   = ',num2str(psi_CH4,'%2.2e'),' mol Pa-1 s-1'])
disp(['psi_CO2   = ',num2str(psi_CO2,'%2.2e'),' mol Pa-1 s-1'])
disp(['N_CO2(0)  = ',num2str(N_CO2_HC_0*1e9,'%2.2f'),' nmol'])
disp(['N_CH4(0)  = ',num2str(N_CH4_HC_0*1e9,'%2.2e'),' nmol'])
disp(['eta_HC(0) = ',num2str(eta_HC_0,'%2.2e'),' ppbv'])
disp(['p_H(0)    = ',num2str(p_HC_0,'%2.2f'),' Pa'])
disp(['t_HC      = ',num2str(t_HC/3600,'%2.2f'),' hours'])

%% Time simulation
for it = 1:ntime
    tt = time(it);
    
    if tt<t_0_s
        % sampling not yet started
        continue
    elseif tt==t_0_s+1
        p_HC_t(1:it-2)     = NaN;
        N_CO2_HC_t(1:it-2) = NaN;
        N_CH4_HC_t(1:it-2) = NaN;
        
        p_HC_t(it-1)       = p_HC_0;
        N_CO2_HC_t(it-1)   = N_CO2_HC_0;
        N_CH4_HC_t(it-1)   = N_CH4_HC_0;
    end
    
    %% Cell pumped out: Parameters
    kk           = time==tF_2;                 %
    delta_p_pump = p_HC_0-p_HC_t(kk);          % pressure difference   [Pa]
    delta_N_pump = V_HC/(rgas*T)*delta_p_pump; % air amount difference [mol]
    dt_pump      = tE_1-tF_2;                  % evacuation duration   [s]
    
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
        Gamma*(p_HC - p_HC_t(it-1))*(SS_heaviside(tt-t_0_s)-SS_heaviside(tt-tF_1)) + ...
        delta_N_pump/dt_pump*(SS_heaviside(tt-tF_2)-SS_heaviside(tt-tE_1));
    
    dN_CO2_dt_diff(it) = psi_CO2*(p_FO-p_HC_t(it-1));
    
    N_CO2_HC_t(it)     = N_CO2_HC_t(it-1)+dN_CO2_dt(it)*dtime;
    
    %% CH4 amount in the Herriott cell N_CH4_t
    eta_CH4        = N_CH4_HC_t(it-1)/(N_CO2_HC_t(it-1)+N_CH4_HC_t(it-1)); % CH4 vmr in cell [/]
    if ( N_CO2_HC_t(it-1)<1e-40 ), eta_CH4 = 1e-40; end
    
    dN_CH4_dt(it)      = psi_CH4*(eta_FO*p_FO - eta_CH4*p_HC_t(it-1)) + ...
        eta_CH4*delta_N_pump/dt_pump*(SS_heaviside(tt-tF_2)-SS_heaviside(tt-tE_1));
    
    dN_CH4_dt_diff(it) = psi_CH4*(eta_FO*p_FO - eta_CH4*p_HC_t(it-1));
    
    N_CH4_HC_t(it)     = N_CH4_HC_t(it-1)+dN_CH4_dt(it)*dtime;
    
    %% Cell pressure p_HC_t
    p_HC_t(it)     = (N_CO2_HC_t(it)+N_CH4_HC_t(it))*rgas*T/V_HC; % [Pa]
end

%% N_CH4_t at first and last times of full-cell runs
N_CH4_HC_tF_1 = N_CH4_HC_t(time==tF_1); % first full-cell run [mol]
N_CH4_HC_tF_2 = N_CH4_HC_t(time==tF_2); % last full-cell run  [mol]

N_CO2_HC_tF_1 = N_CO2_HC_t(time==tF_1); % first full-cell run [mol]
N_CO2_HC_tF_2 = N_CO2_HC_t(time==tF_2); % last full-cell run  [mol]

disp(['N_CH4_HC(tF_1) = ',num2str(N_CH4_HC_tF_1*1e9,'%2.2e'),' nmol'])
disp(['N_CH4_HC(tF_2) = ',num2str(N_CH4_HC_tF_2*1e9,'%2.2e'),' nmol'])

%% CH4 vmr in the Herriott cell eta_HC_t
eta_HC_t = N_CH4_HC_t./(N_CO2_HC_t+N_CH4_HC_t)*1e9; % [ppbv]

%% Statistics
disp('Statistics')
SS_diffusion_model_statistics

%% Make figures
disp('Make figures')
SS_diffusion_model_main_figure

%% Function 'q(t)' == > Figure S3
SS_diffusion_model_q_mod_vs_q_TLS