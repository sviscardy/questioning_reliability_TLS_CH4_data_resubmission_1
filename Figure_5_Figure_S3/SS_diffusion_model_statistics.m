%% Simulation of {eta_F(k)} and {eta_E(k)} + statistics
sig_F_TLS     = std(F_Wefg_CH4);                   % std dev. full-cell runs         [ppbv]
sig_E_TLS     = std(E_Wefg_CH4);                   % std dev. empty-cell runs        [ppbv]

contrib_FO    = N_FO/(V_FO)*l_FO;                  % FO contribution
contrib_HC    = N_CH4_HC_t/V_HC*l_HC;              % cell contribution
q_fct         = contrib_HC./( contrib_FO);         % function 'q(t)'                 [/]

eta_mod       = eta_E*(1+q_fct)*1e9;               % function 'eta_m(t)'             [ppbv]

eta_F_mod_TLS = interp1(time,eta_mod,F_elapsed_t); % 'eta_m' at 'full-cell' times    [ppbv]
eta_E_mod_TLS = interp1(time,eta_mod,E_elapsed_t); % 'eta_m' at 'empty-cell' times   [ppbv]

%% TLS data: linear regression

% 1) Full-cell runs (TLS)
xl_F     = [tF_1 tF_2];
x        = F_elapsed_t;
y        = F_Wefg_CH4;
X        = [ones(length(x),1) x];
bT_F     = X\y;
y_linT_F = bT_F(2)*xl_F+bT_F(1);
a1T_F    = bT_F(2);
b1T_F    = bT_F(1);

% 2) Empty-cell runs (TLS)
xl_E     = [tE_1 tE_2];
x        = E_elapsed_t;
y        = E_Wefg_CH4;
X        = [ones(length(x),1) x];
bT_E     = X\y;
y_linT_E = bT_E(2)*xl_E+bT_E(1);
a1T_E    = bT_E(2);
b1T_E    = bT_E(1);

%% size of the sample
nsim     = 1e6;

%% Load or generate random values
sdpath     = 'data'; mkdir(sdpath)
sdfilename = [t_exp,'_',num2str(sol_index),'_nsim_',num2str(nsim),'_data_',init_cell,'.mat'];
sdflnm     = fullfile(sdpath,sdfilename);

if exist(sdflnm,'file')
    load(sdflnm)
    disp('statistics: data loaded')
else
    %% Generate large number of instances
    a_s_F_list      = zeros(1,nsim);
    a_s_E_list      = zeros(1,nsim);
    b_s_F_list      = zeros(1,nsim);
    b_s_E_list      = zeros(1,nsim);
    y_lins_F_list   = zeros(nsim,2);
    y_lins_E_list   = zeros(nsim,2);
    
    mod_eta_HC_list = zeros(1,nsim);
    mod_sig_HC_list = zeros(1,nsim);
    eta_F_mod_list  = zeros(nsim,nFpts);
    eta_E_mod_list  = zeros(nsim,nEpts);
    
    %% Loop over samples
    for isim = 1:nsim
        eta_F_mod      = eta_F_mod_TLS + sig_F_TLS*randn(nFpts,1);
        eta_E_mod      = eta_E_mod_TLS + sig_E_TLS*randn(nEpts,1);
        
        mean_eta_F_mod = mean(eta_F_mod);
        mean_eta_E_mod = mean(eta_E_mod);
        
        eta_H_mod      = mean_eta_F_mod-mean_eta_E_mod;
        
        sig_F_mod      = std(eta_F_mod);
        sig_E_mod      = std(eta_E_mod);
        
        sig_HC_mod     = sqrt(sig_F_mod^2/nFpts + sig_E_mod^2/nEpts);
        
        % Full-cell runs (sim.)
        xl_F     = [tF_1 tF_2];
        x        = F_elapsed_t;
        y        = eta_F_mod;
        X        = [ones(length(x),1) x];
        b        = X\y;
        y_lins_F = b(2)*xl_F+b(1);
        a_s_F    = b(2);
        b_s_F    = b(1);
        
        % Empty-cell runs (sim.)
        xl_E     = [tE_1 tE_2];
        x        = E_elapsed_t;
        y        = eta_E_mod;
        X        = [ones(length(x),1) x];
        b        = X\y;
        y_lins_E = b(2)*xl_E+b(1);
        a_s_E    = b(2);
        b_s_E    = b(1);
        
        %% Store data
        a_s_F_list(isim)       = a_s_F;
        a_s_E_list(isim)       = a_s_E;
        
        b_s_F_list(isim)       = b_s_F;
        b_s_E_list(isim)       = b_s_E;
        
        y_lins_F_list(isim,:)  = y_lins_F;
        y_lins_E_list(isim,:)  = y_lins_E;
        
        mod_eta_HC_list(isim)  = eta_H_mod;
        mod_sig_HC_list(isim)  = sig_HC_mod;
        
        eta_F_mod_list(isim,:) = eta_F_mod;
        eta_E_mod_list(isim,:) = eta_E_mod;
    end
    
    %% Save data
    save(sdflnm,'a_s_F_list','a_s_E_list','b_s_F_list','b_s_E_list', ...
        'y_lins_F_list','y_lins_E_list', ...
        'mod_eta_HC_list','mod_sig_HC_list', ...
        'eta_F_mod_list','eta_E_mod_list','-v7.3')
end


mu_F_mod  = mean(a_s_F_list)*3600;
sig_F_mod = std(a_s_F_list)*3600;

mu_E_mod  = mean(a_s_E_list)*3600;
sig_E_mod = std(a_s_E_list)*3600;

disp(' ')
disp('Statistics')
disp('----------')
disp('')
disp(['mu_F_mod  = ',num2str(mu_F_mod,'%2.2f'),' ppbv h-1'])
disp(['sig_F_mod = ',num2str(sig_F_mod,'%2.2f'),' ppbv h-1'])
disp(['mu_E_mod  = ',num2str(mu_E_mod,'%2.2f'),' ppbv h-1'])
disp(['sig_E_mod = ',num2str(sig_E_mod,'%2.2f'),' ppbv h-1'])
