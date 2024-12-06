%% eta_sim
sig_F_TLS     = std(F_Wefg_CH4);                   % std dev. full-cell runs         [ppbv]
sig_E_TLS     = std(E_Wefg_CH4);                   % std dev. empty-cell runs        [ppbv]

contrib_FO    = N_FO/(V_FO)*l_FO;                  % FO contribution
contrib_HC    = N_CH4_HC_t/V_HC*l_HC;              % cell contribution
q_fct         = contrib_HC./( contrib_FO);         % function 'q(t)'                 [/]

eta_sim       = eta_E*(1+q_fct)*1e9;               % function 'eta_s(t)'             [ppbv]

%% generating random data points
tE1_pts = tE1_1:dtE/(nEpts-1):tE1_2;
tF1_pts = tF1_1:dtF/(nFpts-1):tF1_2;
tF2_pts = tF2_1:dtF/(nFpts-1):tF2_2;
tE2_pts = tE2_1:dtE/(nEpts-1):tE2_2;

eta_E1_sim_TLS = interp1(time,eta_sim,tE1_pts); % (1) 'eta_s' at 'empty-cell' times   [ppbv]
eta_F1_sim_TLS = interp1(time,eta_sim,tF1_pts); % (1) 'eta_s' at 'full-cell' times    [ppbv]
eta_F2_sim_TLS = interp1(time,eta_sim,tF2_pts); % (2) 'eta_s' at 'full-cell' times    [ppbv]
eta_E2_sim_TLS = interp1(time,eta_sim,tE2_pts); % (2) 'eta_s' at 'empty-cell' times    [ppbv]

eta_E1_sim     = eta_E1_sim_TLS' + sig_E_TLS*randn(nEpts,1);
eta_F1_sim     = eta_F1_sim_TLS' + sig_F_TLS*randn(nFpts,1);
eta_F2_sim     = eta_F2_sim_TLS' + sig_F_TLS*randn(nFpts,1);
eta_E2_sim     = eta_E2_sim_TLS' + sig_E_TLS*randn(nEpts,1);

%% Statistics
mean_eta_E1_sim  = mean(eta_E1_sim);
mean_eta_F1_sim  = mean(eta_F1_sim);
mean_eta_F2_sim  = mean(eta_F2_sim);
mean_eta_E2_sim  = mean(eta_E2_sim);

delta_eta_E1_sim = std(eta_E1_sim);
delta_eta_F1_sim = std(eta_F1_sim);
delta_eta_F2_sim = std(eta_F2_sim);
delta_eta_E2_sim = std(eta_E2_sim);

mean_tE1 = mean(tE1_pts);
mean_tF1 = mean(tF1_pts);
mean_tF2 = mean(tF2_pts);
mean_tE2 = mean(tE2_pts);

disp(['Step1 / empty-cell runs: ',num2str(mean_eta_E1_sim,'%2.2f'), ...
    ' +/- ',num2str(delta_eta_E1_sim,'%2.2f'),' ppbv'])
disp(['Step1 / full-cell runs:  ',num2str(mean_eta_F1_sim,'%2.2f'), ...
    ' +/- ',num2str(delta_eta_F1_sim,'%2.2f'),' ppbv'])
disp(['Step2 / full-cell runs:  ',num2str(mean_eta_F2_sim,'%2.2f'), ...
    ' +/- ',num2str(delta_eta_F2_sim,'%2.2f'),' ppbv'])
disp(['Step2 / empty-cell runs: ',num2str(mean_eta_E2_sim,'%2.2f'), ...
    ' +/- ',num2str(delta_eta_E2_sim,'%2.2f'),' ppbv'])

disp(' ')

%% Step 1
eta_HC_1_sim = mean_eta_F1_sim - mean_eta_E1_sim;
sig_HC_1_sim = sqrt(delta_eta_F1_sim^2/nFpts + delta_eta_E1_sim^2/nEpts);

eta_HC_2_sim = mean_eta_F2_sim - mean_eta_E2_sim;
sig_HC_2_sim = sqrt(delta_eta_F2_sim^2/nFpts + delta_eta_E2_sim^2/nEpts);

disp(['Step1 : ',num2str(eta_HC_1_sim,'%2.2f'),' +/- ',num2str(sig_HC_1_sim,'%2.2f'),' ppbv'])
disp(['Step2 : ',num2str(eta_HC_2_sim,'%2.2f'),' +/- ',num2str(sig_HC_2_sim,'%2.2f'),' ppbv'])

disp(['Step1 : ',num2str(eta_HC_1_sim/n_E,'%2.2f'),' +/- ',num2str(sig_HC_1_sim/n_E,'%2.2f'),' ppbv'])
disp(['Step2 : ',num2str(eta_HC_2_sim/n_E,'%2.2f'),' +/- ',num2str(sig_HC_2_sim/n_E,'%2.2f'),' ppbv'])