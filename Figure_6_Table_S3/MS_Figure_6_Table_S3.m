% ===================
% Figure 6 + Table S3
% ===================
%
% This script produces Figure 6 of the manuscript.
%
% Analysis of the trends of retrieved full-cell and empty-cell CH4 vmr
% during 18 TLS experiment documented in Webster et al. (2015, 2021).
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

%% Figure: main infos
figid      = 911;
sfigtype   = 'png';
ftsz       = 12;
ftnm       = 'times';
savefig    = 1;
savedata   = 1;
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

%% Define variables
a1T_F_list = zeros(1,nsol);
a1T_E_list = zeros(1,nsol);
eta_list   = zeros(1,nsol);
sigma_list = zeros(1,nsol);

%% Loop over TLS experiments
for isol = 1:nsol
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
    
    %% Key times
    tF_1 = F_elapsed_t(1);               % first full-cell run  [s]
    tF_2 = F_elapsed_t(end);             % last full-cell run   [s]
    tE_1 = E_elapsed_t(1);               % first empty-cell run [s]
    tE_2 = E_elapsed_t(end);             % last empty-cell run  [s]
    
    %% Calculation of the slopes
    
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
    
    a1T_F_list(isol) = a1T_F*3600; % full-cell runs   [ppbv h-1]
    a1T_E_list(isol) = a1T_E*3600; % empty-cell runs  [ppbv h-1]
end

%% Make figure
figure(figid)
delete(get(gcf,'children'))
set(gcf,'paperpositionmode','auto')
set(gcf,'position',[200 200 1200 500])
set(gca,'fontsize',ftsz,'fontname',ftnm)
nrow = 1; ncol = 2; isub = 1;

dxbin = 5;

%% All TLS experiments: full-cell runs
subplot(nrow,ncol,isub); isub = isub + 1;
clear legname
ileg    = 1;
legname = cell(1,3);

collist = collistF;
hhF = histogram(a1T_F_list,'BinWidth',dxbin); hold on
set(hhF,'FaceColor',collist,'EdgeColor',collist,'EdgeAlpha',0)
legname{ileg} = 'TLS data (full-cell runs)'; ileg = ileg + 1;

N_TLS     = nsol;
x_a       = -40:dxbin/100:40;

% TLS: mean and standard deviation + statistical distribution 'f_a1T_F'
mu_F_TLS  = mean(a1T_F_list);
sig_F_TLS = std(a1T_F_list);
f_F_TLS   = N_TLS/(sqrt(2*pi)*sig_F_TLS)*exp(-0.5*(mu_F_TLS-x_a).^2/(sig_F_TLS)^2);

plot(x_a,f_F_TLS*dxbin,'color',collist,'linewidth',2); hold on
text_leg = '$$N_{\rm TLS} \; f(a|\mu_{{\rm F, TLS}},\sigma_{{\rm F, TLS}}^2) \; \Delta a $$';
legname{ileg} = text_leg; ileg = ileg + 1;

% Model: mean and standard deviation + statistical distribution 'f_a1T_F_mod'
mu_F_mod  = 3.01;
sig_F_mod = 8.91;
% mu_F_mod  = 3.46;
% sig_F_mod = 8.9;

f_F_mod   = N_TLS/(sqrt(2*pi)*sig_F_mod)*exp(-0.5*(mu_F_mod-x_a).^2/(sig_F_mod)^2);

plot(x_a,f_F_mod*dxbin,'color',collist/2,'LineStyle','-','linewidth',2); hold on
text_leg = '$$N_{\rm TLS} \; f(a|\mu_{{\rm F,m}},\sigma_{{\rm F,m}}^2) \; \Delta a $$';
legname{ileg} = text_leg;

hlF       = legend(legname,'location','north','fontsize',ftsz,'fontname',ftnm,'interpreter','latex');
h_posF    = get(hlF,'position');
h_posF(2) = h_posF(2)-0.03;
h_posF(4) = h_posF(4)+0.03;
set(hlF,'position',h_posF);

% Axes and labels
xlabel('trend {\ita} [ppbv h^{-1}]','fontsize',ftsz,'fontname',ftnm,'Interpreter','tex')
ylabel('count','fontsize',ftsz,'fontname',ftnm)

xlim([-40 40])
ylim([0 8])

frac_pos_F = length(find(a1T_F_list>0))/nsol*100;

set(gca,'XTick',-40:10:40,'XMinorTick','on','YMinorTick','on','box','on')
set(gca,'fontsize',ftsz,'fontname',ftnm)

SS_mk_label_panel(lab_panel{isub-1},xlim,'lin',ylim,'lin',ftsz,ftnm)

disp(['N_TLS      = ',num2str(N_TLS)])
disp(['da         = ',num2str(dxbin),' ppbv h-1'])
disp(['mu(F,TLS)  = ',num2str(mu_F_TLS,'%2.2f'),' ppbv h-1'])
disp(['sig(F,TLS) = ',num2str(sig_F_TLS,'%2.2f'),' ppbv h-1'])
disp(['mu(F,mod)  = ',num2str(mu_F_mod,'%2.2f'),' ppbv h-1'])
disp(['sig(F,mod) = ',num2str(sig_F_mod,'%2.2f'),' ppbv h-1'])

%% All TLS experiments: empty-cell runs
subplot(nrow,ncol,isub); isub = isub + 1;
clear legname
ileg    = 1;
legname = cell(1,3);

collist = collistE;
hhE = histogram(a1T_E_list,'BinWidth',dxbin); hold on
set(hhE,'FaceColor',collist,'EdgeColor',collist,'EdgeAlpha',0)
legname{ileg} = 'TLS data (empty-cell runs)'; ileg = ileg + 1;

mu_E_TLS  = mean(a1T_E_list);
sig_E_TLS = std(a1T_E_list);

f_E_TLS   = N_TLS/(sqrt(2*pi)*sig_E_TLS)*exp(-0.5*(mu_E_TLS-x_a).^2/(sig_E_TLS)^2);

plot(x_a,f_E_TLS*dxbin,'color',collist,'linewidth',2); hold on
text_leg = '$$N_{\rm TLS} \; f(a|\mu_{{\rm E, TLS}},\sigma_{{\rm E, TLS}}^2) \; \Delta a $$';
legname{ileg} = text_leg; ileg = ileg + 1;

mu_E_mod  = 3.01;
sig_E_mod = 6.21;
% mu_E_mod  = 3.45;
% sig_E_mod = 6.2;
f_E_mod   = N_TLS/(sqrt(2*pi)*sig_E_mod)*exp(-0.5*(mu_E_mod-x_a).^2/(sig_E_mod)^2);

plot(x_a,f_E_mod*dxbin,'color',collist/2,'LineStyle','-','linewidth',2); hold on
text_leg = '$$N_{\rm TLS} \; f(a|\mu_{{\rm E,m}},\sigma_{{\rm E,m}}^2) \; \Delta a $$';
legname{ileg} = text_leg;

hlE       = legend(legname,'location','north','fontsize',ftsz,'fontname',ftnm,'interpreter','latex');
h_posE    = get(hlE,'position');
h_posE(2) = h_posE(2)-0.03;
h_posE(4) = h_posE(4)+0.03;
set(hlE,'position',h_posE);

% Axes and labels
xlabel('trend {\ita} [ppbv h^{-1}]','fontsize',ftsz,'fontname',ftnm,'Interpreter','tex')
ylabel('count','fontsize',ftsz,'fontname',ftnm)

xlim([-40 40])
ylim([0 8])

frac_pos_E = length(find(a1T_E_list>0))/nsol*100;

set(gca,'XTick',-40:10:40,'XMinorTick','on','YMinorTick','on','box','on')

set(gca,'fontsize',ftsz,'fontname',ftnm)

SS_mk_label_panel(lab_panel{isub-1},xlim,'lin',ylim,'lin',ftsz,ftnm)

disp(' ')
disp(['mu(E,TLS)  = ',num2str(mu_E_TLS,'%2.2f'),' ppbv h-1'])
disp(['sig(E,TLS) = ',num2str(sig_E_TLS,'%2.2f'),' ppbv h-1'])
disp(['mu(E,mod)  = ',num2str(mu_E_mod,'%2.2f'),' ppbv h-1'])
disp(['sig(E,mod) = ',num2str(sig_E_mod,'%2.2f'),' ppbv h-1'])

%% Save figure
if (savefig == 1)
    sfpath     = '../Figures/';
    sffilename = 'Figure_06.png';
    sfflnm     = fullfile(sfpath,sffilename);
    exportgraphics(gcf,sfflnm,'Resolution',400)
end

%% Save data
if (savedata == 1)
    % Save leak output
    T = table(sol_list',a1T_F_list',a1T_E_list');
    
    sdpath     = '../output_tables/';
    sdfilename = 'table_S3.xlsx';
    sdflnm     = fullfile(sdpath,sdfilename);
    writetable(T,sdflnm,'Sheet',1,'Range','A1')
end

toc