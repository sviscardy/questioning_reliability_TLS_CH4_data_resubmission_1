% =========
% Figure S2
% =========
% 
% This script fits the number of CH4 nanomoles measured in the FO chamber
% through the Curiosity mission.
% 
% Two periods can be defined:
%   - Period 1: Sol 450 - Sol 1000
%   - Period 2: Sol 1000 - Sol 2350
% both starting when the FO chamber is pumped out. The initial conditions
% are:
%   - Period 1: N_FO(0) = 0.75 nmol
%   - Period 2: N_FO(0) = 0 nmol
% 
% Further details provided in the Supplementary Information, Section S4.2.2
% 
% Author: sebastien.viscardy@aeronomie.be
% 
%% Link to functions
addpath('../functions/');

%%
clearvars
disp(' ')

%% Physical constants
marscst

%% Parameters
tau     = 150;       % characteristic time                     [sols]
Ninf    = 2.4;       % steady-state number of CH4 amount in FO [nmol]
F_CH4   = Ninf/tau;  % release rate of CH4 in FO               [nmol sol-1]

disp(' ')
disp('Parameters:')
disp('-----------')
disp(['tau     = ',num2str(tau,'%2.0f'),' sols'])
disp(['N(inf.) = ',num2str(Ninf,'%2.1f'),' nmol'])
disp(['F_CH4   = ',num2str(F_CH4,'%2.1e'),' nmol sol-1'])

%% Estimate of time t_H
V_FO          = 988.2e-6;                    % FO volume                  [m3]
T             = 320;                         % temperature                [K]
tau_s         = tau*sol_s;                   % timescale (150 sols)       [s]

psi_CH4_tot   = V_FO/(tau_s*rgas*T);         % parameter                  [mol Pa-1 s-1]
psi_CH4_FO_HC = 1.01e-14;                    % transport coeff. (FO-cell) [mol Pa-1 s-1]
psi_CH4_FO_MA = psi_CH4_tot - psi_CH4_FO_HC; % transport coeff. (FO-Mars) [mol Pa-1 s-1]

disp(' ')
disp('Results:')
disp('--------')
disp(['psi (tot.)     = ',num2str(psi_CH4_tot,'%2.2e'),' mol Pa-1 s-1'])
disp(['psi (FO-Hcell) = ',num2str(psi_CH4_FO_HC,'%2.2e'),' mol Pa-1 s-1'])
disp(['psi (FO-Mars)  = ',num2str(psi_CH4_FO_MA,'%2.2e'),' mol Pa-1 s-1'])

%% Estimate of the amount of CH4 molecules released in the FO chamber
N_tot = F_CH4*2644;           % total number of nanmoles CH4 [nmol]
M_tot = M_CH4*N_tot*1e-9; % total mass of CH4            [kg]

disp(['N_tot   = ',num2str(N_tot,'%2.1f'),' nmol'])
disp(['M_tot   = ',num2str(M_tot,'%2.2e'),' kg'])
disp(['M_tot   = ',num2str(M_tot*1e12,'%2.0f'),' ng'])

%% Figure: main infos
figid    = 401;
sfigtype = 'png';
ftsz     = 12;
ftnm     = 'times';
savefig  = 1;
savedata = 0;

%% Load FO data (nanomoles of CH4 in FO vs. Sols)
curr_dir = pwd;
cd data/
SS_CH4_nanomoles_in_FO
cd(curr_dir)

FO_CH4_sol_list  = CH4_nanomoles_in_FO(:,1);
FO_CH4_nmol_list = CH4_nanomoles_in_FO(:,2);

%% Time
time  = 0:1500; % [sols]

t_sh0 = 305;    % FO chamber pumped out [sols]
t_sh1 = 450;    % FO chamber pumped out [sols]
t_sh2 = 1000;   % FO chamber pumped out [sols]

%% Period 1: Sol 450 - Sol 1000
N0_1  = 0.75; 
N_1   = N0_1*exp(-time/tau) + tau*F_CH4 * (1-exp(-time/tau));

%% Period 2: Sol 1000 - 2350
N0_2  = 0;
N_2   = N0_2*exp(-(time)/tau) + tau*F_CH4 * (1-exp(-(time)/tau));

%% Make figure
close all
lfflnm = '~/mars/matlab/Figures_from_papers/Webster_Science_2018/Webster_Science_2018_FigS43b.fig';
open(lfflnm)
pause(1)
legname = cell(1,4);
ileg    = 1;

% Colors
collist1   = [0.8 0 0];
collist2   = [0 0.8 0];
collistpts = [0.8 0 0.8];

% Data points calculated in this paper
kk = FO_CH4_sol_list <= 684 | FO_CH4_sol_list >= 2442;
plot(FO_CH4_sol_list(kk),FO_CH4_nmol_list(kk),'Marker','o','MarkerSize',6, ...
    'MarkerFaceColor',collistpts,'MarkerEdgeColor',collistpts,'LineStyle','none'); hold on
legname{ileg} = 'calculated in this study'; ileg = ileg + 1;

% Data points taken from Figure S43B in the SM for Webster et al. (2018)
kk = FO_CH4_sol_list > 684 & FO_CH4_sol_list < 2442;
plot(FO_CH4_sol_list(kk),FO_CH4_nmol_list(kk),'Marker','v','MarkerSize',6, ...
    'MarkerFaceColor',collistpts,'MarkerEdgeColor',collistpts,'LineStyle','none'); hold on
legname{ileg} = 'taken from Webster et al. (2018)'; ileg = ileg + 1;

% production-loss model for Period 1
plot(time+t_sh1,N_1,'color',collist1,'linewidth',2); hold on
txt1 = ['$$N_{\rm FO}(t) = ',num2str(N0_1),' \exp \left (-\frac{t-', ...
    num2str(t_sh1),'}{\tau} \right ) + F_{\rm CH_4} \tau \left[ 1 - \exp \left (-\frac{t-', ...
    num2str(t_sh1),'}{\tau} \right ) \right]$$'];
legname{ileg} = txt1; ileg = ileg + 1;

% production-loss model for Period 2
plot(time+t_sh2,N_2,'color',collist2,'linewidth',2); hold on
txt2 = ['$$N_{\rm FO}(t) = F_{\rm CH_4} \tau \left[ 1 - \exp \left (-\frac{t-', ...
    num2str(t_sh2),'}{\tau} \right ) \right]$$'];
legname{ileg} = txt2;

hl = legend(legname,'location','north','interpreter','latex','fontsize',ftsz+2,'fontname',ftnm);
set(hl,'position',[0.28167 0.81554 0.54166 0.17488])

set(gca,'XTick',[],'YTick',[])

%% Save figure
if (savefig == 1)
    sfpath     = '';
    sffilename = 'FO_CH4_amount_production_loss_model_vs_Webster_Science_2018_Fig_S43B.png';
    sfflnm     = fullfile(sfpath,sffilename);
    exportgraphics(gcf,sfflnm,'Resolution',400)
end

%% Make figure 2: new figure
figure(figid)
delete(get(gcf,'children'))
set(gcf,'paperpositionmode','auto')
set(gcf,'position',[200 200 1200 600])
set(gca,'fontsize',ftsz,'fontname',ftnm)
legname = cell(1,4);
ileg    = 1;

% Data points calculated in this paper
kk = FO_CH4_sol_list <= 684 | FO_CH4_sol_list >= 2442;
plot(FO_CH4_sol_list(kk),FO_CH4_nmol_list(kk),'Marker','o','MarkerSize',6, ...
    'MarkerFaceColor',collistpts,'MarkerEdgeColor',collistpts,'LineStyle','none'); hold on
legname{ileg} = 'calculated in this study'; ileg = ileg + 1;

% Data points taken from Figure S43B in the SM for Webster et al. (2018)
kk = FO_CH4_sol_list > 684 & FO_CH4_sol_list < 2442;
plot(FO_CH4_sol_list(kk),FO_CH4_nmol_list(kk),'Marker','v','MarkerSize',6, ...
    'MarkerFaceColor',collistpts,'MarkerEdgeColor',collistpts,'LineStyle','none'); hold on
legname{ileg} = 'taken from Webster et al. (2018)'; ileg = ileg + 1;

% production-loss model for Period 1
plot(time+t_sh1,N_1,'color',collist1,'linewidth',2); hold on
legname{ileg} = txt1; ileg = ileg + 1;

% production-loss model for Period 2
plot(time+t_sh2,N_2,'color',collist2,'linewidth',2); hold on
legname{ileg} = txt2; ileg = ileg + 1;

xline(t_sh0,'--k','linewidth',2,'handlevisibility','off'); hold on
xline(t_sh1,'--k','linewidth',2,'handlevisibility','off'); hold on
xline(t_sh2,'--k','linewidth',2,'handlevisibility','off'); hold on

legend(legname,'interpreter','latex','fontsize',ftsz+2,'fontname',ftnm)

xlabel('time since Curiosity''s landing [sols]','fontsize',ftsz,'fontname',ftnm)
ylabel('CH_4 amount in the FO chamber {\itN}_{FO} [nmol]','fontsize',ftsz,'fontname',ftnm)

set(gca,'XMinorTick','on','YMinorTick','on','box','on')
set(gca,'fontsize',ftsz,'fontname',ftnm)

title(['$$\tau = ',num2str(tau),'\; {\rm sols} \; ; \; F_{\rm CH_4} = ', ...
    num2str(F_CH4,'%2.3f'),'\; {\rm nmol \; sol^{-1}}$$'], ...
    'fontsize',ftsz+4,'fontname',ftnm,'FontWeight','bold','interpreter','latex')

%% Save figure
if (savefig == 1)
    sfpath     = '../Figures/';
    sffilename = 'Figure_S02.png';
    sfflnm     = fullfile(sfpath,sffilename);
    exportgraphics(gcf,sfflnm,'Resolution',400)
end