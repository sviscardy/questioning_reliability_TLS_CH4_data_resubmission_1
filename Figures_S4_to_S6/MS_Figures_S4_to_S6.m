% ===============
% Figures S4 - S6
% ===============
% 
% Assessment of the stationary-flow assumption ==> Figures S4 - S6
% 
% The solution-diffusion mechanism is considered to describe the transport
% of CO2 and CH4 through the O-ring, using Fick's second law:
% 
%   dn_i(x,t)/dt = D_i d^2 n_i(x,t)/dx^2
% 
% where
% 
%   - t   : time
%   - x   : spatial coordinate (0 <= x <= L) 
%   - L   : thickness of the O-ring
%   - n_i : concentration of species i across the O-ring
%   - D_i : diffusion coefficient
% 
% Further details provided in the Supplementary Information, Section S4.3.8.
% 
% Author: sebastien.viscardy@aeronomie.be
%%
clearvars

tic

%% Physical constants
physcst_ref

%% Select gas species
spnm_list = {'co2' 'ch4'}; nsp = length(spnm_list);
L_list    = [1 10]*1e-3;   nL  = length(L_list);

%% Figure: main infos
figid      = 625;
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
lab_figs   = {'S05' 'SXX' 'S06' 'S04'};
jlab       = 0;

%% Type of experiment
t_exp_list = {'D' 'D' 'D' 'D' 'D' 'D' 'D' 'D' 'D' 'D' 'E' 'E' 'D' 'E' 'E' 'E' 'E' 'E'};
sol_list   = [79 81 106 292 306 313 466 474 504 526 573 684 684 2442 2446 2615 2627 2644];
nsol       = length(sol_list);

%% Select experiments
isol        = 12; % E / Sol 684
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
eta_FO  = p_HC/p_FO*l_HC/l_FO*eta_E; % TLS: mean CH4 vmr in FO      [/]
N_FO    = p_HC/(rgas*T) ...
    *l_HC/l_FO*eta_E*V_FO;           % TLS: CH4 amount in FO        [mol]
N_HC    = eta_HC*p_HC*V_HC/(rgas*T); % TLS: CH4 amount in HC        [mol]

%% Key times
tF_1     = F_elapsed_t(1);           % first full-cell run          [s]
tF_2     = F_elapsed_t(end);         % last full-cell run           [s]
tE_1     = E_elapsed_t(1);           % first empty-cell run         [s]
tE_2     = E_elapsed_t(end);         % last empty-cell run          [s]

t_0_s    = tF_1-dt_samp;             % start of air sampling        [s]
t_0_h    = t_0_s/3600;               % start of air sampling        [h]

time_sel = [t_0_h 1.5 2 3 4 4.7 5];  % time selected (conc. prf.)   [h]
nt_sel   = length(time_sel);

%% Loop over thickness L
for iL = 1:nL
    L = L_list(iL); % thickness [m]
    disp(['L = ',num2str(L),' m'])
    
    %% Parameters
    psi_CH4      = 1.16e-14;          % transport coeff. for CH4   [mol Pa-1 s-1]
    A            = 1000e-6;           % surface area               [m2]
    
    % CH4
    K_CH4        = psi_CH4*L/A;       % CH4 permeability           [mol m-1 Pa-1 s-1]
    S_CH4        = 1e-4;              % CH4 solubility coefficient [mol m-3 Pa-1]
    D_CH4        = K_CH4/S_CH4;       % CH4 diffusion coefficient  [m2 s1]
    
    % CO2
    omega        = 10;                % CO2/CH4 selectivity        [/]
    psi_CO2      = omega*psi_CH4;     % transport coeff. for CO2   [mol Pa-1 s-1]
    K_CO2        = psi_CO2*L/A;       % CO2 permeability           [mol m-1 Pa-1 s-1]
    S_CO2        = S_CH4;             % CO2 solubility coefficient [mol m-3 Pa-1]
    D_CO2        = K_CO2/S_CO2;       % CO2 diffusion coefficient  [m2 s1]
    
    p_CO2_FO     = p_FO;              % CO2 part. prs in FO        [Pa]
    p_CH4_FO     = eta_FO*p_FO;       % CH4 part. prs in FO        [Pa]
    
    p_CO2_2_end  = p_HC;              % final pressure 2           [Pa]
    tau_samp     = 1200;              % sampling timescale         [s]
    
    dtime_ini    = dt_samp + (tF_2 - tF_1)/2;                 % time diff. ini/mid-full [s]
    N_CH4_HC_0   = N_HC - psi_CH4*N_FO*rgas*T/V_FO*dtime_ini; % initial CH4 amount      [mol]
    if (N_CH4_HC_0 < 0), error('initial CH4 amount < 0 nmol'), end
    
    N_CO2_HC_0   = omega*N_CH4_HC_0/eta_FO;                   % Initial CO2 amount      [mol]
    p_HC_0       = (N_CO2_HC_0+N_CH4_HC_0)*rgas*T/V_HC;       % Initial pressure        [Pa]
    p_CO2_HC     = p_HC_0;
    
    eta_CH4_HC_0 = N_CH4_HC_0/(N_CO2_HC_0 + N_CH4_HC_0);      % Init. CH4 vmr           [/]
    
    p_CH4_HC     = eta_CH4_HC_0*p_HC_0;                       % Init. CH4 p.p. in HC    [Pa]
    
    %% Print key information
    disp(['psi_CH4        = ',num2str(psi_CH4),' mol Pa-1 s-1'])
    disp(['K_CH4          = ',num2str(K_CH4),' mol m-1 Pa-1 s-1'])
    disp(['D_CH4          = ',num2str(D_CH4),' m2 s-1'])
    disp(['S_CH4          = ',num2str(S_CH4),' mol m-3 Pa-1'])
    disp(' ')
    disp(['psi_CO2        = ',num2str(psi_CO2),' mol Pa-1 s-1'])
    disp(['K_CO2          = ',num2str(K_CO2),' mol m-1 Pa-1 s-1'])
    disp(['D_CO2          = ',num2str(D_CO2),' m2 s-1'])
    disp(['S_CO2          = ',num2str(S_CO2),' mol m-3 Pa-1'])
    disp(' ')
    disp(['thickness L    = ',num2str(L),' m'])
    disp(['surface area A = ',num2str(A),' m2'])
    disp(['p_HC_0         = ',num2str(p_HC_0),' Pa'])
    disp(' ')
    
    %% discretization parameters
    dt          = 0.1;               % timestep                   [s]
    t_end       = tE_2;              % final time                 [s]
    Nt          = t_end/dt;          % number of timesteps
    
    Nx          = 51;                % number of grid points
    x           = linspace(0,L,Nx);  % spatial grid               [m]
    dx          = x(2) - x(1);       % grid cell size             [m]
    
    disp(['dt = ',num2str(dt),' s'])
    disp(['dx = ',num2str(dx),' m'])
    
    %% Verification of the Courant stability criterion for the explicit method
    alpha_CO2 = D_CO2 * dt / dx^2;
    disp(['alpha_CO2 = ',num2str(alpha_CO2)])
    
    if alpha_CO2 > 0.5
        error(['Reduce timestep for stability (alpha_CO2 = ',num2str(alpha_CO2),' > 0.5)']);
    end
    
    alpha_CH4 = D_CH4 * dt / dx^2;
    disp(['alpha_CH4 = ',num2str(alpha_CH4)])
    if alpha_CH4 > 0.5
        error(['Reduce timestep for stability (alpha_CH4 = ',num2str(alpha_CH4),' > 0.5)']);
    end
    
    %% Define concentration profile
    c_CH4     = zeros(Nx, 1); % initial concentration [mol m-3]
    c_CH4_new = c_CH4;
    
    c_CO2     = zeros(Nx, 1); % initial concentration [mol m-3]
    c_CO2_new = c_CO2;
    
    %% Boundary conditions
    c_CO2(1)  = p_CO2_FO * S_CO2; % Concentration at x = 0 [mol m-3]
    c_CO2(Nx) = p_CO2_HC * S_CO2; % Concentration at x = L [mol m-3]
    
    c_CH4(1)  = p_CH4_FO * S_CH4; % Concentration at x = 0 [mol m-3]
    c_CH4(Nx) = p_CH4_HC * S_CH4; % Concentration at x = L [mol m-3]
    
    %% Initialization of concentration profile
    dc_dx = (c_CO2(Nx) - c_CO2(1))/L;
    c_CO2 = c_CO2(1) + dc_dx * x;
    
    dc_dx = (c_CH4(Nx) - c_CH4(1))/L;
    c_CH4 = c_CH4(1) + dc_dx * x;
    
    %% Define variable
    time_h          = zeros(1,Nt);
    c_CO2_list      = zeros(1,Nx);
    p_CO2_HC_list   = zeros(1,Nt);
    N_CO2_list      = zeros(1,Nt);
    c_CH4_list      = zeros(1,Nx);
    p_CH4_HC_list   = zeros(1,Nt);
    N_CH4_list      = zeros(1,Nt);
    flux_CO2_list   = zeros(1,Nt);
    flux_CO2_s_list = zeros(1,Nt);
    flux_CH4_list   = zeros(1,Nt);
    flux_CH4_s_list = zeros(1,Nt);
    
    %% Loop over time
    grad_c_CO2 = gradient(c_CO2);
    dN_CO2_dt  = -D_CO2*A*grad_c_CO2(Nx)/dx;
    N_CO2_HC   = N_CO2_HC_0;
    grad_c_CH4 = gradient(c_CH4);
    dN_CH4_dt  = -D_CH4*A*grad_c_CH4(Nx)/dx;
    N_CH4_HC   = N_CH4_HC_0;
    jsel       = 1;
    for it = 1:Nt
        time       = dt*it;     % [s]
        time_h(it) = time/3600; % [h]
        if (time<t_0_s), continue, end
        
        %% Cell pumped out: Parameters
        delta_p_pump = p_HC_0-p_CO2_2_end;         % pressure difference   [Pa]
        delta_N_pump = V_HC/(rgas*T)*delta_p_pump; % air amount difference [mol]
        dt_pump      = tE_1-tF_2;                  % evacuation duration   [s]
        
        %% Update concentrations using the explicit finite difference method
        for i = 2:Nx-1
            c_CO2_new(i) = c_CO2(i) + alpha_CO2 * (c_CO2(i+1) - 2*c_CO2(i) + c_CO2(i-1));
            c_CH4_new(i) = c_CH4(i) + alpha_CH4 * (c_CH4(i+1) - 2*c_CH4(i) + c_CH4(i-1));
        end
        
        %% Updating concentrations at the boundaries
        c_CO2_new(1) = c_CO2(1);
        c_CH4_new(1) = c_CH4(1);
        N_CH4_HC     = N_CH4_HC + dN_CH4_dt*dt;
        p_CH4_HC     = N_CH4_HC*rgas*T/V_HC;
        
        eta_CH4  = N_CH4_HC/(N_CO2_HC + N_CH4_HC);
        
        if (time >= t_0_s && time < tF_1)   % Air sampling
            p_CO2_HC = p_CO2_HC + (p_CO2_2_end - p_CO2_HC)/tau_samp*dt;
            N_CO2_HC = p_CO2_HC*V_HC/(rgas*T);
        elseif (time > tF_2 && time < tE_1) % Air evacuated
            N_CO2_HC = N_CO2_HC + delta_N_pump/dt_pump*dt;
            p_CO2_HC = N_CO2_HC*rgas*T/V_HC;
            N_CH4_HC = N_CH4_HC + eta_CH4*delta_N_pump/dt_pump*dt;
        else
            N_CO2_HC = N_CO2_HC + dN_CO2_dt*dt;
            p_CO2_HC = N_CO2_HC*rgas*T/V_HC;
        end
        
        %% Store data
        p_CO2_HC_list(it) = p_CO2_HC;
        N_CO2_list(it)    = N_CO2_HC;
        
        p_CH4_HC_list(it) = p_CH4_HC;
        N_CH4_list(it)    = N_CH4_HC;
        
        c_CH4_new(Nx)     = p_CH4_HC * S_CH4; % CH4 concentration at x = L
        c_CO2_new(Nx)     = p_CO2_HC * S_CO2; % CO2 concentration at x = L
        
        %% Update concentrations
        c_CO2      = c_CO2_new;
        c_CH4      = c_CH4_new;
        
        %% Select concentration profiles
        if time == time_sel(jsel)*3600
            disp(['time selected = ',num2str(time/3600),' h'])
            c_CO2_list(jsel,:) = c_CO2;
            c_CH4_list(jsel,:) = c_CH4;
            jsel = jsel + 1; jsel = min(jsel,nt_sel);
        end
        
        %% Molar fluxes
        grad_c_CO2          = gradient(c_CO2);
        dN_CO2_dt           = -D_CO2*A*grad_c_CO2(Nx)/dx;
        
        grad_c_CH4          = gradient(c_CH4);
        dN_CH4_dt           = -D_CH4*A*grad_c_CH4(Nx)/dx;
        
        flux_CO2_s          = D_CO2*A/L*(c_CO2(1)-c_CO2(Nx));
        flux_CH4_s          = D_CH4*A/L*(c_CH4(1)-c_CH4(Nx));
        
        flux_CO2_list(it)   = dN_CO2_dt;
        flux_CO2_s_list(it) = flux_CO2_s;
        
        flux_CH4_list(it)   = dN_CH4_dt;
        flux_CH4_s_list(it) = flux_CH4_s;
    end
    
    %% Loop over gas species
    for isp = 1:nsp
        jlab = jlab + 1;
        
        spnm = spnm_list{isp};
        
        %% Make figure
        figure(figid)
        delete(get(gcf,'children'))
        set(gcf,'paperpositionmode','auto')
        set(gcf,'position',[400 0 900 900])
        set(gca,'fontsize',ftsz,'fontname',ftnm)
        nrow = 3; ncol = 1; isub = 1;
        switch spnm
            case 'co2'
                c_prf          = c_CO2_list;
                p_part_HC_list = p_CO2_HC_list;
                flux_list      = flux_CO2_list;
                flux_s_list    = flux_CO2_s_list;
                
            case 'ch4'
                c_prf          = c_CH4_list;
                p_part_HC_list = p_CH4_HC_list;
                flux_list      = flux_CH4_list;
                flux_s_list    = flux_CH4_s_list;
        end
        
        kk                 = time_h*3600<t_0_s;
        p_part_HC_list(kk) = NaN;
        flux_list(kk)      = NaN;
        flux_s_list(kk)    = NaN;
        
        %% Concentration across the O-ring seal
        subplot(nrow,ncol,isub); isub = isub + 1;
        legname = cell(1,nt_sel);
        
        if exist('h1', 'var'), delete(h1); end
        if exist('h2', 'var'), delete(h2); end
        
        collist = jet(nt_sel);
        
        for jsel = 1:nt_sel
            h1 = plot(x*1e3,c_prf(jsel,:),'color',collist(jsel,:),'linewidth',2); hold on
            if (time_sel(jsel)*3600 == t_0_s)
                legname{jsel} = '{\itt}_0';
            else
                legname{jsel} = [num2str(time_sel(jsel)),' h'];
            end
        end
        
        if (L == 1e-3)
            xl  = [-0.5 L*1e3+0.5];
            xt1 = [0 0.5]+L*1e3;
            xt2 = [0 L*1e3];
            xt3 = [-0.5 0];
            yl  = [0 50e-3];
        elseif (L == 10*1e-3)
            xl  = [-5 L*1e3+5];
            xt1 = [0 5]+L*1e3;
            xt2 = [0 L*1e3];
            xt3 = [-5 0];
            yl  = [0 50e-3];
        end
        
        switch spnm
            case 'co2'
                yl  = [0 50e-3];
            case 'ch4'
                yl  = [0 6e-7];
        end
        
        xlim(xl)
        ylim(yl)
        
        col_seal = [1 1 1]*0.6;
        ha10 = area([0 0 L L]*1e3,[ylim fliplr(ylim)]); hold on
        set(ha10,'facecolor',col_seal,'edgecolor',col_seal,'FaceAlpha',0.2,'EdgeAlpha',0, ...
            'LineStyle','none','HandleVisibility','off')
        ha10.BaseLine.LineStyle = 'none';
        ha10.BaseValue = yl(1);
        uistack(ha10,'bottom');
        
        legend(legname,'location','southeast','fontsize',ftsz,'fontname',ftnm)
        
        text(mean(xt1),yl(2)*0.9,'Herriott cell', ...
            'fontsize',ftsz,'fontname',ftnm,'FontWeight','bold','HorizontalAlignment','center')
        
        text(mean(xt2),yl(2)*0.9,'O-ring', ...
            'fontsize',ftsz,'fontname',ftnm,'FontWeight','bold','HorizontalAlignment','center')
        
        text(mean(xt3),yl(2)*0.9,'FO chamber', ...
            'fontsize',ftsz,'fontname',ftnm,'FontWeight','bold','HorizontalAlignment','center')
        
        xlabel('Position x [mm]','fontsize',ftsz,'fontname',ftnm);
        ylabel(sprintf([SS_splab(spnm),' concentration\n {\\itn}_{', ...
            SS_splab(spnm),'}({\\itx,t}) [mol m^{-3}]']),'fontsize',ftsz,'fontname',ftnm);
        
        set(gca,'XDir','reverse','XMinorTick','on','YMinorTick','on','box','on')
        set(gca,'fontsize',ftsz,'fontname',ftnm)
        
        x1   = xl(1); x2 = xl(2);
        xlab = x2 - (x2-x1)/20;
        
        yl   = ylim;
        y1   = yl(1); y2 = yl(2);
        ylab = y2+(y2-y1)/6;
        
        ht = text(xlab,ylab,lab_panel{isub-1});
        set(ht,'fontsize',ftsz+2,'FontWeight','bold','fontname',ftnm,'HorizontalAlignment','center')
        
        %% Time evolution of pressure p_HC
        subplot(nrow,ncol,isub); isub = isub + 1;
        if exist('h3', 'var'), delete(h3); end
        if exist('ha1', 'var'), delete(ha1); end
        if exist('ha2', 'var'), delete(ha2); end
        
        h3 = plot(time_h,p_part_HC_list,'color',[0.8 0 0],'linewidth',2); hold on
        xlim([0 t_end/3600])
        
        xlabel('time [h]','fontsize',ftsz,'fontname',ftnm)
        ylabel(sprintf([SS_splab(spnm),' partial pressure\n {\\itp}_{', ...
            SS_splab(spnm),',H}({\\itt}) [Pa]']),'fontsize',ftsz,'fontname',ftnm)
        
        switch spnm
            case 'co2'
                ylim([0 500])
            case 'ch4'
                ylim([0 1.5e-5])
        end
        
        % Full-cell and empty-cell periods
        ha1 = area([tF_1 tF_1 tF_2 tF_2]/3600,[ylim fliplr(ylim)]); hold on
        set(ha1,'facecolor',collistF,'edgecolor',collistF,'FaceAlpha',0.2,'EdgeAlpha',0, ...
            'LineStyle','none','HandleVisibility','off')
        ha1.BaseLine.LineStyle = 'none';
        ha1.BaseValue = -50;
        uistack(ha1,'bottom');
        
        ha2 = area([tE_1 tE_1 tE_2 tE_2]/3600,[ylim fliplr(ylim)]);
        set(ha2,'facecolor',collistE,'edgecolor',collistE,'FaceAlpha',0.2,'EdgeAlpha',0, ...
            'LineStyle','none','HandleVisibility','off')
        ha2.BaseValue = -50;
        ha2.BaseLine.LineStyle = 'none';
        uistack(ha2,'bottom');
        
        set(gca,'XMinorTick','on','YMinorTick','on','box','on')
        set(gca,'fontsize',ftsz,'fontname',ftnm)
        
        xl   = xlim;
        x1   = xl(1); x2 = xl(2);
        xlab = (x2-x1)/20;
        yl   = ylim;
        y1   = yl(1); y2 = yl(2);
        ylab = y2+(y2-y1)/6;
        
        ht = text(xlab,ylab,lab_panel{isub-1});
        set(ht,'fontsize',ftsz+2,'FontWeight','bold','fontname',ftnm,'HorizontalAlignment','center')
        
        %% Time evolution of flux
        subplot(nrow,ncol,isub); isub = isub + 1;
        set(gca,'fontsize',ftsz,'fontname',ftnm)
        clear legname
        ileg = 1;
        
        if exist('h5', 'var'), delete(h5); end
        if exist('h6', 'var'), delete(h6); end
        if exist('ha3', 'var'), delete(ha3); end
        if exist('ha4', 'var'), delete(ha4); end
        
        N2p = rgas*T/V_HC; % conversion: N --> pressure
        
        h5 = plot(time_h,flux_list*N2p,'color',[0.8 0 0],'linewidth',3); hold on
        legname{ileg} = 'non-stationary flow'; ileg = ileg + 1;
        
        h6 = plot(time_h,flux_s_list*N2p,'k','linewidth',2); hold on
        legname{ileg} = 'stationary flow'; ileg = ileg + 1;
        
        legend(legname,'location','north','fontsize',ftsz,'fontname',ftnm)
        
        xlim([0 t_end/3600])
        yl = ylim;
        if (yl(1)<0), yline(0,'--k','handlevisibility','off'); end
        
        switch spnm
            case 'co2'
                if (L == 1e-3)
                    ylim([-2e-4 3e-4])
                elseif (L == 10*1e-3)
                    ylim([-5e-4 10.0001e-4])
                end
            case 'ch4'
                if (L == 1e-3)
                    ylim([3.94e-10 3.97e-10])
                elseif (L == 10*1e-3)
                    ylim([3.9e-10 4.05e-10])
                end
        end
        
        collist = jet(nt_sel);
        yl      = ylim;
        for isel = 1:nt_sel
            plot(time_sel(isel),yl(2),'Marker','v','MarkerSize',6, ...
                'MarkerFaceColor',collist(isel,:),'MarkerEdgeColor',collist(isel,:), ...
                'LineStyle','none','HandleVisibility','off');
        end
        
        % Full-cell and empty-cell periods
        yl = ylim;
        ha3 = area([tF_1 tF_1 tF_2 tF_2]/3600,[ylim fliplr(ylim)]); hold on
        set(ha3,'facecolor',collistF,'edgecolor',collistF,'FaceAlpha',0.2,'EdgeAlpha',0, ...
            'LineStyle','none','HandleVisibility','off')
        ha3.BaseLine.LineStyle = 'none';
        ha3.BaseValue = yl(1);
        uistack(ha3,'bottom');
        
        ha4 = area([tE_1 tE_1 tE_2 tE_2]/3600,[ylim fliplr(ylim)]);
        set(ha4,'facecolor',collistE,'edgecolor',collistE,'FaceAlpha',0.2,'EdgeAlpha',0, ...
            'LineStyle','none','HandleVisibility','off')
        ha4.BaseValue = yl(1);
        ha4.BaseLine.LineStyle = 'none';
        uistack(ha4,'bottom');
        
        xlabel('time {\itt} [h]','fontsize',ftsz,'fontname',ftnm);
        ylabel(sprintf('pressure change\n {\\itdp}({\\itt})/{\\itdt} [Pa s^{-1}]'), ...
            'fontsize',ftsz,'fontname',ftnm);
        
        set(gca,'XMinorTick','on','YMinorTick','on','box','on')
        set(gca,'fontsize',ftsz,'fontname',ftnm)
        
        xl   = xlim;
        x1   = xl(1); x2 = xl(2);
        xlab = (x2-x1)/20;
        yl   = ylim;
        y1   = yl(1); y2 = yl(2);
        ylab = y2+(y2-y1)/6;
        
        ht = text(xlab,ylab,lab_panel{isub-1});
        set(ht,'fontsize',ftsz+2,'FontWeight','bold','fontname',ftnm,'HorizontalAlignment','center')
        
        %% Save figure
        if (savefig == 1)
            sfpath     = '../Figures/';
            mkdir(sfpath)
            sffilename = ['Figure_',lab_figs{jlab},'.png'];
            sfflnm     = fullfile(sfpath,sffilename);
            exportgraphics(gcf,sfflnm,'Resolution',400)
        end
    end
end

toc