% =================
% Figures S13 - S30
% =================
% 
% This script plots the pressure data in FO chamber and Herriott cell
% during the full- and empty-cell runs.
% 
% Linear regression and correlation coefficient ares calculated
%
% Author: sebastien.viscardy@aeronomie.be
% 
%% Link to functions
addpath('../functions/');

%%
clearvars

tic

marscst

%% Type of experiment
t_exp_list = {'D' 'E'}; n_exp = length(t_exp_list);
D_sol_list = [79 81 106 292 306 313 466 474 504 526 684];
E_sol_list = [573 684 2442 2446 2615 2627 2644];

%% Infos on figure
figid     = 421;
sfigtype  = 'png';
ftsz      = 12;
ftnm      = 'times';
savefig   = 1; % = 1: save Figure
savedata  = 0; % = 1: save Data (Table)
collistFF = [0.8 0 0.8];
collistFH = [1 0 0];
collistE  = [0.2 0.8 0.6];
alpha_lev = 0.2;
mksz      = 6;
xsc       = 'lin';
ysc       = 'lin';
lab_panel = {'(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)'};
leg_loc   = {'east' 'east' 'southeast' 'southeast'};
jexp      = 0;

%% Loop over types of experiments
for iexp = 1:n_exp
    t_exp = t_exp_list{iexp};
    switch t_exp
        case 'D', sol_list = D_sol_list; n_enr = 1;
        case 'E', sol_list = E_sol_list; n_enr = 25;
    end
    nsol = length(sol_list);
    
    %% Loop over run sols
    for isol = 1:nsol
        sol_index   = sol_list(isol);
        jexp = jexp + 1;
        
        %% Load full data
        if ( sol_index <  2442 )
            SS_MSL_full_data_Webster_2015
        else
            SS_MSL_full_data_Webster_2021
        end
        
        % Mean Wefg: Empty-cell and Full-cell runs
        E_mean_Wefg_CH4 = mean(E_Wefg_CH4);
        F_mean_Wefg_CH4 = mean(F_Wefg_CH4);
        
        % number of runs: Empty-cell and Full-cell runs
        nE              = length(E_Wefg_CH4);
        nF              = length(F_Wefg_CH4);
        
        % std Wefg: Empty-cell and Full-cell runs
        E_sig_Wefg_CH4  = sqrt( 1/(nE-1)*sum( (E_Wefg_CH4-E_mean_Wefg_CH4).^2 ) );
        F_sig_Wefg_CH4  = sqrt( 1/(nF-1)*sum((F_Wefg_CH4-F_mean_Wefg_CH4).^2) );
        
        % CH4 vmr and error in the Herriott cell
        eta_H           = F_mean_Wefg_CH4 - E_mean_Wefg_CH4;
        sig_H           = sqrt(1/nF * F_sig_Wefg_CH4^2 + 1/nE * E_sig_Wefg_CH4^2);
        
        eta             = eta_H/n_enr;
        sig             = sig_H/n_enr;
        
        data_eta_sig(jexp).eta   = eta;
        data_eta_sig(jexp).sigma = sig;
        
        disp(['Sol ',num2str(sol_index),' : ',num2str(eta,'%2.2f'),' +/- ', ...
            num2str(sig,'%2.2f'),' ppbv'])
        
        %% Make figure
        figure(figid)
        delete(get(gcf,'children'))
        set(gcf,'paperpositionmode','auto')
        set(gcf,'position',[200 200 1200 500])
        set(gca,'fontsize',ftsz,'fontname',ftnm)
        nrow = 1; ncol = 2; isub = 1;
        
        %% 1) FO pressure: Full-cell runs
        hs1 = subplot(nrow,ncol,isub); isub = isub + 1;
        clear legname
        legname = cell(1,2);
        ileg = 1;
        
        %% 1.a) FO pressure: Full-cell runs
        prs_shift  = min(min(F_FO_prs),0);
        F_FO_prs   = F_FO_prs - prs_shift;
        E_FO_prs   = E_FO_prs - prs_shift;
        x          = F_elapsed_t/3600; % elapsed time [h]
        y          = F_FO_prs;         % FO pressure  [Pa]
        F_R_FO     = corrcoef(x,y);    % Correlation coefficient
        
        % Linear regression
        X          = [ones(length(x),1) x];
        b          = X\y;                   % regression coefficients
        F_y_lr     = b(2)*x+b(1);           % linear regression   [Pa]
        F_delta    = F_y_lr(end)-F_y_lr(1); % pressure difference [Pa]
        
        % Store data
        data_FO_F(jexp).t_exp  = t_exp;
        data_FO_F(jexp).sol    = sol_index;
        data_FO_F(jexp).x      = x;
        data_FO_F(jexp).b      = b;
        data_FO_F(jexp).F_R_FO = F_R_FO(2,1);
        
        % plot TLS data
        hp1 = plot(x,y); hold on
        set(hp1,'LineStyle','none','Marker','o','MarkerFaceColor',collistFF, ...
            'MarkerSize',mksz,'MarkerEdgeColor',collistFF,'HandleVisibility','off')
        
        % plot linear fit
        plot(x,F_y_lr,'color',collistFF/1.5,'linewidth',4); hold on
        b2 = b(2);
        if ( b(1) >= 0 )
            b1 = ['+ ',num2str(b(1),'%2.1f')];
        elseif (b(1) < 0 )
            bcst = num2str(abs(b(1)),'%2.1f');
            b1 = ['- ',bcst(1:end)];
        end
        legname{ileg} = ['F: {\itp}_{FO}({\itt}) = ', ...
            num2str(b2,'%2.1f'),' {\itt} ',b1]; ileg = ileg + 1;
        dp_text_full  = ['\Delta{\itp}_{FO} = ',num2str(F_delta,'%2.1f'),' Pa'];
        
        %% 1.b) FO pressure: Empty-cell runs
        % 1) FO
        x       = E_elapsed_t/3600; % elapsed time [h]
        y       = E_FO_prs;         % FO pressure  [Pa]
        E_R_FO  = corrcoef(x,y);    % Correlation coefficient
        
        % Linear regression
        X       = [ones(length(x),1) x];
        b       = X\y;
        E_y_lr  = b(2)*x+b(1);
        E_delta = E_y_lr(end)-E_y_lr(1);
        
        % Store data
        data_FO_E(jexp).t_exp  = t_exp;
        data_FO_E(jexp).sol    = sol_index;
        data_FO_E(jexp).x      = x;
        data_FO_E(jexp).b      = b;
        data_FO_E(jexp).E_R_FO = E_R_FO(2,1);
        
        % plot TLS data
        hp1 = plot(x,y); hold on
        set(hp1,'LineStyle','none','Marker','o','MarkerFaceColor',collistE, ...
            'MarkerSize',mksz,'MarkerEdgeColor',collistE,'HandleVisibility','off')
        
        % plot linear fit
        plot(x,E_y_lr,'color',collistE/1.5,'linewidth',4); hold on
        b2 = b(2);
        if ( b(1) >= 0 )
            b1 = ['+ ',num2str(b(1),'%2.1f')];
        elseif (b(1) < 0 )
            bcst = num2str(abs(b(1)),'%2.1f');
            b1 = ['- ',bcst(1:end)];
        end
        legname{ileg} = ['E: {\itp}_{FO}({\itt}) = ',num2str(b2,'%2.1f'), ...
            ' {\itt} ',b1]; ileg = ileg + 1;
        dp_text_empty = ['\Delta{\itp}_{FO} = ',num2str(E_delta,'%2.1f'),' Pa'];
        
        yl     = ylim;
        ylmin  = floor(mean(yl)/10)*10 - 50;
        ylmax  = ylmin + 100;
        ylim([ylmin ylmax])
        xl     = xlim; x1 = xl(1); x2 = xl(2);
        yl     = ylim; y1 = yl(1); y2 = yl(2);
        xl_sol = xlim;
        
        ht1 = text(x1+(x2-x1)/4,y2-(y2-y1)/20,dp_text_full,'Color',collistFF, ...
            'FontWeight','bold','fontsize',ftsz,'fontname',ftnm, ...
            'HorizontalAlignment','center');
        
        ht2 = text(x1+(x2-x1)*3/4,y2-(y2-y1)/20,dp_text_empty,'Color',collistE, ...
            'FontWeight','bold','fontsize',ftsz,'fontname',ftnm, ...
            'HorizontalAlignment','center');
        
        % Legend
        hl1 = legend(legname,'location','south','fontsize',ftsz-2,'fontname',ftnm,'EdgeColor','w');
        
        % Axes, labels
        set(gca,'XMinorTick','on','YMinorTick','on')
        
        xlabel('elapsed time {\itt} [h]','fontsize',ftsz,'fontname',ftnm)
        ylabel('pressure {\itp}_{FO} [Pa]','fontsize',ftsz,'fontname',ftnm)
        
        % Label of panel
        xlim(xl_sol)
        xl = xlim;
        yl = ylim;
        SS_mk_label_panel(lab_panel{isub-1},xl,xsc,yl,ysc,ftsz,ftnm);
        
        set(gca,'XMinorTick','on','YMinorTick','on','YTick',yl(1):10:yl(2))
        set(gca,'fontsize',ftsz,'fontname',ftnm)
        
        texttitle = ['Sol ',num2str(sol_index),': ', ...
            num2str(eta,'%2.2f'),' \pm ',num2str(sig,'%2.2f'),' ppbv'];
        title(texttitle,'fontsize',ftsz-1,'fontname',ftnm)
        
        %% 2.a) HCell pressure: full-cell runs
        hs2 = subplot(nrow,ncol,isub); isub = isub + 1;
        clear legname
        
        x         = F_elapsed_t/3600; % elapsed time   [h]
        y         = F_HC_prs;         % HCell pressure [Pa]
        F_R_HCell = corrcoef(x,y);    % Correlation coefficient
        
        % Linear regression
        X         = [ones(length(x),1) x];
        b         = X\y;                   % regression coefficients
        F_y_lr    = b(2)*x+b(1);           % linear regression       [Pa]
        F_delta   = F_y_lr(end)-F_y_lr(1); % pressure difference     [Pa]
        
        % Store data
        data_H_F(jexp).t_exp  = t_exp;
        data_H_F(jexp).sol    = sol_index;
        data_H_F(jexp).x      = x;
        data_H_F(jexp).b      = b;
        data_H_F(jexp).E_R_H  = F_R_HCell(2,1);
        
        % plot TLS data
        hp1 = plot(x,y); hold on
        set(hp1,'LineStyle','none','Marker','v','MarkerFaceColor',collistFF, ...
            'MarkerSize',mksz,'MarkerEdgeColor',collistFF,'HandleVisibility','off')
        
        % plot linear fit
        plot(x,F_y_lr,'color',collistFF/1.5,'linewidth',4); hold on
        b2 = b(2);
        if ( b(1) >= 0 )
            b1 = ['+ ',num2str(b(1),'%2.1f')];
        elseif (b(1) < 0 )
            bcst = num2str(abs(b(1)),'%2.1f');
            b1 = ['- ',bcst(1:end)];
        end
        legname = ['F: {\itp}_{H}({\itt}) = ',num2str(b2,'%2.1f'), ...
            ' {\itt} ',b1];
        
        dp_text = ['\Delta{\itp}_{H} = ',num2str(F_delta,'%2.1f'),' Pa'];
        
        xlim(xl_sol)
        yl    = ylim;
        ylmin = floor(mean(yl)/10)*10 - 50;
        ylmax = ylmin + 100;
        ylim([ylmin ylmax])
        xl    = xlim; x1 = xl(1); x2 = xl(2);
        yl    = ylim; y1 = yl(1); y2 = yl(2);
        
        ht3 = text(x1+(x2-x1)/4,y2-(y2-y1)/20,dp_text,'Color',collistFF, ...
            'FontWeight','bold','fontsize',ftsz,'fontname',ftnm, ...
            'HorizontalAlignment','center');
        
        hl2 = legend(legname,'position',[0.65719 0.17 0.17125 0.047143], ...
            'fontsize',ftsz-2,'fontname',ftnm,'EdgeColor','w');
        
        % Axes, labels
        yl = ylim;
        xlim(xl)
        
        set(gca,'XMinorTick','on','YMinorTick','on','YTick',yl(1):10:yl(2))
        
        if (isub >=nrow*ncol)
            xlabel('elapsed time {\itt} [h]','fontsize',ftsz,'fontname',ftnm)
        end
        ylabel('pressure {\itp}_H during full-cell runs [Pa]','fontsize',ftsz,'fontname',ftnm)
        
        % Label of panel
        xl = xlim; yl = ylim;
        SS_mk_label_panel(lab_panel{isub-1},xl,xsc,yl,ysc,ftsz,ftnm);
        
        set(gca,'fontsize',ftsz,'fontname',ftnm)
        title(texttitle,'fontsize',ftsz-1,'fontname',ftnm)
        
        %% 2.b) HCell pressure (empty-cell runs)
        if (sol_index >= 2442)
            set(gca,'YColor',collistFF)
            
            % Second axes
            ax1 = gca;
            ax2 = axes('Position',get(ax1,'Position'),...
                'XAxisLocation','top',...
                'YAxisLocation','right',...
                'Color','none',...
                'XColor','k','YColor',collistE,'yMinorTick','off','fontsize',ftsz);
            set(ax1,'box','off','xminortick','on');
            set(ax2,'box','off','xminortick','on');
            set(ax2,'fontsize',ftsz,'fontname',ftnm)
            set(ax2,'XTickLabel',[])
            xlim(ax2,xl_sol)
            
            clear legname
            
            x         = E_elapsed_t/3600; % elapsed time   [h]
            y         = E_HC_prs;         % HCell pressure [Pa]
            E_R_HCell = corrcoef(x,y);    % Correlation coefficient
            
            % Linear regression
            X         = [ones(length(x),1) x];
            b         = X\y;                   % regression coefficients
            E_y_lr    = b(2)*x+b(1);           % linear regression       [Pa]
            E_delta   = E_y_lr(end)-E_y_lr(1); % pressure difference     [Pa]
            
            % Store data
            data_H_E(jexp).t_exp  = t_exp;
            data_H_E(jexp).sol    = sol_index;
            data_H_E(jexp).x      = x;
            data_H_E(jexp).b      = b;
            data_H_E(jexp).E_R_H  = E_R_HCell(2,1);
            
            % plot TLS data
            hp1 = line(ax2,x,y); hold on
            set(hp1,'LineStyle','none','Marker','v','MarkerFaceColor',collistE, ...
                'MarkerSize',mksz,'MarkerEdgeColor',collistE,'HandleVisibility','off')
            
            % plot linear fit
            line(ax2,x,E_y_lr,'color',collistE/1.5,'linewidth',4); hold on
            b2 = b(2);
            if ( b(1) >= 0 )
                b1 = ['+ ',num2str(b(1),'%2.1f')];
            elseif (b(1) < 0 )
                bcst = num2str(abs(b(1)),'%2.1f');
                b1 = ['- ',bcst(1:end)];
            end
            legname = ['E: {\itp}_{H}({\itt}) = ',num2str(b2,'%2.1f'), ...
                ' {\itt} ',b1];
            
            dp_text = ['\Delta{\itp}_H = ',num2str(E_delta,'%2.1f'),' Pa'];
            
            xlim(xl_sol)
            ylim([-50 50])
            xl    = xlim; x1 = xl(1); x2 = xl(2);
            yl    = ylim; y1 = yl(1); y2 = yl(2);
            
            ht3 = text(x1+(x2-x1)*3/4,y2-(y2-y1)/20,dp_text,'Color',collistE, ...
                'FontWeight','bold','fontsize',ftsz,'fontname',ftnm, ...
                'HorizontalAlignment','center');
            
            hl3 = legend(legname,'position',[0.65719 0.12524 0.17125 0.047143], ...
                'fontsize',ftsz-2,'fontname',ftnm,'color','w','EdgeColor','w');
            
            % Axes, labels
            yl = ylim;
            xlim(xl)
            set(ax2,'XMinorTick','on','YMinorTick','on','YTick',yl(1):10:yl(2))
            
            ylabel('pressure {\itp}_H during empty-cell runs [Pa]','fontsize',ftsz,'fontname',ftnm)
        end
        
        %% Save figure
        if (savefig == 1)
            sfpath     = '../Figures/';
            sffilename = ['Figure_S',num2str(13+jexp-1,'%02i'),'.png'];
            sfflnm     = fullfile(sfpath,sffilename);
            exportgraphics(gcf,sfflnm,'Resolution',400)
        else
            waitforbuttonpress
        end
    end
end

toc