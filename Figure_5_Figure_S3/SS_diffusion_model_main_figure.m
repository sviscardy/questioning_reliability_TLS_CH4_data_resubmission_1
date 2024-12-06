tic

%% Mean cell pressure (empty-cell runs) / (E) Sol 2627
ref_prs_sol     = 2627;
ref_E_p_H       = 2.6;  % [Pa]
ref_E_delta_p_H = 12.6; % [Pa]

%% Loop over randomly generated data points (model)
njsim = 20;
for jsim = 13 % 1:njsim
    disp([num2str(jsim),' / ',num2str(njsim)])
    
    %% Make figure
    figure(figid)
    close(figid)
    figure(figid)
    delete(get(gcf,'children'))
    set(gcf,'paperpositionmode','auto')
    
    set(gcf,'position',[200 0 1200 900])
    set(gca,'fontsize',ftsz,'fontname',ftnm)
    nrow = 3; ncol = 2; isub = 1;
    
    %% 1) p_HC_t vs. time
    hs1 = subplot(nrow,ncol,isub); isub = isub + 1;
    legname = cell(1,3);
    ileg = 1;
    
    plot(F_elapsed_t/3600,F_HC_prs,'marker','o','MarkerSize',4,'MarkerFaceColor',collistF, ...
        'MarkerEdgeColor',collistF,'LineStyle','none'); hold on
    legname{ileg} = ['$$\left \lbrace p_{\rm H,F}(k) ', ...
        '\right \rbrace_{k=1}^{N_{\rm F}} \; {\rm (TLS)}$$']; ileg = ileg + 1;
    
    if (sol_index >= 2442)
        plot(E_elapsed_t/3600,E_HC_prs,'marker','o','MarkerSize',4,'MarkerFaceColor',collistE, ...
            'MarkerEdgeColor',collistE,'LineStyle','none'); hold on
        legname{ileg} = ['$$\left \lbrace p_{\rm H,E}(k)', ...
            '\right \rbrace_{k=1}^{N_{\rm F}} \; {\rm (TLS)}$$']; ileg = ileg + 1;
    end
    
    plot(time_h,p_HC_t,'color','k','linewidth',2); hold on
    legname{ileg} = '$$p_{\rm H}(t) \; {\rm (model)}$$'; ileg = ileg + 1;
    
    if (sol_index < 2442)
        errorbar(mean(E_elapsed_t)/3600,ref_E_p_H,ref_E_delta_p_H,'marker','o','MarkerSize',4, ...
            'MarkerFaceColor',collistE,'MarkerEdgeColor',collistE, ...
            'color',collistE,'HandleVisibility','off')
    end
    
    errorbar(t_0_s/3600,ref_E_p_H,ref_E_delta_p_H,'marker','o', ...
        'MarkerSize',4,'MarkerFaceColor',collistE, ...
        'MarkerEdgeColor',collistE,'color',collistE,'LineStyle','none')
    legname{ileg} = ['$${\overline p}_{\rm H,E} = ', ...
        num2str(ref_E_p_H),' \pm ',num2str(ref_E_delta_p_H), ...
        ' \; {\rm Pa} \; {\rm (TLS, Sol \; ',num2str(ref_prs_sol),')}$$'];
    
    xl = xlim; xlim([0 xl(2)])
    
    % Full-cell and empty-cell periods
    yl = ylim;
    ylim([-50 yl(2)*1.1])
    ha1 = area([tF_1 tF_1 tF_2 tF_2]/3600,[ylim fliplr(ylim)]);
    set(ha1,'facecolor',collistF,'edgecolor',collistF,'FaceAlpha',0.2,'EdgeAlpha',0, ...
        'LineStyle','none','HandleVisibility','off')
    ha1.BaseLine.LineStyle = 'none';
    ha1.BaseValue = -50;
    uistack(ha1,'bottom');
    
    ha2 = area([tE_1 tE_1 tE_2 tE_2]/3600,[ylim fliplr(ylim)]);
    set(ha2,'facecolor',collistE,'edgecolor',collistE,'FaceAlpha',0.2,'EdgeAlpha',0, ...
        'LineStyle','none','HandleVisibility','off')
    ha1.BaseValue = -50;
    ha2.BaseLine.LineStyle = 'none';
    uistack(ha2,'bottom');
    
    % Legend
    hl1 = legend(legname,'location','northwest', ...
        'interpreter','latex','fontsize',ftsz,'fontname',ftnm);
    set(hl1,'position',[0.22854      0.90231       0.1128     0.080611])
    
    % Axes and labels
    xlabel('elapsed time {\itt} [h]','fontsize',ftsz,'fontname',ftnm)
    ylabel(sprintf('Herriott cell pressure\n{\\itp}_H({\\itt}) [Pa]'),'fontsize',ftsz,'fontname',ftnm)
    
    SS_mk_label_panel(lab_panel{isub-1},xlim,'lin',ylim,'lin',ftsz,ftnm)
    
    set(gca,'XMinorTick','on','YMinorTick','on','box','on')
    set(gca,'fontsize',ftsz,'fontname',ftnm)
    
    %% 2) N_CH4_HC_t vs. time
    hs2 = subplot(nrow,ncol,isub); isub = isub + 1;
    clear legname
    legname = cell(1,2);
    ileg = 1;
    
    plot(time_h,N_CH4_HC_t*1e9,'color','k','linewidth',2); hold on
    legname{ileg} = '$$N_{\rm CH_4}(t) \; {\rm (model)}$$'; ileg = ileg + 1;
    
    errorbar(mean(F_elapsed_t)/3600,N_HC*1e9,delta_N_HC*1e9,'Marker','o','MarkerFaceColor', ...
        [1 0 0],'MarkerEdgeColor',[1 0 0],'color',[1 0 0],'LineStyle','none'); hold on
    legname{ileg} = '$$N_{\rm H} \pm \delta N_{\rm H} \; {\rm (TLS)}$$';
    
    % Full-cell and empty cell periods
    xl = xlim; xlim([0 xl(2)])
    yl = ylim; ylim([0 yl(2)*1.1])
    
    ha1 = area([tF_1 tF_1 tF_2 tF_2]/3600,[ylim fliplr(ylim)]);
    set(ha1,'facecolor',collistF,'edgecolor',collistF,'FaceAlpha',0.2, ...
        'EdgeAlpha',0,'HandleVisibility','off')
    uistack(ha1,'bottom');
    
    ha2 = area([tE_1 tE_1 tE_2 tE_2]/3600,[ylim fliplr(ylim)]);
    set(ha2,'facecolor',collistE,'edgecolor',collistE,'FaceAlpha',0.2, ...
        'EdgeAlpha',0,'HandleVisibility','off')
    uistack(ha2,'bottom');
    
    % Legend
    hl2 = legend(legname,'interpreter','latex','fontsize',ftsz,'fontname',ftnm);
    set(hl2,'position',[0.68476      0.91231      0.10607     0.053625])
    
    % Axes and labels
    xlabel('elapsed time {\itt} [h]','fontsize',ftsz,'fontname',ftnm)
    ylabel(sprintf('CH_4 amount in the Herriott cell\n{\\itN}_{CH_4}({\\itt}) [nmol]'), ...
        'fontsize',ftsz,'fontname',ftnm)
    
    SS_mk_label_panel(lab_panel{isub-1},xlim,'lin',ylim,'lin',ftsz,ftnm)
    
    set(gca,'XMinorTick','on','YMinorTick','on','box','on')
    set(gca,'fontsize',ftsz,'fontname',ftnm)
    
    %% 3) eta_HC_t vs. time
    hs3 = subplot(nrow,ncol,isub); isub = isub + 1;
    clear legname
    legname = cell(1,2);
    ileg = 1;
    
    kk = time > t_0_s;
    
    plot(time_h(kk),eta_HC_t(kk),'color','k','linewidth',2); hold on
    legname{ileg} = '$$\eta_{\rm CH_4}(t) \; {\rm (model)}$$'; ileg = ileg + 1;
    
    errorbar(mean(F_elapsed_t)/3600,eta_HC*1e9,delta_eta_HC*1e9,'Marker','o','MarkerFaceColor', ...
        [1 0 0],'MarkerEdgeColor',[1 0 0],'color',[1 0 0],'LineStyle','none'); hold on
    legname{ileg} = '$$\eta_{\rm H} \pm \sigma_{\rm H} \; {\rm (TLS)}$$';
    
    set(gca,'YScale','log')
    
    xl = xlim; xlim([0 xl(2)])
    
    ha1 = area([tF_1 tF_1 tF_2 tF_2]/3600,[ylim fliplr(ylim)]);
    set(ha1,'facecolor',collistF,'edgecolor',collistF,'FaceAlpha',0.2, ...
        'EdgeAlpha',0,'HandleVisibility','off')
    uistack(ha1,'bottom');
    
    ha2 = area([tE_1 tE_1 tE_2 tE_2]/3600,[ylim fliplr(ylim)]);
    set(ha2,'facecolor',collistE,'edgecolor',collistE,'FaceAlpha',0.2, ...
        'EdgeAlpha',0,'HandleVisibility','off')
    uistack(ha2,'bottom');
    
    % Legend
    hl3 = legend(legname,'location','north', ...
        'interpreter','latex','fontsize',ftsz,'fontname',ftnm);
    
    % Axes and labels
    xlabel('elapsed {\itt} time [h]','fontsize',ftsz,'fontname',ftnm)
    ylabel(sprintf('CH_4 vmr in the Herriott cell\n\\eta_{CH_4}({\\itt}) [ppbv]'),'fontsize',ftsz,'fontname',ftnm)
    
    SS_mk_label_panel(lab_panel{isub-1},xlim,'lin',ylim,'log',ftsz,ftnm)
    
    set(gca,'XMinorTick','on','YMinorTick','on','box','on')
    set(gca,'fontsize',ftsz,'fontname',ftnm)
    
    %% 4) Statistics
    hs5 = subplot(nrow,ncol,isub); isub = isub + 1;
    clear legname
    legname = cell(1,4);
    ileg = 1;
    
    dxbin = 0.2;
    
    hhE = histogram(a_s_E_list*3600,'BinWidth',dxbin); hold on
    set(hhE,'FaceColor',collistE,'FaceAlpha',0.4,'EdgeAlpha',0)
    legname{ileg} = 'empty-cell (model)'; ileg = ileg + 1;
    
    xline(a1T_E*3600,'color',collistE,'linewidth',2); hold on
    legname{ileg} = 'empty-cell (TLS)'; ileg = ileg + 1;
    
    hhF = histogram(a_s_F_list*3600,'BinWidth',dxbin); hold on
    set(hhF,'FaceColor',collistF,'FaceAlpha',0.4,'EdgeAlpha',0)
    legname{ileg} = 'full-cell (model)'; ileg = ileg + 1;
    
    xline(a1T_F*3600,'color',collistF,'linewidth',2); hold on
    legname{ileg} = 'full-cell (TLS)';
    
    legend(legname,'location','northwest','fontsize',ftsz-2,'fontname',ftnm)
    
    % Axes and labels
    xlabel('trend {\ita} [ppbv h^{-1}]','fontsize',ftsz,'fontname',ftnm)
    ylabel('count ','fontsize',ftsz,'fontname',ftnm)
    
    xlim([-40 40])
        
    set(gca,'XTick',-100:10:100,'XMinorTick','on','YMinorTick','on','box','on')
    set(gca,'fontsize',ftsz,'fontname',ftnm)
    
    SS_mk_label_panel(lab_panel{isub-1},xlim,'lin',ylim,'lin',ftsz,ftnm)
    
    %% 5) Simulated eta_mod, eta_E, and eta_F
    hs4 = subplot(nrow,ncol,[5 6]); isub = isub + 1;
    legname = cell(1,1);
    ileg    = 1;
    
    % Select data
    a_s_F      = a_s_F_list(jsim);
    a_s_E      = a_s_E_list(jsim);
    y_lins_F   = y_lins_F_list(jsim,:);
    y_lins_E   = y_lins_E_list(jsim,:);
    mod_eta_HC = mod_eta_HC_list(jsim);
    mod_sig_HC = mod_sig_HC_list(jsim);
    eta_F_mod  = eta_F_mod_list(jsim,:);
    eta_E_mod  = eta_E_mod_list(jsim,:);
    
    % full-cell runs: TLS
    plot(F_elapsed_t/3600,F_Wefg_CH4,'Marker','o','MarkerSize',6, ...
        'MarkerFaceColor',collistF,'MarkerEdgeColor',collistF,'color',collistF); hold on
    legname{ileg} = ['$$\left \lbrace \eta_{\rm F}^\prime({\it k})', ...
        '\right \rbrace_{{\it k}=1}^{{\it N}_{\rm F}} ({\rm TLS})$$']; ileg = ileg + 1;
    
    % full-cell runs: simulated
    colsimF = [1 1 1]*0.4;
    plot(F_elapsed_t/3600,eta_F_mod,'Marker','^','MarkerSize',6, ...
        'MarkerFaceColor',colsimF,'MarkerEdgeColor',colsimF,'color',colsimF); hold on
    legname{ileg} = ['$$\left \lbrace \eta_{\rm F,m}^\prime({\it k})', ...
        '\right \rbrace_{{\it k}=1}^{{\it N}_{\rm F}} ({\rm model})$$']; ileg = ileg + 1;
    
    % empty-cell runs: TLS
    plot(E_elapsed_t/3600,E_Wefg_CH4,'Marker','o','MarkerSize',6, ...
        'MarkerFaceColor',collistE,'MarkerEdgeColor',collistE,'color',collistE); hold on
    legname{ileg} = ['$$\left \lbrace \eta_{\rm E}^\prime({\it k})', ...
        '\right \rbrace_{{\it k}=1}^{{\it N}_{\rm E}} ({\rm TLS})$$']; ileg = ileg + 1;
    
    % empty-cell runs: simulated
    colsimE = [1 1 1]*0.7;
    plot(E_elapsed_t/3600,eta_E_mod,'Marker','^','MarkerSize',6, ...
        'MarkerFaceColor',colsimE,'MarkerEdgeColor',colsimE,'color',colsimE); hold on
    legname{ileg} = ['$$\left \lbrace \eta_{\rm E,m}^\prime({\it k})', ...
        '\right \rbrace_{{\it k}=1}^{{\it N}_{\rm E}} ({\rm model})$$']; ileg = ileg + 1;
    
    % simulated \eta_s(t)
    plot(time_h,eta_mod,'color','k','linewidth',2); hold on
    legname{ileg} = '$$\eta_{\rm m}(t) = {\overline \eta}_{\rm E} \left ( 1+q(t) \right )$$'; ileg = ileg + 1;
    
    yline(0,'--k','handlevisibility','off')
    
    xl = xlim;
    yl = ylim;
    xlim([0 xl(2)])
    ylim([min(0,yl(1)) yl(2)])
    
    ha1 = area([tF_1 tF_1 tF_2 tF_2]/3600,[ylim fliplr(ylim)]);
    set(ha1,'facecolor',collistF,'edgecolor',collistF,'FaceAlpha',0.2, ...
        'EdgeAlpha',0,'HandleVisibility','off')
    uistack(ha1,'bottom');
    ha1.BaseValue = -50;
    ha1.BaseLine.LineStyle = 'none';
    
    ha2 = area([tE_1 tE_1 tE_2 tE_2]/3600,[ylim fliplr(ylim)]);
    set(ha2,'facecolor',collistE,'edgecolor',collistE,'FaceAlpha',0.2, ...
        'EdgeAlpha',0,'HandleVisibility','off')
    uistack(ha2,'bottom');
    ha2.BaseValue = -50;
    ha2.BaseLine.LineStyle = 'none';
    
    % eta_H +/- sigma_H: TLS vs. model
    xl = xlim;
    yl = ylim;
    eta_H_TLS_vs_sim = sprintf([num2str(eta_HC*1e9,'%2.2f'), ...
        ' \\pm ',num2str(delta_eta_HC*1e9,'%2.2f'),' ppbv (TLS)','\n', ...
        num2str(mod_eta_HC,'%2.2f'),' \\pm ',num2str(mod_sig_HC,'%2.2f'),' ppbv (model)']);
    ht = text(xl(2)/50,yl(1)+(yl(2)-yl(1))/8, ...
        eta_H_TLS_vs_sim,'fontsize',ftsz,'fontname',ftnm, ...
        'BackgroundColor','white','HorizontalAlignment','left');
    
    % Legend
    hl = legend(legname,'location','northeast', ...
        'interpreter','latex','fontsize',ftsz,'fontname',ftnm);
    set(hl,'position',[0.75352 0.11753 0.15541 0.20125])
    
    % Linear regressions
    plot(xl_F/3600,y_linT_F,'color',collistF,'linewidth',2,'HandleVisibility','off'); hold on
    plot(xl_F/3600,y_lins_F,'color',colsimF,'linewidth',2,'HandleVisibility','off'); hold on
    plot(xl_E/3600,y_linT_E,'color',collistE,'linewidth',2,'HandleVisibility','off'); hold on
    plot(xl_E/3600,y_lins_E,'color',colsimE,'linewidth',2,'HandleVisibility','off'); hold on
    
    xl = xlim;
    yl = ylim;
    
    text1 = ['$$a \; = \; {\rm ',num2str(a1T_F*3600,'%2.2f'),'\; ppbv \; h^{-1}}$$'];
    httr1 = text(xl(2)/50,yl(1)+(yl(2)-yl(1))*0.92,text1,'color',collistF, ...
        'fontsize',ftsz+1,'fontname',ftnm,'FontWeight','bold','BackgroundColor','white', ...
        'interpreter','latex'); hold on
    
    text2 = ['$$a \; = \; {\rm ',num2str(a_s_F*3600,'%2.2f'),'\; ppbv \; h^{-1}}$$'];
    httr2 = text(xl(2)/50,yl(1)+(yl(2)-yl(1))*0.8,text2,'color',colsimF, ...
        'fontsize',ftsz+1,'fontname',ftnm,'FontWeight','bold','BackgroundColor','white', ...
        'interpreter','latex'); hold on
    
    text3 = ['$$a \; = \; {\rm ',num2str(a1T_E*3600,'%2.2f'),'\; ppbv \; h^{-1}}$$'];
    httr3 = text(xl(2)/50,yl(1)+(yl(2)-yl(1))*0.68,text3,'color',collistE, ...
        'fontsize',ftsz+1,'fontname',ftnm,'FontWeight','bold','BackgroundColor','white', ...
        'interpreter','latex'); hold on
    
    text4 = ['$$a \; = \; {\rm ',num2str(a_s_E*3600,'%2.2f'),'\; ppbv \; h^{-1}}$$'];
    httr4 = text(xl(2)/50,yl(1)+(yl(2)-yl(1))*0.56,text4,'color',colsimE, ...
        'fontsize',ftsz+1,'fontname',ftnm,'FontWeight','bold','BackgroundColor','white', ...
        'interpreter','latex'); hold on
    
    % Axes and labels
    xlabel('elapsed time {\itt} [h]','fontsize',ftsz,'fontname',ftnm)
    ylabel('CH_4 vmr \eta_m ({\itt}) [ppbv]','fontsize',ftsz,'fontname',ftnm)
    
    SS_mk_label_panel(lab_panel{isub-1},xlim,'lin',ylim,'lin',ftsz,ftnm)
    
    set(gca,'XMinorTick','on','YMinorTick','on','box','on')
    set(gca,'fontsize',ftsz,'fontname',ftnm)
    
    set(gca,'position',[0.105 0.109 0.64333 0.21574])
    
    %% Save figure
    if (savefig == 1)
        sfpath     = '../Figures/';
        sffilename = 'Figure_05.png';
        sfflnm     = fullfile(sfpath,sffilename);
        exportgraphics(gcf,sfflnm,'Resolution',400)
    else
        waitforbuttonpress
    end
end

toc