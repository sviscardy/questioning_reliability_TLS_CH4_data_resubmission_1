for jsim = 1
    
    %% Make figure
    figure(figid)
    delete(get(gcf,'children'))
    set(gcf,'paperpositionmode','auto')
    set(gcf,'position',[200 0 1600 900])
    set(gca,'fontsize',ftsz,'fontname',ftnm)
    nrow = 2; ncol = 2; isub = 1;
    
    %% 1) p_HC_t vs. time
    hs1 = subplot(nrow,ncol,isub); isub = isub + 1;
    
    t_end_1 = 4.5;
    t_ini_2 = 24;
    
    t_rem   = t_ini_2 - t_end_1 - 0.5;
    
    kx1 = time_h <= t_end_1;
    plot(time_h(kx1),p_HC_t(kx1),'color','k','linewidth',2); hold on
    
    kx2 = time_h >= t_ini_2;
    plot(time_h(kx2)-t_rem,p_HC_t(kx2),'color','k','linewidth',2); hold on
    
    xl = xlim; xlim([0 xl(2)])
    
    % Full-cell and empty-cell periods
    ylim([-50 500])
    ha1 = area([tF1_1 tF1_1 tF1_2 tF1_2]/3600,[ylim fliplr(ylim)]);
    set(ha1,'facecolor',collistF,'edgecolor',collistF,'FaceAlpha',0.2,'EdgeAlpha',0, ...
        'LineStyle','none','HandleVisibility','off')
    ha1.BaseLine.LineStyle = 'none';
    ha1.BaseValue = -50;
    uistack(ha1,'bottom');
    
    ha2 = area([tE1_1 tE1_1 tE1_2 tE1_2]/3600,[ylim fliplr(ylim)]);
    set(ha2,'facecolor',collistE,'edgecolor',collistE,'FaceAlpha',0.2,'EdgeAlpha',0, ...
        'LineStyle','none','HandleVisibility','off')
    ha2.BaseValue = -50;
    ha2.BaseLine.LineStyle = 'none';
    uistack(ha2,'bottom');
    
    ha3 = area([tF2_1 tF2_1 tF2_2 tF2_2]/3600-t_rem,[ylim fliplr(ylim)]);
    set(ha3,'facecolor',collistF,'edgecolor',collistF,'FaceAlpha',0.2,'EdgeAlpha',0, ...
        'LineStyle','none','HandleVisibility','off')
    ha3.BaseLine.LineStyle = 'none';
    ha3.BaseValue = -50;
    uistack(ha3,'bottom');
    
    ha4 = area([tE2_1 tE2_1 tE2_2 tE2_2]/3600-t_rem,[ylim fliplr(ylim)]);
    set(ha4,'facecolor',collistE,'edgecolor',collistE,'FaceAlpha',0.2,'EdgeAlpha',0, ...
        'LineStyle','none','HandleVisibility','off')
    ha4.BaseValue = -50;
    ha4.BaseLine.LineStyle = 'none';
    uistack(ha4,'bottom');
    
    ha5 = area([t_end_1 t_end_1 t_ini_2-t_rem t_ini_2-t_rem],[ylim fliplr(ylim)]);
    set(ha5,'facecolor',[1 1 1]*0,'edgecolor',[1 1 1]*0,'FaceAlpha',1,'EdgeAlpha',1, ...
        'LineStyle','none','HandleVisibility','off')
    ha5.BaseLine.LineStyle = 'none';
    uistack(ha5,'top');
    
    text(t_end_1+0.25,mean(ylim),'air sample preserved until the second night', ...
        'Rotation',90,'color','w','HorizontalAlignment','center','fontsize',ftsz-2,'fontname',ftnm)
    
    % Axes and labels
    xlabel('elapsed time {\itt} [h]','fontsize',ftsz,'fontname',ftnm)
    ylabel(sprintf('Herriott cell pressure\n{\\itp}_H({\\itt}) [Pa]'),'fontsize',ftsz,'fontname',ftnm)
    
    SS_mk_label_panel(lab_panel{isub-1},xlim,'lin',ylim,'lin',ftsz,ftnm)
    
    set(gca,'XMinorTick','on','YMinorTick','on','box','on')
    
    xt  = get(gca,'XTick');
    xtl = get(gca,'XTickLabel');
    kx  = find(xt>=t_ini_2-t_rem); nkx = length(kx);
    for ik = 1:nkx
        xtl(kx(ik)) = {num2str(xt(kx(ik))+t_rem)};
    end
    set(gca,'XTickLabel',xtl)
    
    set(gca,'fontsize',ftsz,'fontname',ftnm)
    
    %% 2) N_CH4_HC_t vs. time
    hs2 = subplot(nrow,ncol,isub); isub = isub + 1;
    
    kx1 = time_h <= t_end_1;
    plot(time_h(kx1),N_CH4_HC_t(kx1)*1e9,'color','k','linewidth',2); hold on
    
    kx2 = time_h >= t_ini_2;
    plot(time_h(kx2)-t_rem,N_CH4_HC_t(kx2)*1e9,'color','k','linewidth',2); hold on
    
    % Full-cell and empty cell periods
    xl = xlim; xlim([0 xl(2)])
    yl = ylim; ylim([0 yl(2)])
    
    ha1 = area([tF1_1 tF1_1 tF1_2 tF1_2]/3600,[ylim fliplr(ylim)]);
    set(ha1,'facecolor',collistF,'edgecolor',collistF,'FaceAlpha',0.2,'EdgeAlpha',0, ...
        'LineStyle','none','HandleVisibility','off')
    ha1.BaseLine.LineStyle = 'none';
    ha1.BaseValue = -50;
    uistack(ha1,'bottom');
    
    ha2 = area([tE1_1 tE1_1 tE1_2 tE1_2]/3600,[ylim fliplr(ylim)]);
    set(ha2,'facecolor',collistE,'edgecolor',collistE,'FaceAlpha',0.2,'EdgeAlpha',0, ...
        'LineStyle','none','HandleVisibility','off')
    ha2.BaseValue = -50;
    ha2.BaseLine.LineStyle = 'none';
    uistack(ha2,'bottom');
    
    ha3 = area([tF2_1 tF2_1 tF2_2 tF2_2]/3600-t_rem,[ylim fliplr(ylim)]);
    set(ha3,'facecolor',collistF,'edgecolor',collistF,'FaceAlpha',0.2,'EdgeAlpha',0, ...
        'LineStyle','none','HandleVisibility','off')
    ha3.BaseLine.LineStyle = 'none';
    ha3.BaseValue = -50;
    uistack(ha3,'bottom');
    
    ha4 = area([tE2_1 tE2_1 tE2_2 tE2_2]/3600-t_rem,[ylim fliplr(ylim)]);
    set(ha4,'facecolor',collistE,'edgecolor',collistE,'FaceAlpha',0.2,'EdgeAlpha',0, ...
        'LineStyle','none','HandleVisibility','off')
    ha4.BaseValue = -50;
    ha4.BaseLine.LineStyle = 'none';
    uistack(ha4,'bottom');
    
    ha5 = area([t_end_1 t_end_1 t_ini_2-t_rem t_ini_2-t_rem],[ylim fliplr(ylim)]);
    set(ha5,'facecolor',[1 1 1]*0,'edgecolor',[1 1 1]*0,'FaceAlpha',1,'EdgeAlpha',1, ...
        'LineStyle','none','HandleVisibility','off')
    ha5.BaseLine.LineStyle = 'none';
    uistack(ha5,'top');
    
    text(t_end_1+0.25,mean(ylim),'air sample preserved until the second night', ...
        'Rotation',90,'color','w','HorizontalAlignment','center','fontsize',ftsz-2,'fontname',ftnm)
    
    % Axes and labels
    xlabel('elapsed time {\itt} [h]','fontsize',ftsz,'fontname',ftnm)
    ylabel(sprintf('CH_4 amount in the Herriott cell\n{\\itN}_{CH_4}({\\itt}) [nmol]'), ...
        'fontsize',ftsz,'fontname',ftnm)
    
    SS_mk_label_panel(lab_panel{isub-1},xlim,'lin',ylim,'lin',ftsz,ftnm)
    
    set(gca,'XMinorTick','on','YMinorTick','on','box','on')
    
    xt  = get(gca,'XTick');
    xtl = get(gca,'XTickLabel');
    kx  = find(xt>=t_ini_2-t_rem); nkx = length(kx);
    for ik = 1:nkx
        xtl(kx(ik)) = {num2str(xt(kx(ik))+t_rem)};
    end
    set(gca,'XTickLabel',xtl)
    
    set(gca,'fontsize',ftsz,'fontname',ftnm)
    
    pos_ax = get(gca,'position'); pos_ax(1) = pos_ax(1)-0.03;
    set(gca,'position',pos_ax)
    
    %% 3) eta_HC_t vs. time
    hs3 = subplot(nrow,ncol,isub); isub = isub + 1;
    
    kx1 = time_h <= t_end_1;
    plot(time_h(kx1),eta_HC_t(kx1),'color','k','linewidth',2); hold on
    
    kx2 = time_h >= t_ini_2;
    plot(time_h(kx2)-t_rem,eta_HC_t(kx2),'color','k','linewidth',2); hold on
    
    set(gca,'YScale','log')
    
    xl = xlim; xlim([0 xl(2)])
    ylim([1 1e3])
    
    ha1 = area([tF1_1 tF1_1 tF1_2 tF1_2]/3600,[ylim fliplr(ylim)]);
    set(ha1,'facecolor',collistF,'edgecolor',collistF,'FaceAlpha',0.2,'EdgeAlpha',0, ...
        'LineStyle','none','HandleVisibility','off')
    ha1.BaseLine.LineStyle = 'none';
    uistack(ha1,'bottom');
    
    ha2 = area([tE1_1 tE1_1 tE1_2 tE1_2]/3600,[ylim fliplr(ylim)]);
    set(ha2,'facecolor',collistE,'edgecolor',collistE,'FaceAlpha',0.2,'EdgeAlpha',0, ...
        'LineStyle','none','HandleVisibility','off')
    ha2.BaseLine.LineStyle = 'none';
    uistack(ha2,'bottom');
    
    ha3 = area([tF2_1 tF2_1 tF2_2 tF2_2]/3600-t_rem,[ylim fliplr(ylim)]);
    set(ha3,'facecolor',collistF,'edgecolor',collistF,'FaceAlpha',0.2,'EdgeAlpha',0, ...
        'LineStyle','none','HandleVisibility','off')
    ha3.BaseLine.LineStyle = 'none';
    uistack(ha3,'bottom');
    
    ha4 = area([tE2_1 tE2_1 tE2_2 tE2_2]/3600-t_rem,[ylim fliplr(ylim)]);
    set(ha4,'facecolor',collistE,'edgecolor',collistE,'FaceAlpha',0.2,'EdgeAlpha',0, ...
        'LineStyle','none','HandleVisibility','off')
    ha4.BaseLine.LineStyle = 'none';
    uistack(ha4,'bottom');
    
    ha5 = area([t_end_1 t_end_1 t_ini_2-t_rem t_ini_2-t_rem],[ylim fliplr(ylim)]);
    set(ha5,'facecolor',[1 1 1]*0,'edgecolor',[1 1 1]*0,'FaceAlpha',1,'EdgeAlpha',1, ...
        'LineStyle','none','HandleVisibility','off')
    ha5.BaseLine.LineStyle = 'none';
    uistack(ha5,'top');
    
    text(t_end_1+0.25,exp(mean(log(ylim))),'air sample preserved until the second night', ...
        'Rotation',90,'color','w','HorizontalAlignment','center','fontsize',ftsz-2,'fontname',ftnm)
    
    % Axes and labels
    xlabel('elapsed {\itt} time [h]','fontsize',ftsz,'fontname',ftnm)
    ylabel(sprintf('CH_4 vmr in the Herriott cell\n\\eta_{CH_4}({\\itt}) [ppbv]'),'fontsize',ftsz,'fontname',ftnm)
    
    SS_mk_label_panel(lab_panel{isub-1},xlim,'lin',ylim,'log',ftsz,ftnm)
    
    set(gca,'XMinorTick','on','YMinorTick','on','box','on')
    
    xt  = get(gca,'XTick');
    xtl = get(gca,'XTickLabel');
    kx  = find(xt>=t_ini_2-t_rem); nkx = length(kx);
    for ik = 1:nkx
        xtl(kx(ik)) = {num2str(xt(kx(ik))+t_rem)};
    end
    set(gca,'XTickLabel',xtl)
    
    set(gca,'fontsize',ftsz,'fontname',ftnm)
    
    %% 4) Simulated eta_s, eta_E, and eta_F
    hs4 = subplot(nrow,ncol,isub); isub = isub + 1;
    legname = cell(1,1);
    ileg    = 1;
    
    yl = [0 160];
    ylim(yl)
    
    SS_sim_data
    
    % eta_sim
    kx1 = time_h <= t_end_1;
    plot(time_h(kx1),eta_sim(kx1),'color','k','linewidth',2); hold on
    
    kx2 = time_h >= t_ini_2;
    plot(time_h(kx2)-t_rem,eta_sim(kx2),'color','k','linewidth',2); hold on
    
    % empty-cell runs 1
    colsimE = [1 1 1]*0.7;
    plot(tE1_pts/3600,eta_E1_sim,'Marker','o','MarkerSize',6, ...
        'MarkerFaceColor',colsimE,'MarkerEdgeColor',colsimE,'color',colsimE); hold on
    
    errorbar(mean_tE1/3600,mean_eta_E1_sim,delta_eta_E1_sim,'Marker','o','MarkerSize',8, ...
        'MarkerFaceColor',collistE,'MarkerEdgeColor',collistE,'color',collistE,'LineWidth',2); hold on
    
    % full-cell runs 1
    colsimF = [1 1 1]*0.4;
    plot(tF1_pts/3600,eta_F1_sim,'Marker','o','MarkerSize',6, ...
        'MarkerFaceColor',colsimF,'MarkerEdgeColor',colsimF,'color',colsimF); hold on
    
    errorbar(mean_tF1/3600,mean_eta_F1_sim,delta_eta_F1_sim,'Marker','o','MarkerSize',8, ...
        'MarkerFaceColor',collistF,'MarkerEdgeColor',collistF,'color',collistF,'LineWidth',2); hold on
    
    % full-cell runs 2
    colsimF = [1 1 1]*0.4;
    plot(tF2_pts/3600-t_rem,eta_F2_sim,'Marker','s','MarkerSize',6, ...
        'MarkerFaceColor',colsimF,'MarkerEdgeColor',colsimF,'color',colsimF); hold on
    
    errorbar(mean_tF2/3600 - t_rem,mean_eta_F2_sim,delta_eta_F2_sim,'Marker','s','MarkerSize',8, ...
        'MarkerFaceColor',collistF,'MarkerEdgeColor',collistF,'color',collistF,'LineWidth',2); hold on
    
    % empty-cell runs 1
    colsimE = [1 1 1]*0.7;
    plot(tE2_pts/3600-t_rem,eta_E2_sim,'Marker','s','MarkerSize',6, ...
        'MarkerFaceColor',colsimE,'MarkerEdgeColor',colsimE,'color',colsimE); hold on
    
    errorbar(mean_tE2/3600 - t_rem,mean_eta_E2_sim,delta_eta_E2_sim,'Marker','s','MarkerSize',8, ...
        'MarkerFaceColor',collistE,'MarkerEdgeColor',collistE,'color',collistE,'LineWidth',2); hold on
    
    xl = xlim;
    xlim([0 xl(2)])
    ylim([min(0,yl(1)) yl(2)])
    
    ha1 = area([tE1_1 tE1_1 tE1_2 tE1_2]/3600,[ylim fliplr(ylim)]);
    set(ha1,'facecolor',collistE,'edgecolor',collistE,'FaceAlpha',0.2,'EdgeAlpha',0, ...
        'LineStyle','none','HandleVisibility','off')
    ha1.BaseLine.LineStyle = 'none';
    uistack(ha1,'bottom');
    
    ha2 = area([tF1_1 tF1_1 tF1_2 tF1_2]/3600,[ylim fliplr(ylim)]);
    set(ha2,'facecolor',collistF,'edgecolor',collistF,'FaceAlpha',0.2,'EdgeAlpha',0, ...
        'LineStyle','none','HandleVisibility','off')
    ha2.BaseLine.LineStyle = 'none';
    uistack(ha2,'bottom');
    
    ha3 = area([tF2_1 tF2_1 tF2_2 tF2_2]/3600-t_rem,[ylim fliplr(ylim)]);
    set(ha3,'facecolor',collistF,'edgecolor',collistF,'FaceAlpha',0.2,'EdgeAlpha',0, ...
        'LineStyle','none','HandleVisibility','off')
    ha3.BaseLine.LineStyle = 'none';
    uistack(ha3,'bottom');
    
    ha4 = area([tE2_1 tE2_1 tE2_2 tE2_2]/3600-t_rem,[ylim fliplr(ylim)]);
    set(ha4,'facecolor',collistE,'edgecolor',collistE,'FaceAlpha',0.2,'EdgeAlpha',0, ...
        'LineStyle','none','HandleVisibility','off')
    ha4.BaseLine.LineStyle = 'none';
    uistack(ha4,'bottom');
    
    ha5 = area([t_end_1 t_end_1 t_ini_2-t_rem t_ini_2-t_rem],[ylim fliplr(ylim)]);
    set(ha5,'facecolor',[1 1 1]*0,'edgecolor',[1 1 1]*0,'FaceAlpha',1,'EdgeAlpha',1, ...
        'LineStyle','none','HandleVisibility','off')
    ha5.BaseLine.LineStyle = 'none';
    uistack(ha5,'top');
    
    text(t_end_1+0.25,mean(ylim),'air sample preserved until the second night', ...
        'Rotation',90,'color','w','HorizontalAlignment','center','fontsize',ftsz-2,'fontname',ftnm)
    
    eta_HC_1_sim = mean_eta_F1_sim - mean_eta_E1_sim;
    sig_HC_1_sim = sqrt(delta_eta_F1_sim^2/nFpts + delta_eta_E1_sim^2/nEpts);
    
    % Axes and labels
    xlabel('elapsed time {\itt} [h]','fontsize',ftsz,'fontname',ftnm)
    ylabel('CH_4 vmr \eta_m ({\itt}) [ppbv]','fontsize',ftsz,'fontname',ftnm)
        
    ht1 = text(mean(tE1_pts/3600),140,sprintf([num2str(mean_eta_E1_sim,'%2.2f'), ...
        '\n\\pm\n',num2str(delta_eta_E1_sim,'%2.2f'),'\nppbv'])); hold on
    set(ht1,'HorizontalAlignment','center','fontsize',ftsz-1,'fontname',ftnm, ...
        'color',collistE,'FontWeight','bold')
    
    ht2 = text(mean(tF1_pts/3600),140,sprintf([num2str(mean_eta_F1_sim,'%2.2f'), ...
        '\n\\pm\n',num2str(delta_eta_F1_sim,'%2.2f'),'\nppbv'])); hold on
    set(ht2,'HorizontalAlignment','center','fontsize',ftsz-1,'fontname',ftnm, ...
        'color',collistF,'FontWeight','bold')
    
    ht3 = text(mean(tF2_pts/3600)-t_rem,20,sprintf([num2str(mean_eta_F2_sim,'%2.2f'), ...
        '\n\\pm\n',num2str(delta_eta_F2_sim,'%2.2f'),'\nppbv'])); hold on
    set(ht3,'HorizontalAlignment','center','fontsize',ftsz-1,'fontname',ftnm, ...
        'color',collistF,'FontWeight','bold')
    
    ht4 = text(mean(tE2_pts/3600)-t_rem,140,sprintf([num2str(mean_eta_E2_sim,'%2.2f'), ...
        '\n\\pm\n',num2str(delta_eta_E2_sim,'%2.2f'),'\nppbv'])); hold on
    set(ht4,'HorizontalAlignment','center','fontsize',ftsz-1,'fontname',ftnm, ...
        'color',collistE,'FontWeight','bold')
    
    SS_mk_label_panel(lab_panel{isub-1},xlim,'lin',ylim,'lin',ftsz,ftnm)
    
    set(gca,'XMinorTick','on','YMinorTick','on','box','on')
    
    xt  = get(gca,'XTick');
    xtl = get(gca,'XTickLabel');
    kx  = find(xt>=t_ini_2-t_rem); nkx = length(kx);
    for ik = 1:nkx
        xtl(kx(ik)) = {num2str(xt(kx(ik))+t_rem)};
    end
    set(gca,'XTickLabel',xtl)
    
    set(gca,'fontsize',ftsz,'fontname',ftnm)
    pos_ax = get(gca,'position'); pos_ax(1) = pos_ax(1)-0.03;
    set(gca,'position',pos_ax)
    
    %% Save figure
    if (savefig == 1)
        sfpath     = '../Figures/';
        sffilename = 'Figure_07.png';
        sfflnm     = fullfile(sfpath,sffilename);
        exportgraphics(gcf,sfflnm,'Resolution',400)
    end
    
end