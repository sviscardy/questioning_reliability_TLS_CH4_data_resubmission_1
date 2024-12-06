contrib_HC_TLS = zeros(1,ntime)*NaN;

contrib_FO_TLS = N_FO/(V_FO)*l_FO;

kk = time >= F_elapsed_t(1) & time <= F_elapsed_t(end);
contrib_HC_TLS(kk) = N_HC/V_HC*l_HC;

kk = time >= E_elapsed_t(1) & time <= E_elapsed_t(end);
contrib_HC_TLS(kk) = 0;

q_fct_TLS  = contrib_HC_TLS./( contrib_FO_TLS);

tE_mean = mean(E_elapsed_t);

%% Make figure
ftsz = ftsz + 1;

figure(figid+10)
delete(get(gcf,'children'))
set(gcf,'paperpositionmode','auto')
set(gcf,'position',[200 200 1200 600])
set(gca,'fontsize',ftsz,'fontname',ftnm)
clear legname
legname = cell(1,2);
ileg = 1;

plot(time_h,q_fct,'color',[0 0 0],'LineWidth',2); hold on
legname{ileg} = 'model'; ileg = ileg + 1;

% TLS: full-cell runs
kk = time<=tF_2;
plot(time_h(kk),q_fct_TLS(kk),'color',collistF,'LineWidth',2); hold on
legname{ileg} = 'TLS (full-cell runs)'; ileg = ileg + 1;

% TLS: empty-cell runs
kk = time>=tE_1;
plot(time_h(kk),q_fct_TLS(kk),'color',collistE,'LineWidth',2); hold on
legname{ileg} = 'TLS (empty-cell runs)'; ileg = ileg + 1;

legend(legname,'fontsize',ftsz,'fontname',ftnm)

xl = xlim; xlim([0 xl(2)])
yl = ylim; ylim([0 1])

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

% Axes and labels
xlabel('time {\itt} [h]','fontsize',ftsz,'fontname',ftnm)
ylabel('{\itq}({\itt}) [/]','fontsize',ftsz,'fontname',ftnm)

set(gca,'XMinorTick','on','YMinorTick','on','box','on')
set(gca,'fontsize',ftsz,'fontname',ftnm)

%% Save figure
sfpath     = '../Figures/';
sffilename = 'Figure_S03.png';
sfflnm     = fullfile(sfpath,sffilename);
exportgraphics(gcf,sfflnm,'Resolution',400)
