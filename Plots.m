clear; clc;
[meta1, output1] = load_simulation_output('simulation_sigma035c015.txt');
[meta2, output2] = load_simulation_output('simulation_sigma07c015.txt');
[meta3, output3] = load_simulation_output('simulation_test.txt');
[meta4, output4] = load_simulation_output('simulation_4.txt');
t_num=floor(meta3.time/(meta3.FrameSkip*meta3.h))+1;
timegrid = (0:t_num-1) * meta3.FrameSkip * meta3.h;
delstar=(meta3.cstar-1)/24;
N = floor(meta3.time*1.2);
%% Visualization

% Use LaTeX fonts globally for consistency
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultLegendInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(groot, 'defaultAxesFontName', 'Times New Roman') % optional: LaTeX-like font
set(groot, 'defaultAxesFontSize', 15)

%%
timegrid = (0:t_num-1) * meta2.FrameSkip * meta2.h;
fig1 = figure;
fig1.Position = [300, 300, 600, 300];
plot((-N:N), output2.r(6,:), 'b.','MarkerSize', 10); 
ylim( [min(output2.r(:))-delstar/20, max(output2.r(:))+delstar/20]);
xlim( [-(N/10)+timegrid(6)*meta2.c0, ...
             (N/20)+timegrid(6)*meta2.c0]);
xlabel('$i$')
exportgraphics(fig1, 'rt=40.pdf', 'ContentType', 'vector')
%%

t_num=floor(meta1.time/(meta1.FrameSkip*meta1.h))+1;
timegrid = (0:t_num-1) * meta1.FrameSkip * meta1.h;
fig2 = figure;
fig2.Position = [300, 300, 600, 300];

% Plot in color (for screen readability)
plot(timegrid, output1.cfit , 'b', 'LineWidth',1.2); hold on
plot(timegrid, output2.cfit , 'r', 'LineWidth',1.2); hold on
plot(timegrid,meta1.c0*ones(size(timegrid)),'k--')

xlabel('$t$','Interpreter','latex')

% --- Inline labels with manual offsets ---
xlab = timegrid(end);


hold on
dx = -210;   % shift left
dy  = -1.6e-04;        % shift down
xpos = xlab + dx;
ypos = output1.cfit(end)+dy;
plot([xpos-90, xpos-65], [ypos ypos], 'b', 'LineWidth', 2)

% Label for sigma=0.035 (blue curve)
text(xlab+dx, output1.cfit(end)+dy, '$\quad \sigma=0.035$', ...
    'Interpreter','latex', 'FontSize',14, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);

% Label for sigma=0.07 (red curve)
dx = -210;  
dy2 = -3.12e-04;        % larger shift down for 0.07 curve
xpos = xlab + dx;
ypos = output2.cfit(end)+dy2;
plot([xpos-90, xpos-65], [ypos ypos], 'r', 'LineWidth', 2)
text(xlab+dx, ypos, ' $\quad \sigma=0.07$ ', ...
    'Interpreter','latex', 'FontSize',14, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);

exportgraphics(fig2, 'sigmas.pdf', 'ContentType', 'vector')
%%
t_num=floor(meta2.time/(meta2.FrameSkip*meta2.h))+1;
timegrid = (0:t_num-1) * meta2.FrameSkip * meta2.h;
fig3a = figure;
fig3a.Position = [300, 300, 600, 300];

% Plot in color
plot(timegrid, output2.cfit, 'b', 'LineWidth',1.2); hold on
plot(timegrid, 1 + output2.delred/24, 'g', 'LineWidth',1.2);
plot(timegrid, 1.015*ones(size(timegrid)), 'k--', 'LineWidth',1.0); % reference
hold off

xlabel('$t$', 'Interpreter','latex')

% --- Inline labels (boxed) ---
xlab = timegrid(end);



% Blue: c_fitxpos = xlab + dx;
dx   = -80;  % left shift
xpos = xlab + dx;
dy_b = 2e-04;   % blue (c_fit)
ypos = output2.cfit(end) + dy_b;
hold on
plot([xpos-50, xpos-25], [ypos ypos], 'b', 'LineWidth', 2.5)

text(xlab+dx, ypos, '\quad $c(t)$', ...
    'Interpreter','latex', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);

% Green: c_1
dx   = -140;  % left shift
xpos = xlab + dx;
dy_g = 2e-04;   % green (c_1)
ypos = 1 + output2.delred(end)/24 + dy_g;
hold on
plot([xpos-115, xpos-90], [ypos ypos], 'g', 'LineWidth', 2)
text(xlab+dx, ypos, '\quad $c_*+\sigma c_1(t)$', ...
    'Interpreter','latex', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);

% Export
exportgraphics(fig3a, 'amplitude0.pdf', 'ContentType', 'vector');

%%
t_num=floor(meta2.time/(meta2.FrameSkip*meta2.h))+1;
timegrid = (0:t_num-1) * meta2.FrameSkip * meta2.h;



fig3b = figure;
fig3b.Position = [300, 300, 600, 300];
dt=timegrid(2)-timegrid(1);
plot(timegrid, output2.lfit- cumsum([zeros(1,1) output2.cfit(1:end-1)])*dt, 'b', 'LineWidth',1.2); hold on
plot(timegrid, output2.lred- timegrid*1.015, 'g', 'LineWidth',1.2);
plot(timegrid,zeros(size(timegrid)),'k--');

xlabel('$t$', 'Interpreter','latex')

% --- Inline labels (boxed) ---
xlab = timegrid(end);

% Blue: c_fitxpos = xlab + dx;
dx   = -130;  % left shift
xpos = xlab + dx;
dy_b = -0.2;   % blue (c_fit)
ypos = output2.lfit(end)- timegrid(end)*1.015 + dy_b;
hold on
plot([xpos-50, xpos-25], [ypos ypos], 'b', 'LineWidth', 2)

text(xlab+dx, ypos, '\quad $\gamma(t)$', ...
    'Interpreter','latex', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);

% Green: c_1
dx   = -130;  % left shift
xpos = xlab + dx;
dy_g = 0.7;   % green (c_1)
ypos = output2.lred(end)- timegrid(end)*1.015 + dy_g;
hold on
plot([xpos-70, xpos-45], [ypos ypos], 'g', 'LineWidth', 2)
text(xlab+dx, ypos, '\quad $\sigma \gamma_{1}(t)$', ...
    'Interpreter','latex', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);

exportgraphics(fig3b, 'position0.pdf', 'ContentType', 'vector')

%%
totalnorm1 = vecnorm(output2.eta1free, 2, 2) + vecnorm(output2.eta1modapprox, 2, 2);  % size: [nt, 1]
diffnorm1 = vecnorm(output2.eta1free - output2.eta1modapprox, 2, 2);                         % size: [nt, 1]
totalnorm2 = vecnorm(output2.eta2free, 2, 2) + vecnorm(output2.eta2modapprox, 2, 2);  % size: [nt, 1]
diffnorm2 = vecnorm(output2.eta2free - output2.eta2modapprox, 2, 2);                         % size: [nt, 1]

t_num=floor(meta2.time/(meta2.FrameSkip*meta2.h))+1;
timegrid = (0:t_num-1) * meta2.FrameSkip * meta2.h;

% Plot
fig4 = figure;
fig4.Position = [300, 300, 600, 300];
plot(timegrid, diffnorm1 ./ totalnorm1,'b','LineWidth', 1.2)
hold on
plot(timegrid, diffnorm2 ./ totalnorm2,'r','LineWidth', 1.2)
markerr=plot(timegrid(11), diffnorm1(11)/ totalnorm1(11), 'bs', 'MarkerSize',8, 'MarkerFaceColor','b');
markerp=plot(timegrid(11), diffnorm2(11)/ totalnorm2(11), 'r^', 'MarkerSize',8, 'MarkerFaceColor','r');
ylim([0,0.25])
legend([markerr, markerp], {'$r$-component', '$p$-component'}, ...
       'Interpreter','latex')
xlabel('$t$')

exportgraphics(fig4, 'relative.pdf', 'ContentType', 'vector')

%%
fig5a = figure;
fig5a.Position = [300, 300, 600, 300];
plot(timegrid, output4.cfit , 'b','LineWidth', 1.2)
hold on
plot(timegrid, 1+output4.dellintest/24, 'g','LineWidth', 1.2)


dx   = -220;  % left shift
xpos = xlab + dx;
dy_b = -3e-04;   % blue (c_fit)
ypos = output4.cfit(end) + dy_b;
hold on
plot([xpos-50, xpos-25], [ypos ypos], 'b', 'LineWidth', 2)

text(xlab+dx, ypos, '\quad $c(t)$', ...
    'Interpreter','latex', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);


xpos = xlab + dx;
dy_b = 0.5e-03;   % blue (c_fit)
ypos = 1+output4.dellintest(end)/24 + dy_b;
hold on
plot([xpos-200, xpos-175], [ypos ypos], 'g', 'LineWidth', 2)

text(xlab+dx, ypos, '\quad $c_*+\sigma c_1(t)+\sigma^2 c_2(t)$', ...
    'Interpreter','latex', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);

%legend('$c_{\mathrm{fit}}(t)$','$\overline{c}(t)$')
xlabel('$t$')

exportgraphics(fig5a, 'amplitude.pdf', 'ContentType', 'vector')
%%
t_num=floor(meta4.time/(meta4.FrameSkip*meta4.h))+1;
timegrid = (0:t_num-1) * meta4.FrameSkip * meta4.h;

fig5b = figure;
fig5b.Position = [300, 300, 600, 300];
dt=timegrid(2)-timegrid(1);
plot(timegrid, output4.lfit- cumsum([zeros(1,1) output4.cfit(1:end-1)])*dt , 'b','LineWidth', 1.2)
hold on
plot(timegrid, output4.llintest- cumsum([zeros(1,1) (1+output4.delexp(1:end-1)/24)])*dt, 'g','LineWidth', 1.2)
plot(timegrid,zeros(size(timegrid)),'k--')

ylim([-5, 1])

dx   = -260;  % left shift
xpos = xlab + dx;
dy_g = -4;   % green 
ypos =  dy_g;
hold on
plot([xpos-160, xpos-135], [ypos ypos], 'g', 'LineWidth', 2)
text(xlab+dx, ypos, '\quad $\sigma\gamma_1(t)+\sigma^2\gamma_2(t)$', ...
    'Interpreter','latex', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);

dx   = -260;  % left shift
xpos = xlab + dx;
dy= -2;   % green 
ypos =  dy;
hold on
plot([xpos-50, xpos-25], [ypos ypos], 'b', 'LineWidth', 2)
text(xlab+dx, ypos, '\quad $\gamma(t)$', ...
    'Interpreter','latex', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);

xlabel('$t$')
exportgraphics(fig5b, 'position.pdf', 'ContentType', 'vector')
%%

fig6a = figure;
fig6a.Position = [300, 300, 600, 300];
plot(timegrid, (output4.dellintest/24-output4.delexp/24)/(meta4.eps.^2), 'g','LineWidth', 1.2)
hold on
plot(timegrid, (output4.delhtest/24-output4.delexp/24)/(meta4.eps.^2), 'm','LineWidth', 1.2)
%plot(timegrid,0*ones(size(timegrid)),'k--')

dx   = -260;  % left shift
xpos = xlab + dx;
dy_b = -0.08;   % blue (c_fit)
ypos = dy_b;
hold on
plot([xpos-55, xpos-30], [ypos ypos], 'g', 'LineWidth', 2)

text(xlab+dx, ypos, '\quad $c_2(t)$', ...
    'Interpreter','latex', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);
dx   = -260;  % left shift
xpos = xlab + dx;
dy_b = -0.17;   % blue (c_fit)
ypos =  dy_b;
hold on
plot([xpos-55, xpos-30], [ypos ypos], 'm', 'LineWidth', 2)

text(xlab+dx, ypos, '\quad $c^{\mathrm{h}}_2(t)$', ...
    'Interpreter','latex', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);


xlabel('$t$')

exportgraphics(fig6a, '2vsap.pdf', 'ContentType', 'vector')

%%
fig6b = figure;
fig6b.Position = [300, 300, 600, 300];
plot(timegrid, (output4.llintest- cumsum([zeros(1,1) (1+output4.delexp(1:end-1)/24)])*dt-output4.lexp+timegrid*1.015)/(meta4.eps.^2), 'g','LineWidth', 1.2)
hold on
plot(timegrid, (output4.lhtest- cumsum([zeros(1,1) (1+output4.delh(1:end-1)/24)])*dt-output4.lexp+timegrid*1.015)/(meta4.eps.^2), 'm','LineWidth', 1.2)

xlabel('$t$')

dx   = -200;  % left shift
xpos = xlab + dx;
dy_b = -550;   % blue (c_fit)
ypos = output4.llintest(end) - timegrid(end)*1.015 + dy_b;
hold on
plot([xpos-60, xpos-35], [ypos ypos], 'g', 'LineWidth', 2)

text(xlab+dx, ypos, '\quad ${\gamma}_{2}(t)$', ...
    'Interpreter','latex', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);

dx   = -200;  % left shift
xpos = xlab + dx;
dy_b = -300;   % blue (c_fit)
ypos = output4.lhtest(end) - timegrid(end)*1.015 + dy_b;
hold on
plot([xpos-60, xpos-35], [ypos ypos], 'm', 'LineWidth', 2)

text(xlab+dx, ypos, '\quad ${\gamma}_{2}^{\mathrm{h}}(t)$', ...
    'Interpreter','latex', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);

exportgraphics(fig6b, '2vsappos.pdf', 'ContentType', 'vector')


%%
rate_inf=-2.23e-4+6.528e-5;
rate_inftrue = -1.403708886574337e-04;
t_num=floor(meta2.time/(meta2.FrameSkip*meta2.h))+1;
timegrid = (0:t_num-1) * meta2.FrameSkip * meta2.h;

fig8a = figure;
fig8a.Position = [300, 300, 600, 300];

plot(timegrid, output2.cfit, 'b','LineWidth', 1.2)
hold on

plot(timegrid, meta2.c0+meta2.eps^2*rate_inf*timegrid,'k.','LineWidth', 1.2)
plot(timegrid, meta2.c0+meta2.eps^2*rate_inftrue*timegrid,'k--','LineWidth', 1.2)
hold off
legend('$c(t)$')
xlabel('$t$')


exportgraphics(fig8a, 'rate.pdf', 'ContentType', 'vector')

%%
t_num=floor(meta3.time/(meta3.FrameSkip*meta3.h))+1;
timegrid = (0:t_num-1) * meta3.FrameSkip * meta3.h;
fig9a = figure;
fig9a.Position = [300, 300, 600, 300];

% Plot in color
plot(timegrid, output3.cfit, 'b', 'LineWidth',1.2); hold on
plot(timegrid, 1 + output3.delmod/24, 'r', 'LineWidth',1.2);
plot(timegrid, 1.015*ones(size(timegrid)), 'k--', 'LineWidth',1.0); % reference
hold off

xlabel('$t$', 'Interpreter','latex')

% --- Inline labels (boxed) ---
xlab = timegrid(end);



% Blue: c_fitxpos = xlab + dx;
dx   = -100;  % left shift
xpos = xlab + dx;
dy_b = 3e-04;   % blue (c_fit)
ypos = output3.cfit(end) + dy_b;
hold on
plot([xpos-65, xpos-40], [ypos ypos], 'b', 'LineWidth', 2.5)

text(xlab+dx, ypos, '\quad $c_{\mathrm{fit}}(t)$', ...
    'Interpreter','latex', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);

%Red: c_mod
dx=-220;
dy_r = -3.2e-04;   % red  (c_mod) - further down
xpos = xlab + dx;
ypos = 1+output3.delmod(end)/24 + dy_r;
hold on
plot([xpos-80, xpos-55], [ypos ypos], 'r', 'LineWidth', 2)
text(xlab+dx, ypos, '\quad $c_{\mathrm{mod}}(t)$', ...
    'Interpreter','latex', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);

exportgraphics(fig9a, 'modvsfit.pdf', 'ContentType', 'vector');

%%
fig9b = figure;
fig9b.Position = [300, 300, 600, 300];
dt=timegrid(2)-timegrid(1);
plot(timegrid, output3.lfit- cumsum([zeros(1,1) output3.cfit(1:end-1)])*dt, 'b', 'LineWidth',1.2); hold on
cmod=(1 + output3.delmod/24);
plot(timegrid, output3.lmod- cumsum([zeros(1,1) cmod(1:end-1)])*dt, 'r', 'LineWidth',1.2);
plot(timegrid,zeros(size(timegrid)),'k--');

xlabel('$t$', 'Interpreter','latex')

% --- Inline labels (boxed) ---
xlab = timegrid(end);

% Blue: c_fitxpos = xlab + dx;
dx   = -130;  % left shift
xpos = xlab + dx;
dy_b = 1.1;   % blue (c_fit)
ypos = output3.lfit(end)- timegrid(end)*1.015 + dy_b;
hold on
plot([xpos-65, xpos-40], [ypos ypos], 'b', 'LineWidth', 2)

text(xlab+dx, ypos, '\quad $\gamma_{\mathrm{fit}}(t)$', ...
    'Interpreter','latex', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);

% % Red: c_mod
dx=-130;
dy_r = 2.5;   % red  (c_mod) - further down
xpos = xlab + dx;
ypos = output3.lmod(end)- timegrid(end)*1.015 + dy_r;
hold on
plot([xpos-75, xpos-50], [ypos ypos], 'r', 'LineWidth', 2)
text(xlab+dx, ypos, '\quad $\gamma_{\mathrm{mod}}(t)$', ...
    'Interpreter','latex', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);

exportgraphics(fig9b, 'modvsfitpos.pdf', 'ContentType', 'vector')

%%
figtest = figure;
figtest.Position = [300, 300, 600, 300];

% Plot in color
plot(timegrid, output3.cfit, 'b', 'LineWidth',1.2); hold on
plot(timegrid, 1+output3.del2/24, 'r', 'LineWidth',1.2); hold on
plot(timegrid, 1.015*ones(size(timegrid)), 'k--', 'LineWidth',1.0); % reference
hold off

xlabel('$t$', 'Interpreter','latex')


xlab = timegrid(end);

% Blue: c_fitxpos = xlab + dx;
dx   = -200;  % left shift
xpos = xlab + dx;
dy_b = -3.5e-04;   % blue (c_fit)
ypos = output3.cfit(end)+ dy_b;
hold on
plot([xpos-50, xpos-25], [ypos ypos], 'b', 'LineWidth', 2)

text(xlab+dx, ypos, '\quad $c(t)$', ...
    'Interpreter','latex', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);

% Blue: c_fitxpos = xlab + dx;
dx   = -200;  % left shift
xpos = xlab + dx;
dy_b = 5e-04;   % blue (c_fit)
ypos = output3.cfit(end)+ dy_b;
hold on
plot([xpos-55, xpos-30], [ypos ypos], 'r', 'LineWidth', 2)

text(xlab+dx, ypos, '\quad $c_2(t)$', ...
    'Interpreter','latex', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);

exportgraphics(figtest, 'c2.pdf', 'ContentType', 'vector')

%%
[metaens01, ensemble01] = load_simulation_output('ensemble.txt');
[metaens02, ensemble02] = load_simulation_output('ensemble02.txt');
[metaens03, ensemble03] = load_simulation_output('ensemble03.txt');
[metaens035, ensemble035] = load_simulation_output('ensemble035.txt');

t_num=floor(metaens01.time/(metaens01.FrameSkip*metaens01.h))+1;
timegrid = (0:t_num-1) * metaens01.FrameSkip * metaens01.h;

t_num03=floor(metaens03.time/(metaens03.FrameSkip*metaens03.h))+1;
timegrid03 = (0:t_num03-1) * metaens03.FrameSkip * metaens03.h;

t_num035=floor(metaens035.time/(metaens035.FrameSkip*metaens035.h))+1;
timegrid035 = (0:t_num035-1) * metaens035.FrameSkip * metaens035.h;

figmean = figure;
figmean.Position = [300, 300, 600, 300];
plot(timegrid*(0.1)^2,mean(ensemble01.cfit_ens), 'b', 'LineWidth',1.2);

[~, idx008] = min(abs(mean(ensemble01.cfit_ens) - 1.008));
[~, idx006] = min(abs(mean(ensemble01.cfit_ens) - 1.006));

% hold on
% plot(timegrid,cfit_ens(1,:))
hold on
rate_inftrue = -1.403708886574337e-04;
rate_inf=-2.23e-4+6.528e-5;
rate_006=-3.744e-5+1.0951e-5;
rate_006true=-2.581216824156394e-5 ;
rate_008=-6.549e-5+1.9297e-5;
rate_008true= -4.444191090522610e-05;

plot(timegrid*(0.1)^2,metaens01.cstar+(0.1)^2*rate_inftrue*timegrid,'k--','LineWidth', 1.2)
plot(timegrid*(0.1)^2,1.006-(0.1)^2*rate_006true*timegrid(idx006)+(0.1)^2*rate_006true*timegrid,'k--','LineWidth', 1.2)
plot(timegrid*(0.1)^2,1.008-(0.1)^2*rate_008true*timegrid(idx008)+(0.1)^2*rate_008true*timegrid,'k--','LineWidth', 1.2)

plot((0.1)^2*timegrid(idx006), mean(ensemble01.cfit_ens(:,idx006)), 'bs', 'MarkerSize',8, 'MarkerFaceColor','b');
plot((0.1)^2*timegrid(idx008), mean(ensemble01.cfit_ens(:,idx008)), 'bs', 'MarkerSize',8, 'MarkerFaceColor','b');
plot(0, mean(ensemble01.cfit_ens(:,1)), 'bs', 'MarkerSize',8, 'MarkerFaceColor','b');


legend('$\mathbf{E}[c(t)]=\mathbf{E}[c(\tau/\sigma^2)]$')
ylim([1, metaens01.cstar])
xlim([0, 195])
xlabel('$\tau$')
exportgraphics(figmean, 'meantangents.pdf', 'ContentType', 'vector');


%%

[metaens01, ensemble01] = load_simulation_output('ensemble.txt');


t_num=floor(metaens01.time/(metaens01.FrameSkip*metaens01.h))+1;
timegrid = (0:t_num-1) * metaens01.FrameSkip * metaens01.h;



figmean = figure;
figmean.Position = [300, 300, 600, 300];
plot(timegrid*(0.1)^2,mean(ensemble01.cfit_ens), 'b', 'LineWidth',1.2);


% Initialize solution
ceuler = zeros(1,numel(timegrid));
ceuler(1) = 1.015;   % initial condition
dt=200/numel(timegrid);

% Euler forward loop
for n = 1:numel(timegrid)-1
    dc = -0.54*(ceuler(n)-1)^1.93;
    ceuler(n+1) = ceuler(n) +dt *dc;
end

% Initialize solution
ceulercor = zeros(1,numel(timegrid));
ceulercor(1) = 1.015;   % initial condition
dt=200/numel(timegrid);

% Euler forward loop
for n = 1:numel(timegrid)-1
    dc = -0.3471*(ceulercor(n)-1)^1.8577;
    ceulercor(n+1) = ceulercor(n) +dt *dc;
end


hold on
plot((0.1)^2*timegrid, ceuler,'k.', 'LineWidth', 1.5);
hold on
plot((0.1)^2*timegrid, ceulercor,'k--', 'LineWidth', 1.5);
legend('$\mathbf{E}[c(t)]=\mathbf{E}[c(\tau/\sigma^2)]$')
ylim([1, metaens01.cstar])
xlim([0, 180])
xlabel('$\tau$')
%ylabel('$\mathbf{E}[c(t)]$','Rotation',0)
exportgraphics(figmean, 'ODE.pdf', 'ContentType', 'vector');

%% Video

tic
Visualization(db, r, p, kappa, cstar, delstar,time, FrameSkip, h, eta1fit, eta2fit, eta1mod, eta2mod, eta1modapprox, eta2modapprox, eta1free, eta2free, cfit,lfit, cmod,c2, lmod,llin);
toc
