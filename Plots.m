clear; clc;
[meta1, output1] = load_simulation_output('shorter_sigma035c015.txt');
[meta2, output2] = load_simulation_output('shorter_sigma07c015.txt');
[meta3, output3] = load_simulation_output('simulation_sigma07c007.txt');
[meta4, output4] = load_simulation_output('validation_sigma07c015.txt');
t_num=floor(meta4.time/(meta4.FrameSkip*meta4.h))+1;
timegrid = (0:t_num-1) * meta4.FrameSkip * meta4.h;
delstar=(meta4.cstar-1)/24;
N = floor(meta4.time*1.2);


output1.cfit(49)=output1.cfit(48);
output2.cfit(49)=output2.cfit(48);
output1.lfit(49)=output1.lfit(48);
output2.lfit(49)=output2.lfit(48);

output4.cfit(4)=output4.cfit(3);
output4.lfit(4)=output4.lfit(3);
%% Visualization

% Use LaTeX fonts globally for consistency
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultLegendInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(groot, 'defaultAxesFontName', 'Times New Roman') % optional: LaTeX-like font
set(groot, 'defaultAxesFontSize', 20)


%%
t_num=floor(meta2.time/(meta2.FrameSkip*meta2.h))+1;
timegrid = (0:t_num-1) * meta2.FrameSkip * meta2.h;
delstar=(meta2.c0-1)/24;
N = floor(meta2.time*1.2);

fig1 = figure;
fig1.Position = [300, 300, 600, 300];
plot((-N:N), output2.r(1,:), 'b.','MarkerSize', 10); 
ylim( [min(output2.r(:))-delstar/20, max(output2.r(:))+delstar/20]);
xlim( [-(N/10)+timegrid(1)*meta2.c0, ...
             (N/20)+timegrid(1)*meta2.c0]);
xlabel('$i$')
exportgraphics(fig1, 'rt=0.pdf', 'ContentType', 'vector')
%%
t_num=floor(meta1.time/(meta1.FrameSkip*meta1.h))+1;
timegrid = (0:t_num-1) * meta1.FrameSkip * meta1.h;
delstar=(meta1.cstar-1)/24;
N = floor(meta1.time*1.2);




fig2 = figure;
fig2.Position = [300, 300, 600, 300];

% Plot in color (for screen readability)
plot(timegrid, output1.cfit , 'b', 'LineWidth',1.2); hold on
plot(timegrid, output2.cfit , 'r', 'LineWidth',1.2); hold on
plot(timegrid,meta1.cstar*ones(size(timegrid)),'k--')

xlabel('$t$','Interpreter','latex')

% --- Inline labels with manual offsets ---
xlab = timegrid(end);


hold on
dx = -210;   % shift left
dy  = 4e-04;        % shift down
xpos = xlab + dx;
ypos = output1.cfit(end)+dy;
plot([xpos-120, xpos-90], [ypos ypos], 'b', 'LineWidth', 2)

% Label for sigma=0.035 (blue curve)
text(xlab+dx, output1.cfit(end)+dy, '$\quad \sigma=0.035 $', ...
    'Interpreter','latex', 'FontSize',20, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',1);

% Label for sigma=0.07 (red curve)
dx = -210;  
dy2 = -2.0e-04;        % larger shift down for 0.07 curve
xpos = xlab + dx;
ypos = output2.cfit(end)+dy2;
plot([xpos-120, xpos-90], [ypos ypos], 'r', 'LineWidth', 2)
text(xlab+dx, ypos, ' $\quad \sigma=0.07$ ', ...
    'Interpreter','latex', 'FontSize',20, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',1);
ylim([1.0138,1.0155])

exportgraphics(fig2, 'sigmas.pdf', 'ContentType', 'vector')

%%
t_num=floor(meta2.time/(meta2.FrameSkip*meta2.h))+1;
timegrid = (0:t_num-1) * meta2.FrameSkip * meta2.h;
fig3a = figure;
fig3a.Position = [300, 300, 600, 300];

% Plot in color
plot(timegrid, output2.cfit, 'b', 'LineWidth',1.2); hold on
%plot(timegrid, 1 + output2.delmod/24, 'r', 'LineWidth',1.2);
plot(timegrid, 1 + output2.delexp/24, 'g', 'LineWidth',1.2);
plot(timegrid, 1.015*ones(size(timegrid)), 'k--', 'LineWidth',1.0); % reference
hold off

xlabel('$t$', 'Interpreter','latex')

% --- Inline labels (boxed) ---
xlab = timegrid(end);



% Blue: c_fitxpos = xlab + dx;
dx   = -80;  % left shift
xpos = xlab + dx;
dy_b = 4e-04;   % blue (c_fit)
ypos = output2.cfit(end) + dy_b;
hold on
plot([xpos-50, xpos-25], [ypos ypos], 'b', 'LineWidth', 2.5)

text(xlab+dx, ypos, '\quad $c(t)$', ...
    'Interpreter','latex', 'FontSize',20, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',1);


% Green: c_1
dx   = -160;  % left shift
xpos = xlab + dx;
dy_g = 4e-04;   % green (c_1)
ypos = 1 + output2.delexp(end)/24 + dy_g;
hold on
plot([xpos-130, xpos-105], [ypos ypos], 'g', 'LineWidth', 2)
text(xlab+dx, ypos, '\quad $c_*+\sigma c_1(t)$', ...
    'Interpreter','latex', 'FontSize',20, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',1);

ylim([1.014,1.0157])
% Export
exportgraphics(fig3a, 'amplitude0.pdf', 'ContentType', 'vector');

%%
t_num = floor(meta2.time/(meta2.FrameSkip*meta2.h)) + 1;
timegrid = (0:t_num-1) * meta2.FrameSkip * meta2.h;

% --- DEFINE ZOOM REGION ---
xZoom = [0, 100];
yZoom = [1.0147, 1.0155];

% --- FIGURE ---
fig3a = figure;
fig3a.Position = [300, 300, 600, 300];

% --- MAIN AXES ---
axMain = axes;
plot(timegrid, output2.cfit, 'b', 'LineWidth',1.2); hold on
plot(timegrid, 1 + output2.delexp/24, 'g', 'LineWidth',1.2);
plot(timegrid, 1.015*ones(size(timegrid)), 'k--', 'LineWidth',1.0);
hold off

xlabel('$t$', 'Interpreter','latex')
ylim([1.0130,1.0159])

xlab=1200;

% Blue: c_fitxpos = xlab + dx;
dx   = -80;  % left shift
xpos = xlab + dx;
dy_b = -4e-04;   % blue (c_fit)
ypos = output2.cfit(end) + dy_b;
hold on
plot([xpos-50, xpos-25], [ypos ypos], 'b', 'LineWidth', 2.5)

text(xlab+dx, ypos, '\quad $c(t)$', ...
    'Interpreter','latex', 'FontSize',20, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',1);


% Green: c_1
dx   = -160;  % left shift
xpos = xlab + dx;
dy_g = 5e-04;   % green (c_1)
ypos = 1 + output2.delexp(end)/24 + dy_g;
hold on
plot([xpos-130, xpos-105], [ypos ypos], 'g', 'LineWidth', 2)
text(xlab+dx, ypos, '\quad $c_*+\sigma c_1(t)$', ...
    'Interpreter','latex', 'FontSize',20, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',1);


% --- ZOOM RECTANGLE ---
hold on
rectangle('Position', ...
    [xZoom(1), yZoom(1), diff(xZoom), diff(yZoom)], ...
    'EdgeColor','k', ...
    'LineWidth',2, ...
    'LineStyle','--');

% --- INSET AXES (ZOOMED VIEW) ---
axInset = axes('Position',[0.15 0.2 0.32 0.35]);
axInset.XTick = [];
axInset.YTick = [];
box on
hold on

plot(timegrid, output2.cfit, 'b', 'LineWidth',1.2);
plot(timegrid, 1 + output2.delexp/24, 'g', 'LineWidth',1.2);
plot(timegrid, 1.015*ones(size(timegrid)), 'k--', 'LineWidth',1.0);

xlim(xZoom)
ylim(yZoom)
set(axInset,'FontSize',10)

% --- CONNECTOR LINES ---
inPos = axInset.Position;  % inset position [x y w h]

% Bottom-left of zoom rectangle → bottom-left of inset
p_bl = data2norm(axMain, [xZoom(1), yZoom(1)]);
annotation('line', [p_bl(1), inPos(1)], [p_bl(2), inPos(2)], 'Color','k');

% Top-right of zoom rectangle → top-right of inset
p_tr = data2norm(axMain, [xZoom(2), yZoom(2)]);
annotation('line', [p_tr(1), inPos(1)+inPos(3)], [p_tr(2), inPos(2)+inPos(4)], 'Color','k');


% --- EXPORT ---
set(fig3a,'Renderer','painters')
exportgraphics(fig3a, 'amplitudezoom.pdf', 'ContentType','vector');

% ==========================================================
function p = data2norm(ax, pt)
    axPos = ax.Position;
    xl = xlim(ax);
    yl = ylim(ax);

    p(1) = axPos(1) + axPos(3)*(pt(1)-xl(1))/diff(xl);
    p(2) = axPos(2) + axPos(4)*(pt(2)-yl(1))/diff(yl);
end


%%
t_num=floor(meta2.time/(meta2.FrameSkip*meta2.h))+1;
timegrid = (0:t_num-1) * meta2.FrameSkip * meta2.h;



fig3b = figure;
fig3b.Position = [300, 300, 600, 300];
dt=timegrid(2)-timegrid(1);
% Plot in color
%plot(timegrid, output2.lfit- timegrid*1.015, 'b', 'LineWidth',1.2); hold on
plot(timegrid, output2.lfit- cumsum([zeros(1,1) output2.cfit(1:end-1)])*dt, 'b', 'LineWidth',1.2); hold on
%plot(timegrid, output2.lmod- timegrid*1.015, 'r', 'LineWidth',1.2);
%cmod=(1 + output2.delmod/24);
%plot(timegrid, output2.lmod- cumsum([zeros(1,1) cmod(1:end-1)])*dt, 'r', 'LineWidth',1.2);
plot(timegrid, output2.lexp- timegrid*1.015, 'g', 'LineWidth',1.2);
plot(timegrid,zeros(size(timegrid)),'k--');

xlabel('$t$', 'Interpreter','latex')

% --- Inline labels (boxed) ---
xlab = timegrid(end);

% Blue: c_fitxpos = xlab + dx;
dx   = -130;  % left shift
xpos = xlab + dx;
dy_b = 1.0;   % blue (c_fit)
ypos = output2.lfit(end)- timegrid(end)*1.015 + dy_b;
hold on
plot([xpos-50, xpos-25], [ypos ypos], 'b', 'LineWidth', 2)

text(xlab+dx, ypos, '\quad $\gamma(t)$', ...
    'Interpreter','latex', 'FontSize',20, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',1);

% % Red: c_mod
% dx=-130;
% dy_r = 1.8;   % red  (c_mod) - further down
% xpos = xlab + dx;
% ypos = output2.lmod(end)- timegrid(end)*1.015 + dy_r;
% hold on
% plot([xpos-75, xpos-50], [ypos ypos], 'r', 'LineWidth', 2)
% text(xlab+dx, ypos, '\quad $\gamma_{\mathrm{mod}}(t)$', ...
%     'Interpreter','latex', 'FontSize',16, ...
%     'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
%     'EdgeColor','k', 'Margin',3);

% Green: c_1
dx   = -130;  % left shift
xpos = xlab + dx;
dy_g = 0.4;   % green (c_1)
ypos = output2.lexp(end)- timegrid(end)*1.015 + dy_g;
hold on
plot([xpos-80, xpos-55], [ypos ypos], 'g', 'LineWidth', 2)
text(xlab+dx, ypos, '\quad $\sigma \gamma_{1}(t)$', ...
    'Interpreter','latex', 'FontSize',20, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',1);

exportgraphics(fig3b, 'position0.pdf', 'ContentType', 'vector')

%%
totalnorm1 = vecnorm(output2.etarh, 2, 2) + vecnorm(output2.etarlin, 2, 2);  % size: [nt, 1]
diffnorm1 = vecnorm(output2.etarh - output2.etarlin, 2, 2);                         % size: [nt, 1]
totalnorm2 = vecnorm(output2.etaph, 2, 2) + vecnorm(output2.etaplin, 2, 2);  % size: [nt, 1]
diffnorm2 = vecnorm(output2.etaph - output2.etaplin, 2, 2);                         % size: [nt, 1]

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
%legend(hMarker, '$c_{\mathrm{fit}}$', 'Interpreter','latex')
legend([markerr, markerp], {'$r$-component', '$p$-component'}, ...
       'Interpreter','latex')
xlabel('$t$')

exportgraphics(fig4, 'relative.pdf', 'ContentType', 'vector')

%%
fig5a = figure;
fig5a.Position = [300, 300, 600, 300];
plot(timegrid, output2.cfit , 'b','LineWidth', 1.2)
hold on
plot(timegrid, 1+output2.dellintest/24, 'g','LineWidth', 1.2)
%plot(timegrid,1.007*ones(size(timegrid)),'k--')


dx   = -220;  % left shift
xpos = 1100 + dx;
dy_b = -1e-04;   % blue (c_fit)
ypos = output2.cfit(end) + dy_b;
hold on
plot([xpos-50, xpos-25], [ypos ypos], 'b', 'LineWidth', 2)

text(xpos, ypos, '\quad $c(t)$', ...
    'Interpreter','latex', 'FontSize',20, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',1);


dy_b = 0.8e-03;   % blue (c_fit)
ypos = 1+output2.dellintest(end)/24 + dy_b;
hold on
plot([xpos-240, xpos-215], [ypos ypos], 'g', 'LineWidth', 2)

text(xpos, ypos, '\quad $c_*+\sigma c_1(t)+\sigma^2 c_2(t)$', ...
    'Interpreter','latex', 'FontSize',20, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',1);

%legend('$c_{\mathrm{fit}}(t)$','$\overline{c}(t)$')
xlabel('$t$')

exportgraphics(fig5a, 'amplitude.pdf', 'ContentType', 'vector')
%%
t_num=floor(meta2.time/(meta2.FrameSkip*meta2.h))+1;
timegrid = (0:t_num-1) * meta2.FrameSkip * meta2.h;

fig5b = figure;
fig5b.Position = [300, 300, 600, 300];
dt=timegrid(2)-timegrid(1);
plot(timegrid, output2.lfit- cumsum([zeros(1,1) output2.cfit(1:end-1)])*dt , 'b','LineWidth', 1.2)
hold on
plot(timegrid, output2.llintest- cumsum([zeros(1,1) (1+output2.delexp(1:end-1)/24)])*dt, 'g','LineWidth', 1.2)
plot(timegrid,zeros(size(timegrid)),'k--')

% figure
% plot(timegrid,outputHJ.lhtest)
% hold on
% plot(timegrid,metaHJ.cstar*timegrid)
% hold on
% plot(timegrid, outputHJ.lfit)
ylim([-3.5, 0.5])

dx   = -260;  % left shift
xpos = 1100 + dx;
dy_g = -3;   % green 
ypos =  dy_g;
hold on
plot([xpos-190, xpos-165], [ypos ypos], 'g', 'LineWidth', 2)
text(xpos, ypos, '\quad $\sigma\gamma_1(t)+\sigma^2\gamma_2(t)$', ...
    'Interpreter','latex', 'FontSize',20, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',1);

dx   = -260;  % left shift

dy= -0.7;   % green 
ypos =  dy;
hold on
plot([xpos-50, xpos-25], [ypos ypos], 'b', 'LineWidth', 2)
text(xpos, ypos, '\quad $\gamma(t)$', ...
    'Interpreter','latex', 'FontSize',20, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',1);

xlabel('$t$')
exportgraphics(fig5b, 'position.pdf', 'ContentType', 'vector')
%%

fig6a = figure;
fig6a.Position = [300, 300, 600, 300];
plot(timegrid, (output2.dellintest/24-output2.delexp/24)/(meta2.eps.^2), 'g','LineWidth', 1.2)
hold on
plot(timegrid, (output2.delhtest/24-output2.delexp/24)/(meta2.eps.^2), 'm','LineWidth', 1.2)
%plot(timegrid,0*ones(size(timegrid)),'k--')

dx   = -260;  % left shift
xpos = 1100 + dx;
dy_b = -0.08;   % blue (c_fit)
ypos = dy_b;
hold on
plot([xpos-55, xpos-30], [ypos ypos], 'g', 'LineWidth', 2)

text(xpos, ypos, '\quad $c_2(t)$', ...
    'Interpreter','latex', 'FontSize',20, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',1);
dx   = -260;  % left shift

dy_b = -0.17;   % blue (c_fit)
ypos =  dy_b;
hold on
plot([xpos-60, xpos-35], [ypos ypos], 'm', 'LineWidth', 2)

text(xpos, ypos, '\quad $c^{\mathrm{h}}_2(t)$', ...
    'Interpreter','latex', 'FontSize',20, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',1);


xlabel('$t$')

exportgraphics(fig6a, '2vsap.pdf', 'ContentType', 'vector')

%%
fig6b = figure;
fig6b.Position = [300, 300, 600, 300];
plot(timegrid, (output2.llintest- cumsum([zeros(1,1) (1+output2.delexp(1:end-1)/24)])*dt-output2.lexp+timegrid*1.015)/(meta2.eps.^2), 'g','LineWidth', 1.2)
hold on
plot(timegrid, (output2.lhtest- cumsum([zeros(1,1) (1+output2.delh(1:end-1)/24)])*dt-output2.lexp+timegrid*1.015)/(meta2.eps.^2), 'm','LineWidth', 1.2)
%plot(timegrid,zeros(size(timegrid)),'k--')

%legend('$\overline{\gamma}(t)$ for $\overline{\eta}=\eta_1$', '$\overline{\gamma}(t)$ for $\overline{\eta}=\eta_{\mathrm{h}}$')
xlabel('$t$')

dx   = -200;  % left shift
xpos = 1100 + dx;
dy_b = -500;   % blue (c_fit)
ypos = output2.llintest(end) - timegrid(end)*1.015 + dy_b;
hold on
plot([xpos-60, xpos-35], [ypos ypos], 'g', 'LineWidth', 2)

text(xpos, ypos, '\quad ${\gamma}_{2}(t)$', ...
    'Interpreter','latex', 'FontSize',20, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',1);

dx   = -200;  % left shift
dy_b = -240;   % blue (c_fit)
ypos = output2.lhtest(end) - timegrid(end)*1.015 + dy_b;
hold on
plot([xpos-60, xpos-35], [ypos ypos], 'm', 'LineWidth', 2)

text(xpos, ypos, '\quad ${\gamma}_{2}^{\mathrm{h}}(t)$', ...
    'Interpreter','latex', 'FontSize',20, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',1);

exportgraphics(fig6b, '2vsappos.pdf', 'ContentType', 'vector')


%%
rate_inf=-2.23e-4+6.528e-5;
rate_inftrue = -1.403708886574337e-04;
%rate_inf2=-0.0258*(meta2.cstar-1);
t_num=floor(meta2.time/(meta2.FrameSkip*meta2.h))+1;
timegrid = (0:t_num-1) * meta2.FrameSkip * meta2.h;
timegriddot = linspace(0,1200,100);

fig8a = figure;
fig8a.Position = [300, 300, 600, 300];

plot(timegrid, output2.cfit, 'b','LineWidth', 1.2)
% hold on
% plot(timegrid, 1+output2.delap/24, 'm','LineWidth', 1.2)
hold on

plot(timegriddot, meta2.cstar+meta2.eps^2*rate_inf*timegriddot,'k.','LineWidth', 1.2)
plot(timegrid, meta2.cstar+meta2.eps^2*rate_inftrue*timegrid,'k--','LineWidth', 1.2)
%plot(timegrid,1.015*ones(size(timegrid)),'k--')
hold off
legend('$c(t)$')
xlabel('$t$')


exportgraphics(fig8a, 'rate.pdf', 'ContentType', 'vector')

%%

t_num=floor(meta4.time/(meta4.FrameSkip*meta4.h))+1;
timegrid = (0:t_num-1) * meta4.FrameSkip * meta4.h;
delstar=(meta4.cstar-1)/24;
N = floor(meta4.time*1.2);

fig9a = figure;
fig9a.Position = [300, 300, 600, 300];

% Plot in color
plot(timegrid, output4.cfit, 'b', 'LineWidth',1.2); hold on
plot(timegrid, 1 + output4.delmod/24, 'r', 'LineWidth',1.2);
%plot(timegrid, 1 + output2.delexp/24, 'g', 'LineWidth',1.2);
plot(timegrid, 1.015*ones(size(timegrid)), 'k--', 'LineWidth',1.0); % reference
hold off

xlabel('$t$', 'Interpreter','latex')

% --- Inline labels (boxed) ---
xlab = timegrid(end);



% Blue: c_fitxpos = xlab + dx;
dx   = -100;  % left shift
xpos = xlab + dx;
dy_b = 4e-04;   % blue (c_fit)
ypos = output4.cfit(end) + dy_b;
hold on
plot([xpos-65, xpos-40], [ypos ypos], 'b', 'LineWidth', 2.5)

text(xlab+dx, ypos, '\quad $c_{\mathrm{fit}}(t)$', ...
    'Interpreter','latex', 'FontSize',20, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',1);

%Red: c_mod
dx=-220;
dy_r = -1.9e-04;   % red  (c_mod) - further down
xpos = xlab + dx;
ypos = 1+output4.delmod(end)/24 + dy_r;
hold on
plot([xpos-80, xpos-55], [ypos ypos], 'r', 'LineWidth', 2)
text(xlab+dx, ypos, '\quad $c_{\mathrm{mod}}(t)$', ...
    'Interpreter','latex', 'FontSize',20, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',1);

% Green: c_1
% dx   = -80;  % left shift
% xpos = xlab + dx;
% dy_g = 2e-04;   % green (c_1)
% ypos = 1 + output2.delexp(end)/24 + dy_g;
% hold on
% plot([xpos-55, xpos-30], [ypos ypos], 'g', 'LineWidth', 2)
% text(xlab+dx, ypos, '\quad $c_{1}(t)$', ...
%     'Interpreter','latex', 'FontSize',16, ...
%     'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
%     'EdgeColor','k', 'Margin',3);
ylim([1.0138,1.0155])
% Export
exportgraphics(fig9a, 'modvsfit.pdf', 'ContentType', 'vector');

%%
lplot=output4.lfit- cumsum([zeros(1,1) output4.cfit(1:end-1)])*dt;
lplot(4)=lplot(3);

fig9b = figure;
fig9b.Position = [300, 300, 600, 300];
dt=timegrid(2)-timegrid(1);
% Plot in color
plot(timegrid, lplot, 'b', 'LineWidth',1.2); hold on
cmod=(1 + output4.delmod/24);
plot(timegrid, output4.lmod- cumsum([zeros(1,1) cmod(1:end-1)])*dt, 'r', 'LineWidth',1.2);

plot(timegrid,zeros(size(timegrid)),'k--');

xlabel('$t$', 'Interpreter','latex')

% --- Inline labels (boxed) ---
xlab = timegrid(end);

% Blue: c_fitxpos = xlab + dx;
dx   = -180;  % left shift
xpos = xlab + dx;
dy_b = 0;   % blue (c_fit)
ypos = output4.lfit(end)- timegrid(end)*1.015 + dy_b;
hold on
plot([xpos-70, xpos-45], [ypos ypos], 'b', 'LineWidth', 2)

text(xlab+dx, ypos, '\quad $\gamma_{\mathrm{fit}}(t)$', ...
    'Interpreter','latex', 'FontSize',20, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',1);

% % Red: c_mod
dx=-130;
dy_r = 1;   % red  (c_mod) - further down
xpos = xlab + dx;
ypos = output4.lmod(end)- timegrid(end)*1.015 + dy_r;
hold on
plot([xpos-90, xpos-65], [ypos ypos], 'r', 'LineWidth', 2)
text(xlab+dx, ypos, '\quad $\gamma_{\mathrm{mod}}(t)$', ...
    'Interpreter','latex', 'FontSize',20, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',1);

% Green: c_1
% dx   = -130;  % left shift
% xpos = xlab + dx;
% dy_g = 0.7;   % green (c_1)
% ypos = output2.lexp(end)- timegrid(end)*1.015 + dy_g;
% hold on
% plot([xpos-55, xpos-30], [ypos ypos], 'g', 'LineWidth', 2)
% text(xlab+dx, ypos, '\quad $\gamma_{1}(t)$', ...
%     'Interpreter','latex', 'FontSize',16, ...
%     'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
%     'EdgeColor','k', 'Margin',3);
ylim([-3,0.5])
exportgraphics(fig9b, 'modvsfitpos.pdf', 'ContentType', 'vector')

%%
figtest = figure;
figtest.Position = [300, 300, 600, 300];

% Plot in color
plot(timegrid, output4.cfit, 'b', 'LineWidth',1.2); hold on
plot(timegrid, 1+output4.del2/24, 'r', 'LineWidth',1.2); hold on
%plot(timegrid, 1+output4.delmod/24, 'r', 'LineWidth',1.2); hold on
%plot(timegrid, 1+output4.dellin/24, 'm', 'LineWidth',1.2); hold on
%legend('fit','2')
%plot(timegrid, 1 + output2.delmod/24, 'r', 'LineWidth',1.2);
%plot(timegrid, 1 + output2.delexp/24, 'g', 'LineWidth',1.2);
plot(timegrid, 1.015*ones(size(timegrid)), 'k--', 'LineWidth',1.0); % reference
hold off

xlabel('$t$', 'Interpreter','latex')


xlab = timegrid(end);

% Blue: c_fitxpos = xlab + dx;
dx   = -200;  % left shift
xpos = xlab + dx;
dy_b = -3.5e-04;   % blue (c_fit)
ypos = output4.cfit(end)+ dy_b;
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
ypos = output4.cfit(end)+ dy_b;
hold on
plot([xpos-55, xpos-30], [ypos ypos], 'r', 'LineWidth', 2)

text(xlab+dx, ypos, '\quad $c_2(t)$', ...
    'Interpreter','latex', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);

exportgraphics(figtest, 'c2.pdf', 'ContentType', 'vector')

%%

figtest = figure;
figtest.Position = [300, 300, 600, 300];

dt=timegrid(2)-timegrid(1);
% Plot in color
plot(timegrid, output4.lfit- cumsum([zeros(1,1) output4.cfit(1:end-1)])*dt, 'b', 'LineWidth',1.2); hold on
c2=(1 + output4.del2/24);
plot(timegrid, output4.l2- cumsum([zeros(1,1) c2(1:end-1)])*dt, 'r', 'LineWidth',1.2); hold on
% cmod=(1 + output4.delmod/24);
% plot(timegrid, output4.lmod- cumsum([zeros(1,1) cmod(1:end-1)])*dt, 'r', 'LineWidth',1.2); hold on
%plot(timegrid, output4.llin - timegrid*1.015, 'm', 'LineWidth',1.2); hold on
legend('fit','2')
%plot(timegrid, 1 + output2.delmod/24, 'r', 'LineWidth',1.2);
%plot(timegrid, 1 + output2.delexp/24, 'g', 'LineWidth',1.2);
%plot(timegrid, zeros(size(timegrid)), 'k--', 'LineWidth',1.0); % reference
hold off

xlabel('$t$', 'Interpreter','latex')

%%
figtest = figure;
figtest.Position = [300, 300, 600, 300];

% Plot in color
plot(timegrid, output4.cfit+0.015-output4.del0/24, 'b', 'LineWidth',1.2); hold on
% plot(timegrid, 1+output4.del2/24, 'g', 'LineWidth',1.2); hold on
% plot(timegrid, 1+output4.dellin/24, 'r', 'LineWidth',1.2); hold on
%plot(timegrid, 1+output4.del2/24, 'r', 'LineWidth',1.2); hold on
%plot(timegrid, 1+output4.delmod/24, 'r', 'LineWidth',1.2); hold on
%plot(timegrid, 1+output4.dellin/24, 'm', 'LineWidth',1.2); hold on
%legend('fit','2')
%plot(timegrid, 1 + output2.delmod/24, 'r', 'LineWidth',1.2);
%plot(timegrid, 1 + output2.delexp/24, 'g', 'LineWidth',1.2);
%plot(timegrid, 1.015*ones(size(timegrid)), 'k--', 'LineWidth',1.0); % reference
hold off

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
plot(timegrid035*(0.35)^2,mean(ensemble035.cfit_ens), 'k', 'LineWidth',1.2);
hold on
plot(timegrid03*(0.3)^2,mean(ensemble03.cfit_ens), 'm', 'LineWidth',1.2);
hold on
plot(timegrid*(0.1)^2,mean(ensemble02.cfit_ens), 'r', 'LineWidth',1.2);
hold on
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

% plot(timegrid*(0.1)^2,metaens01.cstar+(0.1)^2*rate_inftrue*timegrid,'k--','LineWidth', 1.2)
% plot(timegrid*(0.1)^2,1.006-(0.1)^2*rate_006true*timegrid(idx006)+(0.1)^2*rate_006true*timegrid,'k--','LineWidth', 1.2)
% plot(timegrid*(0.1)^2,1.008-(0.1)^2*rate_008true*timegrid(idx008)+(0.1)^2*rate_008true*timegrid,'k--','LineWidth', 1.2)
% 
% plot((0.1)^2*timegrid(idx006), mean(ensemble01.cfit_ens(:,idx006)), 'bs', 'MarkerSize',8, 'MarkerFaceColor','b');
% plot((0.1)^2*timegrid(idx008), mean(ensemble01.cfit_ens(:,idx008)), 'bs', 'MarkerSize',8, 'MarkerFaceColor','b');
% plot(0, mean(ensemble01.cfit_ens(:,1)), 'bs', 'MarkerSize',8, 'MarkerFaceColor','b');


%legend('$\mathbf{E}[c(t)]$')
legend('$\sigma=0.35$','$\sigma=0.3$','$\sigma=0.2$','$\sigma=0.1$')
ylim([1, metaens01.cstar])
xlim([0, 195])
xlabel('$\sigma^2 t$')
ylabel('$\mathbf{E}[c(t)]$','Rotation',0)
exportgraphics(figmean, 'mean.pdf', 'ContentType', 'vector');


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
plot((0.1)^2*timegrid, ceulercor,'k--', 'LineWidth', 1.5);
hold on
plot((0.1)^2*timegrid, ceuler,'k.', 'LineWidth', 1.5);

legend('$\mathbf{E}[c(t)]$','$c_{\mathrm{lim}}(\sigma^2 t)$','$c^{\mathrm{h}}_{\mathrm{lim}}(\sigma^2 t)$')
ylim([1, metaens01.cstar])
xlim([0, 180])
xlabel('$\sigma^2 t$')
%ylabel('$\mathbf{E}[c(t)]$','Rotation',0)
exportgraphics(figmean, 'ODE.pdf', 'ContentType', 'vector');

%%
t_num=floor(metaHJ.time/(metaHJ.FrameSkip*metaHJ.h))+1;
timegrid = (0:t_num-1) * metaHJ.FrameSkip * metaHJ.h;
figtest = figure;
figtest.Position = [300, 300, 600, 300];
plot(timegrid, outputHJ.cfit , 'b','LineWidth', 1.2)
hold on
plot(timegrid, 1+outputHJ.delhtest/24, 'g','LineWidth', 1.2)
hold on
plot(timegrid, 1+outputHJ.dellintest/24, 'm','LineWidth', 1.2)

%%
corrections=zeros(1,10);

[metaHJ, outputHJ] = load_simulation_output('simulation_forHJ.txt');
t_num=floor(metaHJ.time/(metaHJ.FrameSkip*metaHJ.h))+1;
timegrid = (0:t_num-1) * metaHJ.FrameSkip * metaHJ.h;
correctioncurve015=(outputHJ.dellintest-outputHJ.delexp-outputHJ.delhtest+outputHJ.delexp)./(24*metaHJ.eps.^2*timegrid);
corrections(end)=correctioncurve015(end);
fig6a = figure;
fig6a.Position = [300, 300, 600, 300];
plot(timegrid, correctioncurve015, 'g','LineWidth', 1.2)



[metaode, outputode] = load_simulation_output('simulation_for_ode006.txt');
t_num=floor(metaode.time/(metaode.FrameSkip*metaode.h))+1;
timegrid = (0:t_num-1) * metaode.FrameSkip * metaode.h;
correctioncurve006=(outputode.dellintest-outputode.delexp-outputode.delhtest+outputode.delexp)./(24*metaode.eps.^2*timegrid);
corrections(1)=correctioncurve006(end);
fig6a = figure;
fig6a.Position = [300, 300, 600, 300];
plot(timegrid, correctioncurve006, 'g','LineWidth', 1.2)

[metaode, outputode] = load_simulation_output('simulation_for_ode007.txt');
t_num=floor(metaode.time/(metaode.FrameSkip*metaode.h))+1;
timegrid = (0:t_num-1) * metaode.FrameSkip * metaode.h;
correctioncurve007=(outputode.dellintest-outputode.delexp-outputode.delhtest+outputode.delexp)./(24*metaode.eps.^2*timegrid);
corrections(2)=correctioncurve007(end);
fig6a = figure;
fig6a.Position = [300, 300, 600, 300];
plot(timegrid, correctioncurve007, 'g','LineWidth', 1.2)

[metaode, outputode] = load_simulation_output('simulation_for_ode008.txt');
t_num=floor(metaode.time/(metaode.FrameSkip*metaode.h))+1;
timegrid = (0:t_num-1) * metaode.FrameSkip * metaode.h;
correctioncurve008=(outputode.dellintest-outputode.delexp-outputode.delhtest+outputode.delexp)./(24*metaode.eps.^2*timegrid);
corrections(3)=correctioncurve008(end);
fig6a = figure;
fig6a.Position = [300, 300, 600, 300];
plot(timegrid, correctioncurve008, 'g','LineWidth', 1.2)

[metaode, outputode] = load_simulation_output('simulation_for_ode009.txt');
t_num=floor(metaode.time/(metaode.FrameSkip*metaode.h))+1;
timegrid = (0:t_num-1) * metaode.FrameSkip * metaode.h;
correctioncurve009=(outputode.dellintest-outputode.delexp-outputode.delhtest+outputode.delexp)./(24*metaode.eps.^2*timegrid);
corrections(4)=correctioncurve009(end);
fig6a = figure;
fig6a.Position = [300, 300, 600, 300];
plot(timegrid, correctioncurve009, 'g','LineWidth', 1.2)

[metaode, outputode] = load_simulation_output('simulation_for_ode010.txt');
t_num=floor(metaode.time/(metaode.FrameSkip*metaode.h))+1;
timegrid = (0:t_num-1) * metaode.FrameSkip * metaode.h;
correctioncurve010=(outputode.dellintest-outputode.delexp-outputode.delhtest+outputode.delexp)./(24*metaode.eps.^2*timegrid);
corrections(5)=correctioncurve010(end);
fig6a = figure;
fig6a.Position = [300, 300, 600, 300];
plot(timegrid, correctioncurve010, 'g','LineWidth', 1.2)

[metaode, outputode] = load_simulation_output('simulation_for_ode011.txt');
t_num=floor(metaode.time/(metaode.FrameSkip*metaode.h))+1;
timegrid = (0:t_num-1) * metaode.FrameSkip * metaode.h;
correctioncurve011=(outputode.dellintest-outputode.delexp-outputode.delhtest+outputode.delexp)./(24*metaode.eps.^2*timegrid);
corrections(6)=correctioncurve011(end);
fig6a = figure;
fig6a.Position = [300, 300, 600, 300];
plot(timegrid, correctioncurve011, 'g','LineWidth', 1.2)

[metaode, outputode] = load_simulation_output('simulation_for_ode012.txt');
t_num=floor(metaode.time/(metaode.FrameSkip*metaode.h))+1;
timegrid = (0:t_num-1) * metaode.FrameSkip * metaode.h;
correctioncurve012=(outputode.dellintest-outputode.delexp-outputode.delhtest+outputode.delexp)./(24*metaode.eps.^2*timegrid);
corrections(7)=correctioncurve012(end);
fig6a = figure;
fig6a.Position = [300, 300, 600, 300];
plot(timegrid, correctioncurve012, 'g','LineWidth', 1.2)

[metaode, outputode] = load_simulation_output('simulation_for_ode013.txt');
t_num=floor(metaode.time/(metaode.FrameSkip*metaode.h))+1;
timegrid = (0:t_num-1) * metaode.FrameSkip * metaode.h;
correctioncurve013=(outputode.dellintest-outputode.delexp-outputode.delhtest+outputode.delexp)./(24*metaode.eps.^2*timegrid);
corrections(8)=correctioncurve013(end);
fig6a = figure;
fig6a.Position = [300, 300, 600, 300];
plot(timegrid, correctioncurve013, 'g','LineWidth', 1.2)

[metaode, outputode] = load_simulation_output('simulation_for_ode014.txt');
t_num=floor(metaode.time/(metaode.FrameSkip*metaode.h))+1;
timegrid = (0:t_num-1) * metaode.FrameSkip * metaode.h;
correctioncurve014=(outputode.dellintest-outputode.delexp-outputode.delhtest+outputode.delexp)./(24*metaode.eps.^2*timegrid);
corrections(9)=correctioncurve014(end);
fig6a = figure;
fig6a.Position = [300, 300, 600, 300];
plot(timegrid, correctioncurve014, 'g','LineWidth', 1.2)

save('corrections.mat','corrections');
%% Video

tic
Visualization(db, r, p, kappa, cstar, delstar,time, FrameSkip, h, eta1fit, eta2fit, eta1mod, eta2mod, eta1modapprox, eta2modapprox, eta1free, eta2free, cfit,lfit, cmod,c2, lmod,llin);
toc
