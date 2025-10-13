clear all

parentFolder = fileparts(pwd);  % Get parent directory of current folder
addpath(parentFolder);         % Add it to the path

% Load and preprocess waveprofiles
filename = 'all_sols_out_e.txt';
[c_groups, c_values] = load_phi_data(filename);
db = preprocess_phi_data(c_groups, c_values);
%%
N = 80;              % lattice half-width, so j in -N:N
n = 2*N + 1;
j_vals = (-N:N)';
dt = 0.1;
T_end = 60;
num_steps = round(T_end / dt);
t_vals = linspace(0, T_end, num_steps);
sigma = 0.07;
pvals = 0:0.1:(1-0.1);

load('Rinfstored.mat', 'Rinf');
load('ratessaved.mat', 'rates');
load('corrections.mat', 'corrections')
%%
Sinf1_vals = zeros(length(pvals), n);
Sinf2_vals = zeros(length(pvals), n);

for p_idx=1:length(pvals)
    for j_idx = 1:n
        j = j_vals(j_idx);
        Rinf1_slice = squeeze(Rinf(p_idx,1, :,j_idx));  % vector over m
        Rinf2_slice = squeeze(Rinf(p_idx,2, :,j_idx));  % vector over m
        Sinf1 = sqrt(sum(Rinf1_slice.^2));
        Sinf2 = sqrt(sum(Rinf2_slice.^2));
        Sinf1_vals(p_idx,j_idx) = Sinf1;
        Sinf2_vals(p_idx,j_idx) = Sinf2;
    end
end




fig7b = figure;
fig7b.Position = [300, 300, 600, 300];
for p_idx=1:length(pvals)
    p=pvals(p_idx);
    plot(j_vals-p, Sinf1_vals(p_idx,:), 'b.', 'MarkerSize',10, 'MarkerFaceColor','b');
    hold on
    plot(j_vals-p, Sinf2_vals(p_idx,:), 'r.', 'MarkerSize',10, 'MarkerFaceColor','r');
end

xlab = j_vals(end);


hold on

xpos = 9;
ypos =0.03;
plot([xpos-5, xpos-4], [ypos ypos], 'b', 'LineWidth', 2)


text(xpos, ypos, '$\quad r$-component', ...
    'Interpreter','latex', 'FontSize',14, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);

xpos = -7;
ypos =0.005;
plot([xpos-5, xpos-4], [ypos ypos], 'r', 'LineWidth', 2)


text(xpos, ypos, '$\quad p$-component', ...
    'Interpreter','latex', 'FontSize',14, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);

%legend('$r$-component','$p$-component')
xlabel('$x$');
xlim([-35,15])
grid on;
exportgraphics(fig7b, 'sqrtvar.pdf', 'ContentType', 'vector');

%%
% Final summation against kappa
kappa = 2*sqrt(3)*(rand(1, 2*N+1) - 0.5);
eta_inf = zeros(length(pvals), 2, n);


    for idx_j = 1:n
        Rinf_j = reshape(Rinf(1,:, :, idx_j), [2, n]);  % ensure 2 × n shape
        eta_inf(1,:, idx_j) = Rinf_j * kappa.';        % result is 2 × 1
    end


fig7a = figure;
fig7a.Position = [300, 300, 600, 300];
    
    plot(j_vals, squeeze(eta_inf(1,1,:)), 'bs', 'MarkerSize',4, 'MarkerFaceColor','b');
    hold on
    plot(j_vals , squeeze(eta_inf(1,2,:)), 'r^', 'MarkerSize',4, 'MarkerFaceColor','r');
    legend('$r$-component','$p$-component')

%title('eta\_inf')
xlim([-35,15])
xlabel('j')
grid on
exportgraphics(fig7a, 'etarealisation.pdf', 'ContentType', 'vector');

%%
% Compute rate_inf
cdotinf=zeros(1,length(pvals));
c0=1.015;

for p_idx = 1:length(pvals)
    p=pvals(p_idx);
    D_inf = 0;
    LIIL_inf = 0;
    for m = 1:n
        Rinf1m = squeeze(Rinf(p_idx,1, m, :))';
        Rinf2m = squeeze(Rinf(p_idx,2, m, :))';
        delta_m = zeros(1, n); delta_m(m) = 1;
    
        D_inf = D_inf + D(db, c0, p/c0, j_vals, Rinf1m + delta_m, Rinf1m);
        LIIL_inf = LIIL_inf + LII(db, c0, p/c0, j_vals, Rinf1m, Rinf2m) * LI(db, c0, p/c0, j_vals, delta_m);
    end
    rate_inf =  D_inf + LIIL_inf;
    cdotinf(p_idx) = rate_inf(2)/24+6.528e-5;
end

fig9 = figure;
fig9.Position = [300, 300, 600, 300];
plot(pvals,cdotinf, 'kx', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');
hold on
plot(pvals,-1.403708886574337e-4, 'ko', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');
xlabel('$p$')
ylim([-2e-4, 0e-4])
exportgraphics(fig9, 'pdependence.pdf', 'ContentType', 'vector')
%set(get(gca,'ylabel'),'rotation',0)


%%
% Compute rate_inf
dc=0.001;
c_vals = 1.005:dc:1.015;  % your range of c0 values

rate_results=zeros(1,length(c_vals));

for c_idx = 1:length(c_vals)
    c0=c_vals(c_idx);
    % Construct the filename based on c0
    c_str = num2str(round((c0-1)*1000)); % e.g., 1.015 -> '015'
    if c_idx > 5
        filename = ['R_RK0' c_str '.mat'];
    elseif c_idx < 6
        filename = ['R_RK00' c_str '.mat'];
    end
    % Load R_RK from the file
    load(filename, 'R_RK');
    
    % Initialize sums
    D_inf = 0;
    LIIL_inf = 0;
    
    for m = 1:n
        Rinf1m = squeeze(R_RK(1, m, :, end))';
        Rinf2m = squeeze(R_RK(2, m, :, end))';
        delta_m = zeros(1, n); delta_m(m) = 1;
        
        D_inf = D_inf + D(db, c0, T_end, j_vals, Rinf1m + delta_m, Rinf1m);
        LIIL_inf = LIIL_inf + LII(db, c0, T_end, j_vals, Rinf1m, Rinf2m) ...
                             * LI(db, c0, T_end, j_vals, delta_m);
    end
    
    % Compute rate
    rate_inf = (D_inf(2) + LIIL_inf(2))/24 + compute_rate_correction(db, c0, 100, 0.05);
    
    % Optional: store rate_inf for each c0
    rate_results(c_idx) = rate_inf;
end

%%
dc=0.001;
c_vals=1.005:dc:1.015;


x = log(c_vals-1);
y = log(-rates);
p = polyfit(x, y, 1);      % p(1) = alpha, p(2) = log(k)

alpha = p(1);
k = exp(p(2));

xcor = log(c_vals(:)-1);
ycor = log(-rate_results);
pcor = polyfit(xcor, ycor, 1);      % p(1) = alpha, p(2) = log(k)

alphacor = pcor(1);
kcor = exp(pcor(2));

fig8b = figure;
fig8b.Position = [300, 300, 600, 300];
plot(c_vals, rates(:), 'kx', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');
hold on
plot(c_vals, -k*((c_vals-1).^alpha), '-b','LineWidth',1.2)
hold on
plot(c_vals(:), rate_results, 'ko', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');
hold on
plot(c_vals, -kcor*((c_vals-1).^alphacor), '-r','LineWidth',1.2)
xlabel('$c$')
xlim([1.005, 1.015])
set(get(gca,'ylabel'),'rotation',0)
ylabel('$\overline{\mathcal{Q}}_c(c)$')


% Blue: c_fitxpos = xlab + dx;

xpos=1.012;
ypos=-0.5e-4;
hold on
plot([xpos-0.0013, xpos-0.0010], [ypos ypos], 'r', 'LineWidth', 2)

text(xpos, ypos, '\qquad $-0.35(c-1)^{1.86}$', ...
    'Interpreter','latex', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);

xpos=1.010;
ypos=-1.2e-4;
hold on
plot([xpos-0.0013, xpos-0.0010], [ypos ypos], 'b', 'LineWidth', 2)

text(xpos, ypos, '\qquad $-0.54(c-1)^{1.93}$', ...
    'Interpreter','latex', 'FontSize',16, ...
    'VerticalAlignment','middle', 'HorizontalAlignment','center', ...
    'EdgeColor','k', 'Margin',3);

%ylabel('Amplitude parameter')
exportgraphics(fig8b, 'rates.pdf', 'ContentType', 'vector')

