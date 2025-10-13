
clear; clc;
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
c0 = 1.005;         % wave speed
del0 = 24*(c0-1);
kappa = 2*sqrt(3)*(rand(1, 2*N+1) - 0.5);
sigma = 0.07;

% Use LaTeX fonts and settings globally
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultLegendInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(groot, 'defaultAxesFontName', 'Times New Roman')  % optional LaTeX-style font
set(groot, 'defaultAxesFontSize', 15)


%%
dc=0.001;
c_vals=1.005:dc:1.015;
rates = zeros(1,length(c_vals));
for c_idx = 1:length(c_vals)
    rate_inf = compute_rate_inf(db, c_vals(c_idx), N, dt,30);
    disp(c_idx)
    rates(c_idx)=rate_inf(2)/24+ compute_rate_correction(db, c_vals(c_idx), 100, 0.05);
    
end
save('ratessaved.mat', 'rates');

slope = (rates(end) - rates(1)) / (c_vals(end) - c_vals(1));
intercept = rates(1) - slope * c_vals(1);  % y = m*x + b ⇒ b = y - m*x

% Create the line values over c_vals
line_vals = slope * (c_vals + intercept/slope);

x = log(c_vals-1);
y = log(-rates);
p = polyfit(x, y, 1);      % p(1) = alpha, p(2) = log(k)

alpha = p(1);
k = exp(p(2));

plotcs=1:0.001:max(c_vals);


% First figure
fig1 = figure;
fig1.Position = [300, 300, 600, 300];
plot(c_vals, rates(:), 'kx', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');
% hold on
% plot(c_vals, line_vals, '-r')
hold on
plot(plotcs, -k*((plotcs-1).^alpha), '-r')
%legend('$\sigma=0.035$','$\sigma=0.07$','Location','southwest')
xlabel('$c$')
xlim([1.005, 1.015])
set(get(gca,'ylabel'),'rotation',0)
ylabel('$\mathcal{D}^c(c)$')
%ylabel('Amplitude parameter')
exportgraphics(fig1, 'rates.pdf', 'ContentType', 'vector')




%% Computes Rinf


% Assumes j_vals is n×1, dt is timestep, and db, del0 are defined
tau_max = T_end;
tau_vals = 0:dt:tau_max;
tau_steps = length(tau_vals);
pvals = 0:0.1:(1-0.1);
%pvals=[0];
% Initialize output: 2 × n × n
Rinf = zeros(length(pvals),2, n, n);
tic
for p_idx=1:length(pvals)
    p=pvals(p_idx);
    disp(p_idx)

    % Precompute Tdelta for s = -tau
    Tdelta_all = zeros(2, n, n, tau_steps);  % 2 × n × n × num_steps
    for m_tau = 1:tau_steps
        tau = tau_vals(m_tau);
        Tdelta_all(:, :, :, m_tau) = Tdelta(db, del0, j_vals, p/c0-tau);
    end

    for idx_j = 1:n
        j = j_vals(idx_j);
        
    
        for m_tau = 1:tau_steps
            tau = tau_vals(m_tau);
            T_s = Tdelta_all(:, :, :, m_tau);  % 2 × n × n
    
            % Compute Phi(j−k, τ): k ∈ j_vals
            shifts = j - j_vals';  % 1 × n
            arg = 2 * tau;
            j0 = besselj(2 * shifts,   arg);
            jp = besselj(2 * shifts+1, arg);
            jm = besselj(2 * shifts-1, arg);
    
            % Phi_vals: 2 × 2 × n
            Phi_vals = zeros(2, 2, n);
            Phi_vals(1,1,:) = j0;
            Phi_vals(1,2,:) = -jp;
            Phi_vals(2,1,:) = -jm;
            Phi_vals(2,2,:) = j0;
    
            % Reshape for batched multiplication: 2 × 1 × n × n
            T_s_reshaped = reshape(T_s, 2, 1, n, n);  % (2 × 1 × n × n)
    
            % Multiply Phi with T_s over k ∈ j_vals
            PT = pagemtimes(Phi_vals, T_s_reshaped);  % 2 × 1 × n × n
    
            % Sum over k (dim 3): result is 2 × 1 × n
            sum_k = sum(PT, 3);  % 2 × 1 × n
    
            % Reshape to 2 × n
            sum_k = reshape(sum_k, 2, n);
    
            % Accumulate into Rinf(:, m, j)
            Rinf(p_idx,:, :, idx_j) = reshape(Rinf(p_idx,:,:,idx_j), size(sum_k)) + dt * sum_k;
        end
    end
end
toc

save('Rinfstored.mat', 'Rinf');

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


fig1 = figure;
fig1.Position = [300, 300, 600, 300];
for p_idx=1:length(pvals)
    p=pvals(p_idx);
    plot(j_vals-p, Sinf1_vals(p_idx,:), 'b.', 'MarkerSize', 10);
    hold on
    plot(j_vals-p, Sinf2_vals(p_idx,:), 'r.', 'MarkerSize', 10);
end
legend('$r$-component','$p$-component')
xlabel('$x$');
xlim([-35,15])
grid on;
exportgraphics(fig1, 'sqrtvar.pdf', 'ContentType', 'vector');

%%

% Final summation against kappa
eta_inf = zeros(length(pvals), 2, n);


    for idx_j = 1:n
        Rinf_j = reshape(Rinf(1,:, :, idx_j), [2, n]);  % ensure 2 × n shape
        eta_inf(1,:, idx_j) = Rinf_j * kappa.';        % result is 2 × 1
    end


fig2 = figure;
fig2.Position = [300, 300, 600, 300];
    p = pvals(p_idx);
    plot(j_vals, squeeze(eta_inf(1,1,:)), 'bs', 'MarkerSize',4, 'MarkerFaceColor','b');
    hold on
    plot(j_vals , squeeze(eta_inf(1,2,:)), 'r^', 'MarkerSize',4, 'MarkerFaceColor','r');
    legend('$r$-component','$p$-component')

%title('eta\_inf')
xlim([-35,15])
xlabel('j')
grid on
exportgraphics(fig2, 'etarealisation.pdf', 'ContentType', 'vector');



%% Rate


D_inf = 0;
D2_inf = 0;
LIIL_inf = 0;

for m = 1:n
    Rinf1m = squeeze(Rinf(1, m, :))';
    Rinf2m = squeeze(Rinf(2, m, :))';
    delta_m = zeros(1, n); delta_m(m) = 1;
    D_inf = D_inf + D(db, c0, 0, j_vals, Rinf1m+delta_m, Rinf1m);
    LIIL_inf = LIIL_inf + LII(db, c0, 0, j_vals, Rinf1m, Rinf2m) * LI(db, c0, 0, j_vals, delta_m);
end
rate_inf=sigma^2  * (D_inf  + LIIL_inf);


%% Random Forced Equation
tic



% Initialize eta (2 x n), zero initial condition
    eta = zeros(2, n);
    
    % Preallocate to save eta at each time step (optional, memory-heavy for large T_end)
    eta_t = zeros(2, n, num_steps);
    eta_t(:, :, 1) = eta;
    
    forcing = @(j_vals, t) Forcing(db, del0, j_vals, t, kappa);
    
    
    % RK4 integration loop
    for m = 2:num_steps
        disp(m)
        t = t_vals(m-1);
        [r1, r2] = J(eta(1,:), eta(2,:));

       
        k1 = [r1; r2] + forcing(j_vals, t);
        [r1, r2] = J(eta(1,:) + 0.5*dt*k1(1,:), eta(2,:) + 0.5*dt*k1(2,:));
        k2 = [r1; r2] + forcing(j_vals, t + 0.5*dt);
        
        [r1, r2] = J(eta(1,:) + 0.5*dt*k2(1,:), eta(2,:) + 0.5*dt*k2(2,:));
        k3 = [r1; r2] + forcing(j_vals, t + 0.5*dt);
        
        [r1, r2] = J(eta(1,:) + dt*k3(1,:), eta(2,:) + dt*k3(2,:));
        k4 = [r1; r2] + forcing(j_vals, t + dt);
        
        eta = eta + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
        eta_t(:, :, m) = eta;
    end
toc

first_order_all  = zeros(2, num_steps);    % size 2 x T
second_order_all = zeros(2, num_steps);
N = (n - 1)/2;
[alpha0, alpha1] = alpha(db, del0, N);
A = [0, alpha0; -alpha0, alpha1];

for t_idx = 1:num_steps
    disp(t_idx)
    t = t_vals(t_idx);
    eta1 = eta_t(1,:,t_idx);  % 1 × 101
    eta2 = eta_t(2,:,t_idx);  % 1 × 101

    % Compute LI
    LI_m = LI(db, c0, t, j_vals, kappa);             % 2×1
    first_order_all(:, t_idx) = sigma * LI_m;

    % Compute D1 and D2
    D1_m = D(db, c0, t, j_vals, eta1, eta1);           % 2×1
    D2_m = D(db, c0, t, j_vals, kappa, eta1);          % 2×1

    % Compute LII
    LII_m = LII(db, c0, t, j_vals, eta1, eta2);      % 2×2

    % Total second order term
    second_order_all(:, t_idx) = sigma^2 * (D1_m + LII_m * LI_m + D2_m);
    
end

figure
plot(t_vals, c0+cumsum(first_order_all(2,:)+second_order_all(2,:))*dt/24,'b')
hold on
plot(t_vals, c0+rate_inf(2)*t_vals/24,'k')
legend('c(t)','average_inf')
title('amplitude attenuation')

figure
plot(t_vals, cumsum(first_order_all(1,:)+second_order_all(1,:))*dt,'b')
hold on
plot(t_vals, rate_inf(1)*t_vals,'k')

legend('gamma(t)','average_inf')
title('phase slowdown')
%% Fundamental Solutions
% Initial conditions (Kronecker delta at j=0)
r1 = zeros(n, 1); r1(N+1) = 1;   % delta at j=0 in r
p1 = zeros(n, 1);

r2 = zeros(n, 1);
p2 = zeros(n, 1); p2(N+1) = 1;   % delta at j=0 in p

% Storage for solutions
U1 = zeros(n, num_steps); V1 = zeros(n, num_steps);
U2 = zeros(n, num_steps); V2 = zeros(n, num_steps);

U1(:,1) = r1; V1(:,1) = p1;
U2(:,1) = r2; V2(:,1) = p2;

% Runge-Kutta 4 time stepping
for i = 2:num_steps
    t = t_vals(i-1);

    % --- First delta init (in r) ---
    [k1_r, k1_p] = J(r1, p1);
    [k2_r, k2_p] = J(r1 + 0.5*dt*k1_r, p1 + 0.5*dt*k1_p);
    [k3_r, k3_p] = J(r1 + 0.5*dt*k2_r, p1 + 0.5*dt*k2_p);
    [k4_r, k4_p] = J(r1 + dt*k3_r,     p1 + dt*k3_p);

    r1 = r1 + (dt/6)*(k1_r + 2*k2_r + 2*k3_r + k4_r);
    p1 = p1 + (dt/6)*(k1_p + 2*k2_p + 2*k3_p + k4_p);

    U1(:,i) = r1;
    V1(:,i) = p1;

    % --- Second delta init (in p) ---
    [k1_r, k1_p] = J(r2, p2);
    [k2_r, k2_p] = J(r2 + 0.5*dt*k1_r, p2 + 0.5*dt*k1_p);
    [k3_r, k3_p] = J(r2 + 0.5*dt*k2_r, p2 + 0.5*dt*k2_p);
    [k4_r, k4_p] = J(r2 + dt*k3_r,     p2 + dt*k3_p);

    r2 = r2 + (dt/6)*(k1_r + 2*k2_r + 2*k3_r + k4_r);
    p2 = p2 + (dt/6)*(k1_p + 2*k2_p + 2*k3_p + k4_p);

    U2(:,i) = r2;
    V2(:,i) = p2;
end

% Compare with Bessel solutions at final time
j = j_vals;

% For initial condition delta in r
u1_exact = besselj(2*j, 2*T_end);
v1_exact = -besselj(2*j - 1, 2*T_end);

% For initial condition delta in p
u2_exact = -besselj(2*j + 1, 2*T_end);
v2_exact =  besselj(2*j, 2*T_end);

% Plot
figure;
sgtitle('Fundamental Solutions, t='+string(T_end))
subplot(2,2,1);
plot(j, U1(:,end), 'b.', j, u1_exact, 'r.');
title('u(t) from r(0)=\delta');
legend('Numerical', 'Exact');

subplot(2,2,2);
plot(j, V1(:,end), 'b.', j, v1_exact, 'r.');
title('v(t) from r(0)=\delta');
legend('Numerical', 'Exact');

subplot(2,2,3);
plot(j, U2(:,end), 'b.', j, u2_exact, 'r.');
title('u(t) from p(0)=\delta');
legend('Numerical', 'Exact');

subplot(2,2,4);
plot(j, V2(:,end), 'b.', j, v2_exact, 'r.');
title('v(t) from p(0)=\delta');
legend('Numerical', 'Exact');

%% General initial conditions
% sech^2-type initial condition
tic
u0 = sech(0.3 * j_vals).^2;
v0 = zeros(n, 1);


% Numerical integration using RK4
r = u0;
p = v0;

U = zeros(n, num_steps);
V = zeros(n, num_steps);
U(:,1) = r;
V(:,1) = p;

for i = 2:num_steps
    [k1_r, k1_p] = J(r, p);
    [k2_r, k2_p] = J(r + 0.5*dt*k1_r, p + 0.5*dt*k1_p);
    [k3_r, k3_p] = J(r + 0.5*dt*k2_r, p + 0.5*dt*k2_p);
    [k4_r, k4_p] = J(r + dt*k3_r,     p + dt*k3_p);

    r = r + (dt/6)*(k1_r + 2*k2_r + 2*k3_r + k4_r);
    p = p + (dt/6)*(k1_p + 2*k2_p + 2*k3_p + k4_p);

    U(:,i) = r;
    V(:,i) = p;
end

% Construct exact solution via convolution with Green's function
u_exact = zeros(n,1);
v_exact = zeros(n,1);

for idx = 1:n
    j0 = j_vals(idx);
    for kidx = 1:n
        k = j_vals(kidx);
        Phi_jk = Phi(j0 - k, T_end);  % Uses your function
        uk_vec = [u0(kidx); v0(kidx)];
        uv = Phi_jk * uk_vec;

        u_exact(idx) = u_exact(idx) + uv(1);
        v_exact(idx) = v_exact(idx) + uv(2);
    end
end
toc
% Plot
figure;
sgtitle('General Solutions')
subplot(2,1,1);
plot(j, U(:,end), 'b.', j, u_exact, 'r.');
title('u(t) from sech^2 initial condition');
legend('Numerical', 'Exact');

subplot(2,1,2);
plot(j, V(:,end), 'b.', j, v_exact, 'r.');
title('v(t) from sech^2 initial condition');
legend('Numerical', 'Exact');



%% R via RK
tic
% Initialize eta (2 x n), zero initial condition
    R_RK = zeros(2, n, n, num_steps);
    function LR = L_t(db,R, t, c0,j_vals)
        R1 = R(1,:); 
        R2 = R(2,:);
    
        prof = profile(db, 'r', c0, c0*t, j_vals); 
        prof = prof(:)';   % force row orientation
    
        R1_new = R1 + 2 * (prof .* R1);
        R2_new = R2;
    
        LR = [R1_new; R2_new];
    end

    % Preallocate to save eta at each time step (optional, memory-heavy for large T_end)
    
    % RK4 integration loop
    for m_idx = 1:n
        delta_m = zeros(1, n); 
        delta_m(m_idx) = 1;
        disp(m_idx)
        forcing_m = @(j_vals, t) Forcing(db, del0, j_vals, t, delta_m);
        %rhs = @(R, t) J( L_t(db,R, t, c0, j_vals) ) + forcing_m(j_vals, t);
        R_temp = zeros(2, n);
        for t_idx = 2:num_steps
            t = t_vals(t_idx-1);
        
            % --- k1 ---
            L_R = L_t(db, R_temp, t, c0, j_vals);   % apply L_t
            [r1, r2] = J(L_R(1,:), L_R(2,:));       % then J
            k1 = [r1; r2] + forcing_m(j_vals, t);
        
            % --- k2 ---
            L_R = L_t(db, R_temp + 0.5*dt*k1, t + 0.5*dt, c0, j_vals);
            [r1, r2] = J(L_R(1,:), L_R(2,:));
            k2 = [r1; r2] + forcing_m(j_vals, t + 0.5*dt);
        
            % --- k3 ---
            L_R = L_t(db, R_temp + 0.5*dt*k2, t + 0.5*dt, c0, j_vals);
            [r1, r2] = J(L_R(1,:), L_R(2,:));
            k3 = [r1; r2] + forcing_m(j_vals, t + 0.5*dt);
        
            % --- k4 ---
            L_R = L_t(db, R_temp + dt*k3, t + dt, c0, j_vals);
            [r1, r2] = J(L_R(1,:), L_R(2,:));
            k4 = [r1; r2] + forcing_m(j_vals, t + dt);
        
            % --- update solution ---
            R_temp = R_temp + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
        
            R_RK(:, m_idx, :, t_idx) = R_temp;
        end
      
    end
toc

save('R_RK005.mat', 'R_RK');


SRK1_vals = zeros(1, n);
SRK2_vals = zeros(1, n);

for j_idx = 1:n
    j = j_vals(j_idx);
    RRK1_slice = squeeze(R_RK(1, :,j_idx,end));  % vector over m
    RRK2_slice = squeeze(R_RK(2, :,j_idx,end));  % vector over m
    SRK1 = sqrt(sum(RRK1_slice.^2));
    SRK2 = sqrt(sum(RRK2_slice.^2));
    SRK1_vals(j_idx) = SRK1;
    SRK2_vals(j_idx) = SRK2;
end
figure;
plot(j_vals, squeeze(R_RK(1, 65,:,50)), 'b.');
hold on
plot(j_vals, squeeze(R_RK(2, 65,:,50)), 'r.');
xlabel('j');

% Plot
figure;
plot(j_vals, SRK1_vals, 'b.');
hold on
plot(j_vals, SRK2_vals, 'r.');
xlabel('j');
title('\surd{\Sigma_m Rinf1(j, m, T)^2}');
grid on;



%% Exact solution


tic

eta_exact = zeros(2, n);  % Output: 2 × n

% Precompute forcing basis terms: 2 × n × n × num_steps
Tdelta_all = zeros(2, n, n, num_steps);
for m_s = 1:num_steps
    s = t_vals(m_s);
    Tdelta_all(:, :, :, m_s) = Tdelta(db, del0, j_vals, s);
end


% Initialize: 2 × n × n
R = zeros(2, n, n);

for idx_j = 1:n
    j = j_vals(idx_j);

    for m_s = 1:num_steps
        s = t_vals(m_s);
        t_minus_s = T_end - s;

        % Compute Bessel-based Phi: 2 × 2 × n
        shifts = j - j_vals';
        arg = 2 * t_minus_s;
        j0 = besselj(2 * shifts,   arg);
        jp = besselj(2 * shifts+1, arg);
        jm = besselj(2 * shifts-1, arg);

        Phi_vals = zeros(2, 2, n);
        Phi_vals(1,1,:) = j0;
        Phi_vals(1,2,:) = -jp;
        Phi_vals(2,1,:) = -jm;
        Phi_vals(2,2,:) = j0;

        % Extract all F_k and G_k for this step: 2 × n × n
        Tdelta_s = Tdelta_all(:, :, :, m_s);

        % Reshape to: 2 × 1 × n × n
        Tdelta_s = reshape(Tdelta_s, 2, 1, n, []);

        % Broadcast Phi: 2 × 2 × n, apply over batch dim K
        % Use pagemtimes on 2×2×n and 2×1×n×K -> output: 2×1×n×K
        PT = pagemtimes(Phi_vals, Tdelta_s);

        % Sum over k (dim 3): result is 2 × 1 × K
        sum_k = sum(PT, 3);  % 2 × 1 × K

        % Squeeze and permute to 2 × K
        sum_k = reshape(sum_k, 2, []);  % 2 × K

        % Accumulate into R: 2 × K × j
        R(:, :, idx_j) = R(:, :, idx_j) + dt * sum_k;
    end
end

% Final summation against kappa
for idx_j = 1:n
    % Extract accumulated products: 2 × length(kappa)
    R_j = R(:, :, idx_j);
    % Sum weighted by kappa: 2 × 1
    eta_exact(:, idx_j) = sum(R_j .* kappa, 2);
end
toc

figure;
plot(j_vals, eta_t(1,:,end), 'b.', 'LineWidth', 2);
hold on;
plot(j_vals ,eta_exact(1,:), 'r.', 'LineWidth', 2)
legend('Numerical RK4','Random Sum');
title(['Comparison at t = ', num2str(T_end)]);
xlabel('Lattice site j');
ylabel('\eta_1(j,t)');
grid on;

figure;
plot(j_vals, eta_t(2,:,end), 'b.', 'LineWidth', 2);
hold on;
plot(j_vals, eta_exact(2,:), 'r.', 'LineWidth', 2)
legend('Numerical RK4','Random Sum');
title(['Comparison at t = ', num2str(T_end)]);
xlabel('Lattice site j');
ylabel('\eta_2(j,t)');
grid on;


S1_vals = zeros(1, n);
S2_vals = zeros(1, n);

for j_idx = 1:n
    j = j_vals(j_idx);
    R1_slice = squeeze(R(1, :,j_idx));  % vector over m
    R2_slice = squeeze(R(2, :,j_idx));  % vector over m
    S1 = sqrt(sum(R1_slice.^2));
    S2 = sqrt(sum(R2_slice.^2));
    S1_vals(j_idx) = S1;
    S2_vals(j_idx) = S2;
end


% Plot
figure;
plot(j_vals, S1_vals, 'b.');
hold on
plot(j_vals, S2_vals, 'r.');
xlabel('j');
title('\surd{\Sigma_m R1(j, m, T)^2}');
grid on;

