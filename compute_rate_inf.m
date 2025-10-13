function rate_inf = compute_rate_inf(db, c0, N, dt, tau_max)
% Computes the expected rate in the infinite-time limit using precomputed Tdelta and Rinf.
j_vals = (-N:N)';
n = 2*N+1;
del0 = 24*(c0-1);
tau_vals = 0:dt:tau_max;
tau_steps = length(tau_vals);

% Precompute Tdelta for s = -tau
Tdelta_all = zeros(2, n, n, tau_steps);  % 2 × n × n × tau_steps
for m_tau = 1:tau_steps
    tau = tau_vals(m_tau);
    Tdelta_all(:, :, :, m_tau) = Tdelta(db, del0, j_vals, -tau);
end

% Initialize Rinf: 2 × n × n
Rinf = zeros(2, n, n);

for idx_j = 1:n
    j = j_vals(idx_j);
    disp(idx_j)

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
        Rinf(:, :, idx_j) = Rinf(:, :, idx_j) + dt * sum_k;
    end
end

% Compute rate_inf
D_inf = 0;
LIIL_inf = 0;

for m = 1:n
    Rinf1m = squeeze(Rinf(1, m, :))';
    Rinf2m = squeeze(Rinf(2, m, :))';
    delta_m = zeros(1, n); delta_m(m) = 1;

    D_inf = D_inf + D(db, c0, 0, j_vals, Rinf1m + delta_m, Rinf1m);
    LIIL_inf = LIIL_inf + LII(db, c0, 0, j_vals, Rinf1m, Rinf2m) * LI(db, c0, 0, j_vals, delta_m);
end

rate_inf =  D_inf + LIIL_inf;
end
