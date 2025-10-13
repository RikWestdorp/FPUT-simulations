function [c, l] = fit_ortho(db, N, r, p)

    x = (-N:N);

    % Initial guess
    c0 = min(1.3, max(1.01, 8 * max([r, -p])));
    [~, idx_r_max] = max(r);
    [~, idx_p_max] = max(-p);
    l0 = mean([x(idx_r_max), x(idx_p_max)]);
    initial_guess = [c0, l0];

    % Adjusted bounds: margin to avoid profile() errors
    lb = [1.005, min(x)];
    ub = [1.4, max(x)];

    % System of equations
    system = @(params) [
        Omega_xi(N, params(1), params(2), x, r, p, db, 1);
        Omega_xi(N, params(1), params(2), x, r, p, db, 2)
    ];

    % Solver options
    options = optimoptions('lsqnonlin', 'Display', 'off', 'TolFun', 1e-13, 'TolX', 1e-13);

    % Solve
    best_params = lsqnonlin(system, initial_guess, lb, ub, options);

    c = best_params(1);
    l = best_params(2);
end
