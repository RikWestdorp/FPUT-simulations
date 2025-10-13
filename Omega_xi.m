function val = Omega_xi(N, c, l, x, r, p, db, which)
    % Get model profiles
    rmodel = profile(db, 'r', c, l, x);
    residual_r = r - rmodel;

    pmodel = profile(db, 'p', c, l, x);
    residual_p = p - pmodel;

    % Get quadrature weights and test functions
    del=24*(c-1);
    [xi11, xi12, xi21, xi22] = xi(db, N,del,l);

    if which == 1
        val = Omega(xi11, xi12, residual_r, residual_p);
    else
        val = Omega(xi21, xi22, residual_r, residual_p);
    end
end
