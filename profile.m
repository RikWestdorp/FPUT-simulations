function val = profile(db, label, c_input, pos, x)
    % Extract list of available c values
    c_list = [db.c];
    cmin = min(c_list); 
    cmax = max(c_list);

      % --- NEW: analytic branch for c below table range ---
    if c_input < cmin
        % Variables per your spec
        y  = 24*(c_input - 1);
        xi = x - pos;

        % Helpers (use elementwise ops)
        s = 1./cosh(0.5*sqrt(y).*xi);  % sech((1/2)*sqrt(y)*xi)
        t = tanh(0.5*sqrt(y).*xi);

        switch label
            case 'r'
                val = (1/8)*y.*(s.^2);
            case 'p'
                val = -(1/8)*y.*(s.^2);
            case 'dxr'
                val = -(1/8)*y.^(3/2).*(s.^2).*t;
            case 'dxp'
                val = (1/8)*y.^(3/2).*(s.^2).*t;
            case 'dcr'
                val = -1/16*(s.^2).*(-2 + xi.*sqrt(y).*t);
            case 'dcp'
                val = +1/16*(s.^2).*(-2 + xi.*sqrt(y).*t);
            case 'dxdxr'
                val = -1/16*y.^2.*(s.^2).*((s.^2) - 2*(t.^2));
            case 'dxdxp'
                val = +1/16*y.^2.*(s.^2).*((s.^2) - 2*(t.^2));
            case 'dcdxr'
                val = 1/32*( -xi.*y.*(s.^4) + 2*(s.^2).*(-3*sqrt(y).*t + xi.*y.*(t.^2)) );
            case 'dcdxp'
                val = -1/32*( -xi.*y.*(s.^4) + 2*(s.^2).*(-3*sqrt(y).*t + xi.*y.*(t.^2)) );
            case 'dcdcr'
                val = 1/64*xi.*(s.^2).*( -xi.*(s.^2) + 2*t.*(-3./sqrt(y) + xi.*t) );
            case 'dcdcp'
                val = -1/64*xi.*(s.^2).*( -xi.*(s.^2) + 2*t.*(-3./sqrt(y) + xi.*t) );
            otherwise
                error('Unknown label "%s".', label);
        end
        return;
    end
    % --- end analytic branch ---

    % Keep your original bounds behavior for the upper side
    if c_input > cmax
        error('Requested c = %.6f is out of bounds.', c_input);
    end

    % Find indices of two c-values surrounding c_input
    idx_upper = find(c_list >= c_input, 1);
    idx_lower = find(c_list <= c_input, 1, 'last');

    % Handle exact match
    if idx_upper == idx_lower
        entry = db(idx_upper);
        val = interp1(entry.xi, entry.(label), x - pos, 'spline', 0);
        return;
    end

    % Get entries at both c-values
    c_lo = c_list(idx_lower);
    c_hi = c_list(idx_upper);
    e_lo = db(idx_lower);
    e_hi = db(idx_upper);

    % Shift x
    x_shifted = x - pos;

    % Interpolate each value at x for both entries
    y_lo = interp1(e_lo.xi, e_lo.(label), x_shifted, 'spline', 0);
    y_hi = interp1(e_hi.xi, e_hi.(label), x_shifted, 'spline', 0);

    % Linear interpolation in c
    t = (c_input - c_lo) / (c_hi - c_lo);
    val = (1 - t) * y_lo + t * y_hi;
end
