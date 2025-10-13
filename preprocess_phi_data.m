function db = preprocess_phi_data(c_groups, c_values)
    num_blocks = numel(c_values);
    db = struct();

    for i = 1:num_blocks
        block = c_groups{i};
        c = c_values(i);
        tau = block(:, 2);
        kappa = block(:, 1);
        xi = tau/kappa(1);
        r = block(:, 3) * kappa(1)^2;  
        dxr = block(:, 4) * kappa(1)^3;
        dxdxr = gradient(dxr, xi);
        dcr = block(:, 6) * kappa(1)^2/24;
        dcdxr = block(:, 7) * kappa(1)^3/24;
        dcdcr = block(:, 8)*kappa(1)^2/(24^2);

        xidouble = [-flipud(xi(2:end)); xi];
        rdouble = [flipud(r(2:end)); r];
        dxrdouble = [-flipud(dxr(2:end)); dxr];
        dxdxrdouble = [flipud(dxdxr(2:end)); dxdxr];
        dcrdouble = [flipud(dcr(2:end)); dcr];
        dcdxrdouble = [-flipud(dcdxr(2:end)); dcdxr];
        dcdcrdouble = [flipud(dcdcr(2:end)); dcdcr];
       

        dxi = xi(2) - xi(1);
        Vr = rdouble + rdouble.^2;
        xidouble_shifted = xidouble - 1;  % shift backward in x (past values)

        % Interpolate Vr at shifted locations, using zero for out-of-bounds
        Vr_shifted = interp1(xidouble, Vr, xidouble_shifted, 'spline', 0);
        dxp = -(1 / c) * (Vr-Vr_shifted);
        p = cumsum(dxp) * dxi;
        
        dxdxp = gradient(dxp, xidouble);
        Vcr = dcrdouble+2*rdouble.*dcrdouble;
        Vcr_shifted = interp1(xidouble, Vcr, xidouble_shifted, 'spline', 0);
        dcdxp = -1/(24*c)*dxp-(1/c)*(Vcr-Vcr_shifted);
        dcp = cumsum(dcdxp) * dxi;
        
        Vccr = dcdcrdouble+2*dcrdouble.^2+2*rdouble.*dcdcrdouble;
        Vccr_shifted = interp1(xidouble, Vccr, xidouble_shifted, 'spline', 0);
        dcdcdxp= -1/(12*c)*dcdxp-(1/c)*(Vccr-Vccr_shifted);
        dcdcp=cumsum(dcdcdxp) * dxi;

        db(i).c = c_values(i);
        db(i).xi = xidouble;

        db(i).r = rdouble;
        db(i).dxr = dxrdouble;
        db(i).dcr = dcrdouble;
        db(i).dxdxr = dxdxrdouble;
        db(i).dcdxr = dcdxrdouble;
        db(i).dcdcr = dcdcrdouble;

        db(i).p = p;
        db(i).dxp = dxp;
        db(i).dcp = dcp;
        db(i).dxdxp = dxdxp;
        db(i).dcdxp = dcdxp;
        db(i).dcdcp = dcdcp;

    end
end
