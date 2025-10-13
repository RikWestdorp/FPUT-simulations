function basis_tensor = ForcingBasisTermsF(db, del0, j_vals, t)
    c0 = 1 + del0 / 24;
    N = (length(j_vals) - 1) / 2;
    n = length(j_vals);

    % Profiles (row vectors)
    rsol    = profile(db, 'r',    c0, t * c0, j_vals);  % 1 × n
    dlrsol  = -profile(db, 'dxr', c0, t * c0, j_vals);  % 1 × n
    dcrsol  =  profile(db, 'dcr', c0, t * c0, j_vals);  % 1 × n

    % Convert to row vectors
    rsol    = rsol(:).';     % 1 × n
    dlrsol  = dlrsol(:).';
    dcrsol  = dcrsol(:).';

    % Build 2 × n × n tensor of rhs vectors for each unit kappa = δ_m
    rhs = zeros(2, 1, n);  % 2 × 1 × n
    rhs(1,1,:) = dlrsol .* rsol;  % inner product part 1
    rhs(2,1,:) = dcrsol .* rsol;  % inner product part 2

    % Construct A and compute A⁻¹ * rhs for all m using pagemtimes
    [alpha0, alpha1] = alpha(db, del0, N);
    A = [0, alpha0; -alpha0, alpha1];  % 2 × 2
    Ainv = inv(A);                     % 2 × 2
    coeffs = pagemtimes(Ainv, -rhs);   % 2 × 1 × n

    % Evaluate xi terms
    [xi11, xi12, xi21, xi22] = xi(db, N, del0, t * c0);  % each 1 × n

    % Reshape xi terms to 1 × n for broadcasting
    xi11 = xi11(:).'; xi12 = xi12(:).';
    xi21 = xi21(:).'; xi22 = xi22(:).';

    % Initialize output tensor: 2 × n × n
    basis_tensor = zeros(2, n, n);

    % Combine xi terms with coefficients
    for m = 1:n
        gamma_dot = coeffs(1,1,m);  % scalar
        delta_dot = coeffs(2,1,m);  % scalar

        e_m = zeros(1, n); e_m(m) = 1;  % δ_m
        prod = e_m .* rsol;            % 1 × n
        d_kappar = prod - circshift(prod, 1);  % forward difference

        basis_tensor(1,:,m) = -gamma_dot * xi11 - delta_dot * xi21;  % row
        basis_tensor(2,:,m) = -gamma_dot * xi12 - delta_dot * xi22;  % row

    end
end
