function basis_tensor = ForcingBasisTermsG(db, del0, j_vals, t)
    c0 = 1 + del0 / 24;
    N = (length(j_vals)-1)/2;
    n = length(j_vals);

    rsol = profile(db, 'r', c0, t * c0, j_vals);  % 1 × n
    rsol = rsol(:).';  % force row vector 1 × n

    % Compute the effect of kappa = δ_m for each m
    basis_tensor = zeros(2, n, n);  % 2 × n × n

    for m = 1:n
        e_m = zeros(1, n); e_m(m) = 1;  % δ_m
        prod = e_m .* rsol;            % 1 × n
        d_kappar = prod - circshift(prod, 1);  % forward difference

        basis_tensor(1, :, m) = 0;           % doteta1 = 0
        basis_tensor(2, :, m) = d_kappar;    % doteta2
    end
end
