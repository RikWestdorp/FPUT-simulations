function res = Forcing(db, del0,j_vals, t, kappa)
    basis_tensorF = ForcingBasisTermsF(db, del0, j_vals, t);
    basis_tensorG = ForcingBasisTermsG(db, del0, j_vals, t);
    kappanew  = reshape(kappa, 1, 1, []);
    res = sum(basis_tensorF .* kappanew, 3)+sum(basis_tensorG .* kappanew, 3);
    
end