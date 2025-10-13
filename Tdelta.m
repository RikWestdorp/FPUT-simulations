function basis_tensor = Tdelta(db, del0, j_vals, t)
    basis_tensor = ForcingBasisTermsF(db, del0, j_vals, t)+ForcingBasisTermsG(db, del0, j_vals, t);
end
