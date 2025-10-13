
function LII = LII(db,c0,t,j_vals,h1,h2)
    del0=24*(c0-1);
    N = (length(j_vals)-1)/2;
    dldlr  = profile(db, 'dxdxr', c0, c0*t, j_vals');
    dcdlr  = -profile(db, 'dcdxr', c0, c0*t, j_vals');
    dcdcr  =  profile(db, 'dcdcr', c0, c0*t, j_vals');
    dldlp  = profile(db, 'dxdxp', c0, c0*t, j_vals');
    dcdlp  = -profile(db, 'dcdxp', c0, c0*t, j_vals');
    dcdcp  =  profile(db, 'dcdcp', c0, c0*t, j_vals');
    [alpha0, alpha1] = alpha(db, del0, N);
    A = [0, alpha0; -alpha0, alpha1];
    B=[Omega(dldlr, dldlp, h1, h2) , Omega(dcdlr, dcdlp, h1, h2); Omega(dcdlr, dcdlp, h1, h2), Omega(dcdcr, dcdcp, h1, h2)];
    %LII=-A\B;
    LII=A\B;
end