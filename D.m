function D = D(db,c0,t,j_vals,f,g)
    del0=24*(c0-1);
    N = (length(j_vals)-1)/2;
    dlrsol  = -profile(db, 'dxr', c0, c0*t, j_vals');
    dcrsol  =  profile(db, 'dcr', c0, c0*t, j_vals');
    [alpha0, alpha1] = alpha(db, del0, N);
    A = [0, alpha0; -alpha0, alpha1];
    S1=sum(dlrsol.*f.*g);
    S2=sum(dcrsol.*f.*g);
    %D=-0.5*A\[S1;S2] ;
    D=-A\[S1;S2] ;
end
