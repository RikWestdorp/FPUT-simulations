
function LI = LI(db,c0,t,j_vals,f)
    del0=24*(c0-1);
    N = (length(j_vals)-1)/2;
    rsol  = profile(db, 'r', c0, c0*t, j_vals');
    dlrsol  = -profile(db, 'dxr', c0, c0*t, j_vals');
    dcrsol  =  profile(db, 'dcr', c0, c0*t, j_vals');
    [alpha0, alpha1] = alpha(db, del0, N);
    A = [0, alpha0; -alpha0, alpha1];
    S1=sum(dlrsol.*rsol.*f);
    S2=sum(dcrsol.*rsol.*f);
    %LI=-2*A\[S1;S2] ;
    LI=-A\[S1;S2] ;
end
