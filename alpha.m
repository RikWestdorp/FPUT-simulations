function [alpha0, alpha1] = alpha(db, del,N)
    [xi11, xi12, xi21, xi22]= xi(db, N,del,0);
    alpha0=Omega(xi11, xi12,xi21, xi22);
    alpha1=Omega(xi21,xi22,xi21,xi22);
end