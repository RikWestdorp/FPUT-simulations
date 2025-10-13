function [doteta1, doteta2] = Etadotapprox(db, eta1,eta2,ldot,deldot, del,l, kappa,potentialType,N)
    c=1+del/24;
    rsol = profile(db, 'r', c, l, (-N:N));
    psol = profile(db, 'p', c, l, (-N:N));
    [JHprime1, JHprime2] = JHprime(rsol+eta1, psol+eta2, potentialType);
    [JHprimesol1, JHprimesol2] = JHprime(rsol, psol, potentialType);
    zer=zeros(size(rsol));
    [quad1, quad2] = J(eta1.^2,zer);
    [xi11, xi12, xi21, xi22] = xi(db, N,del,l);
   
    [lin1,lin2]=J(eta1+2*rsol.*eta1,eta2);
    
    d_kappar=kappa.*(rsol+eta1)- circshift(kappa .* (rsol+eta1), 1);

    doteta1 = JHprime1-JHprimesol1-quad1+(c-ldot)*xi11-deldot*xi21;
    doteta2 = JHprime2-JHprimesol2-quad2+d_kappar+(c-ldot)*xi12-deldot*xi22;
end
