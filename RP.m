function [R, P] = RP(db, eta1,eta2,N,del,l, potentialType,kappa)
    rsol=profile(db, 'r', 1+del/24, l, (-N:N));
    r=rsol+eta1;
    psol=profile(db, 'p', 1+del/24, l, (-N:N));
    p=psol+eta2;
    [alpha0, alpha1] = alpha(db, del,N);
    [xi11, xi12, xi21, xi22]= xi(db, N,del,l);
    [JHprime1, JHprime2]= JHprime(rsol+eta1, psol+eta2, potentialType);
    [JHprime01, JHprime02]= JHprime(rsol, psol, potentialType);
    R=[Omega(xi11, xi12,JHprime1, JHprime2);Omega(xi21, xi22,JHprime1, JHprime2)]-[Omega(xi11, xi12,JHprime01, JHprime02);Omega(xi21, xi22,JHprime01, JHprime02)]+[0 ; -(1+del/24)*alpha0];
    delta_kappar=kappa.*r -circshift(kappa.*r, 1);
    zer=zeros(size(xi11));
    P=[Omega(xi11, xi12,zer, delta_kappar);Omega(xi21, xi22,zer, delta_kappar)];
end