function [doteta1, doteta2] = EtaApproxDot(db, eta1,eta2,del0,t, kappa,N)
    c0=1+del0/24;
    
    rsol= profile(db, 'r', c0, t*c0, (-N:N));
    dlrsol= -profile(db, 'dxr', c0, t*c0, (-N:N));
    dcrsol= profile(db, 'dcr', c0, t*c0, (-N:N));
    [alpha0, alpha1] = alpha(db, del0,N);
    A=[0, alpha0;-alpha0, alpha1];
    temp = -A\[sum(dlrsol.*rsol.*kappa);sum(dcrsol.*rsol.*kappa)];
    gammadotapprox=temp(1);
    deldotapprox=temp(2);
  
    [JHlin1,JHlin2] = J(eta1+2*rsol.*eta1,eta2);
    [xi11, xi12, xi21, xi22] = xi(db, N,del0,t*c0);
    
    d_kappar=kappa.*(rsol)- circshift(kappa .* (rsol), 1);
    
    doteta1 = JHlin1+(-gammadotapprox)*xi11-deldotapprox*xi21;
    doteta2 = JHlin2+d_kappar+(-gammadotapprox)*xi12-deldotapprox*xi22;
end
