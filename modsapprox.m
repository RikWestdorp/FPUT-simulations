function [ldot, cdot] = modsapprox(db, N, potentialType, eta1, eta2, del, l, kappa)
    %x=(-N:N)-l;
    % dcdcphi = (1/64)*x.*sech(x*sqrt(del)/2).^2.*(-x.*sech(x*sqrt(del)/2).^2 + 2*tanh(x*sqrt(del)/2).*((-3/sqrt(del)) + x.*tanh(x*sqrt(del)/2)));
    % dldlphi=-(1/16)* del^2*sech(sqrt(del)*x/2).^2.*(sech(sqrt(del)*x/2).^2 - 2*tanh(x*(1/2)*sqrt(del)).^2);
    % dcdxphi=(1/16)*(-2*sech((1/2)*x*sqrt(del)).^2.*tanh((1/2)*x*sqrt(del))*(1/2)*sqrt(del).*(2 - x*sqrt(del).*tanh((1/2)*x*sqrt(del))) + sech((1/2)*x*sqrt(del)).^2.*(-sqrt(del).*tanh((1/2)*x*sqrt(del)) - (1/2)*x*del.*sech((1/2)*x*sqrt(del)).^2));
    % dcdlphi=-dcdxphi;
    dldlr = profile(db, 'dxdxr', 1+del/24, l, (-N:N));
    dldlp = profile(db, 'dxdxp', 1+del/24, l, (-N:N));
    dcdlr = -profile(db, 'dcdxr', 1+del/24, l, (-N:N));
    dcdlp = -profile(db, 'dcdxp', 1+del/24, l, (-N:N));
    dcdcr = profile(db, 'dcdcr', 1+del/24, l, (-N:N));
    dcdcp = profile(db, 'dcdcp', 1+del/24, l, (-N:N));

    [alpha0, alpha1]= alpha(db, del,N);
    
    A=[0, alpha0;-alpha0, alpha1];
    Ainv=inv(A);
    B=[Omega(dldlr, dldlp, eta1, eta2) , Omega(dcdlr, dcdlp, eta1, eta2); Omega(dcdlr, dcdlp, eta1, eta2), Omega(dcdcr, dcdcp, eta1, eta2)];
    Kinv0= Ainv; 
    Kinv1=Ainv*B*Ainv;


    [R0, R1, R2, P1, P2] = RPsplit(db, eta1,eta2,N,del,l,kappa,potentialType);

    temp=Kinv0*(R0+R1+R2+P1+P2)+Kinv1*(R0+R1+P1);
    ldot=temp(1);
    cdot=temp(2);
    
end