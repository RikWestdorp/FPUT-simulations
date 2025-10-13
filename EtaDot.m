function [doteta1, doteta2] = EtaDot(db, eta1,eta2,ldot,deldot, del,l, springforce,potentialType,N)
    c=1+del/24;
    % doteta1 = circshift(p, -1) - p;
    % doteta2 = eta1+circshift(q0,-ceil(c0*t)).*(2*eta1+springforce)- circshift(eta1+circshift(q0,-ceil(c0*t)).*(2*eta1+springforce), 1);
    % rsol=(1/8)*del.*sech((1/2)*sqrt(del).*((-N:N)-l)).^2;
    % psol=-(1/8)*del.*sech((1/2)*sqrt(del).*((-N:N)-l)).^2;
    rsol = profile(db, 'r', c, l, (-N:N));
    % psol = profile(db, 'p', c, l, (-N:N));
    % [JHprime1, JHprime2] = JHprime(rsol+eta1, psol+eta2, potentialType);
    % [JHprime01, JHprime02] = JHprime(rsol, psol, potentialType);
    [JH1,JH2] = J(eta1+2*rsol.*eta1+eta1.^2,eta2);
    [xi11, xi12, xi21, xi22] = xi(db, N,del,l);
    
    d_kappar=springforce.*(rsol+eta1)- circshift(springforce .* (rsol+eta1), 1);
    % dpluseta2 = circshift(eta2, -1) - eta2;
    % eta1term= eta1+2*rsol.*eta1+eta1.^2;
    % d_eta1term = eta1term-circshift(eta1term,1);

    %doteta1 = JHprime1-JHprime01+(c-ldot)*xi11-deldot*xi21;
    %doteta2 = JHprime2-JHprime02+d_kappar+(c-ldot)*xi12-deldot*xi22;
    doteta1 = JH1+(c-ldot)*xi11-deldot*xi21;
    doteta2 = JH2+d_kappar+(c-ldot)*xi12-deldot*xi22;
end
