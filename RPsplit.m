function [R0, R1, R2, P1, P2] = RPsplit(db, eta1,eta2,N,del,l,kappa, potentialType)
        rsol=profile(db, 'r', 1+del/24, l, (-N:N));
        psol=profile(db, 'p', 1+del/24, l, (-N:N));
        [alpha0, alpha1]= alpha(db, del,N);
        [xi11, xi12, xi21, xi22]= xi(db, N,del,l);
        zer=zeros(size(xi11));
               
        %[JHprime1sol, JHprime2sol]= JHprime(rsol, psol, potentialType);
        [lin1,lin2]=J(eta1+2*rsol.*eta1,eta2);
        [N1,N2]=J(eta1.^2,zer);
        % [JHprime1, JHprime2]= JHprime(rsol+eta1, psol+eta2, potentialType);
        % [JHprime1eta, JHprime2eta]= JHprime(eta1, eta2, potentialType);
        %R0=[Omega(xi11, xi12,JHprime1sol, JHprime2sol);Omega(xi21, xi22,JHprime1sol, JHprime2sol)];

        R0=[0;-(1+del/24)*alpha0];
        %R1=-[sum(xi11.*(eta1+2*rsol.*eta1))+sum(xi12.*eta2);
            %sum(xi21.*(eta1+2*rsol.*eta1))+sum(xi22.*eta2)];
        %R2= -[sum(xi11.*(eta1.^2));sum(xi21.*(eta1.^2))];
        R1=[Omega(xi11, xi12,lin1, lin2);Omega(xi21, xi22,lin1, lin2)];
        R2=[Omega(xi11, xi12,N1, N2);Omega(xi21, xi22,N1, N2)];
        % P1=-[sum(xi11.*kappa.*(rsol));sum(xi21.*kappa.*(rsol))];
        % P2=-[sum(xi11.*kappa.*(eta1));sum(xi21.*kappa.*(eta1))];
        delta_kapparsol=kappa.*rsol -circshift(kappa.*rsol, 1);
        delta_kappaeta1=kappa.*eta1 -circshift(kappa.*eta1, 1);
        
        P1=[Omega(xi11, xi12,zer, delta_kapparsol);Omega(xi21, xi22,zer, delta_kapparsol)];
        
        P2=[Omega(xi11, xi12,zer, delta_kappaeta1);Omega(xi21, xi22,zer, delta_kappaeta1)];
   end