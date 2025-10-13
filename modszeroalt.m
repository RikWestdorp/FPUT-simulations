function [ldot, cdot] = modszeroalt(db, N, potentialType, del, l, kappa)
    rsol=profile(db, 'r', 1+del/24, l, (-N:N));
    dcr=profile(db, 'dcr', 1+del/24, l, (-N:N));
    dxr=profile(db, 'dxr', 1+del/24, l, (-N:N));

    zer=zeros(size(rsol));

    [alpha0, alpha1]= alpha(db, del,N);
    
    A=[0, alpha0;-alpha0, alpha1];
    Ainv=inv(A);
    S=[-sum(rsol.*dxr.*kappa);sum(rsol.*dcr.*kappa) ];
    temp=-Ainv*S;
    ldot=1+del/24+temp(1);
    cdot=sum(rsol.*dxr.*kappa)/alpha0;%temp(2);
    
end