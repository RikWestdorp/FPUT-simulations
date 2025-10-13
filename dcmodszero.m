function [ldot, cdot] = dcmodszero(db, N, potentialType, del, l, kappa)
    rsol=profile(db, 'r', 1+del/24, l, (-N:N));
    dcr=profile(db, 'dcr', 1+del/24, l, (-N:N));
    dxr=profile(db, 'dxr', 1+del/24, l, (-N:N));
    dxdxr=profile(db, 'dxdxr', 1+del/24, l, (-N:N));
    dxdcr=profile(db, 'dcdxr', 1+del/24, l, (-N:N));
    dcdcr=profile(db, 'dcdcr', 1+del/24, l, (-N:N));

    zer=zeros(size(rsol));

    [alpha0, alpha1]= alpha(db, del,N);
    [dcalpha0, dcalpha1]= dcalpha(db, del,N);
    A=[0, alpha0;-alpha0, alpha1];
    Ainv=inv(A);
    dcA=[0, dcalpha0;-dcalpha0, dcalpha1];
    dcAinv=inv(dcA);
    S=[-sum(rsol.*dxr.*kappa);sum(rsol.*dcr.*kappa) ];
    dcS=[-sum((dcr.*dxr+rsol.*dxdcr).*kappa);sum((dcr.*dcr+rsol.*dcdcr).*kappa) ];
    temp=-Ainv*dcS-dcAinv*S;
    ldot=temp(1);
    cdot=temp(2);
    
end