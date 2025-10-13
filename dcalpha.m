function [alpha0, alpha1] = dcalpha(db, del,N)
    dlr = -profile(db, 'dxr', 1+del/24, 0, (-N:N));
    dlp = -profile(db, 'dxp', 1+del/24, 0, (-N:N));
    dcr = profile(db, 'dcr', 1+del/24, 0, (-N:N));
    dcp = profile(db, 'dcp', 1+del/24, 0, (-N:N));
    dcdcr = profile(db, 'dcdcr', 1+del/24, 0, (-N:N));
    dcdcp = profile(db, 'dcdcp', 1+del/24, 0, (-N:N));
    dcdlr = -profile(db, 'dcdxr', 1+del/24, 0, (-N:N));
    dcdlp = -profile(db, 'dcdxp', 1+del/24, 0, (-N:N));
    alpha0=Omega(dcdlr, dcdlp,dcr, dcp)+Omega(dlr, dlp,dcdcr, dcdcp);
    alpha1=Omega(dcdcr,dcdcp,dcr,dcp)+Omega(dcr,dcp,dcdcr,dcdcp);
end