function [xi11, xi12, xi21, xi22] = xi(db, N,del,l)
    %x=((-N:N)-l);
    xi11 = -profile(db, 'dxr', 1+del/24, l, (-N:N));
    xi12 = -profile(db, 'dxp', 1+del/24, l, (-N:N));
    xi21 = profile(db, 'dcr', 1+del/24, l, (-N:N));
    xi22 = profile(db, 'dcp', 1+del/24, l, (-N:N));
    % xi11=(1/8)*del^(3/2)*sech((1/2)*sqrt(del)*x).^2.*tanh((1/2)*sqrt(del)*x);
    % xi12=-xi11;
    % xi21=(1/16)*sech(x*(1/2)*sqrt(del)).^2.*(2-x.*sqrt(del).*tanh(x*(1/2)*sqrt(del)));
    % xi22=-xi21;
    
end