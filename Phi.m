function P = Phi(j, t)
    P = [besselj(2*j, 2*t), -besselj(2*j+1, 2*t);
        -besselj(2*j-1, 2*t), besselj(2*j, 2*t)];
end




