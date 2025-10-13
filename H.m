function H = H(r,p)
    H=sum(p.^2 / 2 + r.^2 / 2 + r.^3 / 3, 2);  % Sum over columns (N)
end