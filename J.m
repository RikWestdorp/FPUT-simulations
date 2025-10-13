function [res1, res2] = J(r, p)
    res1=circshift(p, -1) - p; 
    res2=r-circshift(r,1);
end