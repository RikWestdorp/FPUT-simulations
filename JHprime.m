function [res1, res2] = JHprime(r, p, potentialType)
    [res1, res2] = J(Potential(r, potentialType),p);
end