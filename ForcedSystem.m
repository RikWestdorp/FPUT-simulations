function [dotr, dotp] = ForcedSystem(r, p, springforce, potentialType)
    % Compute derivatives for the system with varying springforce
    [JHprime1, JHprime2] = JHprime(r, p, potentialType);
    dotr = JHprime1;
    dotp = JHprime2 + springforce .* r - circshift(springforce .* r, 1);

end
