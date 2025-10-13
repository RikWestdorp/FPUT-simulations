function out = Potential(v, potentialType)
    % Uprime computes the derivative of the potential energy U with respect to v.
    % potentialType allows switching between different potential models.
    
    switch potentialType
        case 'quadratic'  % Quadratic potential: U'(v) = v
            out = v;
            
        case 'cubic'  % Cubic potential: U'(v) = v + v^2
            out = v + v.^2;
            
        otherwise
            error('Invalid potential type. Choose from quadratic, cubic.');
    end
end
