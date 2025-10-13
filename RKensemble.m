function [r, p] = RKensemble(db, r0, p0, h, kappa, time, potentialType, b,delstar,N)
    
    iter = time * (1/h);  % Total number of iterations
    cstar=1+delstar/24;
    saveIndices = 1:b:iter+1;  % Indices to save every b-th step
    numSaved = length(saveIndices);  % Number of stored timepoints

    r = zeros(numSaved, 2*N+1);  % Space-time solution matrix for r
    p = zeros(numSaved,2*N+1);  % Space-time solution matrix for p
    
    % del0= zeros(1, numSaved);
    % l0= zeros(1, numSaved);
   
    
    % Set initial conditions
    r(1, :) = r0;
    p(1, :) = p0;
    
    % del0(1) = delstar;
    % l0(1) = 0;

    % Runge-Kutta 4th order explicit scheme
    rCurrent = r0;
    pCurrent = p0;

    % del0Current = delstar;
    % l0Current = 0;
    % 
    % l0dot = cstar;
    % del0dot = 0;
    
    saveIdx = 2;  % Index for storing results

    for i = 1:iter

       
        % Runge Kutta step for r, p
        [r1, p1] = ForcedSystem(rCurrent, pCurrent, kappa, potentialType);
        [r2, p2] = ForcedSystem(rCurrent + h*r1/2, pCurrent + h*p1/2, kappa, potentialType);
        [r3, p3] = ForcedSystem(rCurrent + h*r2/2, pCurrent + h*p2/2, kappa, potentialType);
        [r4, p4] = ForcedSystem(rCurrent + h*r3, pCurrent + h*p3, kappa, potentialType);

        % zer=zeros(1,2*N+1);
        % [l0dot, del0dot]=modsapprox(db, N, potentialType,zer,zer, del0Current, l0Current, kappa);
        
     
        
        % Updating
       
        % del0Current = del0Current+h*del0dot;
        % l0Current = l0Current+h*l0dot;
        
        rCurrent = rCurrent + h*(r1 + 2*r2 + 2*r3 + r4)/6;
        pCurrent = pCurrent + h*(p1 + 2*p2 + 2*p3 + p4)/6;
        

        % Save only every b-th timepoint
        if mod(i, b) == 0
            r(saveIdx, :) = rCurrent;
            p(saveIdx, :) = pCurrent;
          
      
            % del0(saveIdx) = del0Current;
            % l0(saveIdx) = l0Current;

            saveIdx = saveIdx + 1;
        end
    end  
end
