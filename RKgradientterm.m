function [delhtest, lhtest, delhtestnew, lhtestnew,delexp, lexp, xi1, E123, E4] = RKgradientterm(db, h, kappa, time, potentialType, b,delstar,N)
    
    iter = time * (1/h);  % Total number of iterations
    cstar=1+delstar/24;
    saveIndices = 1:b:iter+1;  % Indices to save every b-th step
    numSaved = length(saveIndices);  % Number of stored timepoints
    
    delexp= zeros(1, numSaved);
    lexp= zeros(1, numSaved);
    xi1=zeros(1, numSaved);

    E123=zeros(1, numSaved);
    E4=zeros(1, numSaved);
 
    delhtest= zeros(1, numSaved);
    lhtest= zeros(1, numSaved);

    delhtestnew= zeros(1, numSaved);
    lhtestnew= zeros(1, numSaved);
    
    % Set initial conditions
    
    
    delexp(1) = delstar;
    lexp(1) = 0;
    xi1(1) = 0;
    E123(1)=0;
    E4(1)=0;

    delhtest(1) = delstar;
    lhtest(1) = 0;
    delhtestnew(1) = delstar;
    lhtestnew(1) = 0;

    delexpCurrent = delstar;
    lexpCurrent = 0;
    xi1Current = 0;
    delhtestCurrent = delstar;
    lhtestCurrent = 0;
    delhtestnewCurrent = delstar;
    lhtestnewCurrent = 0;

    E123Current = 0;
    E4Current = 0;

    saveIdx = 2;  % Index for storing results

    for i = 1:iter
        %disp(i);

        % Compute dot of approximations
        
        [lexpdot, delexpdot] = modszeroalt(db, N, potentialType, delstar, i*h*cstar, kappa);

        [lhtestdot, delhtestdot] = modszeroalt(db, N, potentialType, delexpCurrent, xi1Current, kappa);

        [dldxi, ddeldxi] = dximodszero(db, N, potentialType, delstar, i*h*cstar, kappa);
        [dldc, ddeldc] = dcmodszero(db, N, potentialType, delstar, i*h*cstar, kappa); 

        lhtestnewdot = lexpdot+(xi1Current-i*h*cstar)*dldxi+(delexpCurrent-delstar)*dldc; 
        delhtestnewdot = delexpdot+(xi1Current-i*h*cstar)*ddeldxi+(delexpCurrent-delstar)*ddeldc; 

        E123dot = (xi1Current-i*h*cstar)*ddeldxi;
        E4dot = (delexpCurrent-delstar)*ddeldc;

        
        % Updating
        xi1Current = xi1Current+h*((delexpCurrent-delstar)/24+lexpdot); 
        delexpCurrent= delexpCurrent+h*delexpdot;
        lexpCurrent= lexpCurrent+h*lexpdot;    

        E123Current = E123Current+h*E123dot;
        E4Current = E4Current+h*E4dot;
        
        delhtestCurrent= delhtestCurrent+h*delhtestdot;
        lhtestCurrent= lhtestCurrent+h*lhtestdot;
          
        delhtestnewCurrent= delhtestnewCurrent+h*delhtestnewdot;
        lhtestnewCurrent= lhtestnewCurrent+h*lhtestnewdot;

        % Save only every b-th timepoint
        if mod(i, b) == 0
          
            delhtest(saveIdx) = delhtestCurrent;
            lhtest(saveIdx) = lhtestCurrent;
            delhtestnew(saveIdx) = delhtestnewCurrent;
            lhtestnew(saveIdx) = lhtestnewCurrent;
            delexp(saveIdx) = delexpCurrent;
            lexp(saveIdx) = lexpCurrent;
            xi1(saveIdx) = xi1Current;
            E123(saveIdx) = E123Current;
            E4(saveIdx) = E4Current;

            saveIdx = saveIdx + 1;
        end
    end  
end
