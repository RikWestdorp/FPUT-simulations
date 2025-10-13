function [r, p, etarmod, etapmod, etarlin, etaplin,etarh, etaph, delmod, lmod, del2, l2, dellin, llin,delh, lh,delhtest, lhtest,dellintest, llintest, delexp, lexp] = RK(db, r0, p0, h, kappa, time, potentialType, b,delstar,N)
    
    iter = time * (1/h);  % Total number of iterations
    cstar=1+delstar/24;
    saveIndices = 1:b:iter+1;  % Indices to save every b-th step
    numSaved = length(saveIndices);  % Number of stored timepoints

    r = zeros(numSaved, 2*N+1);  % Space-time solution matrix for r
    p = zeros(numSaved,2*N+1);  % Space-time solution matrix for p
    etarmod = zeros(numSaved,2*N+1);
    etapmod = zeros(numSaved,2*N+1);
    etarlin = zeros(numSaved,2*N+1);
    etaplin = zeros(numSaved,2*N+1);
    etarh = zeros(numSaved,2*N+1);
    etaph = zeros(numSaved,2*N+1);
    
    delmod= zeros(1, numSaved);
    lmod= zeros(1, numSaved);
    del2= zeros(1, numSaved);
    l2= zeros(1, numSaved);
    delexp= zeros(1, numSaved);
    lexp= zeros(1, numSaved);
    dellin= zeros(1, numSaved);
    llin= zeros(1, numSaved);
    delh= zeros(1, numSaved);
    lh= zeros(1, numSaved);
    delhtest= zeros(1, numSaved);
    lhtest= zeros(1, numSaved);
    dellintest= zeros(1, numSaved);
    llintest= zeros(1, numSaved);
    
    % Set initial conditions
    r(1, :) = r0;
    p(1, :) = p0;
    etarmod(1, :) = zeros(size(p0));
    etapmod(1, :) = zeros(size(p0));
    etar_1(1, :) = zeros(size(p0));
    etap_1(1, :) = zeros(size(p0));
    etarlin(1, :) = zeros(size(p0));
    etaplin(1, :) = zeros(size(p0));
    etarh(1, :) = zeros(size(p0));
    etaph(1, :) = zeros(size(p0));
    
    delmod(1) = delstar;
    lmod(1) = 0;
    del2(1) = delstar;
    l2(1) = 0;
    delexp(1) = delstar;
    lexp(1) = 0;
    dellin(1) = delstar;
    llin(1) = 0;
    delh(1) = delstar;
    lh(1) = 0;
    delhtest(1) = delstar;
    lhtest(1) = 0;
    dellintest(1) = delstar;
    llintest(1) = 0;


    % Runge-Kutta 4th order explicit scheme
    rCurrent = r0;
    pCurrent = p0;
    etarmodCurrent = zeros(size(r0));
    etapmodCurrent = zeros(size(r0));
    etarlinCurrent = zeros(size(r0));
    etaplinCurrent = zeros(size(r0));
    etarhCurrent = zeros(size(r0));
    etaphCurrent = zeros(size(r0));
    
    delmodCurrent = delstar;
    lmodCurrent = 0;
    del2Current = delstar;
    l2Current = 0;
    delexpCurrent = delstar;
    lexpCurrent = 0;
    dellinCurrent = delstar;
    llinCurrent = 0;
    delhCurrent = delstar;
    lhCurrent = 0;
    delhtestCurrent = delstar;
    lhtestCurrent = 0;
    dellintestCurrent = delstar;
    llintestCurrent = 0;
    xi1Current = 0;
    
    ldot = cstar;
    deldot = 0;
    l2dot = cstar;
    del2dot = 0;

    saveIdx = 2;  % Index for storing results

    for i = 1:iter
        disp(i);

       
        % Runge Kutta step for r, p
        [r1, p1] = ForcedSystem(rCurrent, pCurrent, kappa, potentialType);
        [r2, p2] = ForcedSystem(rCurrent + h*r1/2, pCurrent + h*p1/2, kappa, potentialType);
        [r3, p3] = ForcedSystem(rCurrent + h*r2/2, pCurrent + h*p2/2, kappa, potentialType);
        [r4, p4] = ForcedSystem(rCurrent + h*r3, pCurrent + h*p3, kappa, potentialType);

        % Runge Kutta for etarmod, etapmod
        [etar1, etap1] = EtaDot(db, etarmodCurrent, etapmodCurrent,ldot, deldot, delmodCurrent, lmodCurrent, kappa,potentialType,N);
        [etar2, etap2] = EtaDot(db, etarmodCurrent + h*etar1/2, etapmodCurrent + h*etap1/2,ldot, deldot, delmodCurrent, lmodCurrent, kappa,potentialType,N);
        [etar3, etap3] = EtaDot(db, etarmodCurrent + h*etar2/2, etapmodCurrent + h*etap2/2,ldot, deldot, delmodCurrent, lmodCurrent, kappa,potentialType,N);
        [etar4, etap4] = EtaDot(db, etarmodCurrent + h*etar3, etapmodCurrent + h*etap3,ldot, deldot, delmodCurrent, lmodCurrent, kappa,potentialType,N);

        % Runge Kutta for etarlin, etaplin 
     
        [etarlin1, etaplin1] = EtaApproxDot(db, etarlinCurrent, etaplinCurrent,delstar,i*h, kappa,N);
        [etarlin2, etaplin2] = EtaApproxDot(db, etarlinCurrent + h*etarlin1/2, etaplinCurrent + h*etaplin1/2,delstar,i*h, kappa,N);
        [etarlin3, etaplin3] = EtaApproxDot(db, etarlinCurrent + h*etarlin2/2, etaplinCurrent + h*etaplin2/2,delstar,i*h, kappa,N);
        [etarlin4, etaplin4] = EtaApproxDot(db, etarlinCurrent + h*etarlin3, etaplinCurrent + h*etaplin3,delstar,i*h, kappa,N);

         % Runge Kutta for etarh, etaph 

        [etarh1, etaph1] = EtaFreeDot(db, etarhCurrent, etaphCurrent,delstar,i*h, kappa,N);
        [etarh2, etaph2] = EtaFreeDot(db, etarhCurrent + h*etarh1/2, etaphCurrent + h*etaph1/2,delstar,i*h, kappa,N);
        [etarh3, etaph3] = EtaFreeDot(db, etarhCurrent + h*etarh2/2, etaphCurrent + h*etaph2/2,delstar,i*h, kappa,N);
        [etarh4, etaph4] = EtaFreeDot(db, etarhCurrent + h*etarh3, etaphCurrent + h*etaph3,delstar,i*h, kappa,N);

        % Compute cdot, ldot via modulation system 
        [ldot, deldot] = mods(db, N, potentialType,etarmodCurrent,etapmodCurrent, delmodCurrent, lmodCurrent, kappa);

        % Compute dot of approximations
        
        zer=zeros(1,2*N+1);

        [l2dot, del2dot]=modsapprox(db, N, potentialType,etarmodCurrent,etapmodCurrent, delmodCurrent, lmodCurrent, kappa);

        [lexpdot, delexpdot] = modsapprox(db, N, potentialType,zer,zer, delstar, i*h*cstar, kappa);

        [llindot, dellindot] = modsapprox(db, N, potentialType,etarlinCurrent,etaplinCurrent, delstar, i*h*cstar, kappa);

        [lhdot, delhdot] = modsapprox(db, N, potentialType,etarhCurrent,etaphCurrent, delstar, i*h*cstar, kappa);
            
        [helpl1, helpdel1] = modsapprox(db, N, potentialType,etarhCurrent,etaphCurrent, delstar, i*h*cstar, kappa);
        [helpl2, helpdel2] = modsapprox(db, N, potentialType,zer,zer, delstar, i*h*cstar, kappa);
     
        [dldxi, ddeldxi] = dximodszero(db, N, potentialType, delstar, i*h*cstar, kappa);
        [dldc, ddeldc] = dcmodszero(db, N, potentialType, delstar, i*h*cstar, kappa); 
         
        helpl3 = lexpdot+(xi1Current-i*h*cstar)*dldxi+(delexpCurrent-delstar)*dldc; 
        helpdel3 = delexpdot+(xi1Current-i*h*cstar)*ddeldxi+(delexpCurrent-delstar)*ddeldc; 

        lhtestdot = helpl1 - helpl2 +helpl3;
        delhtestdot = helpdel1 - helpdel2 +helpdel3;

        [helpl1, helpdel1] = modsapprox(db, N, potentialType,etarlinCurrent,etaplinCurrent, delstar, i*h*cstar, kappa);

        llintestdot = helpl1 - helpl2 +helpl3;
        dellintestdot = helpdel1 - helpdel2 +helpdel3;

        % Updating
        delmodCurrent = delmodCurrent+h*deldot;
        lmodCurrent = lmodCurrent+h*ldot;

        xi1Current = xi1Current+h*((delexpCurrent-delstar)/24+lexpdot); 

        del2Current = del2Current+h*del2dot;
        l2Current = l2Current+h*l2dot;
       
        delexpCurrent= delexpCurrent+h*delexpdot;
        lexpCurrent= lexpCurrent+h*lexpdot;    

        dellinCurrent= dellinCurrent+h*dellindot;
        llinCurrent= llinCurrent+h*llindot;

        delhCurrent= delhCurrent+h*delhdot;
        lhCurrent= lhCurrent+h*lhdot;
       
        delhtestCurrent= delhtestCurrent+h*delhtestdot;
        lhtestCurrent= lhtestCurrent+h*lhtestdot;

        dellintestCurrent= dellintestCurrent+h*dellintestdot;
        llintestCurrent= llintestCurrent+h*llintestdot;

        
        etarmodCurrent = etarmodCurrent + h*(etar1 + 2*etar2 + 2*etar3 + etar4)/6;
        etapmodCurrent = etapmodCurrent + h*(etap1 + 2*etap2 + 2*etap3 + etap4)/6;

        etarlinCurrent = etarlinCurrent + h*(etarlin1 + 2*etarlin2 + 2*etarlin3 + etarlin4)/6;
        etaplinCurrent = etaplinCurrent + h*(etaplin1 + 2*etaplin2 + 2*etaplin3 + etaplin4)/6;

        etarhCurrent = etarhCurrent + h*(etarh1 + 2*etarh2 + 2*etarh3 + etarh4)/6;
        etaphCurrent = etaphCurrent + h*(etaph1 + 2*etaph2 + 2*etaph3 + etaph4)/6;
        
        rCurrent = rCurrent + h*(r1 + 2*r2 + 2*r3 + r4)/6;
        pCurrent = pCurrent + h*(p1 + 2*p2 + 2*p3 + p4)/6;
        
        
        % Save only every b-th timepoint
        if mod(i, b) == 0
            r(saveIdx, :) = rCurrent;
            p(saveIdx, :) = pCurrent;
            etarmod(saveIdx, :) = etarmodCurrent;
            etapmod(saveIdx, :) = etapmodCurrent;
            etarlin(saveIdx, :) = etarlinCurrent;
            etaplin(saveIdx, :) = etaplinCurrent;
            etarh(saveIdx, :) = etarhCurrent;
            etaph(saveIdx, :) = etaphCurrent;
            delmod(saveIdx) = delmodCurrent;
            lmod(saveIdx) = lmodCurrent;  
            del2(saveIdx) = del2Current;
            l2(saveIdx) = l2Current;
            dellin(saveIdx) = dellinCurrent;
            llin(saveIdx) = llinCurrent;
            delh(saveIdx) = delhCurrent;
            lh(saveIdx) = lhCurrent;
            delhtest(saveIdx) = delhtestCurrent;
            lhtest(saveIdx) = lhtestCurrent;
            dellintest(saveIdx) = dellintestCurrent;
            llintest(saveIdx) = llintestCurrent;
            delexp(saveIdx) = delexpCurrent;
            lexp(saveIdx) = lexpCurrent;

            saveIdx = saveIdx + 1;
        end
    end  
end
