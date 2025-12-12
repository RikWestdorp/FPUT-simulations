
clear; clc;
load('small.mat','db')

%% Configuration Section
time = 600; %Total timesteps
N = floor(time*1.2); % Number of lattice points (half)
cstar=1.014; % 
delstar=24*(cstar-1); %eps^2 in Friesecke Pego
h = 0.05; % Numerical time-step 
eps = 0.07;  % Forcing strength (take in [0,1/sqrt(3)))
springforce_variation_type = 'iid'; % Type of variation for spring force (constant, iid, sine, periodic)
noise_type = 'iid'; % Type of noise (uniform, iid)
potentialType = 'cubic'; % Choose potential type (quadratic, cubic)
FrameSkip = 160; % Shortens runtime of visualisation. Take low for accurate plots
t_num=floor(time/(FrameSkip*h))+1;
timegrid = (0:t_num-1) * FrameSkip * h;


%% Initializing

r0 = profile(db, 'r', cstar, 0, (-N:N));
p0 = profile(db, 'p', cstar, 0, (-N:N));


% Generate spring force variation ~ Unif(-\sqrt(3), sqrt(3))
kappa = eps*2*sqrt(3)*(rand(1, 2*N+1) - 0.5);


%% Run new term
M = 10;                % number of realisations (change as needed)
rng(12345);             % seed for reproducibility of kappas
KAPPA = eps*2*sqrt(3) * ( rand(M, 2*N+1) - 0.5 );  % pre-generate all kappas

% Preallocate storage
% rens and pens returned by RKensemble are (t_num x (2*N+1))

delexpens = zeros(M, t_num);   
lexpens   = zeros(M, t_num);
xi1ens = zeros(M, t_num);
delhtestens = zeros(M, t_num);
lhtestens = zeros(M, t_num);
delhtestnewens = zeros(M, t_num);
lhtestnewens = zeros(M, t_num);
E123ens = zeros(M, t_num);
E4ens = zeros(M, t_num);

%parpool('local');  
tic
for i = 1:M
    %try
        kappa_i = KAPPA(i, :);
        [delhtest, lhtest, delhtestnew, lhtestnew,delexp, lexp, xi1, E123, E4] =  RKgradientterm(db, h, kappa_i, time, potentialType, FrameSkip,delstar,N);
        delexpens(i, :) = delexp;
        lexpens(i, :) = lexp;
        delhtestens(i, :) = delhtest;
        lhtestens(i, :) = lhtest;
        delhtestnewens(i, :) = delhtestnew;
        lhtestnewens(i, :) = lhtestnew;
        xi1ens(i,:) = xi1;
        E123ens(i,:) = E123;
        E4ens(i,:) = E4;

    
    % catch ME
    %         fprintf('Iteration %d failed: %s\n', i, ME.message);
    % end
end
toc

%%

[alpha0, alpha1]= alpha(db, delstar,N);
[dcalpha0, dcalpha1]= dcalpha(db, delstar,N);
A=[0, alpha0;-alpha0, alpha1];
Ainv=inv(A);
dcA=[0, dcalpha0;-dcalpha0, dcalpha1];
dcAinv=inv(dcA);

t_num=floor(time/(FrameSkip*h))+1;
timegrid = (0:t_num-1) * FrameSkip * h;

cexp=zeros(1, t_num);

Y1=zeros(1, t_num);
Y2=zeros(1, t_num);

E1direct=zeros(1, t_num);
E2direct=zeros(1, t_num);
E3direct=zeros(1, t_num);
E4direct=zeros(1, t_num);

Int = zeros(2*N+1, t_num);
Int2 = zeros(2*N+1, t_num);
Int3 = zeros(2*N+1, t_num);


js=(-N:N);

% Loop over j
for jidx = 1:(2*N+1)
    j = js(jidx);
    r0=profile(db, 'r', cstar,0 , j+zeros(1,length(timegrid)));
    rtime=profile(db, 'r', cstar,0 , j-cstar*timegrid);
    dcrtime=profile(db, 'dcr', cstar,0 , j-cstar*timegrid);
    dxrtime=profile(db, 'dxr', cstar,0, j-cstar*timegrid);
    % Numerical integration
    Int(jidx,:) = cumsum(rtime.*dcrtime) * h;
    Int2(jidx,:) = 2*Int(jidx,:);
    Int3(jidx,:) = cumsum(rtime.^2-r0.^2) * h;    
end



r0=profile(db, 'r', cstar, 0, (-N:N));
for tidx=1:length(timegrid)
    t=timegrid(tidx);

    r=profile(db, 'r', cstar, cstar*t, (-N:N));
    dxr=profile(db, 'dxr', cstar, cstar*t, (-N:N));
    dxdxr=profile(db, 'dxdxr', cstar, cstar*t, (-N:N));
    dcr=profile(db, 'dcr', cstar, cstar*t, (-N:N));
    dcdcr=profile(db, 'dcdcr', cstar, cstar*t, (-N:N));
    dcdxr=profile(db, 'dcdxr', cstar, cstar*t, (-N:N));
   
    cexp(tidx)=cstar-1/(24*2*cstar*alpha0)*sum(KAPPA(1, :).*(r.^2-r0.^2));
    Y1(tidx)=-alpha1/(2*cstar*alpha0.^2)*sum(KAPPA(1, :).*(r.^2-r0.^2));
    Y2(tidx)=1/(alpha0)*sum(KAPPA(1, :).*Int(:,tidx)');
    
    S1 = 2*[ dxr.^2+r.*dxdxr; -dcr.*dxr-r.*dcdxr];
    S4 = [-2*dxr.*r; 2*dcr.*r];
    dcS4 = 2*[ -dcr.*dxr-r.*dcdxr; dcr.^2+r.*dcdcr];
    temp1 = Ainv*S1;
    temp4 = Ainv*dcS4+dcAinv*S4;

    E1direct(tidx)=alpha1/(24*4*cstar*alpha0^2)*sum(temp1(2,:).*(r.^2-r0.^2));
    E2direct(tidx)=-1/(24*4*alpha0)*sum(temp1(2,:).*Int2(:,tidx)');
    E3direct(tidx)=1/(24*4*cstar*alpha0)*sum(temp1(2,:).*Int3(:,tidx)');
    E4direct(tidx)=1/(24*4*cstar*alpha0)*sum(temp4(2,:).*(r.^2-r0.^2));
end


%% Figures gradient term

set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultLegendInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(groot, 'defaultAxesFontName', 'Times New Roman') % optional: LaTeX-like font
set(groot, 'defaultAxesFontSize', 15)

figmean = figure;
figmean.Position = [300, 300, 600, 300];
plot(timegrid,mean(1+delhtestens/24), 'b', 'LineWidth',1.2);
hold on
plot(timegrid,mean(1+delhtestnewens/24), 'r', 'LineWidth',1.2);
hold on
plot(timegrid,cstar+mean(E123ens+E4ens)/24,'g')
hold on
plot(timegrid,cstar+eps^2*(E1direct(end)+E2direct(end)+E3direct(end))*timegrid,'r')
xlabel('$t$')
title('mean of correction term')


figmean = figure;
figmean.Position = [300, 300, 600, 300];
plot(timegrid,cstar+mean(E4ens)/24,'g')
hold on
plot(timegrid,cstar+eps^2*E4direct(end)*timegrid,'r')
title('E4')

figmean = figure;
figmean.Position = [300, 300, 600, 300];
plot(timegrid,cstar+mean(E123ens)/24,'g')
hold on
plot(timegrid,cstar+eps^2*(E1direct(end)+E2direct(end)+E3direct(end))*timegrid,'r')
title('E123')

figmean = figure;
figmean.Position = [300, 300, 600, 300];
plot(timegrid,1+delexpens(1,:)/24, 'b', 'LineWidth',1.2);
hold on
plot(timegrid,cexp, 'r', 'LineWidth',1.2);
xlabel('$t$')
title('$c_*+\sigma c_{1}$')

figmean = figure;
figmean.Position = [300, 300, 600, 300];
plot(timegrid,lexpens(1,:)-cstar*timegrid, 'b', 'LineWidth',1.2);
hold on
plot(timegrid,Y1+Y2, 'r', 'LineWidth',1.2);
xlabel('$t$')
title('$\sigma\gamma_{1}$')

figmean = figure;
figmean.Position = [300, 300, 600, 300];
plot(timegrid,xi1(1,:)-cstar*timegrid, 'b', 'LineWidth',1.2);
xlabel('$t$')
title('$\sigma\xi_{1}$')
