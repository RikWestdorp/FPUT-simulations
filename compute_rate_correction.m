function rate_correction = compute_rate_correction(db, cstar, time, h)

N = floor(time*1.2); % Number of lattice points (half)
delstar=24*(cstar-1); %eps^2 in Friesecke Pego
eps = 0.07;  % Forcing strength (take in [0,1/sqrt(3)))
springforce_variation_type = 'iid'; % Type of variation for spring force (constant, iid, sine, periodic)
noise_type = 'iid'; % Type of noise (uniform, iid)
potentialType = 'cubic'; % Choose potential type (quadratic, cubic)
FrameSkip = 1; % Shortens runtime of visualisation. Take low for accurate plots
t_num=floor(time/(FrameSkip*h))+1;
timegrid = (0:t_num-1) * FrameSkip * h;

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
    %Int(jidx,:) = cumsum(rtime.*dcrtime) * h;
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

rate_correction = E1direct(end)+E2direct(end)+E3direct(end)+E3direct(end); 
end