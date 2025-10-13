% main.m
clear; clc;
parentFolder = fileparts(pwd);  % Get parent directory of current folder
addpath(parentFolder);         % Add it to the path
% Load and preprocess waveprofiles
filename = 'all_sols_out_e.txt';
[c_groups, c_values] = load_phi_data(filename);
db = preprocess_phi_data(c_groups, c_values);

%% Configuration Section
time = 600; %Total timesteps
N = floor(time*1.2); % Number of lattice points (half)
cstar=1.014; % 
delstar=24*(cstar-1); %eps^2 in Friesecke Pego
h = 0.5; % Numerical time-step 
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

%% Run ensemble
M = 2000;                % number of realisations (change as needed)
rng(12345);             % seed for reproducibility of kappas
KAPPA = eps*2*sqrt(3) * ( rand(M, 2*N+1) - 0.5 );  % pre-generate all kappas

% Preallocate storage
% rens and pens returned by RKensemble are (t_num x (2*N+1))
r_ens    = zeros(M, t_num, 2*N+1);
p_ens    = zeros(M, t_num, 2*N+1);
del0_ens = zeros(M, t_num);   
l0_ens   = zeros(M, t_num);
cfit_ens = zeros(M, t_num);
lfit_ens = zeros(M, t_num);

%parpool('local');  
tic
parfor i = 1:M
    try
        kappa_i = KAPPA(i, :);
        [r, p] = RKensemble(db, r0, p0, h, kappa_i, time, potentialType, FrameSkip, delstar, N);
        r_ens(i, :, :) = r;
        p_ens(i, :, :) = p;
       
    
        % compute fits for each time slice
        censfit_local = zeros(1, t_num);
        lensfit_local = zeros(1, t_num);
    
        for t = 1:t_num
            r_t = r(t, :);
            p_t = p(t, :);
            [c_val, l_val] = fit_ortho(db, N, r_t, p_t);
            censfit_local(t) = c_val;
            lensfit_local(t) = l_val;
        end
    
        % write local results into sliced arrays
        cfit_ens(i, :) = censfit_local;
        lfit_ens(i, :) = lensfit_local;
    catch ME
            fprintf('Iteration %d failed: %s\n', i, ME.message);
    end
end
toc



filename = 'ensemble035.txt';
fid = fopen(filename, 'w');

% Write metadata
fprintf(fid, 'Metadata:\n');
fprintf(fid, 'time = %.0f\n', time);
fprintf(fid, 'cstar = %.4f\n', cstar);
fprintf(fid, 'h = %.5f\n', h);
fprintf(fid, 'eps = %.4f\n', eps);
fprintf(fid, 'FrameSkip = %d\n', FrameSkip);
fprintf(fid, '\n');

write_matrix(fid,'cfit_ens', cfit_ens);
write_matrix(fid,'lfit_ens', lfit_ens);




%% Pathwise Simulation


tic
[r, p, etarmod, etapmod, etarlin, etaplin,etarh, etaph, delmod, lmod, del2, l2, dellin, llin,delh, lh,delhtest, lhtest,dellintest, llintest, delexp, lexp] = RK(db, r0, p0, h, kappa, time, potentialType, FrameSkip, delstar,N);
toc

cfit = zeros(1,t_num);
lfit = zeros(1,t_num);
etarfit = zeros(t_num, 2*N+1);
etapfit = zeros(t_num, 2*N+1);

for t = 1:t_num
    r_t = r(t, :);
    p_t = p(t, :);

    [cfit(t), lfit(t)] = fit_ortho(db, N, r_t, p_t);

    etarfit(t, :) = r_t - profile(db, 'r', cfit(t), lfit(t), -N:N);
    etapfit(t, :) = p_t - profile(db, 'p', cfit(t), lfit(t), -N:N);
end

filename = 'simulation_for_ode014.txt';
fid = fopen(filename, 'w');

% Write metadata
fprintf(fid, 'Metadata:\n');
fprintf(fid, 'time = %.0f\n', time);
fprintf(fid, 'cstar = %.4f\n', cstar);
fprintf(fid, 'h = %.5f\n', h);
fprintf(fid, 'eps = %.4f\n', eps);
fprintf(fid, 'FrameSkip = %d\n', FrameSkip);
fprintf(fid, '\n');

function write_matrix(fid, name, mat)
    fprintf(fid, '%s:\n', name);
    fmt = [repmat('% .6e ', 1, size(mat,2)) '\n'];
    for i = 1:size(mat,1)
        fprintf(fid, fmt, mat(i,:));
    end
    fprintf(fid, '\n');
end

% Write main variables
write_matrix(fid,'r', r);
write_matrix(fid,'p', p);
write_matrix(fid,'etarmod', etarmod);
write_matrix(fid,'etapmod', etapmod);
write_matrix(fid,'etarlin', etarlin);
write_matrix(fid,'etaplin', etaplin);
write_matrix(fid,'etarh', etarh);
write_matrix(fid,'etaph', etaph);
write_matrix(fid,'delmod', delmod);
write_matrix(fid,'lmod', lmod);
write_matrix(fid,'del2', del2);
write_matrix(fid,'l2', l2);
write_matrix(fid,'dellin', dellin);
write_matrix(fid,'llin', llin);
write_matrix(fid,'delh', delh);
write_matrix(fid,'lh', lh);
write_matrix(fid,'delhtest', delhtest);
write_matrix(fid,'lhtest', lhtest);
write_matrix(fid,'dellintest', dellintest);
write_matrix(fid,'llintest', llintest);
write_matrix(fid,'delexp', delexp);
write_matrix(fid,'lexp', lexp);

% Write fitted data
write_matrix(fid,'cfit', cfit);
write_matrix(fid,'lfit', lfit);
write_matrix(fid,'etarfit', etarfit);
write_matrix(fid,'etapfit', etapfit);

fclose(fid);