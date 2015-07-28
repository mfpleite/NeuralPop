%% Test Parameters - From NeuroElectro.com
%% Priors
% --------------------
% --- Expectations ---
% --------------------
close all; clear all;
P.ang    = pi;

% --- Conectivity matrices ---

P.A(1).M = [ -32  log(.01); % excitatory
             -32  -32];
    
P.A(2).M = [-32 -32;       % inhibitory
            log(.5) -32];
        
P.A(3).M = [ -32  -32;     % excitatory
             -32  log(.0005)];
    
P.A(4).M = [log(1) -32;  % inhibitory
            -32 -32];
        
P.C(1).M = [log(4) log(4)]'; % excitatory exogenous input
P.C(2).M = [log(.1) log(.1)]'; % inhibitory exogenous input

P.C(3).M = [-32 -32]';

P.C(4).M = [-32 -32]';

% --- CA1 inhibitory interneurons ---
P.P(1).C   = log(.084);                 % [nF] Cell capacitance
P.P(1).Vl  =    -62;                    % [mV] Leak reversal portential
P.P(1).gl  = P.P(1).C -log(.012);       % [nF/s] = [nS] Leak conductance
P.P(1).Vt  = -55;                       % [mV] Spike threshold voltage
P.P(1).Vr  = -70;                       % [mV] Reset voltage
P.P(1).Rt  = log(1/1000);               % [s]  Refractory time
P.P(1).Vg  = [0,-70,0,-70];             % [mV] Reversal portential of the families of synaptic channels
P.P(1).T   = log([16,5,16,5]/1000);   % [s] Decay time constants for the synaptic conductances
P.P(1).D   = log(2000);                 % (1/2*sigma^2) [(mV^2)/ms] Membrane noise term
P.P(1).d   = log(1/1000);

% --- CA1 excitatory pyramidal cells ---
P.P(2).C   = log(.150);                 % [nF] Cell capacitance 
P.P(2).Vl  =    -65;                    % [mV] Leak reversal portential
P.P(2).gl  = P.P(2).C - log(0.030);     % [nF/s] = [nS] Leak conductance
P.P(2).Vt  = -48;                       % [mV] Spike threshold voltage
P.P(2).Vr  = -70;                       % [mV] Reset voltage
P.P(2).Rt  = log(8/1000);               % [s]  Refractory time
P.P(2).Vg  = [0,-70,0,-70];             % [mV] Reversal portential of the families of synaptic channels 
P.P(2).T   = log([16,5,16,5]/1000);   % [s] Decay time constants for the synaptic conductances
P.P(2).D   = log(4000);                 % (1/2*sigma^2) [(mV^2)/s] Membrane noise term
P.P(2).d   = log(1/1000);

% -------------------
% --- Covariances ---
% -------------------

CP.ang    = 10*pi;

CP.A(1).M = exp([ -32 4;
              -32 -32]);
    
CP.A(2).M = exp([-32 -32;
             4 -32]);

CP.A(3).M = exp([ -32 -32;
              -32 4]);
    
CP.A(4).M = exp([4 -32;
             -32 -32]);

CP.C(1).M = exp([2 4]');
CP.C(2).M = exp([2 2]');

CP.C(3).M = exp([-32 -32]');

CP.C(4).M = exp([-32 -32]');

% --- inhibitory interneurons ---
CP.P(1).C   = log(8);           % [nF] Cell capacitance
CP.P(1).Vl  =    3;             % [mV] Leak reversal portential
CP.P(1).gl  = log(8);           % [pF/ms] = [nS] Leak conductance
CP.P(1).Vt  = 100;              % [mV] Spike threshold voltage
CP.P(1).Vr  = exp(-32);         % [mV] Reset voltage
CP.P(1).Rt  = log(4);           % [s]  Refractory time
CP.P(1).Vg  = [100,100,100,100];% [mV] Reversal portential of the families of synaptic channels
CP.P(1).T   = log([8,8,8,8]);   % [s] Decay time constants for the synaptic conductances
CP.P(1).D   = log(16);           % (1/2*sigma^2) [(mV^2)/ms] Membrane noise term

% --- excitatory pyramidal cells ---
CP.P(2).C   = log(8);           % [nF] Cell capacitance
CP.P(2).Vl  =    4;           % [mV] Leak reversal portential
CP.P(2).gl  = log(8);           % [pF/ms] = [nS] Leak conductance
CP.P(2).Vt  = 50;              % [mV] Spike threshold voltage
CP.P(2).Vr  = exp(-32);              % [mV] Reset voltage
CP.P(2).Rt  = log(4);           % [s]  Refractory time
CP.P(2).Vg  = [100,100,100,100];% [mV] Reversal portential of the families of synaptic channels
CP.P(2).T   = log([8,8,8,8]);   % [s] Decay time constants for the synaptic conductances
CP.P(2).D   = log(16);           % (1/2*sigma^2) [(mV^2)/ms] Membrane noise term

%%
M.opt.N              = 4000; % Maximum time points for time integration
M.opt.dttol          = 1e-3;
M.opt.LCtol          = 1e-3;
M.opt.r              = 1;
M.opt.popCurrents    = [2 10; 2 -70; 1 10; 1 -70]; % Pop and Voltage clamp values [mV]
M.opt.popFiringRates = [2; 1]; % Pop for fring rates
M.opt.dpsize         = 10^-6;  % For DCM 
M.opt.svdtol         = 10^-6;  % For LC perturbation
M.opt.Drawplots      = 0;      % Progress plots of pdfs and conductances
M.opt.vmin           = -16;    % For DCM
M.opt.vmax           = 8;      % For DCM
M.opt.vini           = -4;     % For DCM
M.opt.Nmax           = 32;     % For DCM
M.opt.sigma          = 0;      % For DCM

Vp = 0*spm_vec(P); Vp = [Vp Vp]; Vp(4,1) = 1; Vp(4,2) = -1; % For DCM tests
% Vp     = spm_svd(diag(spm_vec(CP)),exp(-32));
%% Data 

M.pE = P;
M.pC = CP;
M.hE = 0*[3 9 -3 -3 3 3 3 9]' + 15;
M.hC = eye(length(M.hE))*.01;
M.Ep = P;

M.opt.Ncoefs = 100;     % For DCM
M.converged = 0;        % For DCM
M.F = -Inf;             % For DCM
M.opt.LCtol = .1;       % Convergence tolerance for LC
M.opt.Drawprogress = 1; % For DCM

% if ~isfield(M.opt,'dttol'); M.opt.dttol = 10^-3; end;
% if ~isfield(M.opt,'dt'); M.opt.dt = .001; end;
% if ~isfield(M.opt,'T0max'); M.opt.T0max = .12; end; 
% if ~isfield(M.opt,'svdtol'); M.opt.svdtol = 10^-6; end; 
% if ~isfield(M.opt,'dpsize'); M.opt.dpsize = exp(-3); end;
% if ~isfield(M.opt,'Ncoefs'); M.opt.Ncoefs = 1000; end;
% if ~isfield(M.opt,'LCtol'); M.opt.LCtol = 10^-3; end;
% M.Gv0 = zeros(M.np,M.np); 
% M.LFP      = zeros(M.np,M.opt.N);
% M.Currents = zeros(size(M.opt.popCurrents,1),M.opt.N);
% 
% P   = spm_ExpP(Ep);
% M   = spm_lifpopsys_LC_prepare(P,M);
% M.opt.Drawplots = 1;
% M = spm_DCM_lifpopsys_LC_int(P,M);
% 
% clear Ang I n phaseshift y y0 y2 yf yf0 ans Vp CP P

%% Parameter Initialization
% 
% M.P.C(1).M = [log(.1) log(8)]';
% M.P.C(2).M = [log(.1) log(.1)]';
% M.P.C(3).M = [-32 -32]';
% M.P.C(4).M = [-32 -32]';
% 

P0  = P;                % Save original parameters
%%
P   = spm_ExpP(P0);     % Adapt parameters for integration
M   = spm_lifpopsys_LC_prepare(P,M); % Complete M structure
cflag = 0;              % convergence flag for Fixed point search

%% --- Coment this section if modeling more than 2 populations ---
[~, ~, M, cflag] = fx_LIFpopMEJpar(P,M); % Try to find a fixed point
% --- end of comented section ---

if cflag
    M = spm_perturb_fp(P,M); % Kick the system out of the fixed point
end

% M.opt.Drawplots = 1; % Option for drawing progress plots
% 
% % This will integrate the model until a limit cycle (or Nmax) is reached 
% if ~cflag || max(real(M.J.S))>10e-4
%     
%     M = spm_DCM_lifpopsys_LC_int(P,M);
%     
%     ts  = timeseries( squeeze(M.Currents(1,(1:M.opt.Kend))), M.opt.t(1:M.opt.Kend));
%     ts1 = resample(ts,0:.001:M.opt.t(M.opt.Kend)); % resample to 1KHz
%     figure
%     plot(ts1)
%     fprintf('LC frequency: %f Hz \n', 1/M.opt.T0)
% else
%     fprintf('Stable Fixed Point')
% end

% The structure M will have the results of the integration
% figure
% TFanalysis(squeeze(ts1.data)',40,200,16,500,1000,.000001);
%%

M.opt.Drawplots = 0;
M.opt.dttol          = 1e-4;
% M.opt.LCconvTol = .1*M.opt.dttol;
tic
M = spm_DCM_lifpopsys_LC_int(P,M);
fprintf('LC frequency: %f Hz \n', 1/M.opt.T0)
toc

%%

try
    k = M.opt.Kend; 
catch
    k = 2;
end
k = 2;
nstates = 0;
sumS = zeros(M.np,1); 
clear S
for l = 1:M.np   % cycle through populations
    P1      = M.P(l);
    SSS2    = [reshape(M.SS(l).SS(:,k-1),[],1,1); squeeze(M.GV(l,:,k-1))'];
    S(nstates+1 : nstates + M.P(l).LVV + M.nc) = SSS2;
    nstates = nstates + M.P(l).LVV + M.nc;
    sumS(l) = 1/M.P(l).Vres; % sum constrain on pdf's
end

y = S';

odefun = @(a,b) fx_LIFpopMEJparOdeFun(0,b,P,M);
opt=[];
% opt.RelTol = M.opt.dttol*1;
% opt.RelTol = 10^(-0);
tic
[T,Y] = ode15s(odefun,[0 M.opt.t(M.opt.Kend)],y,opt);
% [T,Y] = ode45(odefun,[0 M.opt.t(M.opt.Kend)],y,opt);
% [T,Y] = ode23tb(odefun,[0 M.opt.t(M.opt.Kend)],y,opt);
toc


%%
figure;
plot(T,Y(:,M.P(1).LVV+1),'c.')
hold on
plot(M.opt.t(1:M.opt.Kend),squeeze(M.GV(1,1,1:M.opt.Kend)),'.')
%
% plot(T1,Y1(:,M1.P(1).LVV+1),'y.')
% plot(M1.opt.t(1:M1.opt.Kend),squeeze(M1.GV(1,1,1:M1.opt.Kend)),'k.')



