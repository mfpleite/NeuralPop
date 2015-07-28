function M = spm_DCM_lifpopsys_LC_int(P,M,InlineProgFlag)

if nargin<3
    InlineProgFlag = 1;
end

if InlineProgFlag
    fprintf('\n Burn in:')
end

if ~isfield(M.opt,'LCconvTol'); M.opt.LCconvTol = 10*M.opt.dttol; end;
       
tol = M.opt.dttol;% tol for burn in
k  =1;
Gout       = zeros(M.nc,M.np);
t0 = 0;
r0 = Inf;

for l = 1:M.np   % cycle through populations
    M.SS(l).SS(:,1) = M.SS(l).SS(:,k);
    M.GV(l,:,1) = M.GV(l,:,k);
end

r =  Inf*ones(M.opt.N,1);
dt = M.opt.dt;
M.opt.t = zeros(M.opt.N,1);

while k < M.opt.N && r(k)>M.opt.LCconvTol % integrate through time
    
    % --- Prediction step ---
    
    k=k+1;
    for l = 1:M.nc % Compute recieved conductances
        Gout(l,:) = P.A(l).M*squeeze(M.GV(:,l,k-1)) + P.C(l).M;
    end
    
    for l = 1:M.np   % cycle through populations
        
        P1      = M.P(l);
        SSS2    = [reshape(M.SS(l).SS(:,k-1),[],1,1); squeeze(M.GV(l,:,k-1))'];
        [Qdpdt] = fx_LIFpopME(Gout(:,l),P1);
        M1      = expm(.5*dt*Qdpdt);
        SSS     = M1*SSS2; % Midpoint
        SSS1    = M1*SSS;  % Endpoint
        
        M.SS(l).SS(:,k) = SSS(1:P1.LVV);
        M.SS0(l).SS = SSS1(1:P1.LVV);
        
        M.GV(l,:,k)     = SSS(P1.LVV+1:end);
        M.GV0(l,:)      = SSS1(P1.LVV+1:end);
        
    end
    
    % --- Correction step ---
  
    for l = 1:M.nc % Compute recieved conductances
        Gout(l,:) = P.A(l).M*(squeeze(M.GV(:,l,k))) + P.C(l).M; % Midpoint
    end
    
    for l = 1:M.np   % cycle through population sexp(-6)
        P1      = M.P(l);
        SSS2    = [reshape(M.SS(l).SS(:,k-1),[],1,1); squeeze(M.GV(l,:,k-1))'];
        [Qdpdt] = fx_LIFpopME(Gout(:,l),P1);
        SSS     = expm(dt*Qdpdt)*SSS2;
        
        M.SS(l).SS(:,k) = SSS(1:P1.LVV);
        M.GV(l,:,k)     = SSS(P1.LVV+1:end);
    end
    
    M.opt.t(k) = M.opt.t(k-1) + dt;
    
    D = (M.GV0 - squeeze(M.GV(:,:,k))).^2;
    D = sum(D(:));
    
    for l=1:M.np
        D = D + sum((M.SS0(l).SS - squeeze(M.SS(l).SS(:,k))).^2);
    end
    
    
    
    if D >tol
        dt = 2^(-1/2)*dt;
        k  = k - 1;
        %             disp(['timestep -: ' num2str(dt)])
%         fprintf('-')
    elseif D < tol/8 %&& dt<M.opt.dt*2
        dt = 2^(1/4)*dt;
        if InlineProgFlag, inlineprogress(k,M.opt.N); end
        %             disp(['timestep +: ' num2str(dt)])
%         fprintf('+')
    elseif InlineProgFlag
%         dt = .9*dt*(tol/D);
        inlineprogress(k,M.opt.N)
%         fprintf('*')
    end
    
    if ~(sum(D(:)) > tol)
        if isfield(M,'J')
            lambmax  = real(M.J.S(M.J.IX(end)));
            expectr0 = lambmax*exp(lambmax*M.opt.t(k)); %avoid early stoping
        else
            expectr0 = Inf;
        end
        if M.opt.t(k) > .005 && ~mod(k,10)
            [r0, t0] = convergence_check(M.GV,M.opt.t,k,M.opt.T0max);
            if r0 > .1*expectr0, 
                r0 = Inf; 
            end
            r(k) = r0;
        end
        % -----------------------
        % --- LFP computation ---
        % -----------------------
        
        if ~isfield(M.opt,'wLFP'); M.opt.wLFP = ones(size(M.P(1).T)); end;
        if ~isfield(M.opt,'LocalElectrode'); M.opt.LocalElectrode = 0; end;
        
        for l = 1:M.np
            P1 = M.P(l);
            FV = P1.FvarV;
            FV(P1.Tr+1:end,:) = FV(repmat(P1.R,P1.LQ,1),:); % heuristic on refractory period
            M.LFP(l,k) = M.opt.wLFP*(Gout(:,l).*(FV'*squeeze(M.SS(l).SS(:,k)*P1.Vres))); % [nS*mV] = [pA] in pico Ampere (per neuron)
            if M.opt.LocalElectrode
                M.LFP(l,k) = M.LFP(l,k) + squeeze(M.SS(l).SS(P1.Tr+1,k)*P1.Vres)*P1.C*(P1.Vt-P1.Vr)/.001; % Current for repolarization in 1 ms (per neuron)
            end
        end
        
        % --- Current computation ---
        
        for m= 1:size(M.opt.popCurrents,1)
            P1              = M.P(M.opt.popCurrents(m,1));
            M.Currents(m,k) = - P1.gl*(P1.Vl - M.opt.popCurrents(m,2)); % [nS*mV] = [pA] in pico Ampere (per neuron)
            for l = 1:M.nc
                M.Currents(m,k) = M.Currents(m,k) - Gout(l,M.opt.popCurrents(m,1))*(P1.Vg(l) - M.opt.popCurrents(m,2)); % [nS*mV] = [pA] in pico Ampere (per neuron)
            end
        end
        
        
        if ~isfield(M.opt,'pop'); M.opt.pop = 2; end;
        pop = M.opt.pop;
        
        % --- Plots ---        
        
        if M.opt.Drawplots
            
            for l = 1:M.np   % cycle through populationsexp(-6)
                P1 = M.P(l);
                subplot(M.np+2,1,l)
                plot(P1.VV,(M.SS(l).SS(:,k)));axis([P1.VV(1) P1.VV(end) -.001 .2]);
                xlabel('Voltage space (mV)');
                ylabel('Probability density (mv^-^1)');
                drawnow;
            end
            
            
            subplot(M.np+2,1,M.np+1)
            plot(M.opt.t(1:k),M.LFP(pop,1:k)','k');
            hold on
            plot(M.opt.t,M.Currents','.','MarkerSize',1)
            hold off
            xlabel('Time (s)');
            ylabel('Currents (pA)')
            
            subplot(M.np+2,1,M.np+2)
            plot(- M.Currents(2,2:k), M.Currents(1,2:k),'.','MarkerSize',1);
            xlabel('-E (pA)'); %xlim([for l = 1:M.np
            ylabel('I (pA)'); %ylim([-100 100]+EIcurrents(k,2));
            drawnow;
        end
    end
    
    
end

M.opt.Kend = k;
M.opt.T0 = M.opt.t(k)-t0;



end
