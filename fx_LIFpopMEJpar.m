function [f, J, M, cflag] = fx_LIFpopMEJpar(P,M)
% Full flow for a system of interconnected LIFpop
% G(:,:,1) = GE; G(:,:,2) = GI; U =  U(:,1)*C'; S = [SSoverT(:,:,kini-1); squeeze(GV(:,:,kini-1))];

cflag =0;
try
    k = M.opt.Kend; 
catch
    k = 2;
end

nstates = 0;
for l = 1:M.np
    nstates = nstates + M.P(l).LVV + M.nc;
end

S  = zeros(nstates,1);
f  = zeros(nstates,1);
J  = zeros(nstates);
Gout = zeros(M.nc,M.np);


for l = 1:M.nc % Compute recieved conductances
    Gout(l,:) = P.A(l).M*squeeze(M.GV(:,l,k-1)) + P.C(l).M;
end
    

nstates = 0;
sumS = zeros(M.np,1); 
% - Linear part
for l = 1:M.np   % cycle through populations
    P1      = M.P(l);
    SSS2    = [reshape(M.SS(l).SS(:,k-1),[],1,1); squeeze(M.GV(l,:,k-1))'];
    S(nstates+1 : nstates + M.P(l).LVV + M.nc) = SSS2;
    nstates = nstates + M.P(l).LVV + M.nc;
    sumS(l) = 1/M.P(l).Vres; % sum constrain on pdf's
end


nstates = 0;
for k = 1:M.np
    ns = M.P(k).LVV + M.nc;
    J(nstates+1:nstates+ns, nstates+1:nstates+ns) = fx_LIFpopME(Gout(:,k),M.P(k));
    nstates = nstates + M.P(l).LVV + M.nc;
end

f = J*S(:);
nstatesk = 0;
pvars = zeros(M.np,length(f));
% - Nonlinear part 
for k = 1:M.np
    nsk = M.P(l).LVV; 
    nstates = 0;
    for l = 1:M.np
        ns = M.P(l).LVV;
        for g = 1:M.nc
            
            dJ0  = (P.A(g).M(l,k)/M.P(l).Vres)*M.P(l).FvarV(:,g).*S(nstates+1:nstates+ns);
         
            J(nstates+1:nstates+ns ,nstatesk+nsk+g) = ...
                        + [dJ0(2:end); 0] - [-dJ0(1); dJ0(1:end-1)];
        
        end
        nstates = nstates + ns + M.nc;
    end
    pvars(k,nstatesk+1:nstatesk+nsk) = 1;
    nstatesk = nstatesk + nsk + M.nc;
end


if nargout<3, return, end

%%
iterN = 0;
f0 = Inf;
Sreset = S;
S0 = S;
dS = Inf;
fdiverge = zeros(1,10);
if ~isfield(M.opt,'wRoot'); M.opt.wRoot = 1; end;
%
while (~cflag) && (iterN<150)  %norm(dS) > 10^-6 ||
%     disp([norm(dS) norm(f)]); plot(S0)
%     pause
    iterN = iterN +1;
    fprintf('.'); if ~mod(iterN,25), fprintf('\n'); end
    [f, J] = fx_LIFpopMEJpar(P,M);
    % constrain x sum to remain constant! otherwise J is singular
    f1 = [J*S0 - f;sumS]; % review this line
    J1 = [J;pvars];
    S = J1\f1; %lsqnonneg(J1,f1); % positively constrained implementation to eschew singular Js
    dS = S - S0;
%     dS = - J\f;
    S0 = S0(:) + M.opt.wRoot*dS;
    

    nstates = 0;
    for l = 1:M.np   % cycle through populations
        M.SS(l).SS(:,1) = S0(nstates+1:nstates+M.P(l).LVV);
        M.GV(l,:,1)     = S0(nstates+M.P(l).LVV+1:nstates+M.P(l).LVV+M.nc);
        nstates = nstates + M.P(l).LVV + M.nc;
    end
    
    fdiverge = 0*[fdiverge(2:end) norm(f)<norm(f0)]; %is the method diverging?
    f0 =f;
    
    if prod(fdiverge)
        M.opt.wRoot  = .5*M.opt.wRoot;
         nstates = 0;
         for l = 1:M.np   % cycle through populations
             M.SS(l).SS(:,1) = Sreset(nstates+1:nstates+M.P(l).LVV);
             M.GV(l,:,1)     = Sreset(nstates+M.P(l).LVV+1:nstates+M.P(l).LVV+M.nc);
             nstates = nstates + M.P(l).LVV + M.nc;
         end
        [f, J, M] = fx_LIFpopMEJpar(P,M);
        M.opt.dttol  = M.opt.dttol*10;
        break
    end
    
%     plot(S,'r')
%     hold on
%     plot(S0)
%     hold off
%     drawnow
    if ~all(isfinite(S)), break;end % stop if the method is badly conditioned
    cflag = (norm(f) < 10^-6);
end

fprintf('%i ',iterN);%fprintf('%d',norm(f))

% Check eigen values

if cflag
    s = size(J);
    [V S1] = eig(full(J(1:s(2),1:s(2))));
    [B,IX] = sort(real(diag(S1)));
    M.J.V = V;
    M.J.S = diag(S1);
    M.J.IX = IX;
else % reset M to original state
    nstates = 0;
    for l = 1:M.np   % cycle through populations
        M.SS(l).SS(:,1) = Sreset(nstates+1:nstates+M.P(l).LVV);
        M.GV(l,:,1)     = Sreset(nstates+M.P(l).LVV+1:nstates+M.P(l).LVV+M.nc);
        nstates = nstates + M.P(l).LVV + M.nc;
    end
end

end
