function [f, J] = fx_LIFpopMEJparOdeFun(t,y,P,M)
% Full flow for a system of interconnected LIFpop
% G(:,:,1) = GE; G(:,:,2) = GI; U =  U(:,1)*C'; S = [SSoverT(:,:,kini-1); squeeze(GV(:,:,kini-1))];

nstates = 0;
for l = 1:M.np
    nstates = nstates + M.P(l).LVV + M.nc;
end

% S  = zeros(nstates,1);
S = y';
f  = zeros(nstates,1);
J  = zeros(nstates);
Gout = zeros(M.nc,M.np);


for l = 1:M.nc % Compute recieved conductances
    
%     Gout(l,:) = P.A(l).M*squeeze(M.GV(:,l,k-1)) + P.C(l).M;
    Gout(l,:) = P.A(l).M*squeeze(S(M.P(1).LVV+l:( M.P(1).LVV+M.nc):end))' + P.C(l).M;

    
end
    
nstates = 0;
for k = 1:M.np
    ns = M.P(k).LVV + M.nc;
    J(nstates+1:nstates+ns, nstates+1:nstates+ns) = fx_LIFpopME(Gout(:,k),M.P(k));
    nstates = nstates + M.P(k).LVV + M.nc;
end

f = J*S(:);

nstatesk = 0;

% - Nonlinear part (implemented for 2 pop system for now)
for k = 1:M.np
    nsk = M.P(k).LVV; 
    nstates = 0;
    for l = 1:M.np
        ns = M.P(l).LVV;
        for g = 1:M.nc
            
            dJ0  = (P.A(g).M(l,k)/M.P(l).Vres)*M.P(l).FvarV(:,g).*S(nstates+1:nstates+ns)';
         
            J(nstates+1:nstates+ns ,nstatesk+nsk+g) = ...
                        + [dJ0(2:end); 0] - [-dJ0(1); dJ0(1:end-1)];
        
        end
        nstates = nstates + ns + M.nc;
    end
    nstatesk = nstatesk + nsk + M.nc;
end



end
