function [Qdpdt] = fx_LIFpopME(Go,P)

% optimization is possible by avoiding the sum of all matrix entries
% Matlab sparse implementation takes longer for typical P.LVV values

m = P.LVV;
n = length(P.G);
T = P.Tr;
R = P.R;
Fsteady = P.Fsteady;
FvarV   = P.FvarV;

FvarG = zeros(size(Fsteady)); 
for k = 1:n
    FvarG(:) = FvarG(:) - (Go(k)*P.G(k) + P.DGL(k)*P.DL(k))*FvarV(:,k);
end

F = (Fsteady + FvarG);

% FS1  = diag(F(:));% sparse(1:m,1:m,F(:),m,m);
F    = F/P.Vres;
D1   = diag(F(2:end),1); % spdiags(ones(m,1),1,m,m); 
D1(P.Tr+1:end,P.Tr+1:end) = 0;
D0   = diag(F); % spdiags(ones(m,1),0,m,m); 
D0(1:P.Tr,1:P.Tr) = 0;
D_1  = diag(F(1:end-1),-1); % spdiags(ones(m,1),-1,m,m); 
D    = (D1 + D0 - D_1);
D(1,1)    = F(1); % Reflecting barrier
D(P.Tr,P.Tr+1) = 0; % Absorbing barrier

% FS1(1:T,1:T)  =  FS1(1:T,1:T)/P.Vres;
Flow          =  D;

% D2              = diag(ones(m-1,1),-1) - diag(2*ones(m,1)) + diag(ones(m-1,1),1); % spdiags([ones(m,1) -2*ones(m,1) ones(m,1)],[-1 0 1],m,m);%LIFDiff21M(m,1);
D2              = zeros(m);
D2(1:m+1:m*m)     = -2;
D2(2:m+1:m*m)     = 1;
D2(m+1:m+1:m*m)   = 1;
D2(T,T-2:T+2)   = 2*[0 1/2 -1 0 0];   %2*[1/2 -1 1/2 0 0];
D2(T-1,T-3:T+1) = 2*[0 1/2 -1 1/2 0]; %2*[0 1/2 -1 1/2 0];
D2(1,1)         = -1;
% D2(T-2,T-4:T)   = 2*[-1/24 2/3 -5/4 2/3 0];
D2(T+1:end,:)   = 0;
D2(T+1,T)       = 1;

% noise for all channel families

% DS1 = zeros(m,m);
% 
% for k = 1:length(Go)
%     DS1 = DS1 + sparse(1:m,1:m, P.DGL(k)*.5*(P.DL(k)^2)*abs(FvarV(:,k)) ,m,m);
% end

Diffusion = (1/2/(P.Vres^2)*P.D)*D2; % 1/2*D2*(P.D*speye(m) + DS1)/(P.Vres^2);

Qdpdt         = -Flow + Diffusion;
Qdpdt(R,end)    = -Qdpdt(end,end);
% Qdpdt(R,end)    = Qdpdt(R,end)   + F(end);

% SrateOperator and Gi flow
Qdpdt(m+1:m+n,P.delay) = P.Vres*Qdpdt(P.delay+1,P.delay).*(1./P.T); % updated Version
Qdpdt(m+1:m+n,m+1:m+n) = -diag((1./P.T)); 

if isfield(P,'optControl')
    
    for c = 1:length(P.optControl)
        chan      = P.optControl(c).chan;
        condgains = P.optControl(c).condgains;
        LdL       = P.optControl(c).LdLgains;
        
        Qdpdt(m+chan,P.delay) = condgains*(1/P.T(chan))*Qdpdt(m+1:m+n,P.delay)*LdL(2); % updated Version
        Qdpdt(m+chan,m+1:m+n) = (1/P.T(chan))*(-condgains*diag((1./P.T))*LdL(2) + condgains*LdL(1));
        Qdpdt(m+chan,m+chan)  = -1/P.T(chan);
    end
end

end

