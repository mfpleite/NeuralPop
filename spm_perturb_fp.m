function M = spm_perturb_fp(P,M,r)


M0=M; P0=P;

S1 = M.J.V(:,M.J.IX(end));
k=2;
delta = .1;



nstates = 0;
for l = 1:M.np   % cycle through populations
    P1      = M.P(l);
    SSS2    = [reshape(M.SS(l).SS(:,k-1),[],1,1); squeeze(M.GV(l,:,k-1))'];
    S(nstates+1 : nstates + M.P(l).LVV + M.nc) = SSS2;
    nstates = nstates + M.P(l).LVV + M.nc;
end

f0=1;f1=-1;
 
while norm((f0 + f1)/norm(f0-f1))<.001
    S0 = S' + delta*real(S1);
    
    nstates = 0;
    for l = 1:M.np   % cycle through populations
        M.SS(l).SS(:,1) = S0(nstates+1:nstates+M.P(l).LVV);
        M.GV(l,:,1)     = S0(nstates+M.P(l).LVV+1:nstates+M.P(l).LVV+M.nc);
        nstates = nstates + M.P(l).LVV + M.nc;
    end
    
    [f0, J] = fx_LIFpopMEJpar(P,M);
    
    S01 = S' - delta*real(S1);
    
    nstates = 0;
    for l = 1:M.np   % cycle through populations
        M.SS(l).SS(:,1) = S01(nstates+1:nstates+M.P(l).LVV);
        M.GV(l,:,1)     = S01(nstates+M.P(l).LVV+1:nstates+M.P(l).LVV+M.nc);
        nstates = nstates + M.P(l).LVV + M.nc;
    end
    
    [f1, J] = fx_LIFpopMEJpar(P,M);
    delta = delta*2;

end
    
if nargin>2
    M.opt.dttol = M.opt.dttol*10^r;
    M.opt.LCconvTol = M.opt.dttol*10^-1;
    M = spm_DCM_lifpopsys_LC_int(P,M);
    M = spm_lifpopsys_LC_prepare(P,M);
    M.opt.dttol = M.opt.dttol*10^-r;
    M.opt.LCconvTol = M.opt.dttol*10;
end
    
end