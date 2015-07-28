function [ P ] = spm_ExpP(P0)

P = P0;
P.ang = P0.ang;

for k = 1:length(P0.A)
    P.A(k).M = exp(P0.A(k).M);
end

for k = 1:length(P0.C)
    P.C(k).M = exp(P0.C(k).M);
end

for k = 1:length(P0.P)
    P.P(k).C   = exp(P0.P(k).C);   % [nF] Cell capacitance
    P.P(k).Vl  =     P0.P(k).Vl;   % [mV] Leak reversal portential
    P.P(k).gl  = exp(P0.P(k).gl);  % [pF/ms] = [nS] Leak conductance
    P.P(k).Vt  =     P0.P(k).Vt;   % [mV] Leak reversal portential
    P.P(k).Vr  =     P0.P(k).Vr;   % [mV] Leak reversal portential
    P.P(k).Rt  = exp(P0.P(k).Rt);  % [s]  Refractory time
    P.P(k).Vg  =     P0.P(k).Vg;   % [mV] Reversal portential of the families of synaptic channels
    P.P(k).T   = exp(P0.P(k).T);   % [s] Decay time constants for the synaptic conductances
    P.P(k).D   = exp(P0.P(k).D);   % (1/2*sigma^2) [(mV^2)/ms] Membrane noise term
    
    if isfield(P0.P(k),'d'); 
        P.P(k).d   = exp(P0.P(k).d); % [s]  Afferent axonal delay
    end
    
end
