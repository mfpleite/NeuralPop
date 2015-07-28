function [r, to]= convergence_check(GV,t,k,T0)
% GV = M.GV; t = M.opt.t; T0=M.opt.T0;
r=Inf;
to = 0;

GV(:,:,k+1:end) = [];
t(k+1:end) = [];

if sum(~(t<t(end)-T0))>10
    GV(:,:,t<t(end)-T0)=[];
    t(t<t(end)-T0)=[];
end

s = size(GV);
VecGV = reshape(GV,[],s(end));

D = sum(bsxfun(@minus,VecGV(:,end),VecGV(:,1:end-4)).^2);
% length(D)
[pks,locs] = findpeaks(-D);

[C, I] =  max(pks); % find nearest recurrences

delta = max(abs(D(locs(I)) - D(locs(I)-1:locs(I)+1))); % remove contigous points

locs(pks<-delta+C) = [];

try
    I = locs(end);
    
    curvexy = VecGV(:,I-3:I+3)';
    mapxy   = VecGV(:,end)';
    [xy,d,t_a]=distance2curve(curvexy,mapxy,'spline');

    
    r = norm(d);
    if nargout>1
        to = interp1(0:1/6:1,t(I-3:I+3),t_a,'spline');
end




end



