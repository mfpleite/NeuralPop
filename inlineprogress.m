function [] = inlineprogress(r,t)

d = floor([0:.05:1]*t);

if r == d(21)
    fprintf(' 100%% \n')
elseif r == d(16)
    fprintf(' 75%% ')
elseif r == d(11)
    fprintf(' 50%% ')
elseif r == d(6)
    fprintf(' 25%% ')
elseif sum(r==d)
    fprintf('.')
end

end