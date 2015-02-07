function [N,dN] = nodefun(u,diffflag)
if nargin<2 || isempty(diffflag)
% compute univariate hat functions at all points u in (-1,1)
    for i=length(u):-1:1
        N(:,i)  = .5*[1 - u(i); 1 + u(i)];
        dN(:,i) = [-0.5; 0.5];
    end
else
    [~,dN] = nodefun(u);
    N = dN;
end

