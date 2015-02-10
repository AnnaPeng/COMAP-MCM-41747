function M = edgefun(u)
% compute univariate constant function at u in (-1,1)
    for i=length(u):-1:1
        M(:,i) = 0.5;
    end
end

