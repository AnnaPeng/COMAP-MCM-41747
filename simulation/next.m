function [P,q] = next(S,P,Pmove)

Nmove = numel(Pmove);
Ppadded = zeros(S.dim+Nmove-1);
% Ppadded(Nmove:(end-Nmove+1),Nmove:(end-Nmove+1)) = P;

for i = 1:size(P,1)
    for j = 1:size(P,2)
        Pcell = P(i,j);
%         for prop = 1:Nmove
            Ppadded(i+[0 2],j+(0:2)) = Ppadded(i+[0 2],j+(0:2)) + Pmove(2)*Pcell;
            Ppadded(i+1,j+[0 2]) = Ppadded(i+1,j+[0 2]) + Pmove(2)*Pcell;
            Ppadded(i+1,j+1) = Ppadded(i+1,j+1) + Pmove(1)*Pcell;
%         end
    end
end

P = Ppadded(Nmove:(end-Nmove+1),Nmove:(end-Nmove+1));
% renormalization
P = P./sum(P(:));
q = sum(sum(Ppadded([1:Nmove-1 (end-Nmove+1):end],:))) + ...
    sum(sum(Ppadded(:,[1:Nmove-1 (end-Nmove+1):end])));

end