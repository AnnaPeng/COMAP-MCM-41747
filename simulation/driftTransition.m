function [P,Pescape] = driftTransition(S,dt,dV,P)

% assumes knowledge of local drift direction.

Ppadded = zeros(S.dim+1);
Ppadded(2:(end-1),2:(end-1)) = P;

% Td = zeros(numel(S.dim+1));

% average debris drift speed
ds = dV*dt; %m
% standard deviation
sigma_s = .5*ds; %m
% grid size
dx1 = mean(diff(unique(S.xnodes))); %m
dx2 = mean(diff(unique(S.ynodes))); %m
% vertical 
dv = dx2/2;
% horizontal
dh = dx1/2;
% diagonal
dd = sqrt(dv^2+dh^2);

sSim = randn(S.numelements)*sigma_s + ds;
Tdv = sSim > dv;
Tdh = sSim > dh;
Tdd = sSim > dd;

Tmove = Tdv*2+Tdh*2+Tdd*4;
Tstay = Tmove == 0;
Tdv = Tdv ./ Tmove;
Tdh = Tdh ./ Tmove;
Tdd = Tdd ./ Tmove;

Tdh(isnan(Tdh)) = 0;
Tdv(isnan(Tdv)) = 0;
Tdd(isnan(Tdd)) = 0;

for i = 1:size(sSim,1)
    for j = 1:size(sSim,2)
        Pcell = P(i,j);
        Ppadded(i+[0 2],j) = Ppadded(i+[0 2],j) + Tdh(i,j)*Pcell;
        Ppadded(i,j+[0 2]) = Ppadded(i,j+[0 2]) + Tdv(i,j)*Pcell;
        Ppadded(i+2,j+[0 2]) = Ppadded(i+2,j+[0 2]) + Tdd(i,j)*Pcell;
        Ppadded(i,j+[0 2]) = Ppadded(i,j+[0 2]) + Tdd(i,j)*Pcell;
        Ppadded(i+1,j+1) = Ppadded(i+1,j+1) + Tstay(i,j)*Pcell;
    end
end

P = Ppadded(2:(end-1),2:(end-1));
% renormalization
P = P./sum(P(:));
Pescape = sum(sum(Ppadded([1 end],:))) + sum(sum(Ppadded(:,[1 end])));

end