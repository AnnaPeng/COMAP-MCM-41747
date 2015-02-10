function [P,Pescape] = driftTransitionOLD(S,dt,dV,P)

% assumes perfect knowledge of local drift direction but perfect knowledge
% of speed. (WEIRD!)

Ppadded = zeros(S.dim+1);
Ppadded(2:(end-1),2:(end-1)) = P;

% average debris drift speed
ds = dV*dt; %m
% standard deviation
sigma_s = .5*ds; %m
% grid size
dx = mean(diff(unique(S.xnodes))); %m
% neighbor
dn = dx/2;
% diagonal
dd = sqrt(2)*dn;

sSim = randn(S.numelements)*sigma_s + ds;
Tdn = sSim > dn;
Tdd = sSim > dd;

Tmove = Tdn*4+Tdd*4;
Tstay = Tmove == 0;
Tdn = Tdn ./ Tmove;
Tdd = Tdd ./ Tmove;

Tdn(isnan(Tdn)) = 0;
Tdd(isnan(Tdd)) = 0;

for i = 1:size(sSim,1)
    for j = 1:size(sSim,2)
        Pcell = P(i,j);
        Ppadded(i+[0 2],j) = Ppadded(i+[0 2],j) + Tdn(i,j)*Pcell;
        Ppadded(i,j+[0 2]) = Ppadded(i,j+[0 2]) + Tdn(i,j)*Pcell;
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