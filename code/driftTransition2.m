function [P,Pescape,Pstay] = driftTransition2(S,dt,dV,P)

% assumes zero knowledge of local drift direction and speed.

% average debris drift speed
ds = dV*dt; %m
% standard deviation
sigma_s = .5*ds; %m
% grid size
dx = mean(diff(unique(S.xnodes))); %m
% nextcell distance (approximate by diam)
dd = sqrt(2)*(dx/2);

% range of motion in cells (1sigma)
Nmove = ceil((ds+sigma_s)/dd);
Pmove = zeros(Nmove,1);
Pmove(1) = normcdf(dd,ds,sigma_s);
for i=2:Nmove-1
    Pmove(i) = normcdf(dd*i,ds,sigma_s) - Pmove(i-1);
end
Pmove(end) = 1 - sum(Pmove(1:end-1));
Pdiff = 8*(1:Nmove)-8;
Pdiff(1) = 1;
Pmove = Pmove(:)./Pdiff(:);

Ppadded = zeros(S.dim+Nmove-1);
Ppadded(Nmove:(end-Nmove+1),Nmove:(end-Nmove+1)) = P;

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
Pescape = sum(sum(Ppadded([1:Nmove-1 (end-Nmove+1):end],:))) + ...
    sum(sum(Ppadded(:,[1:Nmove-1 (end-Nmove+1):end])));

Pstay = Ppadded(
end