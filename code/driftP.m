function [Pmove,Ncell] = driftP(S,dt,dV)

%% assumes zero knowledge of local drift direction and speed.

% average debris drift speed
ds = dV*dt; %m
% standard deviation (GUESS!)
sigma_s = .5*ds; %m
% grid size
dx1 = mean(diff(unique(S.xnodes))); %m
dx2 = mean(diff(unique(S.ynodes))); %m
% nextcell distance (approximate by diam)
dd = sqrt(dx1^2+dx2^2)/2;

% range of motion in cells (1sigma)
Nmove = ceil((ds+sigma_s)/dd);
Pmove = zeros(Nmove,1);
Pmove(1) = normcdf(dd,ds,sigma_s);
for i=2:Nmove-1
    Pmove(i) = normcdf(dd*i,ds,sigma_s) - Pmove(i-1);
end
Pmove(end) = 1 - sum(Pmove(1:end-1));
Ncell = 8*(1:Nmove)-8;
Ncell(1) = 1;
Pmove = Pmove(:)./Ncell(:);

end