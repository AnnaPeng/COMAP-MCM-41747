function [P,Pescape] = driftTransition(S,dt,dV,P)

% assumes perfect knowledge of local drift direction or speed. (BAD!)

Ppadded = zeros(S.dim+1);
Ppadded(2:(end-1),2:(end-1)) = P;

% average debris drift speed
ds = dV*dt; %m
% standard deviation
sigma_s = .5*ds; %m
% grid size
dx = mean(diff(unique(S.xnodes))); %m
% nextcell distance (approximate by diam)
dd = sqrt(2)*(dx/2);

sSim = randn(S.numelements)*sigma_s + ds;
dirSim = randi(8,S.numelements);
Tmove = (sSim > dd).*dirSim;

Pnext = zeros(size(Ppadded));
for ii = 1:size(sSim,1)
    for jj = 1:size(sSim,2)
        i = ii+1; j=jj+1;
        switch Tmove(ii,jj)
            case 0
                Pnext(i,j) = Pnext(i,j) + Ppadded(i,j);
            case 1
                Pnext(i+1,j) = Pnext(i+1,j) + Ppadded(i,j);
            case 2
                Pnext(i+1,j+1) = Pnext(i+1,j+1) + Ppadded(i,j);
            case 3
                Pnext(i,j+1) = Pnext(i,j+1) + Ppadded(i,j);
            case 4
                Pnext(i-1,j+1) = Pnext(i-1,j+1) + Ppadded(i,j);
            case 5
                Pnext(i-1,j) = Pnext(i-1,j) + Ppadded(i,j);
            case 6
                Pnext(i-1,j-1) = Pnext(i-1,j-1) + Ppadded(i,j);
            case 7
                Pnext(i,j-1) = Pnext(i,j-1) + Ppadded(i,j);
            case 8
                Pnext(i+1,j-1) = Pnext(i+1,j-1) + Ppadded(i,j);
        end
    end
end

P = Pnext(2:(end-1),2:(end-1));
% renormalization
P = P./sum(P(:));
Pescape = sum(sum(Pnext([1 end],:))) + sum(sum(Pnext(:,[1 end])));

end