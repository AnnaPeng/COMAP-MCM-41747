clc; clear; close all;
% constants
gE = 9.81; %m/s2
%% plane specs (B737-900ER; G280; A380)
% Cruise speed (m/s)
B737.Vc = 243;
G280.Vc = 250;
A380.Vc = 262;
% crash distance
% fire, collision, glide, other
A380.R = [24 50 160 80]*1e3; %m
B737.R = [15 30 140 63]*1e3; %m
G280.R = [10 20 135 56]*1e3; %m

%% Assumptions/parameters

% target aircraft make
ACcase = 1;
if ACcase == 1
    AC = B737; acname = 'B737-900ER';
elseif ACcase == 2
    AC = A380; acname = 'Airbus380';
else
    AC = G280; acname = 'G280';
end
% remaining fuel ratio at incident
Fr = .5;
% communication interval/distance
interval = 15*60; %s
rint = AC.Vc*interval; %m
% continuous probability (rtil) riemann sum resolution
Nrtil = 50;
% grid resolution
GRIDcase = 2;
if GRIDcase == 1
    N1grid = 100;
    N2grid = 60;
else
    N1grid = 25;
    N2grid = 20;
end
% # of 1D quadrature points (use even number for symm)
Nquad = 2;
% search domain bounds [-x1 +x1 -x2 +x2]
bdry = [-150 350 -150 150]*1e3; %m
% reversing probability param
s = .05; q = 1/pi-2*s;

% average debris drift speed
dV = 10; %m/s

%% incident to crash range pdf

% fire, collision, glide, other
Ncrash = [1406, 1901, 901, 498];

% stdev of the statistics (guess)
Rstdev = [.2 .2 .2 .2];

Rhist = [];
for i=1:numel(Ncrash)
    Rhist = [Rhist; randn(Ncrash(i),1)*Rstdev(i)*AC.R(i) + AC.R(i)];
end

% smoothed profile + visuialization
[PR,R]=ksdensity(Rhist,[0:1e3:50e3 55e3:5e3:300e3 310e3:10e3:500e3]);
save([acname '_CrashRadius.mat'],'R','PR');

%% continuous probability at x=(x1,x2)
    rtil = linspace(0,rint,Nrtil);
    Rtil = @(x) sqrt(sum((repmat(x,1,Nrtil)-[rtil;zeros(size(rtil))]).^2));
    theta = @(x) abs(atan2(x(2),x(1)-rtil));
    ptheta = @(x) q + theta(x)*(s-q)/pi;
    pdfx = @(x) interp1(R,PR,Rtil(x)).*ptheta(x);
    %% discritized probability at cell (eu,ev)

    % pre-compute gauss-legendre points and weights
    [u,wu] = gaussquad(Nquad);
    [v,wv] = gaussquad(Nquad);
    Nu = nodefun(u);
    Nv = nodefun(v);

    % total domain area
    Agrid = (bdry(2)-bdry(1))*(bdry(4)-bdry(3)); %m2
    % grid pt construction
    S = Surface(N1grid,N2grid,bdry);
    P = zeros(S.numelements);

    %% loops
    % note the vertical symmetry!!!
    wv = wv*2;
    for ev=1:S.numelements(2)/2
        for eu=1:S.numelements(1)
            for l=1:Nquad/2
                for k=1:Nquad
                    % pullback quadrature pts coordinate 
                    [x,y] = S.coords(eu,ev,Nu(:,k),Nv(:,l));
                    Ptemp = sum(pdfx([x;y]))*rint/Nrtil;
                    P(eu,ev) = P(eu,ev) + wu(k)*wv(l)*Ptemp;
                end
            end
        end
        P(:,S.numelements(2)-ev+1) = P(:,ev);
    end
    %% renormalize and save data
    P = P / sum(P(:));
    save(num2str(ACcase,'prior%d.mat'),'P','S');
    % load(num2str(casenum,'prior%d.mat'))

    %% crash probability distribution graph

    Splot = Surface(N1grid,N2grid,bdry/1e3);
    figure(); hold all; grid on;
    plottwoform(Splot,P,3); colorbar;
    xlabel('Tangent Direction [km]'); ylabel('Lateral Direction [km]');
    title('Aircraft Debris Location Density at t=0 hr');
    hold all;
    traj = plot3([-1e6 0 rint 1e6]/1e3, [0 0 0 0], [1 1 1 1],'rx--');
    set(traj,'linewidth',2,'markersize',15)
    saveas(gcf,[acname '_PriorDistribution.png']);

%% drift/diffusion simulation
    Tsim = 96*3600; %s

    PP = P;
    % update interval that is appropriate
    % i.e. allow only single cell diffusion given grid resolution
    if GRIDcase == 1
        dt = .1*60*60; %s
        Nsim = Tsim/dt; %steps
    else
        dt = .4*60*60; %s
        Nsim = Tsim/dt; %steps
    end
    [Pmove,Ncell] = driftP(S,dt,dV);

    % propogation steps
    
    % Probability of escape at t
    qt = zeros(Nsim+1,1);
    
    for t = 1:Nsim
        [PP,qt(t+1)] = next(S,PP,Pmove);
    end
    %% Location density if no search initiates
    figure(); hold all; grid on;
    plottwoform(Splot,PP,3); colorbar;
    xlabel('Tangent Direction [km]'); ylabel('Lateral Direction [km]');
    title(num2str(Tsim/3600, 'Aircraft Debris Location Density at t=%d hr'));
    saveas(gcf,[acname '_NoSearchDistribution.png']);
    %% Graph of escape probability over time
    tVec = (0:dt:Tsim)/3600; %hr
    figure(); hold all; grid on;
    plot(tVec,qt,'k-');
    xlabel('Time [hr]'); ylabel('Probability');
    title('Probability of Aircraft Debris Escaping the Search Domain');
    saveas(gcf,[acname '_NoSearchEscape.png']);
%% Search Agent Data

% given 99% detection range find sigma
ncdf = @(sig,Rd) normcdf(Rd,0,sig)-normcdf(-Rd,0,sig);
cdf2sig = @(Rd) fminsearch(@(sig) abs(ncdf(sig,Rd)-.99),Rd/2);

% Marine Vessel (Damen SAR vessel 1816/1906 range= 600km+)
MV.Vs = 15.5; %m/s
MV.Rdetect = 1e3; %m
MV.FA = 0;
MV.sig = cdf2sig(MV.Rdetect); %m

% UAV (Hermes)
UAV.Vs = 49; %m/s
UAV.Rdetect = 5e3; %m
UAV.alt = 5e3; %m
UAV.FA = .05;
UAV.sig = cdf2sig(UAV.Rdetect); %m

% Helicoptor (USCG Dolphin MH65C range= 650km+)
heli.Vs = 90; %m/s
heli.Rdetect = 5e3; %m
heli.alt = 5.5e3; %m
heli.FA = .05;
heli.sig = cdf2sig(heli.Rdetect); %m

%% search agent initial states
% number of agents on one side
Nagent = 100;

% Initial position

% Ps0 = [S.xnodes(S.getglobalboundarynodes) ...
%     S.ynodes(S.getglobalboundarynodes)];

% Detection Probability are assumed to be normal
mvncdf(x1r,x2r,xs,sig);
%% Distributed Search Plan (edge first and chase the highest cell)
    Tsim = 96*3600; %s

    PP = P;
    % update interval that is appropriate
    % i.e. allow only single cell diffusion given grid resolution
    if GRIDcase == 1
        dt = .1*60*60; %s
        Nsim = Tsim/dt; %steps
    else
        dt = .4*60*60; %s
        Nsim = Tsim/dt; %steps
    end
    [Pmove,Ncell] = driftP(S,dt,dV);

    % propogation steps
    
    % Probability of escape at t
    qt = zeros(Nsim+1,1);
    
    for t = 1:Nsim
        [PP,qt(t+1)] = next(S,PP,Pmove);
    end
    %% Location density
    figure(); hold all; grid on;
    plottwoform(Splot,PP,3); colorbar;
    xlabel('Tangent Direction [km]'); ylabel('Lateral Direction [km]');
    title(num2str(Tsim/3600, 'Aircraft Debris Location Density at t=%d hr'));
    saveas(gcf,[acname '_NoSearchDistribution.png']);
    %% Graph of escape probability over time
    tVec = (0:dt:Tsim)/3600; %hr
    figure(); hold all; grid on;
    plot(tVec,qt,'k-');
    xlabel('Time [hr]'); ylabel('Probability');
    title('Probability of Aircraft Debris Escaping the Search Domain');
    saveas(gcf,[acname '_NoSearchEscape.png']);
%% Distributed Search Plan (center first)


%% concentrated "single" agent Search Plan (steepest descent)


%% Monte Carlo evaluation of time of capture


% "particle filter" ?
