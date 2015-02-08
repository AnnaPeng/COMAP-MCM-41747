clc; clear; close all;
% constants
gE = 9.81; %m/s2
%% plane specs (B737-900ER; G280; A380)
% % Specific Fuel Consumption (kg/N/s)
% B737.SFC = 20E-6;
% G280.SFC = 20E-6;
% A380.SFC = 20E-6;
% % Lift to Drag Ratio
% B737.L2D = 20;
% G280.L2D = 20;
% A380.L2D = 20;
% % Dry Weight = Operation Empty Weight + Total Payload (kg)
% B737.Wdry = 40E3+25E3;
% G280.Wdry = 40E3+25E3;
% A380.Wdry = 40E3+25E3;
% % Total fuel load = Maximum Take Off Weight - Dry Weight (kg)
% B737.Wfuel = 85E3 - B737.Wdry;
% G280.Wfuel = 85E3 - G280.Wdry;
% A380.Wfuel = 85E3 - A380.Wdry;
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
casenum = 1;
if casenum == 1
    AC = B737; acname = 'B737-900ER';
elseif casenum == 2
    AC = A380; acname = 'Airbus 380';
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
% N1grid = 100;
% N2grid = 60;
N1grid = 25;
N2grid = 20;
% # of 1D quadrature points (use even number for symm)
Nquad = 2;
% boundary distance [-x1 +x1 -x2 +x2]
bdry = [-150 350 -150 150]*1e3; %m
% reverse probability
s = .05; q = 1/pi-2*s;

%% incident to crash range pdf

AC.Rc = AC.Vc/gE/AC.SFC*AC.L2D*log(1+Fr*AC.Wfuel/AC.Wdry);
AC.Tc = AC.Rc/AC.Vc;
% fire, collision, glide, other
Ncrash = [1406, 1901, 901, 498];

% stdev of the statistics (guess)
Rstdev = [.2 .2 .2 .2];

Rhist = [];
for i=1:numel(Ncrash)
    Rhist = [Rhist; randn(Ncrash(i),1)*Rstdev(i)*AC.R(i) + AC.R(i)];
end

% smoothed profile + visuialization
[PR,R]=ksdensity(Rhist,[0:1e3:50e3 53e3:3e3:300e3]);
figure(); hold all; grid on;
plot(R/1e3,PR,'rx--');
xlim([0 300]); xlabel('Crash Radius [km]'); ylabel('Probability')
title(['Estimated Crash Distance for ' acname])
saveas(gcf,[acname '_CrashDistance.png']);

%% continuous probability at x=(x1,x2)
rtil = linspace(0,rint,Nrtil);
Rtil = @(rtil,x) norm(x-[rtil;0]);
theta = @(rtil,x) abs(atan2(x(2),x(1)-rtil));
ptheta = @(rtil,x) q+theta(rtil,x)*(s-q)/pi;
pdfx = @(rtil,x) interp1q(R,PR,Rtil(rtil,x))*ptheta(rtil,x);
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

%% heavy loops
% note the vertical symmetry!!!
wv = wv*2;
for ev=1:S.numelements(2)/2
    for eu=1:S.numelements(1)
        for l=1:Nquad/2
            for k=1:Nquad
                % pullback quadrature pts coordinate 
                [x,y] = S.coords(eu,ev,Nu(:,k),Nv(:,l));
                Ptemp = 0;
                for i = 1:Nrtil
                    Ptemp = Ptemp + pdfx(rtil(i),[x;y])*rint/Nrtil;
                end
                P(eu,ev) = P(eu,ev) + wu(k)*wv(l)*Ptemp;
            end
        end
    end
    P(:,S.numelements(2)-ev+1) = P(:,ev);
end
%% renormalize and save data
P = P / sum(P(:));
save(num2str(casenum,'prior%d.mat'),'P','S');
% load(num2str(casenum,'prior%d.mat'))

%% crash probability distribution graph

Splot = Surface(N1grid,N2grid,bdry/1e3);
figure(); hold all; grid on;
plottwoform(Splot,P,3); colorbar;
xlabel('Tangent Direction [km]'); ylabel('Lateral Direction [km]');
title('Probability of Aircraft Landed in Cell');
hold all;
traj = plot3([-1e6 0 rint 1e6]/1e3, [0 0 0 0], [1 1 1 1],'rx--');
set(traj,'linewidth',2,'markersize',15)
saveas(gcf,[acname 'PriorDistribution.png']);

%% drift diffusion transition

% decimate P to 25x20 for this process
N1drift = 25; N2drift = 20;
PP = zeros(N1drift,N2drift);
Sdrift = Surface(N1drift,N2drift,bdry);
SdriftP = Surface(N1drift,N2drift,bdry/1e3);
% scaling factors
sc1 = N1grid/N1drift;
sc2 = N2grid/N2drift;
for i = 1:25
    for j = 1:20
        PP(i,j) = sum(sum(P(sc1*i-sc1+(1:sc1),sc2*j-sc2+(1:sc2))));
    end
end
% update interval
dt = .5*60*60; %s
% average debris drift speed
dV = 10; %m/s
% escape/out of range probability
Pescape = 0;

figure(); hold all; grid on;
% traj = plot3([-1e6 0 rint 1e6]/1e3, [0 0 0 0], [1 1 1 1],'rx--');
plottwoform(SdriftP,PP,3); colorbar;
xlabel('Tangent Direction [km]'); ylabel('Lateral Direction [km]');
title('Probability of Aircraft Landed in Cell');
hold all;
set(traj,'linewidth',2,'markersize',15)
for tstep = 1:100
    pause(.1);
    [PP,temp] = driftTransition(Sdrift,dt,dV,PP);
    Pescape = Pescape*temp;
    plottwoform(SdriftP,PP,3);
end

%% add in (wind) random field disturbance

% distribution


%% lateral range curve/function

% magnetometer (gaussian)

% visual (tapered linear)

% 


%% Search planning

% search effort

% search density

% detection probability (transitions)

% optimal search plan


%% Monte Carlo, evaluation

% "particle filter" sample prior, propagate both search party and lost a/c

% calculate time till capture
