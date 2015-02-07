function q = glq1d(f,Np,varargin)
% N-point Guassian-Legendre Quadrature on unit interval [0,1]
% q = glq1d(f,Np,varargin)
%
% INPUT
% f        - function handle of quadrature target
% Np       - number of quadrature points
% varargin - leftover arguments for function handle f
% 
% OUTPUT
% q        - quadrature result
% 
if nargin < 2 || isempty(Np),Np = 3;end
if Np > 5, error('order 6+ not available yet'); end
switch Np
    case 1
        x = .5; w = 1;
    case 2
        t = sqrt(3)/3;
        x = ([-t t]+1)/2;
        w = [1 1]/2;
    case 3
        t = sqrt(3/5);
        x = ([-t 0 t]+1)/2;
        w = [5/9 8/9 5/9]/2;        
    case 4
        t1 = sqrt((3-2*sqrt(6/5))/7);
        t2 = sqrt((3+2*sqrt(6/5))/7);
        x = ([-t2 -t1 t1 t2]+1)/2;
        t = sqrt(30)/36;
        w = [.5-t .5+t .5+t .5-t]/2;        
    case 5
        t1 = sqrt(5-2*sqrt(10/7))/3;
        t2 = sqrt(5+2*sqrt(10/7))/3;
        x = ([-t2 -t1 0 t1 t2]+1)/2;
        t1 = (322+13*sqrt(70))/900;
        t2 = (322-13*sqrt(70))/900;
        w = [-t2 -t1 128/225 t1 t2]/2;
end
q = 0;
for ii=1:numel(x)
    q = q + w(ii)*feval(f,x(ii),varargin{:});
end
end
