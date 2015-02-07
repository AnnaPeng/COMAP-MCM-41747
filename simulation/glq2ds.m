function q = glq2ds(f,Np,varargin)
% N-point Guassian-Legendre Quadrature on unit square [0,1]^2
% q = glq2ds(f,Np,varargin)
%
% INPUT
% f        - function handle of quadrature target
% Np       - number of quadrature points
% varargin - leftover arguments for function handle f
% 
% OUTPUT
% q        - quadrature result
%
%% 
% if nargin < 2 || isempty(Np),Np = 4;end
% if Np > 4, error('order 5+ not available yet'); end
% switch Np
%     case 1
%         x = 1/3; y = 1/3; w = 1/2;
%     case 2
%         error('order 2 not available');
%     case 3
%         x = [.5 .5 0];
%         y = [0 .5 .5];
%         w = [1 1 1]/6;
%     case 4
%         x = [1/3 2/15 2/15 11/15];
%         y = [1/3 11/15 2/15 2/15];
%         w = [-27 25 25 25]/96;
% end
% q = 0;
% for ii=1:numel(x)
%     q = q + w(ii)*feval(f,x(ii),y(ii),varargin{:});
% end
% end
%%
if nargin < 2 || isempty(Np),Np = 3;end
if Np > 25, error('point 26+ not available yet'); end
switch Np
    case 1
        x = .5; w = 1;
    case 4
        t = sqrt(3)/3;
        x = ([-t t]+1)/2;
        w = [1 1]/2;
    case 9
        t = sqrt(3/5);
        x = ([-t 0 t]+1)/2;
        w = [5/9 8/9 5/9]/2;        
    case 16
        t1 = sqrt((3-2*sqrt(6/5))/7);
        t2 = sqrt((3+2*sqrt(6/5))/7);
        x = ([-t2 -t1 t1 t2]+1)/2;
        t = sqrt(30)/36;
        w = [.5-t .5+t .5+t .5-t]/2;        
    case 25
        t1 = sqrt(5-2*sqrt(10/7))/3;
        t2 = sqrt(5+2*sqrt(10/7))/3;
        x = ([-t2 -t1 0 t1 t2]+1)/2;
        t1 = (322+13*sqrt(70))/900;
        t2 = (322-13*sqrt(70))/900;
        w = [-t2 -t1 128/225 t1 t2]/2;
    otherwise
        error('not a square number');
end
q = 0;
c = @(w_x,w_y) w_x*w_y;
for ii=1:numel(x)
    for jj=1:numel(x)
        q = q + c(w(ii),w(jj))*...
            feval(f,x(ii),x(jj),varargin{:});
    end
end
end
