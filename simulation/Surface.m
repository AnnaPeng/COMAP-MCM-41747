classdef Surface
    %SURFACE Summary of this class goes here
    %   Detailed explanation goes here
    
    % x-,y- coordinates of nodes
    properties
        xnodes
        ynodes
    end
    
    % number of mesh-primitives
    properties
        dim
        numelements
        numnodes
        numedges
        numfaces
    end
    
    methods 
        
        % constructor
        function S = Surface(varargin)
           
            switch nargin
                case 0
                    n = 1; m = 1;
                    [S.ynodes,S.xnodes] = meshgrid(linspace(0,1,n+1),linspace(0,1,m+1)); 
                    
                case 2
                    n = varargin{1}; 
                    m = varargin{2};
                    [X,Y] = meshgrid(linspace(0,1,n+1),linspace(0,1,m+1));
                    S.xnodes = X'; S.ynodes = Y';
                    
                case 3
                    n = varargin{1}; 
                    m = varargin{2};
                    bdry = varargin{3};
                    [X,Y] = meshgrid(linspace(bdry(1),bdry(2),n+1),...
                                     linspace(bdry(3),bdry(4),m+1));
                    S.xnodes = X'; S.ynodes = Y';
                case 4
                    n = varargin{1}; 
                    m = varargin{2};
                    S.xnodes = varargin{3};
                    S.ynodes = varargin{4};
            
                otherwise
                    error('To many input arguments');
            end
            
            % set dependent values
            S.dim = [n+1,m+1];
            S.numelements = [n,m];
            S.numnodes    = prod([n+1,m+1]);
            S.numedges(2) = prod([n,m+1]);
            S.numedges(1) = prod([n+1,m]);
            S.numfaces    = prod([n,m]);
        end
        
        % compute jacobian 2 by 2 matrix
        function J = jacobian(S,eu,ev,Nu,Nv,dNu,dNv)
            % input: Surface S, element (eu,ev) and
            % basis functions Nu, Nv and derivatives dNu,dNv
            
            % active nodes in element (eu,ev)
            X = S.xnodes(eu:eu+1,ev:ev+1);
            Y = S.ynodes(eu:eu+1,ev:ev+1);
            
            % compute Jacobian
            J = [dNu' * X * Nv, Nu' * X * dNv;
                 dNu' * Y * Nv, Nu' * Y * dNv];
        end
        
        % compute coordinates of mapping
        function [x,y] = coords(S,eu,ev,Nu,Nv)
            % input: Surface S, element (eu,ev) and
            % basis functions Nu, Nv
            
            % active nodes in element (eu,ev)
            X = S.xnodes(eu:eu+1,ev:ev+1);
            Y = S.ynodes(eu:eu+1,ev:ev+1);
            
            % compute Jacobian
            x = Nu' * X * Nv;
            y = Nu' * Y * Nv;
        end
        
        % get global nodenumbers on the boundary
        function gnn = getglobalboundarynodes(S)
            n = S.dim(1); m = S.dim(2);
            
            temp = reshape([1:S.numnodes],n,m);
            gnn = setdiff(temp,temp(2:end-1,2:end-1));
        end
        
        % get global edgenumbers on the boundary
        function gnn = getglobalboundaryedges(S)
            n = S.dim(1); m = S.dim(2);
            
            temp = reshape([1:S.numedges(1)],n,m-1);
            gnn{1} = setdiff(temp,temp(2:end-1,:));
            
            temp = reshape([1:S.numedges(2)],n-1,m) + S.numedges(1);
            gnn{2} = setdiff(temp,temp(:,2:end-1));
            gnn = [gnn{1}(:);gnn{2}(:)];
        end
   
        % get global nodenumbers for element (eu,ev)
        function gnn = getglobalnodenumber(S,eu,ev)
            % global numbering in subs
            [nnu,nnv] = meshgrid(eu:eu+1,ev:ev+1);

            % global numbering in linear indexing
            gnn = sub2ind(S.dim, nnu', nnv');
        end
        
        % get global edgenumber for element (eu,ev) in z-dir 
        function gnn = getglobaledgenumber(S,eu,ev,z)
            
            % direction 1
            if z==2
                [nnu,nnv] = meshgrid(eu,[ev:ev+1]');       % global numbering in subs-format
                gnn = sub2ind([S.dim(1)-1,S.dim(2)], nnu', nnv') + S.numedges(1); % global numbering in linear indexing format
            
            % direction 2
            elseif z==1
                [nnu,nnv] = meshgrid(eu:eu+1,ev);       % global numbering in subs-format
                gnn = sub2ind([S.dim(1),S.dim(2)-1], nnu', nnv');     % global numbering in linear indexing format
            else
                error('Wrong input arguments');
            end
        end

        % get global facenumbers for element (eu,ev)
        function gnn = getglobalfacenumber(S,eu,ev)
            % global numbering in linear indexing
            gnn = sub2ind(S.dim-1, eu, ev);
        end
        
        % plot surface
        function plot(S)  
            
            % Create figure
            figure1 = figure;

            % Create axes
            axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[434 342.3 684.6],...
            'DataAspectRatio',[1 1 1]);
            grid(axes1,'on');
            hold(axes1,'on');

            % Create surf
            surf(S.xnodes,S.ynodes,zeros(S.dim),'Parent',axes1,'FaceLighting','none',...
                'EdgeLighting','flat',...
                'MarkerSize',20,...
                'Marker','.',...
                'LineWidth',2,...
                'FaceColor','none');
            hold off;
        end
        
        
    end
    
    
    
end

