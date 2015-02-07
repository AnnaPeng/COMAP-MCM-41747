function  plottwoform(S,P,gran)
    
    % dimension of space
    n = S.dim(1); m = S.dim(2);
    
    % partition
    u = linspace(-1,1,gran); v = u;
    
    % get basis functions on (-1,1)
    Nu = nodefun(u); Nv = Nu;
    Mu = edgefun(u); Mv = Mu;
    
    % plot mesh
    hold all; view(2);
    for ev=1:1:S.numelements(2)
        for eu=1:1:S.numelements(1)
            [x,y] = S.coords(eu,ev,Nu,Nv);
            surf(x,y,Mu'* P(eu,ev) * Mv,'EdgeColor','none');
        end
    end
    % plot mesh on top
    xlim([S.xnodes(1,1) S.xnodes(end,1)]);
    ylim([S.ynodes(1,1) S.ynodes(1,end)]);
    hold off;
end

