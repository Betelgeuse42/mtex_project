function [ot1, ci, work] = recon(ot1, ci, X, Y, h, g, SGN)
% Reconstruct 3+

% dx = zeros(size(X));
% dy = zeros(size(X));
k = length(X);              % number of points

bad_h_counter = 0;
work = 0;
ind_work = [];

for i = 1:k
    if ci(i) <= 0.1
        dx(1:k) = (X(1:k)-X(i));
        dy(1:k) = (Y(1:k)-Y(i));
        Radius = (dy.^2+dx.^2).^0.5; 
        % Just in case
        Radius(i) = 30; % больше чем в области
        
        % Search for the minimum
        z = min(Radius);
        
        % Tolerance 20%
        if h == 12
            eps = 1.8;
        elseif h == 6
            eps = 1.3;
        else
            error('Bad h!');
        end;
        
        % Indices of neighours
        ind = find(Radius < z*eps);
        
        % Is boundary?
        if numel(ind) ~= h
            bad_h_counter = bad_h_counter+1;
            continue;
        end;
        
        renumbering(ind);
        
        % ot1(i) = orientation('Euler',0*degree,0*degree,0*degree,symmetry('m-3m'));
        
        % Number of good neigbours
        b = 0;
        n = [];
        for j = 1:h
            if ci(ind(j)) > 0.1
               b = b+1;
               n = [n ind(j)]; %#ok<AGROW> % last good orientation
            end
        end;
        
        % For good point ??? IF H = 12
        if b >= g
            work = work + 1;
            ind_work = [ind_work i]; %#ok<AGROW>
            ot1(i) = ot1(selectGoodNeighbour(n, SGN));
        else
            continue;
        end;
    end;
    
    if (mod(i, fix(k/15)) == 0)
        fprintf('0');
    end;
end;
fprintf(' done \n');

ci(ind_work) = 1.55;
end

function i = selectGoodNeighbour(n, SGN)
if SGN == 0
    i = n(end);
else
    i = n(randi(length(n),1));
end
end
