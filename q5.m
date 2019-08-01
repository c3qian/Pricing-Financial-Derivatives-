
s = 0.3; K =140;S = 150;r = 0.05;T = 1;n = 5;h = T/n;u = exp(s*sqrt(h));d = exp(-s*sqrt(h));
% RNP
p = (exp(r*h)-d)/(u-d);
q = 1 - p;
Max = -1;

% Setup
V = zeros(1,n);
V(1) = S;
height = 1;
V = Tree(S, n, height, u, d, V);

% All trajectories stores in All matrix,each column is one path.
% Have two strategies find Max: 
% 1. Keep tracking "Max" value when construct the tree.
% 2.Take the maximum of the column,which is the max price of the path.
% I use method to double check the payoff for each path is correct.
All = zeros(n,2^(n-1));
k = 0;
maximums = zeros(1,2^(n)-1 - 2^(n-1));
for i = 2^(n-1):2^n-1
    ind = i;
    k = k + 1;
    path = zeros(n,1);
    cont = 1;
    while  ind>=1
        path(cont,1) = V(1, ind);
        Max = max(Max,V(1, ind));
        ind = floor(ind / 2);
        cont = cont +1;
    end
    All(:, k) = path;
    maximums(1,k) = max(All(:, k));
    
end

payoff = maximums - All(1,:);
 
% Premium matrix "Ma"
Ma = payoff;
for l= 1:n-1  
          for Q=2^(n-1):-2^(n-1)/2:1 
                for O=1:round(Q/2)
           
            Ma(l +1,O)= exp(-r*h)*(q * Ma(l,2*O-1) + p * Ma(l,2*O));
                 
                end
         end
end
Premium = Ma(n,1);

% Construct binomial tree from Andrew
function V = Tree(S_0,n, parent, u, d, V)
    if (n==1)
    else
       up = S_0*u;
       down = S_0*d;
       V(2*parent) = down;
       V(2*parent+1) = up;
       V = Tree(down,n-1, 2*parent, u,d, V);
       V = Tree(up,n-1, 2*parent + 1, u, d, V);
    end

end
