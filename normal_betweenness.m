% Test Graph
% Adjacency matrix
A = [0 1 0 1 0 0 0 0;
     1 0 0 1 1 1 0 0;
     0 0 0 1 0 0 0 0;
     1 1 1 0 0 0 1 0;
     0 1 0 0 0 1 0 0;
     0 1 0 0 1 0 1 0;
     0 0 0 1 0 1 0 1;
     0 0 0 0 0 0 1 0];
G = graph(A);
 % return the length of shortest path of each pair of nodes
d = distances(G);

% ********************* Betweenness *********************
% The number of shortest path matrix Q
Q = zeros(8,8);
for i = 1:8
    for j = 1:8
        if d(i,j) == 1
            Q(i,j) = 1;
        else
            for k = 1:8
                if d(i,j) > 1 && d(i,j) < 9 && d(k,j) == 1 && d(i,k) == d(i,j) - 1
                    Q(i,j) = Q(i,j) + 1;
                    %B(k) = B(k) + 1; 
                end
            end
        end
    end 
end

B = zeros(8,8);
BC = zeros(1,8);
for i = 1:8
    for j = 1:8
        for k = 1:8 
            if k~=i && k~=j && d(i,j) == d(i,k) + d(k,j) && d(i,j) < 101 && d(i,j) > 1
                B(i,j) = B(i,j) + 1;
                BC(k) = BC(k) +  B(i,j) / (2 * Q(i,j));
            end
        end
    end
end
plot(G);
disp(BC);