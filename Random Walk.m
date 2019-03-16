% ################# Random Walk #################
% Create a zero matrix as Adjacency matrix
A = zeros(10,10);
% Probability of edge exists between any two nodes
p = 0.8;
% For every pair of nodes check if r > p
for i = 1:10
    for j = 1:10
        r = rand;
        if r > p && i ~= j
            A(i,j) = 1;
            A(j,i) = 1;
        else
            A(i,j) = 0;
            A(j,i) = 0;
        end
    end
end
% The adjacency matrix
disp(A);
% The diagonal matrix
v = zeros(1,10);
for i = 1:10
    for j = 1:10
        if A(i,j) == 1
            v(1,i) = v(1,i) + 1;
        end
    end
end
D = diag(v);
disp(D);
% Laplacian matrix L = D - A
L = D - A;
disp(L);
% The graph
G = graph(A);
plot(G)
