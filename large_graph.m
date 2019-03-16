% ################# Random Graph #################
% Create a zero matrix as Adjacency matrix
A = zeros(50,50);
% Probability of edge exists between any two nodes
p = 0.9;
% For every pair of nodes check if r > p
for i = 1:50
    for j = 1:50
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
% The adjacency matrix of the random graph
% disp(A);
G = graph(A);
% Diagonal matrix 
v = zeros(1,50);
for i = 1:50
    for j = 1:50
        if A(i,j) == 1
            v(1,i) = v(1,i) + 1;
        end
    end
end
D = diag(v);
% disp(D);
% Construct matrix L = D - A
L = D - A;
% disp(L);
% Remove the 5th row and 5th column of D - A (t1 = Dv - Av)
t1 = zeros(49,49);
for i = 1:4
    for j = 1:4
        t1(i,j) = L(i,j);
    end
end
for i = 1:4
    for j = 5:49
        t1(i,j) = L(i,j+1);
    end
end
for i = 5:49
    for j = 1:4
        t1(i,j) = L(i+1,j);
    end
end   
for i = 5:49
    for j = 5:49
        t1(i,j) = L(i+1,j+1);
    end
end
% disp(t1);        
% Compute the inverse of D - A (t2 = inverse(Dv - Av))
t2 = inv(t1);
% disp(t2);
% add zeros to 5th row and 5th column of t2 ( T )
T = zeros(49,49);
for i = 1:4
    for j = 1:4
        T(i,j) = t2(i,j);
    end
end
for i = 1:4
    for j = 5:49
        T(i,j+1) = t2(i,j);
    end
end
for i = 5:49
    for j = 1:4
        T(i+1,j) = t2(i,j);
    end
end
for i = 5:49
    for j = 5:49
        T(i+1,j+1) = t2(i,j);
    end
end
% disp(T);
% The current flows through each vertex denoted as I
B = zeros(1,50);
for s = 1:50
    for t = 1:50
        if s < t 
            I = zeros(1,50);
            for i = 1:50 
                for j = 1:50 
                    if i~= s && i~= t
                    I(i) = 0.5 * A(i,j) * abs( T(i,s) - T(i,t) - T(j,s) + T(j,t) );
                    I(s) = 1;
                    I(t) = 1;
                    % disp(I);
                    B(i) = B(i) + I(i);
                    end
                end
            end
        end
    end
end
% disp(B);
% The betweenness centrality of each vertex
b = zeros(1,50);
for i = 1:50
    b(i) = B(i) / (0.5 * 50 * 49);
end
disp(b);
plot(G);