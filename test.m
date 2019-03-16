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
% disp(A);
% Diagonal matrix 
v = zeros(1,8);
for i = 1:8
    for j = 1:8
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
t1 = zeros(7,7);
for i = 1:4
    for j = 1:4
        t1(i,j) = L(i,j);
    end
end
for i = 1:4
    for j = 5:7
        t1(i,j) = L(i,j+1);
    end
end
for i = 5:7
    for j = 1:4
        t1(i,j) = L(i+1,j);
    end
end   
for i = 5:7
    for j = 5:7
        t1(i,j) = L(i+1,j+1);
    end
end
disp(t1);        
% Compute the inverse of D - A (t2 = inverse(Dv - Av))
t2 = inv(t1);
% disp(t2);
% add zeros to 5th row and 5th column of t2 ( T )
T = zeros(8,8);
for i = 1:4
    for j = 1:4
        T(i,j) = t2(i,j);
    end
end
for i = 1:4
    for j = 5:7
        T(i,j+1) = t2(i,j);
    end
end
for i = 5:7
    for j = 1:4
        T(i+1,j) = t2(i,j);
    end
end
for i = 5:7
    for j = 5:7
        T(i+1,j+1) = t2(i,j);
    end
end
disp(T);
% Assume node1 for starting point ( s ), node 8 for ending point ( t )
% V = zeros(8,8);
% for s = 1:7
%     for t = 2:8
%         if s < t
%             for i = 1:8 && i~= s && i~= t
%                 for j = 1:8 && i~= s && i~= t
%                     V(i,j) = T(i,s) - T(i,t);
%                 end
%             end
%         end
%     end
% end
% disp(V);
% The current flows through each vertex denoted as I
B = zeros(1,8);
for s = 1:8
    for t = 1:8
        if s < t 
            I = zeros(1,8);
            for i = 1:8 
                for j = 1:8 
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
b = zeros(1,8);
for i = 1:8
    b(i) = B(i) / (0.5 * 8 * 7);
end
% disp(b);
G = graph(A);
plot(G);
