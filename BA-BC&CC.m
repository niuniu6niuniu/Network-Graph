% ################# BA Graph #################
% Adjacency matrix
A = zeros(100,100);
% The boundary for two sets
k = randi([1 100]);
% Probability of edges
p = 0.8;
for i = 1:k
    for j = k+1:100
        r = rand;
        if r > p
            A(i,j) = 1;
            A(j,i) = 1;
        else
            A(i,j) = 0;
            A(j,i) = 0;
        end
    end
end
% disp(r);
% disp(C);
G = graph(A);
plot(G);

% ********************* Betweenness *********************
% The number of shortest path matrix Q
Q = zeros(100,100);
for i = 1:100
    for j = 1:100
        if d(i,j) == 1
            Q(i,j) = 1;
        else
            for k = 1:100
                if d(i,j) > 1 && d(i,j) < 101 && d(k,j) == 1 && d(i,k) == d(i,j) - 1
                    Q(i,j) = Q(i,j) + 1;
                    %B(k) = B(k) + 1; 
                end
            end
        end
    end 
end
% disp(Q);

B = zeros(100,100);
BC = zeros(1,100);
for i = 1:100
    for j = 1:100
        for k = 1:100  
            if k~=i && k~=j && d(i,j) == d(i,k) + d(k,j) && d(i,j) < 101 && d(i,j) > 1
                B(i,j) = B(i,j) + 1;
                BC(k) = BC(k) +  B(i,j) / (2 * Q(i,j));
            end
        end
    end
end

% ********************* Cluster Coefficient *********************
C = zeros(1,100);
% The number of neighbors each node has
n1 = zeros(1,100);
for i = 1:100
     for j = 1:100
         if d(i,j) == 1 
             n1(i) = n1(i) + 1;
         end
     end
end
% The number of edges between each node's neighbors
n2 = zeros(1,100);
for i = 1:100
    for j =1:100
        for k = 1:100
            if d(i,j) == 1 && d(i,k) == 1 && d(j,k) == 1
                n2(i) = n2(i) + 1;
            end
        end
    end
end
%disp(n1);
%disp(n2);
for i  = 1:100
    C(i) = n2(i) / n1(i) * (n1(i) - 1);
end

% The top 20 nodes with highest betweenness centrality 
B = maxk(BC,20);
% The top 20 nodes with highest cluster coefficient
Cm = maxk(C,20);
% Compute the overlaps
count = 0;
for i =1:100
    for j =1:20
        for k = 1:20
            if BC(i) == B(j) && C(i) == Cm(k)
                disp(i);
                count = count + 1;
            end
        end
    end
end
fraction = count / 20;
disp(fraction);