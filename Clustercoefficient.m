% Cluster coefficient
C = zeros(1,5);
% The number of neighbors each node has
n1 = zeros(1,5);
for i = 1:5
     for j = 1:5
         if d(i,j) == 1 
             n1(i) = n1(i) + 1;
         end
     end
end
% The number of edges between each node's neighbors
n2 = zeros(1,5);
for i = 1:5
    for j =1:5
        for k = 1:5
            if d(i,j) == 1 && d(i,k) == 1 && d(j,k) == 1
                n2(i) = n2(i) + 1;
            end
        end
    end
end
%disp(n1);
%disp(n2);
for i  = 1:5
    C(i) = n2(i) / n1(i) * (n1(i) - 1);
end
disp(C);
