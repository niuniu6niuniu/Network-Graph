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
% Graph A
 G = graph(A);
 plot(G);
% Diagonal matrix 
% v = zeros(1,8);
% for i = 1:8
%     for j = 1:8
%         if A(i,j) == 1
%             v(1,i) = v(1,i) + 1;
%         end
%     end
% end
% D = diag(v);
% Build a path 
path = zeros(1,8);
% Randomly select a starting point
s = randi([1 8],1);
path(1,1) = s;
disp(s);
% Decide the next stop for starting point s
for j = 1:8
    for i = 1:8
        if A(path(1,i),j) == 1 
            temp = zeros(1,8);
            temp(1,j) = rand;
            for k = 1:8
                if temp(1,k) == max(temp) && k ~= path(1,i)
                    path(1,i+1) = j;
                else
                    path(1,i+1) = 0;
                end
            end
        end
    end
end

disp(path);
        
    
