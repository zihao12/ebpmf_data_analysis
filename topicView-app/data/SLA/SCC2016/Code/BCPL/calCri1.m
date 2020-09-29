function [value, M] = calCri1(A, X, K)

   % Caculate the criteria value
   
   % INPUT
   % A:  the adjacency matrix
   % X:  the membership vector
   % K:  the number of blocks
   
   % OUTPUT
   % value:  the criteria value
   % M:  the matrix of total edges between blocks
 
   M = zeros(K);
   for m = 1:K
       for l= 1:K
           M(m,l)= sum(sum(A(X==m,X==l))); 
       end
   end
   value =0;
   for m = 1:K
       for l= 1:K
           if (M(m,l)> 0) 
               value = value + M(m,l)*log( M(m,l)/(sum(M(m,:))*sum(M(:,l))) ); 
           end
       end
   end

end
