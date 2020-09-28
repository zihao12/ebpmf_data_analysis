function [value, M] = calCri2(A, X, K, I, new, M)

   % Caculate the criteria value
   
   % INPUT
   % A:  the adjacency matrix
   % X:  the membership vector
   % K:  the number of blocks
   % I:  the node to be moved
   % new:    index of the block to which node I is moved
   % M:  the matrix of total edges between blocks before the new block is
   % formed
   
   % OUTPUT
   % value:  the criteria value
   % M1:  the matrix of total edges between blocks after the new block is
   % formed
 
   
   old = X(I);   % index of the block that node I currently belongs to 
   
   for k = 1:K
       if   (k~=old)  &&  (k~=new)
           s= sum(sum(A(I,X == k))); 
           M(k,new) = M(k,new)+s;
           M(new,k) = M(k,new);
           M(k,old) = M(k,old)-s;
           M(old,k) = M(k,old);
       end
   end
   s = sum(A(I,(X==old)))- A(I,I);
   t = sum(A(I,(X==new)));
   M(old,new) = M(old,new)-t+s;
   M(new,old) = M(old,new);
   M(old,old) = M(old,old)- 2*s - A(I,I);
   M(new,new) = M(new,new)+ 2*t + A(I,I);
   X(I) = new;
   
   value = 0;
   for m = 1:K
    for l = 1:K
       if (M(m,l)> 0) 
           value = value + log(M(m,l)/(sum(M(m,:))*sum(M(:,l))))*M(m,l);
       end 
    end
   end


end