function [XOptimal, VOptimal] = mutiExp(A, init, K) 

   % Coordinate descent search for the optimal membership
   
   % INPUT
   % A:  the adjacency matrix
   % init:  the initial membership vector
   % K:  the number of blocks
   
   % OUTPUT
   % XOptimal:  optimal membership vector
   % VOptimal:  optimal criterial value
   
   
   n = size(A, 1);
   LTenure = min(20,ceil(n/4));
   m = 1000*n;  % times of iteration
   
   % initialization
   X = init;
   XOptimal = X;
   [VOptimal, M] = calCri1(A, X, K);
   
   L = zeros(n, K);
   localTime = 0;   
   t=0;
   
   % "Coordinate descent algorithm": 
   % At each iteration, given the current membership X, 
   % consider all possible changes of moving only one node i to block k. 
   % Search for the optimal change and let the corresponding criteria=VLocal, i=I and k=iterchange. 
   
   while (t<m) && (localTime<1000)  % If the total number of iterations >= m, terminate. 
                                    % At each iteration, if the number of possible pairs (i,k) exceeds 1000, only search 1000 of them. 
      VLocal = -Inf;   
      flagMove = 0;

      for  i = 1:n 
          % For node i, search for optimal k. 
          VK = -Inf;
          change = 0;
          for  k = 1:K
              if (L(i,k)==0) && (X(i)~=k) && (sum(X==X(i))>2)    
                  flagMove = 1;
                  [V, ~] = calCri2(A, X, K, i, k, M);
                  if (V>VK) 
                      VK = V; 
                      change = k;
                  end
              end
          end 
          % If VK > VOptimal, no need to consider other node i. Clear the localTime. 
          if (VK>VOptimal)    
              I = i;
              iterChange = change;
              VOptimal = VK;
              XOptimal = X;
              XOptimal(I) = iterChange;
              localTime = 0;
              break 
          end
          % If VK > VLocal, update VLocal. 
          if (VK>VLocal)   
              I = i;
              iterChange = change;
              VLocal = VK;
          end
      end
      
      % If it is impossible to improve any more, terminate. 
      if (flagMove==0) 
          break
      end
      
      % Update the membership X.   
      localTime = localTime + 1;
      [~, M] = calCri2(A,X,K,I,iterChange,M);
      X(I) = iterChange;
      
      % If the choice "i movs to k" is selected at this iteration, the the algorithm will not select it in next LTenure iterations. 
      L = L-1;
      L(L<0) = 0;
      L(I,iterChange) = LTenure;
      
   %%   t = t+1;
   end