function estlabel = ProfileLike(A, K, flag) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Based on Zhao et al's algorithm. 
%%%  It needs a guess to start. 
%%%  flag = 0, 1, 2, 3 correspond to a random start, or the initial label vector
%%%             computed from oPCA, nPCA, and SCORE
%%% Reference: 
%%%      Bickel and Chen (2009). A nonparametric view of network models and 
%%%           Newman–Girvan and other modularities. PNAS 106 (50), 21068-21073
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nr, nc] =  size(A);
d = sum(A);
n = length(d);
D = diag(d.^(-1/2));
B = D*A*D;

[vb,db] = eigs(A, K);
%vb1 = vb(:,1);
% vb2 = vb(:,2);
rb = vb(:, 2:K) ./ (vb(:, 1)*ones(1, K-1)); % element-wise ratios


[vc,dc] = eigs(B, K);

% compute initial labels
initRandom = unidrnd(K, [nr,1]); % random start
inita = kmeans(vb,K, 'replicates', 100 ); % oPCA   start (worst?)
initb = kmeans(vc,K, 'replicates', 100 );  % nPCA   start (best?)
initc = kmeans(rb,K, 'replicates', 100 ); %  SCORE start

init = ones(n,1);
if flag == 0
    init = initRandom;
else if flag == 1
        init = inita;
    else if flag == 2
            init = initb;
        else if flag ==3
                init = initc;
            end
        end
    end
end

[XOptimal, VOptimal] = mutiExp(A, init, K);  % maximize the profile likelihood
estlabel = XOptimal; 

