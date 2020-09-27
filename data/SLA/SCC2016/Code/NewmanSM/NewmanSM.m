function label = NewmanSM(A, K)

% By Pengsheng Ji

% Run Newman's Spectral Method recursively up to K communities
% Works for directed and undirected networks.  For directed network, A(i, j)=1 iff an edge from j to i  
% Stop if the extra modularity is equal or less than 0, or the number of
%      communities reaches K; for each iteration split the one with largest
%      positive extra modularity; the actual number of communities may be
%      smaller than K
% MEJ Newman (06): Finding community structure using the eigenvectors of matrices, Physical Review E 74 036104, 2006. 
% MEJ Newman (08): directed network
  
    [nr, nc]=size(A);
    kin = sum(A, 2 );  % row sum  as a col
    kout = sum(A); % col sum as a row
    m = sum(kin); 
    b = kin * kout/m;
    B = A - b;

    lab = ones(nr,K); % labels for each step
    Qg = zeros(1,K); % extra modularity if this subgraph is divided
    for i = 1:(K-1)
    % at this moment, there are already  i communities 
    for j = 1:i 
        % calculate the extra modularity if we split community j
        tmp = (lab(:,i)==j);
        Bg = B(tmp, tmp) - diag(sum(B(tmp,tmp), 2 ));  % eq(9) in Newman(08); correction only on the diagonal
        [xi, ev]=eigs(Bg+Bg',1, 'la'); % use the (positive) largest eigenvalue
        tmpLabel= sign(xi);
        Qg(j) = tmpLabel' * (Bg+Bg') * tmpLabel; %  extra modularity 
        % no split yet
    end
    % find the community with largest extra modularity
    [maxQg, idx] = max( Qg(1:i)); 
    if(maxQg>0)   
         % indeed split the community currently labeled as idx
        tmp = (lab(:,i)==idx);
        Bg = B(tmp, tmp) - diag(sum(B(tmp,tmp), 2 ));  % eq(9) in Newman(08); correction only on the diagonal
        [xi, ev]=eigs(Bg+Bg',1, 'la'); % use the (positive) largest eigenvalue
        tmpLabel= sign(xi);
            if(sum(tmpLabel==-1)==0 || sum(tmpLabel==1)==0 )
                  % only + or - is returned; stop
                   label=lab(:, i);
                    return 
            end

            tmpLabel(tmpLabel==1) = idx; % relabel
            tmpLabel(tmpLabel==-1) = i+1; % relabel
            % make new label
            lab(:,i+1) = lab(:,i);
            lab(tmp, i+1)= tmpLabel;
            % at this moment, there are i+1 communities 
    else
       % no community to split; stop with i communities
       label=lab(:, i);
       return 
    end

    end

    label  = lab(:, K);

   
end

