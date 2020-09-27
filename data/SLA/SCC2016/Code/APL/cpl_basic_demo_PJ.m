%% Creating the block model
n = 100; % number of nodes
K = 3; % number of communities 

X12 = sparse(rand(n/2) < 0.05); 
X11 = sparse(rand(n/2) < 0.1); X11 = X11 | X11';
X22 = sparse(rand(n/2) < 0.1); X22 = X22 | X22';
X = [X11, X12; X12', X22]; % Adjacency matrix

figure(1), clf
spy(X)




%%
% options for the init method and cpl/upl
rng(0)
init_opts = struct('verbose',false);
T = 20;
cpl_opts = struct('verbose',false,'delta_max',0.1, ...
                  'itr_num',T,'em_max',80,'track_err',false);

%% initilize with SCP, Can use 'bi_deg' instead of 'scp' for degree clustering
e = initLabel5b(X, K, 'scp', init_opts);     


%% Apply CPl with initial labels e
SLLabel = cpl4c(X, K, e, [], 'cpl', cpl_opts);


%% sequential:  no result from the second step
K=2
e = initLabel5b(X, K, 'scp', init_opts);     
Label = cpl4c(X, K, e, [], 'cpl', cpl_opts);
tmp = cpl4c(X(Label==2, Label==2), K, e, [], 'cpl', cpl_opts);
Label(Label==2)  =  1+ tmp;


%%
set(figure(2),'position',[200 300 800 200]), clf, hold on
stem(e,'bo')
stem(Label,'r.')


