
% make sure the folders APL, BCPL and NewmanSM are added to the search path


networknames = {'coauthorThresh2Giant'; 'coauthorGiant'; 'citGiant'};
NCluster = [2;3;3];

for(i = 1:3)
  K = NCluster(i);
  inputfile = strcat(networknames{i},'Adj.txt');
  outputfile = strcat(networknames{i},'CommLabelK', num2str(K), 'Matlab.txt');
  if ((exist(inputfile, 'file') == 2)&&(~(exist(outputfile,'file'))))

        adj  = load(inputfile);
        [nauthor, tmp] = size(adj);


        %  Newman Spectral
        rng(0)
        newman = NewmanSM(adj, K) ;

        %
        BCPL = -ones(nauthor,1);
        APL  = -ones(nauthor,1);

        if(i < 3)
            % BCPL  ProfileLike: may take 20 minutes for the network (B)
            rng(10) 
            BCPL  =  ProfileLike(adj, K, 0); % heavily depends on the random start/seed

            % APL pseudo likelihood by Amini
            X = adj;
            % options for the init method and cpl/upl
            rng(0)
            init_opts = struct('verbose',false);
            T = 20;
            cpl_opts = struct('verbose',false,'delta_max',0.1, ...
                              'itr_num',T,'em_max',80,'track_err',false);
            % initilize with SCP, Can use 'bi_deg' instead of 'scp' for degree clustering
            e = initLabel5b(X, K, 'scp', init_opts);     
            % Apply CPl with initial labels e
            APL = cpl4c(X, K, e, [], 'cpl', cpl_opts);

        end
        %  save all results
        dlmwrite(outputfile, [newman, BCPL, APL]) 
        disp(strcat('Success: A file has been created.'))

    end
end

