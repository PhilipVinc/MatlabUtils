function res = AverageExtractData( obj, data, params )
%AVERAGEEXTRACTDATA Compute all the observables for a chunk of trajectories
%   This virtual function takes as input a chunk of TWMC Trajectories and
%   computes all average and accumulated (Quntities) observables for those
%   trajectories, storing the result of the computation in a structure 
%   named res. the parameters of the data are stored inside the structure
%   params.
    
    % Store the version of this analysis script inside params.
    res.params = params;
    res.params.VERSION = obj.VERSION;
    %%%-----------------------------------------------------------------%%%
    %%%                           Shortcuts                             %%%
    %%%-----------------------------------------------------------------%%%
    nx = res.params.nx;
    ny = res.params.ny;
    nxy = nx*ny;
    
    %%%-----------------------------------------------------------------%%%
    %%%                         Hardcoded Stuff                         %%%
    %%%-----------------------------------------------------------------%%%
    % this is the point in time where I consider the cut.
    beginCutFrac = 1/5;
    cutFrac = 4/5;
    n_workers = 2;
    computeG2 = true;
    
    %%%-----------------------------------------------------------------%%%
    %%%                   Selecting the right dataset                   %%%
    %%%-----------------------------------------------------------------%%%
    % dataset it of the trajectories themselves
    trajId = obj.GetVarId('traj');
    % dataset ids of the Disordered variables
    randomVarIds = 1:length(obj.varNames);
    randomVarIds(trajId) = [];
    randomVarNames = obj.varNames(randomVarIds);
    randomVarIds = arrayfun(@obj.GetVarId, randomVarNames);

    
    % Get the length of the traces and number of traces, and get the 
    % cut point.
    t_length = size(data{trajId},2);
    n_traces = size(data{trajId},3);
    t_cut = floor(cutFrac*t_length);
    cutted_frames = t_length-t_cut+1;
    res.params.t_cut = t_cut;
    
    %%%-----------------------------------------------------------------%%%
    %%%     Compute the average and std of the Disordered Matrices      %%%
    %%%-----------------------------------------------------------------%%%
    for i=1:length(randomVarNames)
        nameAvg = strrep(char(randomVarNames(i)), '_Realizations', '_avg');
        nameStd = strrep(char(randomVarNames(i)), '_Realizations', '_std');
        nameAvesq = strrep(char(randomVarNames(i)), '_Realizations', '_avesq');
        varData = data{randomVarIds(i)};
        
        if ~(size(varData,2) == 1)
            fprintf('Error, the matrix has more than 1 element');
        end
        
        res.ave.(nameAvg) = mean(varData,3);
        res.ave.(nameStd) = std(varData,0,3);
        % Workaround to no license of signal processing toolbox
        % res.ave.(nameAvesq) = (rssq(varData,3).^2)./n_traces;
        res.ave.(nameAvesq) = sum(varData.^2, 3)./n_traces;
    end

    
    %%%-----------------------------------------------------------------%%%
    %%%                    Compute n(t) and similar                     %%%
    %%%-----------------------------------------------------------------%%%
    % Compute  n(t) for each trace.
    n_t = data{trajId}.*conj(data{trajId}) -1/2;
    % Average across traces and cavities.
    n_i_t = squeeze(mean(n_t,3));
    n_i_t_ave = squeeze(mean(n_t,3));
    n_i_t_avesq = squeeze(mean(n_t.^2, 3));

    n_a_t=squeeze(mean(mean(n_t,3),1));
    n_a_t_ave =squeeze(mean(mean(n_t,3),1));
    n_a_t_avesq =squeeze(mean(mean(n_t.^2,3),1));
    n_a_t_normAve = size(n_t, 1);
    
    % Compute the standard deviation (Spread) of n(t) for every cavity.
    % n_i_t_std = std(n_t,0,3);
    % std_n_i_t=n_i_t_std;
    % n_a_t_std = 
    % n_t_avesq = (rssq(n_t,3).^2)./n_traces;
    % Error for the averaged cavities.
    % std_n_a_t = sqrt(sum(std_n_i_t.^2,1))/nx;
    
    % G_2(0)
    tmp=mean(mean(abs(data{trajId}).^4,3),1);
    g2_t=(tmp -2 * n_a_t-0.25)./n_a_t.^2;
    g2_tErr = linspace(0,0,length(g2_t)); %TODO
    
    % N_{k=0}
    beta_k0=sum(data{trajId},1);
    nk0_t=(mean(abs(beta_k0).^2/(nxy)-0.5, 3))./(n_a_t*nxy);
    nk0_tErr = linspace(0,0,length(nk0_t)); %TODO
    
    % Correlation function if it is 2D
    if (params.nx ~= 1 && params.ny ~= 1)
        fprintf('Computing for 2D: ');
        % Our beta_ij_tcut
        beta_ij_t = reshape(data{trajId}(:,t_cut:end,:), [nx, ny, cutted_frames, n_traces]);
        n_ij_t = reshape(n_t(:,t_cut:end,:), [nx, ny, cutted_frames, n_traces]);
        nat_cut = n_a_t(t_cut:end);
        nat_cut = reshape(nat_cut,[1,1,length(nat_cut)]);
        
        fprintf('nk..');
        % Compute populations of nk neq 0
        beta_kxy_t = fft2(beta_ij_t);
        n_kxy_t = (mean(abs(beta_kxy_t).^2,4)/nxy-0.5)./(nxy*nat_cut);
        
        beta_k0_t = beta_kxy_t(1,1,:,:);
        
    elseif (params.nx == 1 | params.ny == 1)
        fprintf('Computing for 1D: ');
        nat_cut = n_a_t(t_cut:end);
        
       	% Compute populations of nk neq 0
        fprintf('nk..');
        beta_kxy_t = fft(data{trajId}(:,t_cut:end,:));
        n_kxy_t = (mean(abs(beta_kxy_t).^2,3)/nxy-0.5)./(nxy*nat_cut);
        beta_k0_t = beta_kxy_t(1,:,:);
    end

    
    % Compute G2ij
    if ( computeG2 == true && (params.nx ~= 1 && params.ny ~= 1))
        
        fprintf('g1/2..');
        % Compute G1 and G2
        g1ij_t = zeros(params.nx,params.ny,cutted_frames);
        g2ij_t = zeros(params.nx,params.ny,cutted_frames);
        corr1ij_t = zeros(params.nx,params.ny,cutted_frames);
        corr2ij_t = zeros(params.nx,params.ny,cutted_frames);
        % Set the know edges
        g1ij_t(1,:,:) = 1;
        g1ij_t(:,1,:) = 1;
        
        avg_beta_ij_t = mean(beta_ij_t, 4);
        avg_betaStar_ij_t = conj(avg_beta_ij_t);
        avg_beta_ij = mean(avg_beta_ij_t, 3);
        avg_betaStar_ij = mean(avg_betaStar_ij_t,3);

        avg_n_t =mean(n_ij_t, 4);

        corr1ij_t(1,1,:) = ones(size(corr1ij_t(1,1,:))).*nat_cut./(mean(mean(avg_betaStar_ij_t,1),2).*mean(mean(avg_beta_ij_t,1),2));

        g2ij_t(1,1,:) = g2_t(t_cut:end);
        beta_abssq_ij_t = abs(beta_ij_t).^2;
        avg_beta_abssq_ij_t = mean(beta_abssq_ij_t, 4);
        beta_ij_conj_t = conj(beta_ij_t);
        nij_t = beta_abssq_ij_t -1/2;
        sq_avg_nij_t = sqrt(mean(beta_abssq_ij_t,4) -1/2);

        n_ij_t = beta_abssq_ij_t - 0.5;
        avg_n_ij_t = mean(n_ij_t, 4);

        for dx=0:nx-1   % for PBC: dx = nx = 0
            for dy=0:ny-1   % for PBC: dx = nx = 0
                % First element is alredy set.
                % r'=r+d
                % g1ij_t = \sum_r <a_r^\dag a_r'>/sqrt(<n_r><n_r'>)
                                        %   y    x   
                g1ij_t(dx+1,dy+1,:) = mean(mean( ...%traj TR[ 
                                                        mean(mean(...
                                                        (beta_ij_conj_t.*    ...% numerator: a_r^\dag          % avg ]Tr
                                                                circshift(circshift(beta_ij_t, dx, 1), dy, 2) ) ,4),3)./ ...% a_r' % end numerator
                                                            (mean(sq_avg_nij_t,3).* ...
                                                                circshift(circshift(mean(sq_avg_nij_t,3),dx,1),dy,2) ) ...% denom sqrt(<n_r><n_r'>)
                                                           ,1),2); % avg across xy and trajs
                                                        %   x  y   real

                % corr1ij_t = \sum_r <a_r^\dag a_r'> - <a_r><(a_r')^*>) (Alberto's quantity)
                                        %       y   x        
                corr1ij_t(dx+1,dy+1,:) = mean(mean(... % Tr[
                                                            mean(mean(...
                                                             (beta_ij_conj_t.*    ...% numerator: a_r^\dag
                                                                circshift(circshift(beta_ij_t, dx, 1), dy, 2) ) ,4),3),1),2) -...%a_r'
                                              mean(mean(...
                                                             (avg_beta_ij.*...
                                                                circshift(circshift(avg_betaStar_ij,dx,1),dy,2) ) ...
                                                            ,1),2);% avg across xy and trajs
                                                        %    x  y   real

                % corr1ij_t = \sum_r <a_r^\dag a_r a_r'^\dag a_r'> - <n_r><n_r'>) (Alberto's g2)
                                        %       y   x        
                corr2ij_t(dx+1,dy+1,:) = mean(mean(... % Tr[
                                                            mean(mean(...
                                                             ((n_ij_t).*    ...% numerator: a_r^\dag
                                                                circshift(circshift(n_ij_t, dx, 1), dy, 2) ) ,4),3),1),2) -...%a_r'
                                              mean(mean(...
                                                             (avg_n_ij_t.*...
                                                                circshift(circshift(avg_n_ij_t,dx,1),dy,2) ) ...
                                                            ,1),2);% avg across xy and trajs
                                                        %    x  y   real

                g2ij_t(dx+1,dy+1,:) = real(mean(mean(mean(beta_abssq_ij_t.*...
                    circshift(circshift(beta_abssq_ij_t,dx,1),dy,2),4),1),2))./(nat_cut.^2);
            end
        end       

        % dx = 0 dy = 0
        g1ij_t(1,1,:) = g1ij_t(1,1,:)-1/2;
        g2ij_t(1,1,:) = g2ij_t(1,1,:)-mean(mean(mean(mean(n_ij_t,4),3),1),2)./(nat_cut.^2);
        corr1ij_t(1,1,:) = corr1ij_t(1,1,:) - 1/2;
        corr2ij_t(1,1,:) = corr2ij_t(1,1,:) - mean(mean(mean(mean(n_ij_t,4),3),1),2);


        corr1ij_a = mean(corr1ij_t, 3);
        corr2ij_a = mean(corr2ij_t, 3);
        distX = [0:floor(nx/2) , floor((nx-1)/2):-1:1];
        distY = [0:floor(ny/2) , floor((ny-1)/2):-1:1];
        NormFac =sum(sum(corr1ij_a));
        NormFac_g2 =sum(sum(corr2ij_a));
        V = 0;
        V_g2 = 0;
        for i=1:nx
            for j=1:ny
                V = V + corr1ij_a(i,j)*(distX(i)^2+distY(j)^2);
                V_g2 = V_g2 + corr2ij_a(i,j)*(distX(i)^2+distY(j)^2);
            end
        end
        V_nxny = V/(nx*ny); 
        V_g2_nxny = V_g2/(nx*ny); 
        V = V/NormFac;
        V_g2 = V_g2/NormFac_g2;

    elseif ( computeG2 == true && (params.nx == 1 | params.ny == 1))
        fprintf('g1/2..');
        g1ij_t = zeros(nx+1,cutted_frames);
        corr1ij_t = zeros(nx+1,cutted_frames);
        g2ij_t = zeros(nx+1,cutted_frames);
        
        beta_ij_t = data{trajId}(:,t_cut:end,:);
        beta_abssq_ij_t = abs(beta_ij_t).^2;
        nat_cut = n_a_t(t_cut:end);
        
        avg_beta_ij_t = mean(beta_ij_t, 3);
        avg_betaStar_ij_t = conj(avg_beta_ij_t);

        g1ij_t(1,:) = 1;
        g1ij_t(end,:) = 1;
        corr1ij_t(1,:) = (nat_cut)./(mean(avg_betaStar_ij_t,1).*mean(avg_beta_ij_t,1));
        corr1ij_t(end,:) = (nat_cut)./(mean(avg_betaStar_ij_t,1).*mean(avg_beta_ij_t,1));
        
        g2ij_t(1,:) = g2_t(t_cut:end);
        g2ij_t(end,:) = g2_t(t_cut:end);
        
        for r=1:nx-1
            g1ij_t(r+1,:) = real(mean(mean(conj(beta_ij_t).*circshift(beta_ij_t,r,1),1),3)./nat_cut);
            corr1ij_t(r+1,:) = real(mean(mean(conj(beta_ij_t).*circshift(beta_ij_t,r,1),3)./...
                (avg_betaStar_ij_t.*circshift(avg_beta_ij_t, r, 1)),1));
            g2ij_t(r+1,:) = real(mean(mean(beta_abssq_ij_t.*circshift(beta_abssq_ij_t,r,1),1),3))./(nat_cut.^2);
        end       
        
        corr1ij_a = mean(corr1ij_t, 2);
        distX = [0:floor(nx/2) , floor((nx-1)/2):1];
        NormFac =sum(sum(corr1ij_t));
        V = 0;
        for i=1:nx
            V = V + corr1ij_a(i)*distX(i)^2;
        end
        V_nxny = V/(nx*ny); 
        V = V/NormFac

    else
        n_kxy_t = zeros(1,1);
        g1ij_t = zeros(1,1);
        g2ij_t = zeros(1,1);
        corr1ij_t = zeros(1,1);
        corr2ij_t = zeros(1,1);
        V_nxny = zeros(1,1);
        V = zeros(1,1);
        V_g2_nxny = zeros(1,1);
        V_g2 = zeros(1,1);
    end    
    
    %%%-----------------------------------------------------------------%%%
    %%%                           Save back                             %%%
    %%%-----------------------------------------------------------------%%%
    res.ave.n_i_t = n_i_t;  res.ave.n_i_t_avesq = n_i_t_avesq; 
    res.ave.n_t = n_a_t; res.ave.n_t_avesq = n_a_t_avesq; res.ave.n_t_normAve = n_a_t_normAve;
    res.ave.n_end = mean(n_a_t(t_cut:end));

    res.ave.g2_t = g2_t;  %res.ave.g2_t_std = g2_tErr;
    res.ave.nk0_t = nk0_t;  %res.ave.nk0_t_std = nk0_tErr;
    res.ave.nk_xy_t = n_kxy_t;

    % And the cut at the last 3/4 of times, averaged, for the steady state.
    %res.ave.std_n_end = sqrt(sum(std_n_a_t(t_cut:end).^2))/( ...
    %    t_length-t_cut);
    res.ave.g2_end = mean(g2_t(t_cut:end));
    res.ave.g2_endErr = sqrt(sum(g2_tErr(t_cut:end).^2))/( ...
        t_length-t_cut);
    res.ave.nk0_end = mean(nk0_t(t_cut:end));
    res.ave.nk0_endErr = sqrt(sum(nk0_tErr(t_cut:end).^2))/( ...
        t_length-t_cut);

    res.ave.g1ij_t = g1ij_t;
    res.ave.g2ij_t = g2ij_t;
    res.ave.corr1ij_t = corr1ij_t;
    res.ave.corr2ij_t = corr2ij_t;
    
    res.ave.V_nxny = V_nxny;
    res.ave.V = V;
    res.ave.V_g2_nxny = V_g2_nxny;
    res.ave.V_g2 = V_g2;

    
    %%%-----------------------------------------------------------------%%%
    %%%                           Quantities                            %%%
    %%%-----------------------------------------------------------------%%%
    fprintf('quans..');
    res.quan.n_hist_vals = n_t(:,end,:);

    t_cut = floor(t_length*0.95);
    res.quan.n_hist_tave_vals = mean(n_t(:,t_cut:end,:),2);
    
    %% Bistability Coefficient
    fprintf('bist..');
    t_begin = ceil(t_length*1/20);
 %   res.quan.bist_coeff_nt = moment(n_t(:,t_begin:end,:),2,2).^2./moment(n_t(:,t_begin:end,:),4,2);
 %   res.quan.bist_coeff_alphat = abs(moment(data{trajId}(:,t_begin:end,:),2,2)).^2./abs(moment(data{trajId}(:,t_begin:end,:),4,2));
    
    fprintf('s..');
    n_t = smoothDirDumb(n_t, 2, 30);
    alpha_t = smoothDirDumb(data{trajId}, 2, 30);
    fprintf('bs..');
    res.quan.bist_coeff_nt_smoothed = moment(n_t(:,t_begin:end,:),2,2).^2./moment(n_t(:,t_begin:end,:),4,2);
    res.quan.bist_coeff_alphat_smoothed = abs(moment(alpha_t(:,t_begin:end,:),2,2)).^2./abs(moment(alpha_t(:,t_begin:end,:),4,2));
 
end

