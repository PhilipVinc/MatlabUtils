function res = AverageExtractData( obj, data, params )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    res.params = params;
    res.params.VERSION = obj.VERSION;
    %%%-----------------------------------------------------------------%%%
    %%%                           Shortcuts                             %%%
    %%%-----------------------------------------------------------------%%%
    nx = res.params.nx;
    ny = res.params.ny;
    nxy = nx*ny;
    
    %%%-----------------------------------------------------------------%%%
    %%%                           Fix Params                            %%%
    %%%-----------------------------------------------------------------%%%
    
    if isfield(params, 'F')
    
        res.params.F = obj.ConvertInterleavedToComplex(params.F, length(params.F) >= 2*3*nxy);
        
        res.params.F_A = params.F(1:3:end);
        res.params.F_B = params.F(2:3:end);
        res.params.F_C = params.F(3:3:end);
        
        res.params.F_tot=sum(abs(params.F_A))+sum(abs(params.F_B))+sum(abs(params.F_C));
        res.params.I_tot=sum(abs(params.F_A).^2)+sum(abs(params.F_B).^2)+sum(abs(params.F_C).^2);

    elseif (isfield(params, 'F_A')&& isfield(params, 'F_B') && isfield(params, 'F_C'))
        res.params.F_A = obj.ConvertInterleavedToComplex(params.F_A, length(params.F_A) >= 2*nxy);
        res.params.F_B = obj.ConvertInterleavedToComplex(params.F_B, length(params.F_B) >= 2*nxy);
        res.params.F_C = obj.ConvertInterleavedToComplex(params.F_C, length(params.F_C) >= 2*nxy);
        
        res.params.F_tot=sum(abs(params.F_A))+sum(abs(params.F_B))+sum(abs(params.F_C));
        res.params.I_tot=sum(abs(params.F_A).^2)+sum(abs(params.F_B).^2)+sum(abs(params.F_C).^2);

    elseif  isfield(params, 'F_t')
        
        res.params.Ftot_t = sum(abs(params.F_t),2);
        res.params.Itot_t = sum(abs(params.F_t).^2,2);
    end
    

    
    %%%-----------------------------------------------------------------%%%
    %%%                         Hardcoded Stuff                         %%%
    %%%-----------------------------------------------------------------%%%
    % this is the point in time where I consider the cut.
    cutFrac = 4/5;
    n_workers = 2;
    computeG2 = true;
    %%%-----------------------------------------------------------------%%%
    
    % Noise stuff
    trajId = obj.GetVarId('traj');
    randomVarNames = obj.varNames(2:end);
    randomVarIds = arrayfun(@obj.GetVarId, randomVarNames);

    
    % Get the length of the traces and number of traces, and get the 
    % cut point.
    t_length = size(data{trajId},2);
    n_traces = size(data{trajId},3);
    t_cut = floor(cutFrac*t_length);
    cutted_frames = t_length-t_cut+1;
    res.params.t_cut = t_cut;
    
    % Compute the average and std of the Random (Disordered) Matrices
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
        res.ave.(nameAvesq) = (rssq(varData,3).^2)./n_traces;
    end

    
    %%%%%%%%%%%%%%%% ---- N(t) ---- %%%%%%%%%%%%%%%%%%%%
    % Compute  n(t) for each trace.
    n_t = data{trajId}.*conj(data{trajId}) -1/2;
    % Average across traces and cavities.
    n_i_t = squeeze(mean(n_t,3));
    n_a_t=squeeze(mean(mean(n_t,3),1));
    % Compute the standard deviation (Spread) of n(t) for every cavity.
    std_n_t = std(n_t,0,3);
    %n_t_avesq = (rssq(n_t,3).^2)./n_traces;
    % Error for the averaged cavities.
    std_n_a_t = sqrt(sum(std_n_t.^2,1))/nx;
    
    % Phase
    phase_i_t = mean(angle(data{trajId}),3);
    
    % G_2(0)
    tmp=mean(mean(abs(data{trajId}).^4,3),1);
    g2_t=(tmp -2 * n_a_t-0.5)./n_a_t.^2;
    g2_tErr = linspace(0,0,length(g2_t)); %TODO
    
    % N_{k=0}
    beta_k0=sum(data{trajId},1);
    nk0_t=(mean(abs(beta_k0).^2/(nxy)-0.5, 3))./(n_a_t*nxy);
    nk0_tErr = linspace(0,0,length(nk0_t)); %TODO
    
    % Population of bands and of k modes
    if (params.nx ~= 1 && params.ny ~= 1)
        fprintf('Computing for 2D: ');
        fprintf('ERROROROROR!!!! ');
        % Our beta_ij_tcut
        %beta_ij_t = reshape(data{trajId}(:,t_cut:end,:), [nx, ny, cutted_frames, n_traces]);
        %n_ij_t = reshape(n_t(:,t_cut:end,:), [nx, ny, cutted_frames, n_traces]);
        %nat_cut = n_a_t(t_cut:end);
        %nat_cut = reshape(nat_cut,[1,1,length(nat_cut)]);
        
        %fprintf('nk..');
        % Compute populations of nk neq 0
        %beta_kxy_t = fft2(beta_ij_t);
        %n_kxy_t = (mean(abs(beta_kxy_t).^2,4)/nxy-0.5)./(nxy*nat_cut);
        
        %beta_k0_t = beta_kxy_t(1,1,:,:);
        
    elseif (params.nx == 1 | params.ny == 1)
        fprintf('Computing for 1D: ');
        
        t_cut_k = 1;
        cutted_frames_k = length(params.t);
        
        nat_cut = n_a_t(t_cut_k:end);

        kx = (0:(nx-1))*pi/nx;
        ky = (0:(ny-1))*pi/ny;

       	% Compute populations of nk neq 0
        fprintf('nk..');
        a_k=fft(data{trajId}(1:3:end,t_cut_k:end,:))*1/sqrt(nx);
        b_k=fft(data{trajId}(2:3:end,t_cut_k:end,:))*1/sqrt(nx);
        c_k=fft(data{trajId}(3:3:end,t_cut_k:end,:))*1/sqrt(nx);
        
        [A,B,C] = obj.Lieb1PartCoeffs(params.J_AB, params.J_BC,nx,1);
        
        phi_b_k_t = zeros(3, nx, cutted_frames_k, n_traces);
        N_b_k_t = zeros(3, nx, cutted_frames_k);
        % For the three bands
        for band=1:3
            phi_b_k_t(band, :,:, :) = bsxfun(@times,A(:,band),a_k) +...
              bsxfun(@times,B(:,band),b_k) + ...
              bsxfun(@times,C(:,band),c_k);
        end
        
        N_b_k_t = mean( (phi_b_k_t.*conj(phi_b_k_t)-1/2) ,4);
        
        N_b_t = squeeze(sum(N_b_k_t,2));
        Nfrac_b_t = N_b_t./sum(N_b_t,1);
        res.params.kx = (0:(nx-1))*pi/nx;
        res.params.ky = (0:(ny-1))*pi/ny;
        %beta_kxy_t = fft(data{trajId}(:,t_cut:end,:));
        %n_kxy_t = (mean(abs(beta_kxy_t).^2,3)/nxy-0.5)./(nxy*nat_cut);
        %beta_k0_t = beta_kxy_t(1,:,:);
    end

    
    % Compute G2ij
    if ( computeG2 == true && (params.nx ~= 1 && params.ny ~= 1))
        
        fprintf('g1/2..');
        % Compute G1 and G2
        g1ij_t = zeros(params.nx+1,params.ny+1,cutted_frames);
        g2ij_t = zeros(params.nx+1,params.ny+1,cutted_frames);
        corr1ij_t = zeros(params.nx+1,params.ny+1,cutted_frames);
        % Set the know edges
        g1ij_t(1,:,:) = 1;
        g1ij_t(:,1,:) = 1;
        g1ij_t(end,:,:) = 1;
        g1ij_t(:,end,:) = 1;
        
        avg_beta_ij_t = mean(beta_ij_t, 4);
        avg_betaStar_ij_t = conj(avg_beta_ij_t);
        avg_n_t =mean(n_ij_t, 4);

        corr1ij_t(1,:,:) = ones(size(corr1ij_t(1,:,1))).*nat_cut./(mean(mean(avg_betaStar_ij_t,1),2).*mean(mean(avg_beta_ij_t,1),2));
        corr1ij_t(:,1,:) = ones(size(corr1ij_t(:,1,1))).*nat_cut./(mean(mean(avg_betaStar_ij_t,1),2).*mean(mean(avg_beta_ij_t,1),2));
        corr1ij_t(end,:,:) = ones(size(corr1ij_t(end,:,1))).*nat_cut./(mean(mean(avg_betaStar_ij_t,1),2).*mean(mean(avg_beta_ij_t,1),2));
        corr1ij_t(:,end,:) = ones(size(corr1ij_t(:,end,1))).*nat_cut./(mean(mean(avg_betaStar_ij_t,1),2).*mean(mean(avg_beta_ij_t,1),2));


        g2ij_t(1,1,:) = g2_t(t_cut:end);
        beta_abs_xy_t = abs(beta_ij_t).^2;
        beta_ij_conj_t = conj(beta_ij_t);

        ny = params.ny;
        
        for rx=2:nx
            for ry=1:ny
                rx_l = rx-1; ry_l = ry-1;
                
                g1ij_t(rx,ry,:) = real(mean(mean(mean(beta_ij_conj_t.*...
                    circshift(circshift(beta_ij_t,rx_l,1),ry_l,2),4),1),2)./nat_cut);
                corr1ij_t(rx,ry,:) = real(mean(mean(mean(beta_ij_conj_t.*...
                    circshift(circshift(beta_ij_t,rx_l,1),ry_l,2),4)./...
                    (avg_betaStar_ij_t.*circshift(circshift(avg_beta_ij_t,rx_l,1),ry_l,2)),1),2));
                g2ij_t(rx,ry,:) = real(mean(mean(mean(beta_abs_xy_t.*...
                    circshift(circshift(beta_abs_xy_t,rx_l,1),ry_l,2),4),1),2))./(nat_cut.^2);
            end
        end       
        for ry=1:nx-1
            rx = 0;
            g1ij_t(1,ry+1,:) = real(mean(mean(mean(beta_ij_conj_t.*...
                circshift(circshift(beta_ij_t,rx,1),ry,2),4),1),2)./nat_cut);
            corr1ij_t(1,ry+1,:) = real(mean(mean(mean(beta_ij_conj_t.*...
                circshift(circshift(beta_ij_t,rx,1),ry,2),4)./...
                (avg_betaStar_ij_t.*circshift(circshift(avg_beta_ij_t,rx,1),ry,2)),1),2));
            g2ij_t(1,ry+1,:) = real(mean(mean(mean(beta_abs_xy_t.*...
                circshift(circshift(beta_abs_xy_t,rx,1),ry,2),4),1),2))./(nat_cut.^2);
        end

    elseif ( computeG2 == true && (params.nx == 1 | params.ny == 1))
        fprintf('g1/2..');
        g1ij_t = zeros(nx+1,cutted_frames);
        corr1ij_t = zeros(nx+1,cutted_frames);
        g2ij_t = zeros(nx+1,cutted_frames);
        
        beta_ij_t = data{trajId}(:,t_cut:end,:);
        beta_abs_xy_t = abs(beta_ij_t).^2;
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
            g2ij_t(r+1,:) = real(mean(mean(beta_abs_xy_t.*circshift(beta_abs_xy_t,r,1),1),3))./(nat_cut.^2);
        end       
        
    else
        n_kxy_t = zeros(1,1);
        g1ij_t = zeros(1,1);
        g2ij_t = zeros(1,1);
        corr1ij_t = zeros(1,1);
    end    
    
    %%%-----------------------------------------------------------------%%%
    %%%                           Save back                             %%%
    %%%-----------------------------------------------------------------%%%
    res.ave.n_i_t = n_i_t;
    res.ave.phase_i_t = phase_i_t;
    res.ave.n_t = n_a_t;  %res.ave.n_t_std = std_n_a_t;
    res.ave.g2_t = g2_t;  %res.ave.g2_t_std = g2_tErr;
    res.ave.nk0_t = nk0_t;  %res.ave.nk0_t_std = nk0_tErr;
    %res.ave.nk_xy_t = n_kxy_t;
    res.ave.N_b_k_t = N_b_k_t;
    res.ave.N_b_t = N_b_t;
    res.ave.Nfrac_b_t = Nfrac_b_t;
        
    % And the cut at the last 3/4 of times, averaged, for the steady state.
    res.ave.n_end = mean(n_a_t(t_cut:end));
    res.ave.std_n_end = sqrt(sum(std_n_a_t(t_cut:end).^2))/( ...
        t_length-t_cut);
    res.ave.g2_end = mean(g2_t(t_cut:end));
    res.ave.g2_endErr = sqrt(sum(g2_tErr(t_cut:end).^2))/( ...
        t_length-t_cut);
    res.ave.nk0_end = mean(nk0_t(t_cut:end));
    res.ave.nk0_endErr = sqrt(sum(nk0_tErr(t_cut:end).^2))/( ...
        t_length-t_cut);
    res.ave.g1ij_t = g1ij_t;
    res.ave.corr1ij_t = corr1ij_t;
    res.ave.g2ij_t = g2ij_t;

    res.params.kx = kx;
    res.params.ky = ky;
    
    
    %%%-----------------------------------------------------------------%%%
    %%%                           Quantities                            %%%
    %%%-----------------------------------------------------------------%%%
    res.quan.n_hist_vals = data{trajId}(:,end,:);


    %%%-----------------------------------------------------------------%%%
    %%%                Special additional analysis Params               %%%
    %%%-----------------------------------------------------------------%%%

    if params.analysisParameters.active == true
        fprintf('Computing additional magic analysis...');
        anPar = params.analysisParameters;
        %pars
        smoothWidth = anPar.smoothWidth;
        tSkip = anPar.tSkip;
        nMinA = anPar.nMinA;  nMinB = anPar.nMinB;  nMinC = anPar.nMinC;
        nMaxA = anPar.nMaxA;  nMaxB = anPar.nMaxB;  nMaxC = anPar.nMaxC;
        nBins = anPar.nBins;

        % do the binning
        nEdgesA = linspace(nMinA, nMaxA, nBins+1);
        dNA = nEdgesA(2)-nEdgesA(1);
        NvalsA = nEdgesA(1:end-1) + diff(nEdgesA);

        nEdgesB = linspace(nMinB, nMaxB, nBins+1);
        dNB = nEdgesB(2)-nEdgesB(1);
        NvalsB = nEdgesB(1:end-1) + diff(nEdgesB);

        nEdgesC = linspace(nMinC, nMaxC, nBins+1);
        dNC = nEdgesC(2)-nEdgesC(1);
        NvalsC = nEdgesC(1:end-1) + diff(nEdgesC);

        % Average across al A,B,C
        nA_t = squeeze(mean(n_t(1:3:end,:,:),1));
        nB_t = squeeze(mean(n_t(2:3:end,:,:),1));
        nC_t = squeeze(mean(n_t(3:3:end,:,:),1));
        
        % smooth the data
        for i=1:n_traces
            nA_t(:,i) = smoothFirst(nA_t(:,i), smoothWidth);
            nB_t(:,i) = smoothFirst(nB_t(:,i), smoothWidth);
            nC_t(:,i) = smoothFirst(nC_t(:,i), smoothWidth);
        end

        % Get only a subset of the data along time
        nA_t = nA_t(1:tSkip:end,:);
        nB_t = nB_t(1:tSkip:end,:);
        nC_t = nC_t(1:tSkip:end,:);

        % fix outliers
        nA_t(nA_t<nMinA) = nMinA;  nA_t(nA_t>nMaxA) = nMaxA;
        nB_t(nB_t<nMinB) = nMinB;  nB_t(nB_t>nMaxB) = nMaxB;
        nC_t(nC_t<nMinC) = nMinC;  nC_t(nC_t>nMaxC) = nMaxC;
        
        %
        tVals = params.t(1:tSkip:end);
        FVals_t = res.params.Ftot_t(1:tSkip:end);

        % Fix for the fact that Fvalues might be slightly asymmetric when going up and down the pulse.
        % Here I take the values going up, and 'interpolate' the values going down to those going up.
        [maxF, maxFId]=max(FVals_t);
        if (maxFId < length(tVals))
            Fgro = FVals_t(1:maxFId);
            Fdes = FVals_t(maxFId+1:end);
            Fdes_corr = interp1(Fgro, Fgro, Fdes, 'nearest');
            FVals_t_corrected = [Fgro; Fdes_corr];
            FVals_t = FVals_t_corrected;
        end

        [Funique, ~, uniqueIds] = unique(FVals_t);
        nFvals = length(Funique);
        PnA = zeros(nFvals,nBins); PnB = PnA; PnC = PnA;


        % iterate over time
        for i=1:size(nA_t, 1)
            [nbinvals_aveA, ~] = histcounts(nA_t(i,:), nEdgesA, 'Normalization', 'count');
            [nbinvals_aveB, ~] = histcounts(nB_t(i,:), nEdgesB, 'Normalization', 'count');
            [nbinvals_aveC, ~] = histcounts(nC_t(i,:), nEdgesC, 'Normalization', 'count');

            Fid = uniqueIds(i);
            PnA(Fid,:) = PnA(Fid,:) + nbinvals_aveA;
            PnB(Fid,:) = PnB(Fid,:) + nbinvals_aveB;
            PnC(Fid,:) = PnC(Fid,:) + nbinvals_aveC;
        end
        
        % normalize
        n = sum(PnA,2)*dNA; % Compute norms of columns
        PnA = bsxfun(@rdivide,PnA,n); % Normalize M
        n = sum(PnB,2)*dNB; % Compute norms of columns
        PnB = bsxfun(@rdivide,PnB,n); % Normalize M
        n = sum(PnC,2)*dNC; % Compute norms of columns
        PnC = bsxfun(@rdivide,PnC,n); % Normalize M

        res.params.NvalsA = NvalsA;
        res.params.NvalsB = NvalsB;
        res.params.NvalsC = NvalsC;
        res.params.Fvals = Funique;
        res.params.Ivals = res.params.Itot_t(uniqueIds);
        res.ave.PnA = PnA;
        res.ave.PnB = PnB;
        res.ave.PnC = PnC;
    end
end

