function output = QuantitiesMergeChunks( obj, n_traj_tot )
%QUANTITIESMERGECHUNKS Summary of this function goes here
%   Detailed explanation goes here

    nChunks = length(obj.chunkNumbers);
    chunks = cell(1, nChunks);

    fprintf(['Merging ', num2str(nChunks), ' quantities - value chunks.']);
    fprintf(['\tReading First chunk...']);
        
    aveCnkPath = obj.QuantitiesChunkPath(obj.chunkNumbers(1));       
    tmp = load(aveCnkPath);
    chunks{1} = tmp;

    fprintf(['\tDone!\n\tSetting up arrays...']);
    
    quanNames = fieldnames(chunks{1}.quantities);
    quanNames(strncmp(quanNames, 'params',6)) = []; 
    params = chunks{1}.params;

    for jj=1:length(quanNames)
        val = quanNames{jj};
        quantities.(val) = 0;
    end
    tmpCnk = chunks{1}.quantities;

    for jj=1:length(quanNames)
        val = quanNames{jj};
        if (length(size(tmpCnk.(val))) == 2)
            quantities.(val) = zeros(size(tmpCnk.(val),1),n_traj_tot);
        elseif (length(size(tmpCnk.(val))) == 3)
            quantities.(val) = zeros(size(tmpCnk.(val),1),size(tmpCnk.(val),2),n_traj_tot);
        elseif (length(size(tmpCnk.(val))) == 4)
            quantities.(val) = zeros(size(tmpCnk.(val),1),size(tmpCnk.(val),2), size(tmpCnk.(val),3),n_traj_tot);
        end
    end

    traj_i = 1;

    fprintf(['Done!\n\tReadSumming chunks...\t[ 00%% ]']); percentCompletionOld = 0;
    for i=1:nChunks
        if i ~= 1
            aveCnkPath = obj.QuantitiesChunkPath(obj.chunkNumbers(i));       
            tmp = load(aveCnkPath);
            chunks{i} = tmp;
        end
        data = chunks{i};

        traj_i_end = traj_i+data.params.n_traj-1;
        for jj=1:length(quanNames)
            val = quanNames{jj};
            qty = quantities.(val);
            
            if (length(size(qty)) == 2)
                qty(:,traj_i:traj_i_end) = data.quantities.(val);
            elseif (length(size(qty)) == 3)
                qty(:,:,traj_i:traj_i_end) = data.quantities.(val);
            elseif (length(size(qty)) == 4)
                qty(:,:,:,traj_i:traj_i_end) = data.quantities.(val);
            end
            
            quantities.(val) = qty;
        end
        traj_i = traj_i + data.params.n_traj;

        chunks{i} = [];

        percentCompletion = floor(i/nChunks*10);
        if percentCompletion > percentCompletionOld
            fprintf(['\b\b\b\b\b',num2str(percentCompletion*10), '%% ]']);
        end
    end
   
    aveDataPath = fullfile(obj.averagedPath, obj.fileQuantitiesFileName);
    fprintf(['Done!\n\tSaving...']);
    save(aveDataPath,'quantities', 'params');
    fprintf(['Done!\n']);
    
    obj.quan = quantities;
end

