function output = AverageMergeChunks( obj )
%AVERAGEMERGECHUNKS Summary of this function goes here
%   Detailed explanation goes here
    
    nChunks = length(obj.chunkNumbers);
    chunks = cell(1, nChunks);

    fprintf(['Averaging ', num2str(nChunks), ' average - value chunks.']);
    fprintf(['\tReading First chunk...']);

    aveCnkPath = obj.AveragedChunkPath(obj.chunkNumbers(1));       
    chunks{1} = load(aveCnkPath);

    fprintf(['\tDone!\n\tSetting up arrays...']);
    
    quanNames = fieldnames(chunks{1}.averaged);
    quanNames(strncmp(quanNames, 'params',6)) = []; 
    params = chunks{1}.params;

    %Remove the std quantities from the averaged ones
    stdQuanNames = quanNames(endsWith(quanNames, '_std'));
    endQuanNames = quanNames(endsWith(quanNames, '_end'));
    quanNames(endsWith(quanNames, '_std')) = [];

    aveSqValNames = quanNames(endsWith(quanNames, '_avesq'));

    for i=1:length(aveSqValNames)
        val = aveSqValNames{i};
        tmp=strrep(val, '_avesq', '_std');
        aveSqValname{i} = tmp;
    end

    for jj=1:length(quanNames)
        val = quanNames{jj};
        averaged.(val) = 0;
    end

    stdQuanErrNames=cell(size(stdQuanNames));
    for jj=1:length(stdQuanNames)
        val = stdQuanNames{jj};
        averaged.(val) = 0;
    end

    n_traj_tot = 0;

    fprintf(['Done!\n\tReadAveraging chunks...\t[ 00%% ]']); percentCompletionOld = 0;
    for i=1:nChunks
        if i ~= 1
            aveCnkPath = obj.AveragedChunkPath(obj.chunkNumbers(i));       
            tmp = load(aveCnkPath);
            chunks{i} = tmp;
        end

        n_traj_tot = n_traj_tot + chunks{i}.params.n_traj;
        for jj=1:length(quanNames)
            val = quanNames{jj};
            averaged.(val) = averaged.(val) + chunks{i}.averaged.(val).* ...
                                                chunks{i}.params.n_traj;
        end
    
        for jj=1:length(stdQuanNames)
            val = stdQuanNames{jj};

            averaged.(val) = averaged.(val) + chunks{i}.averaged.(val).* ...
                                    chunks{i}.params.n_traj;
        end

        chunks{i} = [];
        percentCompletion = floor(i/nChunks*10);
        if percentCompletion > percentCompletionOld
            percentCompletionOld = percentCompletion;
            fprintf(['\b\b\b\b\b',num2str(percentCompletion*10), '%% ]']);
        end
    end

    fprintf(['\n\tFinishing...']);
    params.n_traj = n_traj_tot;

    for jj=1:length(quanNames)
        averaged.(quanNames{jj}) = averaged.(quanNames{jj})./n_traj_tot;
    end

    for jj=1:length(stdQuanNames)
        val = stdQuanNames{jj};
        valErr = strrep(val, '_std', '_err');
        averaged.(val) = averaged.(val)./n_traj_tot;
        averaged.(val) = sqrt(averaged.(val));
        averaged.(valErr) = averaged.(val)/sqrt(n_traj_tot);
    end


    % Add the std of different quantities

    
    % Compute the std of quantities for which we got avesq and ave
    for jj=1:length(aveSqValNames)
        val = aveSqValNames{jj};
        val = strrep(val, '_avesq', '_std');
        averaged.(val) = 0;
    end

    for jj=1:length(aveSqValNames)
        val = aveSqValNames{jj};
        val = strrep(val, '_avesq', '_std');
        valErr = strrep(val, '_std', '_err');
        valNormAve = strrep(val, '_std', '_normAve');
        aveSqValname = strrep(val, '_std', '_avesq');
        aveValName = strrep(val, '_std', '_avg');

        if ~isfield(averaged, aveValName)
            aveValName = strrep(val, '_std', '');
        end

        if ~isfield(averaged, valNormAve)
            averaged.(valNormAve) = 1;
        end
        
        if (isfield(averaged, aveSqValname) && isfield(averaged, aveValName))            
            averaged.(val) = sqrt((averaged.(aveSqValname) - (averaged.(aveValName)).^2));
            averaged.(valErr) = averaged.(val)/sqrt(n_traj_tot*averaged.(valNormAve));
        else
            fprintf(['did not compute std for ', aveSqValname, '\n']);
        end % TODO: Must fix for different number of elements!!!
    end

    for jj=1:length(aveSqValNames)
        val = aveSqValNames{jj};
        rmfield(averaged, val);
        rmfield(averaged, strrep(val, '_avesq', '_normAve'));
    end    

    %% calc errors for _end values
    for jj=1:length(endQuanNames)
        val = endQuanNames{jj};
        tVal = strrep(val, '_end', '_t');
        valErr = strrep(val, '_end', '_end_err');
        if isfield(averaged, tVal)
            data = averaged.(tVal);
            sz = size(data);
            if length(sz) == 1
                averaged.(val) = mean(data(params.t_cut:end));
                averaged.(valErr) = std(data(params.t_cut:end));
            elseif length(sz) == 2
                averaged.(val) = mean(data(:,params.t_cut:end));
                averaged.(valErr) = std(data(:,params.t_cut:end));
            elseif length(sz) == 3
                averaged.(val) = mean(data(:,:,params.t_cut:end));
                averaged.(valErr) = std(data(:,:,params.t_cut:end));
            end
        end
    end

    aveDataPath = fullfile(obj.averagedPath, obj.fileAveragedFileName);
    fprintf(['Done!\n\tSaving...']);
    save(aveDataPath,'averaged', 'params');
    fprintf(['Done!\n']);
    
    obj.ave = averaged;
    obj.params = params;
end

