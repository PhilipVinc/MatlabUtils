function output = ReadAnalizeStoreData( obj )
%READANALIZESTOREDATA Reads and Analizes data, after the path has been correctly 
%set up in the class instance.
%   Detailed explanation goes here
    
    dataPath = fullfile(obj.simPath, obj.simDataFolderName);
    anPath = obj.averagedPath;
    aveDataPath = fullfile(anPath, obj.fileAveragedFileName);
    quanDataPath = fullfile(anPath, obj.fileQuantitiesFileName);
    
    loaded = obj.CheckLoadAnalyzedData();

    if ~loaded
        fprintf(['Analized data is not up to date. Reloading.\n']);
        obj.ReadIniFile();
        obj.ReadDatFiles();
        obj.params.analysisParameters = obj.analysisParameters;
    elseif (loaded && ~obj.keepInMemory)
        fprintf(['Skipping chunk logic.\n']);
        output = 'fast';
        return;
    end
    
    obj.ReadRegisterFast();
    
    obj.chunkData = cell(1, length(obj.chunkNumbers));
    
    % Check if I need to create the folder for the output data
    if ~exist(fullfile(obj.averagedPath, obj.simDataFolderName))
        mkdir(fullfile(obj.averagedPath, obj.simDataFolderName))
    end

    if obj.keepInMemory
        nTrajsToLoad = obj.trajsToKeep;
    else
        nTrajsToLoad = -1;
    end
    nTrajsLoaded = 0;

    % Flag to avoid re-averaging chunks if we did not really recompute new observables.
    analizedSomethingNew = false;


    % Iterate through all chunk numbers
    for i=1:length(obj.chunkNumbers)
        
        chunkId = obj.chunkNumbers(i);
        cIPath = obj.ChunkIndexPath(chunkId);
        aveCnkPath = obj.AveragedChunkPath(chunkId);
        
        % Check if analysis of the chunk is up to date, but we want some trajectories.
        % In this case, just load the few trajectories we want.
        if (obj.IsDataUpToDate(aveCnkPath, cIPath) && ~obj.forceAnalysis && nTrajsToLoad <= 0)
            fprintf(['Skipping chunk #', num2str(chunkId), ': alredy analyzed.\n']);
            continue;
        % if we have to reload because we must
        elseif ((~obj.IsDataUpToDate(aveCnkPath, cIPath) || obj.forceAnalysis) && ~obj.ForceNoAnalysis)
            fprintf( ['Reloading all chunk # ', num2str(chunkId), ' :\n']);
            fprintf( ['\tReading...']);
            nLoaded = obj.ReadAllChunk(i);
            fprintf(['\t Done!\n']);
        % else, I just load because I want the trajectories
        else
            fprintf( ['Reloading partially chunk # ', num2str(chunkId), ' :\n']);
            fprintf( ['\tReading...']);
            nLoaded = obj.ReadAllChunk(i, nTrajsToLoad);
            fprintf(['\t Done!\n']);
        end

        % Check if analysis of the chunk is up to date. In this case, just load the data
        if ((obj.IsDataUpToDate(aveCnkPath, cIPath) && ~obj.forceAnalysis) || obj.ForceNoAnalysis)
            fprintf(['\tSkipping Analysis (alredy done).\n']);
            %continue;
        else 
            fprintf( ['\tAnalizing...']);
            analizedSomethingNew = true;
            res = obj.AverageExtractData(obj.chunkData{i}, obj.params);
            
            fprintf( ['\t Done!\n\tSaving Data...']);
            obj.params = res.params;
            params = obj.params;
            params.n_traj = obj.chunkTrajN(i);
            averaged = res.ave;
            quantities = res.quan;
            if i==1
                save(aveCnkPath, 'averaged', 'params');
            else
                clear params;
                params.n_traj = obj.chunkTrajN(i);
                save(aveCnkPath, 'averaged', 'params');
            end            
            % This fix is needed because thtop has a future date.
            obj.FixEditDate(aveCnkPath, cIPath);
            quanCnkPath = obj.QuantitiesChunkPath(chunkId);
            save(quanCnkPath, 'params', 'quantities');
            obj.FixEditDate(quanCnkPath, cIPath);
            fprintf( ['\t Done!\n']);
            clear params; clear averaged; clear quantities;        
            % params.n_traj = kn
        end      

        % If we do not need to store those trajectories, then remove them from RAM
        if (~obj.keepInMemory || nTrajsToLoad <= 0)
            obj.chunkData{i} = [];
        end
        
        nTrajsLoaded = nTrajsLoaded + nLoaded;
        nTrajsToLoad = nTrajsToLoad - nLoaded;

    end
    
    if (analizedSomethingNew || ~loaded)
        obj.AverageMergeChunks();
        obj.QuantitiesMergeChunks(obj.params.n_traj);
    end

    if obj.keepInMemory
        fprintf(['Reordering trajectories to store...']);
        obj.MergeChunksData();
        fprintf(['\tDone!\n'])
    end
    
    output = 'good';
end

