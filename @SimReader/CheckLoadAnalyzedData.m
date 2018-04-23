function loaded = CheckLoadAnalyzedData( obj )
%CHECKLOADANALYZEDDATA Reads and Analizes data, after the path has been correctly
%set up in the class instance.
%   Detailed explanation goes here

    dataPath = fullfile(obj.simPath, obj.simDataFolderName);
    anPath = obj.averagedPath;
    aveDataPath = fullfile(anPath, obj.fileAveragedFileName);
    quanDataPath = fullfile(anPath, obj.fileQuantitiesFileName);

    loaded = false;

    if ((obj.IsDataUpToDate(aveDataPath, dataPath) && ~obj.forceAnalysis) || obj.ForceNoAnalysis)
        fprintf(['Averaged and Quantities data is up', ...
            ' to date. Loading .mat file.\n']);
        aveData = load(aveDataPath);
        
        if (~obj.ForceNoAnalysis)
            if isfield(aveData.params, 'VERSION')
                if (aveData.params.VERSION < obj.VERSION)
                    fprintf('Ah, no. Analysis VERSION has been bumped. Forcing reanalysis.\n');
                    obj.forceAnalysis = true;
                    return
                end
            else
                if (obj.VERSION ~= 0)
                    fprintf('Ah, no. Analysis VERSION has been bumped. Forcing reanalysis.\n');
                    obj.forceAnalysis = true;
                    return;
                end
            end
            
            if obj.analysisParameters.active == true
                if isfield(aveData.params, 'analysisParameters')
                    if ~isequaln(aveData.params.analysisParameters, obj.analysisParameters)
                        fprintf('Ah, no. Asking special Additional Analysis.\n');
                        obj.forceAnalysis = true;
                        return
                    end
                else
                    fprintf('Ah, no. Asking special Additional Analysis.\n');
                    obj.forceAnalysis = true;
                    return
                end
                fprintf('Special Analysis with those parameters alredy performed.\n');
            end
        end
        
        obj.ave = aveData.averaged;
        obj.params = aveData.params;
        clear aveData;
        quanData = load(quanDataPath);
        obj.quan = quanData.quantities;
        clear quanData;
        
        loaded = true;
    else
        if exist(aveDataPath)
            aveData = load(aveDataPath);
            if (~obj.ForceNoAnalysis)
                if isfield(aveData.params, 'VERSION')
                    if (aveData.params.VERSION < obj.VERSION)
                        fprintf('Ah, no. Analysis VERSION has been bumped. Forcing reanalysis.\n');
                        obj.forceAnalysis = true;
                        return
                    end
                else
                    if (obj.VERSION ~= 0)
                        fprintf('Ah, no. Analysis VERSION has been bumped. Forcing reanalysis.\n');
                        obj.forceAnalysis = true;
                        return;
                    end
                end
            end
            
            if obj.analysisParameters.active == true
                if isfield(aveData.params, 'analysisParameters')
                    if ~isequaln(aveData.params.analysisParameters, obj.analysisParameters)
                        fprintf('Ah, no. Asking special Additional Analysis.\n');
                        obj.forceAnalysis = true;
                        return
                    end
                else
                    fprintf('Ah, no. Asking special Additional Analysis.\n');
                    obj.forceAnalysis = true;
                    return
                end
                fprintf('Special Analysis with those parameters alredy performed.\n');
            end

        else
            obj.forceAnalysis = true;
        end
    end
end