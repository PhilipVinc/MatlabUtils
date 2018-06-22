function isUpdated = IsDataUpToDate( obj, analizedMatPath, dataPath)
%ISDATAUPTODATE Summary of this function goes here
%   Detailed explanation goes here

    isUpdated = false;
    if ~(exist(analizedMatPath, 'file') == 2)
        return
    end

    objToLoad = dir(analizedMatPath);
    try
        dateObj = datetime(objToLoad(1).date, 'Locale',  get(0, 'Language'));
    catch
        dateObj = datetime(objToLoad(1).date, 'Locale',  'en');
    end
    
    if (exist(dataPath, 'file') == 7)
        % get the Folder Name
        folderName =  dir(dataPath);
        folderName = folderName(1).folder;
        folderName = regexp(folderName,'/','split');
        folderName = folderName{end};

        folderToLoad = dir(fullfile(dataPath, '..'));

        for i=1:length(folderToLoad)
            if strcmp(folderToLoad(i).name, folderName)
                try
                    dateFolder = datetime(folderToLoad(i).date, 'Locale',  get(0, 'Language'));
                catch
                    dateFolder = datetime(folderToLoad(i).date, 'Locale',  'en');
                end
                break;
            end
        end
    else
        folderToLoad = dir(dataPath);
        try
            dateFolder = datetime(folderToLoad(1).date, 'Locale',  get(0, 'Language'));
        catch
            dateFolder = datetime(folderToLoad(1).date, 'Locale',  'en');
        end
    end
    
    if dateObj > dateFolder 
        isUpdated = true;
    end
end

