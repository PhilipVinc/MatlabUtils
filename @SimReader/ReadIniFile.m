function  output = ReadIniFile(obj)
%READINIFILE Summary of this function goes here
%   Detailed explanation goes here
    
    tomlFile = dir(fullfile(obj.simPath,'*.toml'));
    
    if (length(tomlFile) == 1)
        iniFiles = tomlFile;

        if (~length(iniFiles) == 1)
            output = ['Error: Found', num2str(length(iniFiles)),'ini Files. 1 Expected.'];
            fprintf(output);
            return;
        end

        iniFilePath = fullfile(obj.simPath, iniFiles(1).name);
        fprintf(['Loading ini file: ', iniFilePath,'\n']);    

        % Read the ini file with ini2struct.m
        params = toml2struct(iniFilePath);
        fields = fieldnames(params);

        % Convert string with numbers to numbers; save the rest as strings
        for ii=1:length(fields)
            field = fields{ii};
            tmp = params.(field);
            
            if isstruct(tmp) 
                if isfield(tmp, 'file')
                    dataFilePath = fullfile(obj.simPath, tmp.file);
                    try
                        varData = load(dataFilePath);
                        params.(field) = varData;
                        success = true;
                        if isfield(tmp, 'complex')
                            if tmp.complex == true
                                tmp = varData(1:2:end)+1j*varData(2:2:end);
                                params.(field) = tmp;
                            end
                        end
                    catch exception
                    end

                elseif isfield (tmp, 'value')
                    params.(field) = tmp.value;
                elseif isfield (tmp, 'real')
                    params.(field) = tmp.real;
                    if isfield(tmp, 'imag')
                        params.(field)  = params.(field) +1j*tmp.imag;
                    end
                end
            else
                params.(field) = tmp;
            end
        end

        % Check if the local basis size is specified
        if ~isfield(obj.params,'local_cell_size')
            params.local_cell_size=1;
        end

        % Store the results
        obj.params = params;

        output = 'success';

    else
        % Find ini File:
        iniFiles = dir(fullfile(obj.simPath,'*.ini'));

        if (~length(iniFiles) == 1)
            output = ['Error: Found', num2str(length(iniFiles)),'ini Files. 1 Expected.'];
            fprintf(output);
            return;
        end

        iniFilePath = fullfile(obj.simPath, iniFiles(1).name);
        fprintf(['Loading ini file: ', iniFilePath,'\n']);    

        % Read the ini file with ini2struct.m
        params = ini2struct(iniFilePath);
        fields = fieldnames(params);

        % Convert string with numbers to numbers; save the rest as strings
        for ii=1:length(fields)
            field = fields{ii};
            tmp = params.(field);
            
            [val, success] = str2num(tmp);
            if (success && (~isstring(val) || isstr(val) ))
                params.(field) = val;
            elseif ((length(tmp)>3 && strcmp(tmp(end-3:end), '.dat')) || (length(tmp)>4 && strcmp(tmp(end-4:end), '.dat"')))
                if isstring(val)
                    tmp = char(erase(val, '"'));
                end
                dataFilePath = fullfile(obj.simPath, tmp);

                success = false;
                if exist(dataFilePath, 'file')
                    try
                        varData = load(dataFilePath);
                        params.(field) = varData;
                        success = true;
                    catch exception
                    end
                end
                if ~success
                    params.(field) = tmp;
                    fprintf(['INI: Could not find file ', dataFilePath, ...
                        ' for key ', field, '\n']);
                end
            else
                params.(field) = tmp;
            end
        end

        % Check if the local basis size is specified
        if ~isfield(obj.params,'local_cell_size')
            params.local_cell_size=1;
        end

        % Store the results
        obj.params = params;

        output = 'success';
    end
end

