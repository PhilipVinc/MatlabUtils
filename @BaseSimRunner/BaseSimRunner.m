classdef BaseSimRunner < handle
    %SIMRUNNER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        params = containers.Map('KeyType', 'char', 'ValueType', 'any');
        
        simName = '';
        simPath = '';
        parentPath = '.';

        fileFormat = 'ini';
    end
    
    properties (Access = protected)
        programName = 'sim';
        programRelativePath = '$HOME/bin/';
        iniFileName = '_sim.ini';
        tomlFileName = '_sim.toml';
        nCores = feature('numcores');
        hasName=false;
    end
    
    methods(Abstract, Access = protected)
        res = GenerateSimulationName( obj );
        PreProcessParams(obj);
        PostProcessParams(obj);
        name = GetManagerName(obj);
    end

    methods
        function obj = BaseSimRunner()
            obj.GenerateName();
            %obj.SetTOML();
        end
        
        CreateFolder(obj, folder);
        paramsText = Params2CellTextToml( obj )
        paramsText = Params2CellTextINI( obj )
        SaveTextFileByLine( obj, filePath, data);
        pulseData = CreateTimeDependentPulse( obj);
        cellText = Pulse2CellText( obj, pulseData );
         
        function SetNCores(obj, nCores)
            obj.nCores = nCores;
        end
        
        function CreateSimData(obj)
            if ~isKey(obj.params, 't_start')
                obj.params('t_start') = 0;
            end
            if (~obj.hasName)
                obj.GenerateRandomName()
            end
            
            obj.simPath = fullfile(obj.parentPath, obj.simName);
            obj.CreateFolder(obj.simPath);
            
            obj.params('Manager') = obj.GetManagerName();
            obj.PreProcessParams();
            
            if (~isKey(obj.params, 'n_frames') ||  ~isKey(obj.params, 't_end'))
                if ~isKey(obj.params, 'frames_freq')
                    obj.params('n_frames') = floor(obj.params('t_end') - obj.params('t_start'));
                else
                    obj.params('n_frames') = floor((obj.params('t_end') - obj.params('t_start'))*obj.params('frames_freq'));
                end
            end
            
            if (isKey(obj.params, 'PBC'))
                if ~islogical(obj.params('PBC'))
                    if isnumeric(obj.params('PBC'))
                        obj.params('PBC') = logical(obj.params('PBC'));
                    elseif isstr(obj.params('PBC')) || isstring(obj.params('PBC'))
                        tmp = str2num(obj.params('PBC'))
                        obj.params('PBC') = logical(str2num(obj.params('PBC')));
                    end
                end
            end
            
            obj.PostProcessParams();
            if (strcmp(obj.fileFormat,'ini'))
                paramsText = obj.Params2CellTextINI();
            else
                paramsText = obj.Params2CellTextToml();
            end
            obj.SaveTextFileByLine(obj.IniFilePath(), paramsText);
        end
        
        function Execute(obj)
            commandStr = [obj.programRelativePath, obj.programName, ...
                ' -i ', obj.simPath,'/', obj.iniFileName];
            commandStr = [commandStr, ' max_processes ', num2str(obj.nCores)];
            tic;
            system(commandStr);
            tt = toc;
            
            fprintf(['ELAPSED TIME: ', num2str(tt), ' s.\n']);
        end
        
        function fName = IniFilePath(obj)
            if (strcmp(obj.fileFormat,'ini'))
                fName = obj.iniFileName;
            else
                fName = obj.tomlFileName;
            end
            fName = fullfile(obj.simPath, fName);
        end
        
        function GenerateName(obj)
            try
                obj.simName = obj.GenerateSimulationName();
            catch
                return
            end
            obj.simName = [obj.simName, '_', datestr(now, 'yy-mm-dd_HH-MM-SS')];
            obj.hasName = true;
        end
        
        function SetName(obj, name)
            obj.simName = [name, '_', datestr(now, 'yy-mm-dd_HH-MM-SS')];
            obj.hasName = true;
        end

        function SetTOML(obj)
            obj.fileFormat = 'toml';
        end
        
        function SetProgram(obj, newProgram)
            obj.programName = newProgram;
        end
    end
end

