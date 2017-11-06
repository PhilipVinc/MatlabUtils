classdef SimRunner < handle
    %SIMRUNNER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        params = containers.Map('KeyType', 'char', 'ValueType', 'any');
        pulseParams = containers.Map('KeyType', 'char', 'ValueType', 'any');
        
        simName = '';
        simPath = '';
        parentPath = '.';
    end
    
    properties (Access = protected)
        programName = 'sim';
        programRelativePath = '$HOME/bin/';
        iniFileName = '_sim.ini';
        pulseFileName = 'F_t.dat'
        nCores = feature('numcores');
        hasName=false;
    end
    
    methods
        function obj = SimRunner()
            obj.GenerateRandomName();
        end
        
        CreateFolder(obj, folder);
        paramsText = Params2CellText( obj )
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
            
            
            if isKey(obj.params, 'F')
                if (prod(size(obj.params('F'))) > 1)
                    dlmwrite(fullfile(obj.simPath,'_F.dat'), obj.params('F'), '\t');
                    obj.params('F') = '_F.dat';
                end
            else
                pulseData = obj.CreateTimeDependentPulse();
                pulseText = obj.Pulse2CellText(pulseData);
                obj.SaveTextFileByLine(fullfile(obj.simPath, ...
                                        obj.pulseFileName), pulseText);
                obj.params('F_t') = obj.pulseFileName;
            end
            
            if ~isKey(obj.params, 't_end')
                obj.params('t_end') = obj.params('t_start') + pulseData.times(end);
            end
            
            if (~isKey(obj.params, 'n_frames') ||  ~isKey(obj.params, 't_end'))
                if ~isKey(obj.params, 'frames_freq')
                    obj.params('n_frames') = floor(obj.params('t_end') - obj.params('t_start'));
                else
                    obj.params('n_frames') = floor((obj.params('t_end') - obj.params('t_start'))*obj.params('frames_freq'));
                end
            end
            
            paramsText = obj.Params2CellText();
            obj.SaveTextFileByLine(obj.IniFilePath(), paramsText);
        end
        
        function Execute(obj)
            commandStr = [obj.programRelativePath, obj.programName, ...
                ' -i ', obj.simPath];
            commandStr = [commandStr, ' --processes ', num2str(obj.nCores)];
            tic;
            system(commandStr);
            tt = toc;
            
            fprintf(['ELAPSED TIME: ', num2str(tt), ' s.\n']);
        end
        
        function fName = IniFilePath(obj)
            fName = fullfile(obj.simPath, obj.iniFileName);
        end
        
        function GenerateRandomName(obj)
            try
                try
                    Fstr = ['_F-', num2str(num2str(obj.params('F')))];
                catch
                    Fstr='';
                end
                
                obj.simName = ['TWMC_',num2str(obj.params('nx')), 'x', num2str(obj.params('ny')), ...
                    Fstr,...
                    '_ts-',num2str(num2str(obj.params('timestep'))),...
                    '_', datestr(now, 'yy-mm-dd_HH-MM-SS')];
                obj.hasName = true;
            catch
                
                obj.hasName = false;
            end
        end
    end
    
end

