function output = ReadSimulationRegister(obj)
%READMATLAB Summary of this function goes here
%   Detailed explanation goes here

    registerFile = fopen(fullfile(obj.simPath, 'data', '_register.bin'), 'r');
    
    % Start reading. Read the Magic Bytes 
    magic = fread(registerFile, 8, 'uchar');
    obj.registerVersion = fread(registerFile, 1, 'uchar');
    
    % Check that file is valid
    if (~isequal(magic, obj.magicSignature))
        output = 'Register file is corrupted';
        fprintf(output)
        return
    end
    
    % If file is valid, start reading in.
    % Number of datasets
    obj.varN = fread(registerFile, 1, 'uchar');
    
    % Variable names
    for i=1:obj.varN
        strl = fread(registerFile, 1, 'uchar');
        obj.varNames(i) = string(char(fread(registerFile, strl, 'char')'));
    end
    
    % Dataset Dimensions
    obj.varDimN = zeros(1,obj.varN);
    obj.varDimensions = cell(1, obj.varN);
    for i=1:obj.varN
        tmp = fread(registerFile, 1, 'uint64');
        obj.varDimN(i) = tmp;
        obj.varDimensions{i} = zeros(1,obj.varDimN(i));
        obj.varDimensions{i} = fread(registerFile, obj.varDimN(i), 'uint64')';
    end
    
    % TypeOfVars
    obj.varTypes = fread(registerFile, obj.varN, 'uchar');
    
    % AdditionalDataSize per Traj
    obj.additionalTrajDataByteSize = fread(registerFile, 1, 'uint64');
    obj.trajDataBegin = fread(registerFile, 1, 'uint64');
    obj.additionalTrajDataItems = obj.additionalTrajDataByteSize/8;
    
    obj.registerFile=registerFile;
    output = 'success';
end

