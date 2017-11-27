function paramsText = Params2CellText( obj )
%CREATEINIFILE Summary of this function goes here
%   Detailed explanation goes here
    
    paramsText = cell(1,length(obj.params));
    
    allKeys = keys(obj.params);
    for i=1:length(obj.params)
        tmp = obj.params(allKeys{i});
        if isnumeric(tmp)          
            if (numel(tmp) ~=1)
                fname = ['_',allKeys{i},'.dat'];
                dlmwrite(fullfile(obj.simPath,fname), tmp, '\t');
                tmp = fname;
            else
                if imag(tmp) == 0
                    tmp = num2str(tmp);
                else % How to write complex numbers
                    tmp = ['( ', num2str(real(tmp)), ', ',...
                                num2str(imag(tmp)), ' )'];
                end
            end
        end
        paramsText{i} = [allKeys{i}, ' = ', tmp];
    end
    
end

