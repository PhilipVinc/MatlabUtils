function paramsText = Params2CellTextINI( obj )
%CREATEINIFILE Summary of this function goes here
%   Detailed explanation goes here
    
    paramsText = {};
    
    allKeys = keys(obj.params);
    for i=1:length(obj.params)
        tmp = obj.params(allKeys{i});
        if isnumeric(tmp)          
            if (numel(tmp) ~=1)
                % If it's complex, interleave real and complex part
                fname = ['_',allKeys{i},'.dat'];
                %paramsText{end+1} = ['"',allKeys{i}, '" = "', fname,'"'];
                paramsText{end+1} = [allKeys{i}, ' = "', fname, '"'];
                if ~isreal(tmp)
                    tmpr=real(tmp);
                    tmpi=imag(tmp);
                    ss=size(tmp); ss(end)=ss(end)*2;
                    tmpEx=zeros(ss);
                    tmpEx(:,1:2:end)=tmpr;
                    tmpEx(:,2:2:end)=tmpi;
                    tmp=tmpEx;
                    tmpStr = sprintf('%i, ', size(tmpr));
                    tmpStr = tmpStr(1:end-2);
                else
                    tmpStr = sprintf('%i, ', size(tmp));
                    tmpStr = tmpStr(1:end-2);
                end

                dlmwrite(fullfile(obj.simPath,fname), tmp, '\t');
            else
                if imag(tmp) == 0
                    tmp = num2str(tmp);
                    paramsText{end+1} = [allKeys{i}, ' = ', tmp];
                else % How to write complex numbers
                    paramsText{end+1} = ['"',allKeys{i}, '.real" = ', num2str(real(tmp))];
                    paramsText{end+1} = ['"',allKeys{i}, '.imag" = ', num2str(imag(tmp))];
                end
            end
        elseif islogical(tmp)
            if (tmp == true)
                tmpStr = 'true';
            else
                tmpStr = 'false';
            end
            paramsText{end+1} = [allKeys{i}, ' = ', tmpStr];
        elseif (isstr(tmp) || isstring(tmp))
            paramsText{end+1} = [allKeys{i}, ' = "', tmp,'"'];
        else
            paramsText{end+1} = [allKeys{i}, ' = ', tmpStr];
        end
    end
end

