function paramsText = Params2CellTextToml( obj )
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
                paramsText{end+1} = ['"',allKeys{i}, '.file" = "', fname,'"'];
                if ~isreal(tmp)
                    tmpr=real(tmp);
                    tmpi=imag(tmp);
                    ss=size(tmp); ss(end)=ss(end)*2;
                    tmpEx=zeros(ss);
                    tmpEx(:,1:2:end)=tmpr;
                    tmpEx(:,2:2:end)=tmpi;
                    tmp=tmpEx;
                    paramsText{end+1} = ['"',allKeys{i}, '.complex" = true'];
                    tmpStr = sprintf('%i, ', size(tmpr));
                    tmpStr = tmpStr(1:end-2);
                else
                    tmpStr = sprintf('%i, ', size(tmp));
                    tmpStr = tmpStr(1:end-2);
                end

                paramsText{end+1} = ['"',allKeys{i}, '.format" = "text"'];
                dlmwrite(fullfile(obj.simPath,fname), tmp, '\t');
                paramsText{end+1} = ['"',allKeys{i}, '.dimensions" = [ ',tmpStr,' ]'];
            else
                if imag(tmp) == 0
                    tmp = num2str(tmp);
                    paramsText{end+1} = ['"',allKeys{i}, '" = ', tmp];
                    %paramsText{end+1} = ['"',allKeys{i}, '"."value" = ', tmp];
                else % How to write complex numbers
                    %tmp = ['( ', num2str(real(tmp)), ', ',...
                    %            num2str(imag(tmp)), ' )'];
                    paramsText{end+1} = ['"',allKeys{i}, '.complex" = true'];
                    paramsText{end+1} = ['"',allKeys{i}, '.real" = ', num2str(real(tmp))];
                    paramsText{end+1} = ['"',allKeys{i}, '.imag" = ', num2str(imag(tmp))];
                end
            end
        else
            paramsText{end+1} = ['"',allKeys{i}, '" = "', tmp,'"'];
        end
    end

    % fix the toml
    paramsText = sort(paramsText);
end

