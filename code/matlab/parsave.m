function parsave(fname, varargin)
    % varargin contains alternating variable names and values
    % Extract names and values
    for i = 1:2:length(varargin)
        varname = inputname(i+1);
        if isempty(varname)
            varname = varargin{i};
        end
        eval([varname ' = varargin{i+1};']);
    end
    
    % Create the save command with all variables
    save_cmd = 'save(fname';
    for i = 1:2:length(varargin)
        varname = inputname(i+1);
        if isempty(varname)
            varname = varargin{i};
        end
        save_cmd = [save_cmd ', ''' varname ''''];
    end
    save_cmd = [save_cmd ', ''-v7.3'');'];
    
    % Execute the save command
    eval(save_cmd);
end
