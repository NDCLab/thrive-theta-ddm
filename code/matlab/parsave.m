function parsave(fname, varargin)
% PARSAVE Save variables to a file inside a parfor loop.
%
% Usage:
%   parsave(fname, 'var1', val1, 'var2', val2, ...)
%
% Inputs:
%   fname    - Name of the file to save (string/char).
%   varargin - Alternating variable names (strings) and their values.
%
% Description:
%   Matlab's 'save' function cannot be called directly inside a parfor loop
%   because it accesses the workspace. This function serves as a wrapper.

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
