% printLogMessage Print a formatted message with timestamp to the screen
%
%   printLogMessage(formatSpec, A1, ..., An) prints a log message from given
%       inputs. The arguments are of the same form as for the printf-family,
%       i.e., a format specifier (string) with placeholders followed by an
%       arbitrary number of arguments replacing those placeholders. The
%       resulting string is printed to the screen prefixed by a time stamp for
%       logging purposes.
% 
%   See also: fprintf, sprintf

function printLogMessage(formatSpec, varargin)

arguments
    formatSpec (1,1) string
end
arguments (Repeating)
    varargin
end

fprintf('%s - %s\n', datestr(datetime('now'), 'HH:MM:SS'), sprintf(formatSpec, varargin{:}));

end