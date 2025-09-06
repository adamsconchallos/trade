classdef root
    %cfg.root  returns the top-level project folder as a constant
    properties (Constant)
        % two .. climb out of .../Programs/+cfg
        path = fileparts(fileparts(mfilename('fullpath')));
    end
end