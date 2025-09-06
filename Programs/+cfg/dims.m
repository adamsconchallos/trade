% +cfg/dims.m
classdef dims
    methods (Static)
        function [N,S] = loadDims(projectRoot, dataSub)
            % loadDims  infer N (#countries) and S (#sectors)
            %   [N,S] = cfg.dims.loadDims(projectRoot, dataSub)
            %
            %   projectRoot : your top folder (string)
            %   dataSub     : subfolder under Data (e.g. 'toy CSV' or 'CSV files')
            
            % build full path to tariff_ijs.csv
            fn = fullfile(projectRoot, 'Data', dataSub, 'tariff_ijs.csv');
            assert(exist(fn,'file')==2, 'dims.loadDims: file not found: %s', fn);
            
            raw = csvread(fn);
            N   = size(raw,2);
            S   = size(raw,1)/N;
            assert(abs(S-round(S))<eps, ...
              'dims.loadDims: rows not an integer multiple of columns.');
        end
    end
end
