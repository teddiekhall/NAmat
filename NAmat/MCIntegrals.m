classdef MCIntegrals < handle
    % Class for accumulating samples to calculate MC integrals in a single loop
    
    properties
        nd 
        save_samples = false
    end

    properties %(Access = private)
        mi
        mi2
        mimj
        mi2mj
        mi2mj2
        N
        samples
    end

    methods
        function obj = MCIntegrals(nd, save_samples)
            % Constructor for the MCIntegrals class

            if  nargin > 0
                obj.nd = nd; 
                obj.save_samples = save_samples;
                obj.mi = zeros(1,nd); % row vector
                obj.mi2 = zeros(1,nd); % row vector
                obj.mimj = zeros(nd); % square matrix
                obj.mi2mj = zeros(nd); % square matrix
                obj.mi2mj2 = zeros(nd); % square matrix
                obj.N = 0;
                if save_samples
                    obj.samples = []; % empty vector of samples; concatenate to add samples
                end
            end
        end

        function accumulate(obj, x)
            % Determines whether input is MCIntegrals or array object, 
            % then accumulates samples depending on input type

            if isa(x, 'MCIntegrals')
                obj.accumulate_mcintegrals(x);
            else
                obj.accumulate_arraylike(x);
            end
        end

        function obj = accumulate_arraylike(obj, x)
            % Accumulates samples when input is array object
            
            obj.mi = obj.mi + x;
            obj.mi2 = obj.mi2 + x.^2;
            obj.mimj = obj.mimj + (x.' * x);
            obj.mi2mj = obj.mi2mj + ((x.^2).' * x);
            obj.mi2mj2 = obj.mi2mj2 + ((x.^2).' * x.^2);
            obj.N = obj.N + 1;
            if obj.save_samples
                % Concatenates samples and x as column vector
                obj.samples = [obj.samples; x];
            end
        end

        function obj = accumulate_mcintegrals(obj, x)
            % Accumulates samples when input is MCIntegrals object
            
            obj.mi = obj.mi + x.mi;
            obj.mi2 = obj.mi2 + x.mi2;
            obj.mimj = obj.mimj + x.mimj;
            obj.mi2mj = obj.mi2mj + x.mi2mj;
            obj.mi2mj2 = obj.mi2mj2 + x.mi2mj2;
            obj.N = obj.N + x.N;
            if ~isempty(obj.samples) & ~isempty(x.samples)
                % Concatenates samples and x.samples as column vector
                obj.samples = [obj.samples; x.samples];
            end
        end

        function output = mean(obj)
            % Calculates mean of the mi vector
            
            output = obj.mi ./ obj.N;
        end

        function output = sample_mean_error(obj)
            % Calculates sample mean error
            
            new_mi = obj.mi ./ obj.N;
            new_mi2 = obj.mi2 ./ obj.N;
            output = sqrt((new_mi2 - new_mi.^2) ./ obj.N);
        end

        function output = covariance(obj)
            % Calculates covariance
            
            new_mi = obj.mi ./ obj.N;
            new_mimj = obj.mimj ./ obj.N;
            output = new_mimj - new_mi.' * new_mi;
        end

        function output = sample_covariance_error(obj)
            % Calculates sample covariance error
            
            new_mi = obj.mi ./ obj.N;
            new_mi2 = obj.mi2 ./ obj.N;
            new_mimj = obj.mimj ./ obj.N;
            new_mi2mj = obj.mi2mj ./ obj.N;
            new_mi2mj2 = obj.mi2mj2 ./ obj.N;

            output = sqrt( ...
                ( ...
                    new_mi2mj2 ...
                    + (new_mi2.' * new_mi.^2) ...
                    + ((new_mi.^2).' * new_mi2) ...
                    - 2 .* new_mi2mj .* new_mi ...
                    - 2 .* new_mi.' .* new_mi2mj.' ...
                    - 4 .* ((new_mi.^2).' * new_mi.^2) ...
                    + 6 .* new_mimj .* (new_mi.' * new_mi) ...
                    - new_mimj.^2 ...
                )...
                ./ obj.N ...
            );
        end
    end
end
