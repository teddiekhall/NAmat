classdef NASearcher < handle
    % Class for the search phase of the neighborhood algorithm
    % 
    % Properties: 
    % objective: objective function to minimize. Must be defined before 
    %       initializing the NASearcher object and take a single argument
    %       of type array and return a float.
    % ns: int, number of samples generated at each iteration.
    % nr: int, number of cells to resample.
    % ni: int, number of samples from initial random search.
    % n: int, number of iterations.
    % bounds: array, set of tuples representing the bounds of the search
    %       space. Column 1 contains the lower bound and column 2 contains 
    %       the upper bound for each parameter.
    % seed: int, optional seed for the random number generator.
    % nspnr: int, number of samples per cell to generate.
    % nt: int, total number of samples.
    % np: int, running total of number of samples
    % nd: int, number of dimensions 
    % lower: column vector, lower bounds of the search space
    % upper: column vector, upper bounds of the search space
    % Cm: column vector, diagonal prior covariance matrix
    % samples: array, samples generated within the search space
    % objective_args: array, additional arguments required by the objective
    %       function, optional.
    % objectives: array, output of the objective function for each 
    % current_best_ind: int, index of the current minimum of the objectives
    %       array.
    % rngs: RandStream object, random number generator with nr-many streams

    properties
        objective 
        ns 
        nr
        ni
        n
        bounds
        seed
        nspnr
        nt
        np
        nd
        lower
        upper
        Cm
        samples
        objective_args
        objectives 
        current_best_ind
        rngs
    end

    methods
        function obj = NASearcher(objective, ns, nr, ni, n, bounds, args, seed)
            % Constructor for the NASearcher class

            if  nargin > 0
                obj.objective = objective;
                obj.objective_args = args;

                obj.ns = ns;
                obj.nr = nr; 
                obj.nspnr = floorDiv(ns, nr);
                obj.ni = ni;
                obj.n = n;
                obj.nt = ni + n * ns;
                obj.np = 0;

                obj.bounds = bounds; % (2 x nd) array 
                obj.nd = size(bounds, 1);
                obj.lower = bounds(:, 1);
                obj.upper = bounds(:, 2);
                obj.Cm = ((obj.upper - obj.lower).^(-2));

                obj.samples = double.empty(0, obj.nd);
                obj.objectives = Inf(obj.np,1);
                obj.current_best_ind = 0;
                
                if nargin < 8
                    seed = 42;
                end
                obj.rngs = RandStream.create("mrg32k3a", ...
                    "NumStreams", obj.nt, "Seed", seed);
            end
        end

        function run(obj, par) 
            % Runs the search phase

            if nargin < 2
                parallel = true;
            else
                parallel = par;
            end
            
            % Display progress bar
            %wb = waitbar(0, "NAI - Initial Random Search");
                        
            % Initial search
            new_samples = obj.initial_random_search();
            obj.update_ensemble(new_samples);
            
            % Update progress bar
            %waitbar((1/(obj.n+1)), wb, "NAI - Optimisation Loop")

            if canUseParallelPool && (parallel == true)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initialize
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parallel
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pool
            end

            for iter = 1:obj.n
                inds = obj.get_best_indices();
                obj.current_best_ind = inds(1);
                cells_to_resample = obj.samples(inds, :);
                
                if parallel == true
                    parfor i = 1:obj.nr
                        stream_num = (iter - 1) * obj.nr + i;
                        new_samples = obj.random_walk_in_voronoi(cells_to_resample(i, :), inds(i), stream_num);
                        obj.update_ensemble(new_samples);
                    end
                else 
                    for i = 1:obj.nr
                        stream_num = (iter - 1) * obj.nr + i;
                        new_samples = obj.random_walk_in_voronoi(cells_to_resample(i, :), inds(i), stream_num);
                        obj.update_ensemble(new_samples);
                    end
                end
                % Update progress bar
                %waitbar(((iter+1)/(obj.n+1)), wb, "NAI - Optimisation Loop")
            end
                
                
            %close(wb);
        end

        function output = objective_(obj, samples)
            % Executes the objective function
            
            output = obj.objective(samples, obj.objective_args);
        end

        function output = initial_random_search(obj)
            % Executes the initial random search by randomly generating
            % samples within the search space
            
            obj.rngs.Substream = 1; % Set to substream 1

            % Scaling terms for the uniformly distributed random numbers 
            scale1 = repmat((obj.upper - obj.lower).', obj.ni, 1);
            scale2 = repmat(obj.lower.', obj.ni, 1);
            
            output = rand(obj.rngs, obj.ni, obj.nd) .* scale1 + scale2;
        end

        function new_samples = random_walk_in_voronoi(obj, vk, k, rng_stream)
            % Completes a random walk in the Voronoi cells created by the
            % samples
            
            % vk is the current voronoi cell
            % k is the index of the current voronoi cell
            obj.rngs.Substream = rng_stream;
            
            old_samples = obj.samples;
            new_samples = double.empty(0, obj.nd); % Concatenate new samples vertically
            walk_length = obj.nspnr;
            if k == obj.current_best_ind
                walk_length = walk_length + mod(obj.ns, obj.nr);
            end

            % Distance to all other cells (sum across rows)
            d2 = sum(obj.Cm.' .* (vk - old_samples).^2 , 2);
            
            % Distance to previous axis
            d2_previous_axis = 0;
            
            % iterate from 1 to walklength
            for step = 1:walk_length
                xA = vk; % start walk at cell center
                
                for i = 1:obj.nd
                    d2_current_axis = obj.Cm(i) * (xA(i) - old_samples(:,i)).^2;
                    d2 = d2 + d2_previous_axis - d2_current_axis;
                    dk2 = d2(k); % distance of cell center to axis
                    
                    % eqn (19) Sambridge 1999
                    vji = old_samples(:, i);
                    vki = vk(i);
                    a = dk2 - d2;
                    b = vki - vji;

                    adb = a ./ b;
                    adb(b==0) = 0; % accounts for division by zero
                    xji = 0.5 * (vki + vji + adb);
                    
                    % eqns (20, 21) Sambridge 199
                    li = max([obj.lower(i); xji(xji < xA(i))], [], "all");
                    ui = min([obj.upper(i); xji(xji > xA(i))], [], "all");
                    % uniformly distributed random number between li and lu
                    xA(i) = rand(obj.rngs) * (ui - li) + li;

                    d2_previous_axis = d2_current_axis;
                end
            
                new_samples = [new_samples; xA]; 
            end
        end
        
        function output = get_best_indices(obj)
            % Gets the best nr-many indices based on the objectives

            [~, sorted_indices] = sort(obj.objectives);
            output = sorted_indices(1:obj.nr);
        end

        function update_ensemble(obj, new_samples)
            % Adds the new samples generated by the random walk to the ensemble

            num_samples = size(new_samples,1);
            obj.samples = [obj.samples; new_samples];
            
            % Compute objectives on each set of new samples
            new_objectives = arrayfun(@(n) objective_(obj, new_samples(n,:)), 1:size(new_samples, 1))'; 
            
            obj.objectives = [obj.objectives; new_objectives];            
            obj.np = obj.np + num_samples;
        end
    end
end