% Computation Time Testing 
% Teddie Hall, NCSU
% April 2026

clear all; close all; clc;
addpath(genpath('..\NA_Matlab_Code'));

%% Simulation Settings
num_sim = 1;
total_param = 16;
samp_sizes = 500:500:5000;

% Define RMSE as the objective function
function output = obj_fn(in, args) 
    pred = args.x * in';
    output = rmse(pred, args.y);
end

%% Function to run MC
function [mean, min_rmse] = run_mc(num_samps, num_params, args, seed)
    rng = RandStream.create("mrg32k3a", "Seed", seed);
    samples = rand(rng, num_samps, num_params) .* 5;

    misfits = zeros(1,num_samps);
    % Begin the model sampling.
    for i = 1:num_samps
        misfits(i) = obj_fn(samples(i, :), args);
    end
    
    [min_rmse, idx] = min(misfits);
    mean = samples(idx,:);
end

%% Function to run NA
function mean = run_na(num_samps, num_params, args, seed)
    bounds = repmat([0 5], num_params, 1);
    ns = 10;                    % number of samples generated per iteration
    nr = floor(0.8 * ns);       % number of cells to resample
    ni = 0.3 * num_samps;       % number of samples in initial search
    N = (num_samps - ni) / ns;  % number of iterations
    
    nas = NASearcher(@obj_fn, ns, nr, ni, N, bounds, args, seed);
    nas.run(false);

    % Initialize and run appraiser 
    num_resample = 500; 
    num_walkers = 1;
    verbose = false;
    naa = NAAppraiser(num_resample, num_walkers, nas, verbose, seed);
    naa.run(false); 
    mean = naa.mean;
end

%% Loops to collect runtimes
method = [];
sample_size = [];
num_params = [];
comp_time = [];
rmse_data = [];
rmse_beta = [];

wb = waitbar(0, "Progress");
for num_param = 2:total_param
    for num_samp = 500:500:5000
        for j = 1:num_sim
            seed = (num_samp/500)*(num_param-2)+(num_samp/500); % seed 1 through 150 for each trial
            
            % Generate data compute y and add noise
            rng = RandStream.create("mrg32k3a", "Seed", seed+150); % seed for generating data in each trial is 150 + seed
            current_x = [ones(50, 1) rand(rng, 50, num_param-1)]; % randomly generate dataset
            current_beta = 5 * rand(rng, 1, num_param); % randomly generate true parameters
            y = current_x * current_beta';
            y_noise = y + 0.3 * randn(size(y));
        
            % Format args
            args.x = current_x;
            args.y = y_noise;

            % run MC and collect data
            tic
            [est_mc, rmse_mc] = run_mc(num_samp, num_param, args, seed);
            comp_time = [comp_time; toc];
            method = [method; "mc"];
            sample_size = [sample_size; num_samp];
            num_params = [num_params; num_param];
            rmse_data = [rmse_data; rmse_mc];
            rmse_beta = [rmse_beta; rmse(current_beta, est_mc)];
    
            % run NA and collect data
            tic
            est_na = run_na(num_samp, num_param, args, seed);
            comp_time = [comp_time; toc];
            method = [method; "na"];
            sample_size = [sample_size; num_samp];
            num_params = [num_params; num_param];
            rmse_data = [rmse_data; obj_fn(est_na, args)];
            rmse_beta = [rmse_beta; rmse(current_beta, est_na)];
    
            waitbar((num_sim*(num_param-2)+j)/(num_sim*(total_param-1)), wb, "Progress");
        end
    end
end

% assemble into table
dat = table(method, sample_size, num_params, comp_time, rmse_data, rmse_beta);

% export to csv
writetable(dat, 'final_comp_times.csv', 'WriteVariableNames', true);
