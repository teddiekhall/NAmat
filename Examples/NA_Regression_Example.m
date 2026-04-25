clear all; close all; clc;

% Sine Wave Example Function
beta1 = 1.3;
beta2 = 2.7;

x = -5:0.05:5;
y = beta1 + beta2*sin(x); 

% Add noise to the line
y_noise = y + 0.7 * randn(size(y));

% Plot 
figure;
p = plot(x, y);
p.LineWidth = 2;
hold on;
scatter(x, y_noise);
title('Noisy data with the actual function')
hold off;


%%%% Search Phase %%%%

% Define bounds for the two parameters
bounds = [0 5; 0 5];

% Define struct for additional objective function arguments
args.x = x;
args.y = y_noise;

% Define RMSE as the objective function
function output = obj_fn(in, args) 
    pred = in(1) + in(2) * sin(args.x);
    output = rmse(pred, args.y);
end

ns = 10;    % number of samples generated per iteration
nr = 8;     % number of cells to resample
ni = 25;    % number of samples in initial search
N = 25;     % number of iterations
seed = 42;

% Initialize NASearcher object
nas = NASearcher(@obj_fn, ns, nr, ni, N, bounds, args, seed);
nas.run(false);

% Retrieve and plot the best parameters
inds = nas.get_best_indices();
best_params = nas.samples(inds, :);

% Trace plots of both parameters
trace_x = 1:size(nas.samples, 1);

figure;
p = plot(trace_x, nas.samples(:, 1));
p.LineWidth = 1;
hold on;
plot(trace_x, ones(size(trace_x))*beta1, '--')
xlabel('Iteration'); ylabel('Beta 1');
title('Trace Plot: Beta 1')
hold off;

figure;
p = plot(trace_x, nas.samples(:, 2));
p.LineWidth = 1;
hold on;
plot(trace_x, ones(size(trace_x))*beta2, '--')
xlabel('Iteration'); ylabel('Beta 2');
title('Trace Plot: Beta 2')
hold off;

% Voronoi Diagram
figure('Color','white'); hold on;
voronoi(nas.samples(:,1), nas.samples(:,2));
xlabel('Beta 1'); ylabel('Beta 2'); title('NA Voronoi Diagram');
hold off;

%%%% Appraisal Phase %%%%

n_resample = 1000;
n_walkers = 1;
verbose = false;
starting_frac = 0.5;
save = true;

naa = NAAppraiser(n_resample, n_walkers, nas, verbose, seed);
naa.run(save, starting_frac);
mean = naa.mean
mean_err = naa.sample_mean_error
cov = naa.covariance
cov_err = naa.sample_covariance_error

% Plot estimated Multivariate Normal distribution of the parameter space
x1 = 0:.2:5;
x2 = 0:.2:5;
[X1,X2] = meshgrid(x1,x2);
X = [X1(:) X2(:)];

y = mvnpdf(X,mean,cov);
y = reshape(y,length(x2),length(x1));

figure;
scatter(nas.samples(:,1),nas.samples(:,2),[size(nas.objectives,1)],nas.objectives,'filled'); 
hold on;
colorbar;xlabel('beta1'); ylabel('beta2'); title('Misfit NA');
contour(x1,x2,y,[0.0001 0.001 0.01 0.05 0.15 0.25 0.35]);
hold off;
