# NAmat

``NAmat`` implements the Neighborhood Algorithm (NA), introduced in two papers by Malcolm Sambridge in 1999, in MATLAB for use in geophysical inversion problems. The MATLAB code is based on the [neighpy Python implementation](https://github.com/auggiemarignier/neighpy) developed by Auggie Marignier. The NA searches for the optimal parameters of complex processes in two phases: the search phase and the appraisal phase. The phases are implemented with the classes ``NASearcher.m`` and ``NAAppraiser.m``.

## Installation

Instructions for package installation to be added.

## Basic Usage

```matlab
%%%% Search Phase %%%%

% Define bounds for the parameters
bounds = [0 5; 0 5];

% Define struct for additional objective function arguments
args.x = x;
args.y = y;

% Define objective function
function output = obj_fn(in, args) 
    pred = args.x * in';
    output = rmse(pred, args.y);
end

ns = 10;    % number of samples generated per iteration
nr = 8;     % number of cells to resample
ni = 25;    % number of samples in initial search
N = 25;     % number of iterations
seed = 42;

% Initialize and run NASearcher object
nas = NASearcher(@obj_fn, ns, nr, ni, N, bounds, args, seed);
nas.run(false);

%%%% Appraisal Phase %%%%

% Appraiser arguments
n_resample = 1000;
n_walkers = 1;
verbose = false;
starting_frac = 0.5;
save = true;

% Initialize and run NAAppraiser object
naa = NAAppraiser(n_resample, n_walkers, nas, verbose, seed);
naa.run(save, starting_frac);

% Retrieve estimates
mean = naa.mean;
mean_err = naa.sample_mean_error;
cov = naa.covariance;
cov_err = naa.sample_covariance_error;
```

More detailed instructions can be found in the NAmat User Manual.

## Licence

This code is distributed under a [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.en.html).

## Contributing

If you have any questions, please to open an issue in this repository.