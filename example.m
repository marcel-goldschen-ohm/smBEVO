%% load data series with baseline drift
load('sample_data_series.mat');

%% plot the data
datafig = figure;
plot(data);

%% smBEVO required parameters

% !!! sigmaX and sigmaY are the most critical parameters and they are
% entirely dependent on the data being analyzed. If you don't like the
% results, start by changing these before you worry about anything else.

% sigmaX:
% The standard deviation of the image filter in x. Since data is a single
% column of y values, the x values default to sample indices. If we had
% supplied a column of x values, we would use those units. A good estimate
% for an appropriate value of sigmaX is the longest period where we can be
% reasonably sure that the baseline won't change much within this period
% wherever we look in the data. Here we choose 100 samples.
sigmaX = 100;

% sigmaY:
% The standard deviation of the image filter in y. A good estimate for an
% appropriate value of sigmaY is ~1/3 of the smallest separtation in y
% between levels that we wish to detect. This will filter noise and not
% levels. Note that good baseline estimation does not always require
% detecting every last minute level. Here the smallest level separation is
% ~0.2, so we choose something like 0.08.
sigmaY = 0.08;

%% smBEVO optional parameters

% We'll use the default of 4 pixels/standard deviation for the image
% representation of the data.

% Smoothing can help to avoid erroneous level jumps do to short bits of
% noisy or hard to decipher data.
% 0 => no smoothing
% >0 => more smoothing

% We'll use the default minimum level separation.

% We'll use active contour snakes to refine the estimated levels based on
% the image. We'll use the default parameters for alpha and beta (i.e.
% constraints on elasticity and curvature of the levels) and gamma (scaling
% for snake movement during each iteration) with the default maximum number
% of iterations.

%% Run smBEVO on data

% The above parameters are very intuitive to estimate, but you may have to
% play around a bit to get the result that you like best. You can of course
% do this via a script as shown below, but the app UI interface makes a
% quick parameter scan simple by providing immediate visual feedback on the
% results, and thus is the preferred method for exploratory analysis.

result = smBEVO(data, sigmaX, sigmaY, 'smoothing', 3, 'snakeLevelRefinement', true);

% Plot data overlaid with baseline (bold) and all identified levels
% !!! Note that good baseline estimation does not necessarily require
% identification of all levels, especially if they are very close to one
% another.
hold off;
plot(data);
hold on;
plot(result.ybaseline + result.ylevelOffsets, 'k--');
plot(result.ybaseline, 'k--', 'linewidth', 2);

%% We can also inspect the image representation that smBEVO used

% Plot smBEVO's image representation of the data and identified levels in
% image coordinates.
imfig = figure;
imagesc(result.im);
colormap(gray(256));
axis xy;
hold on;
ncols = size(result.im, 2);
plot((1:ncols)', (result.imbaseline + result.imlevelOffsets)', 'g-', 'linewidth', 1);






