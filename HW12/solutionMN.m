clear; close all; clc; 
solution()
function solution()
% function solution()
% 
% BMI 500 Solution
% Clinical Motion Analysis 2023
% 
% Name: Nateghi, Masoud

% Declare some global plotting variables, take small joy from Matlab
% complaining. (You should really use an input parser structure or a class
% structure for this, but globals are ok for initial development)
global COL_MAX T_MAX F_MAX Y_MAX FC_LO FC_HI F_FEAT

% Maximum number of variables to plot
COL_MAX = 9;

% Maximum time to plot
T_MAX = 15;

% Maximum frequency to plot
F_MAX = 15;

% Maximum Y limit
Y_MAX = 5;

% Hipass filter cutoff
FC_HI = 2;

% Lowpass filter cutoff
FC_LO = 20;

% Feature vector frequencies
N_FEAT = 50;
F_FEAT = linspace(1,12,N_FEAT);

% Load the labels
labels = readtable('tremor-present-absent/train-labels.csv','TextType','string');

% Prepend a path to the labels
% labels.file = "tremor-present-absent/" + labels.file;

% Plot three cases of tremor absent, mild, and severe
inds = [1442 1718 1723];
for ind = inds
    plot_record(filter_butterworth(readtable(labels.file(ind))));    
end

% Compare the severe tremor case before and after Butterworth filtering
plot_record(filter_fft(readtable(labels.file(ind))));    

% Compare the feature vectors of the three records
for ind = inds
    plot_feature_vec(filter_butterworth(readtable(labels.file(ind))));
end

% Select only the L_Hand records

% Loop over all of the samples and calculate the feature vec
X = nan(size(labels,1),N_FEAT);
for obs = labels.id'
    X(obs,:) = get_feature_vec(filter_butterworth(readtable(labels.file(obs))));
end

% Isolate the labels of tremor present vs. absent
y = labels.tremor;

% Plot the feature vectors vs. the labels
plot_features_vs_labels(X,y);
a = 1;
% Fit an SVM with defaults
SVMModel1 = fitcsvm(X,y);

% or

% Fit an SVM with optimization of slack variable etc.
rng default
SVMModel2 = fitcsvm(X,y,'OptimizeHyperparameters','auto', ...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName', ...
    'expected-improvement-plus'));

label1 = predict(SVMModel1,X);
label2 = predict(SVMModel2,X);

F1(y,label1)
% ans = 0.7679

F1(y,label2)
% ans = 0.8640

end

function str_out = names(tbl_in)
% function str_out = names(tbl_in)
% 
% Helper function to extract variable names

str_out = string(tbl_in.Properties.VariableNames);

end

function out = F1(truth,labels)
% function out = F1(truth,labels)
% 
% Design a function to calculate F1 score
% 
% 1 POINT

%%%%%%%%%%%%%%%%
TP = sum(truth == 1 & labels == 1); 
FP = sum(truth == 0 & labels == 1); 
FN = sum(truth == 1 & labels == 0); 
out = TP / (TP + 0.5 * (FP + FN)); 
%%%%%%%%%%%%%%%%

end

function fig = plot_record(data)
% function fig = plot_record(data)
% 
% Design a function that plots the timeseries data and FFT data of an individual record
% 
% 2 POINTS

global COL_MAX T_MAX F_MAX

% Extract the data and variable names
t = data.Time;
xyz = data{:,3:end};
nms = names(data);

% Discard frame number and time vector
nms = nms(3:end);

%%%%%%%%%%%%%%%%
fs = 1 / (t(2) - t(1)); 
x = data{:, 3:3+COL_MAX-1};
N = length(x); 
X = fft(x); 
f = 0:fs/N:fs/2;
X_abs = abs(X) / N; 
figure; 
subplot(121)
stackedplot(t(t <= T_MAX), x(t <= T_MAX, :))
xlabel('Seconds')
subplot(122)
stackedplot(f(f <= F_MAX), X_abs(f <= F_MAX, :))
xlabel('Hz')
%%%%%%%%%%%%%%%%

end

function vec = get_feature_vec(data)
% function vec = get_feature_vec(data)
%
% Design a function that interpolates the single-sided FFT to from
% forwardfft to the specific frequency components at 2 Hz, 3 Hz, etc. up to
% 20 Hz. Average across all x, y, z indices of all markers. You may use
% interp1; this allows for "nearest" interpoloation, "linear," and other
% methods. Use "nearest" for this purpose.
% 
% Make sure that the data are preprocessed prior to this function.
%  
% 2 POINTS

global F_FEAT

%%%%%%%%%%%%%%%%
[P1,f] = forwardfft(data); 
P1_avg = mean(P1, 2); 
vec = interp1(f, P1_avg, F_FEAT, "nearest"); 
%%%%%%%%%%%%%%%%

end

function fig = plot_feature_vec(data)
% function fig = plot_feature_vec(data)
% 
% Design a function to plot the feature vector of a given record.
%  
% 1 POINT

global F_MAX F_FEAT

%%%%%%%%%%%%%%%%
[P1,f] = forwardfft(data); 
vec = get_feature_vec(data); 
n = size(P1, 2); 
figure;
for i = 1:n
    plot(f(f <= F_MAX), P1(f <= F_MAX, i))
    hold on
end
scatter(F_FEAT, vec, 'black', 'filled')
%%%%%%%%%%%%%%%%

end

function fig = plot_features_vs_labels(X,y)
% function fig = plot_features_vs_labels(X,y)
% 
% Design a function that plots the feature vectors vs. the labels
%  
% 1 POINT

global F_FEAT F_MAX Y_MAX

%%%%%%%%%%%%%%%%
X_0 = X(y == 0, :); 
X_1 = X(y == 1, :); 
figure; 
subplot(211); 
for i = 1:size(X_0, 1)
    scatter(F_FEAT, X_0(i, :), 'MarkerEdgeColor', [0.5, 0.5, 0.5], 'MarkerFaceColor', [0.5, 0.5, 0.5])
    hold on
end
scatter(F_FEAT, mean(X_0), 'black', 'filled')
xlim([0, F_MAX])
ylim([0 Y_MAX])
title('Y == 0')
ylabel('$$mm^{2}_{RMSE}$$', 'Interpreter', 'latex')


subplot(212); 
for i = 1:size(X_1, 1)
    scatter(F_FEAT, X_1(i, :), 'MarkerEdgeColor', [0.5, 0.5, 0.5], 'MarkerFaceColor', [0.5, 0.5, 0.5])
    hold on
end
scatter(F_FEAT, mean(X_1), 'black', 'filled')
xlim([0, F_MAX])
ylim([0, Y_MAX])
title('Y == 1')
ylabel('$$mm^{2}_{RMSE}$$', 'Interpreter', 'latex')
%%%%%%%%%%%%%%%%

end

function data_out = filter_butterworth(data)
% function data_out = filter_butterworth(data)
% 
% Design a function that preprocesses data with highpass and lowpass
% 6th-order butterworth filters at 2 and 20 Hz, respectively
%  
% 2 POINTS

global FC_LO FC_HI

%%%%%%%%%%%%%%%%
fs = 1 / (data.Time(2) - data.Time(1)); 
[b,a] = butter(6, [FC_HI, FC_LO] / (125 / 2), 'bandpass');
data_out = data; 
for i = 3:length(names(data))
    data_out{:, i} = filtfilt(b, a, data{:, i});
end
%%%%%%%%%%%%%%%%

end

function data_out = filter_fft(data)
% function data_out = filter_fft(data)
% 
% Design a function that preprocesses data by setting components below 2
% and higher than 20 Hz to zero in the frequency domain
% note that this uses the syntax from https://www.mathworks.com/help/matlab/ref/fft.html
%  
% 2 POINTS

global FC_LO FC_HI

%%%%%%%%%%%%%%%%
t = data.Time;
fs = 1 / (t(2) - t(1)); 
x = data{:, 3:end};
N = length(x); 
X = fft(x); 
f = 0:fs/N:fs-fs/N;
X((f >= 0 & f <= FC_HI) | (f >= FC_LO & f <= fs-FC_LO) | (f >= fs-FC_HI & f <= fs), :) = 0;
x_filt = ifft(X); 
data_out = data; 
data_out{:, 3:end} = x_filt;
%%%%%%%%%%%%%%%%

end

function [P1,f] = forwardfft(data)
% function [P1,f] = forwardfft(data)
% 
% Design a function that implements the single-sided FFT based on the
% Matlab FFT function. 
% note that this uses the syntax from https://www.mathworks.com/help/matlab/ref/fft.html
%  
% 1 POINT

%%%%%%%%%%%%%%%%
t = data.Time; 
fs = 1/(t(2)-t(1));
x = data{:, 3:end}; 
N = length(x); 
f = 0:fs/N:fs/2;
X = fft(x); 
P2 = abs(X) / N; 
P1 = P2(1:N/2+1, :); 
P1(2:end-1, :) = 2*P1(2:end-1, :); 
%%%%%%%%%%%%%%%%

end

