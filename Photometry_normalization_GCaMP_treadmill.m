%% Fiber Photometry Epoch Averaging Example
%
% <html>
% This example goes through fiber photometry analysis using techniques <br>
% such as data smoothing, bleach detrending, and z-score analysis. <br>
% The epoch averaging was done using TDTfilter. <br><br>
% Author Contributions: <br>
% TDT, David Root, and the Morales Lab contributed to the writing and/or conceptualization of the code. <br>
% The signal processing pipeline was inspired by the workflow developed by <a href="https://doi.org/10.1016/j.celrep.2017.10.066">David Barker et al. (2017)</a> for the Morales Lab. <br>
% The data used in the example were provided by David Root. <br><br>


%% Housekeeping
% Clear workspace and close existing figures. Add SDK directories to Matlab
% path.
close all; clear all; clc;
[MAINEXAMPLEPATH,name,ext] = fileparts(cd); % \TDTMatlabSDK\TDTExamples\TDTExampleData
% DATAPATH = fullfile(MAINEXAMPLEPATH, 'TDTExampleData'); % \TDTMatlabSDK\TDTExamples\ExampleData
DATAPATH = fullfile(MAINEXAMPLEPATH, 'FED_k'); % \TDTMatlabSDK\TDTExamples\ExampleData
[SDKPATH,name,ext] = fileparts(MAINEXAMPLEPATH); % \TDTMatlabSDK
addpath(genpath(SDKPATH));

%% Importing the Data


BLOCKPATH = fullfile(DATAPATH,'R825-SNR-241127-105706');

% Now read the specified data from our block into a Matlab structure.
data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs', 'scalars', 'streams'});

%% check raw demodulated signals, red DA3 with 560 channel
% also plot with TTL pulses timing, for simple direct plotting
start1 = data.streams.x560D.startTime; % session start time
tim1=start1:1/(data.streams.x560D.fs):(start1 + (length(data.streams.x560D.data)-1)/(data.streams.x560D.fs)); % time for streamed data
r560 = double(data.streams.x560D.data); % demodulated 465 raw data;
r405r = double(data.streams.x405C.data); % demodulated 405 raw data;
figure
hold on
plot(tim1,r560) % plot blue signals, gCamp
hold on;
plot(tim1,r405r) % plot uv signals

%% check raw demodulated signals, GCaMP6 with 465 channel
% also plot with TTL pulses timing, for simple direct plotting
start1 = data.streams.x465C.startTime; % session start time
tim1=start1:1/(data.streams.x465C.fs):(start1 + (length(data.streams.x465C.data)-1)/(data.streams.x465C.fs)); % time for streamed data
r465 = double(data.streams.x465C.data); % demodulated 465 raw data;
r405 = double(data.streams.x405C.data); % demodulated 405 raw data;
figure
plot(tim1,r465) % plot blue signals, gCamp
hold on;
plot(tim1,r405) % plot uv signals

%% for 465 channel signal, calculating dF/F using raw 10th percenptile with 5s sliding window, -1.28 based on z score table for 10th percentile
% z-score table based on percentile can be found here: https://www.statology.org/z-table/
i_win = 5;
i_er = -1.28; % F0 percentile, if using 10th percentile, using value -1.28, based on z score table
kk = 1000; % typically exclude first 1000 samples, ~1sec, to get rid of large ringing artifacts when recording started
F0_1 = movmean(r465(kk:end),i_win*data.streams.x465C.fs) - 1.28*movstd(r465(kk:end),i_win*data.streams.x465C.fs);
dF_r465_1 = (r465(kk:end) - F0_1)./F0_1;% initial post-processing of dF/F
F0_2 = movmean(r405(kk:end),i_win*data.streams.x465C.fs) - 1.28*movstd(r405(kk:end),i_win*data.streams.x465C.fs);
dF_r405_1 = (r405(kk:end) - F0_2)./F0_2;
figure
kk1=kk;
plot(tim1(kk1:end),dF_r465_1) % plot blue signals, gCamp
hold on;
plot(tim1(kk1:end),dF_r405_1) % plot uv signals
[r,p] = corrcoef(dF_r465_1,dF_r405_1); % calculate cross correlation between uv and gcamp signals, only r<0.6 will be used

%% for 560 channel signal, calculating dF/F using raw 10th percenptile with 5s sliding window, -1.28 based on z score table for 10th percentile
% z-score table based on percentile can be found here: https://www.statology.org/z-table/
i_win = 5;
kk = 1000; % typically use 1000, use rid of 1st second for large ringing effects
r560 = r560(kk:end);r405r = r405r(kk:end);tim2 = tim2(kk:end);
F0_1 = movmean(r560,i_win*data.streams.x560D.fs) - 1.28*movstd(r560,i_win*data.streams.x560D.fs);
dF_r560_1 = (r560 - F0_1)./F0_1;
F0_2 = movmean(r405r,i_win*data.streams.x560D.fs) - 1.28*movstd(r405r,i_win*data.streams.x560D.fs);
dF_r405_1 = (r405r - F0_2)./F0_2;
figure
plot(tim2,dF_r560_1) % plot blue signals, gCamp
hold on;
plot(tim2,dF_r405_1) % plot uv signals
[r,p] = corrcoef(dF_r560_1,dF_r405_1); % calculate cross correlation between uv and gcamp signals, only r<0.6 will be used
%% a post-processing step of applying some smooth filter for the reference signal, 3rd order butterknive filter
F = 3; %corner frequency
Fs = data.streams.x405C.fs; %sampling rate
[y, x] = butter(3, F/(Fs/2)); % create how pass filter, 3rd order, if other filter time, define Ftype
inputSignal = dF_r405_1; % create input signal
dF_r405_1_f = filter(y, x, inputSignal); % created filtered signals
% [Passing the input signal as an input to the butterworth filter created]
figure
plot(dF_r405_1_f) % plot filtered signals
%% Signal normalization, Method#0, fit 465 with 405 data, using RANSAC least square fit
sampleSize = 2; % number of points to sample per trial
maxDistance = 2; % max allowable distance for inliers
tim2 = tim1(kk1:end);
points = []; 
points(:,1) = dF_r405_1_f; points(:,2) = dF_r465_1;
% points(:,1) = r405; points(:,2) = r465;
% points(:,1) = r405(1:1172*data.streams.x405C.fs); points(:,2) = r465(1:1172*data.streams.x405C.fs);tim3=tim2(1:1172*data.streams.x405C.fs);

fitLineFcn = @(points) polyfit(points(:,1),points(:,2),1); % fit function using polyfit
evalLineFcn = ...   % distance evaluation function
  @(model, points) sum((points(:, 2) - polyval(model, points(:,1))).^2,2);

[modelRANSAC, inlierIdx] = ransac(points,fitLineFcn,evalLineFcn, ...
  sampleSize,maxDistance);

% refit the model
modelInliers = polyfit(points(inlierIdx,1),points(inlierIdx,2),1);

inlierPts = points(inlierIdx,:);
x = [min(inlierPts(:,1)) max(inlierPts(:,1))];
y = modelInliers(1)*points(:,1) + modelInliers(2);% UV transformed 560 signal
figure
plot(tim2, y, 'g-')
% legend('Noisy points','Least squares fit','Robust fit');
hold off



%% final post-processed dF/F signal, red channels
dF_560_f = dF_r560_1-y'; %normalized 560 signal
dF_405_f = y'; % normalized uv signal


%% final post-processed dF/F signal, green channels
dF_465_f = dF_r465_1-y'; %normalized 465 signal
dF_405_f = y'; %normalized uv signal


%% find locomotion, acceleration onset or offset based on video tracking, if use DLC output, use resampled locomotion data
% velocity is always positive, all behavioral data (i.e.velocity) stored in the array called spdR1 

spdR1 = [tim2;data.streams.Wav1.data(kk1:end)];% first column timestamp, 2nd column analog velocity output from treadmill, kk1=1000;
spdR1(:,2) = spdR1(:,2)*13; % multiply scaling factor according to manufactuer instruction to convert unit to cm/s
spdR1(:,3) = [nan;diff(spdR1(:,2))./diff(spdR1(:,1))];% acceleration
spdR1 = double(spdR1); % somehow the original data format was singles  
