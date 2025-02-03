%% Fiber Photometry Epoch Averaging Example
%
% <html>
% This example goes through fiber photometry analysis using techniques <br>
% such as data smoothing, bleach detrending, and z-score analysis. <br>
% The epoch averaging was done using TDTfilter. <br><br>
% Author Contributions: <br>
% TDT, David Root, and the Morales Lab contributed to the writing and/or conceptualization of the SDK code. <br>
% The signal processing pipeline was inspired by the workflow developed by <a href="https://doi.org/10.1016/j.celrep.2017.10.066">David Barker et al. (2017)</a> for the Morales Lab. <br>


%% Housekeeping
% Clear workspace and close existing figures. Add SDK directories to Matlab
% path.
close all; clear all; clc;
[MAINEXAMPLEPATH,name,ext] = fileparts(cd); % \TDTMatlabSDK\TDTExamples\T≥DTExampleData
% DATAPATH = fullfile(MAINEXAMPLEPATH, 'TDTExampleData'); % \TDTMatlabSDK\,TDTExamples\ExampleData
DATAPATH = fullfile(MAINEXAMPLEPATH, 'GABABKO_CHR2_RDA3M'); % \TDTMatlabSDK\TDTExamples\ExampleData
[SDKPATH,name,ext] = fileparts(MAINEXAMPLEPATH); % \TDTMatlabSDK
addpath(genpath(SDKPATH));

% Importing the Data
% This example assumes you downloaded our example data sets
% (<https://www.tdt.com/files/examples/TDTExampleData.zip link>) and extracted
% it into the \TDTMatlabSDK\Examples\ directory. To import your own data, replace
% |BLOCKPATH| with the path to your own data block.
%
% In Synapse, you can find the block path in the database. Go to Menu > History. 
% Find your block, then Right-Click > Copy path to clipboard.

BLOCKPATH = fullfile(DATAPATH,'H156-L-241109-114236');

% Now read the specified data from our block into a Matlab structure.
data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs', 'scalars', 'streams'});


%% check raw demodulated signals, red DA3
% also plot with TTL pulses timing, for simple direct plotting
start1 = data.streams.x560D.startTime; % session start time
tim2=start1:1/(data.streams.x560D.fs):(start1 + (length(data.streams.x560D.data)-1)/(data.streams.x560D.fs)); % time for streamed data
r560 = double(data.streams.x560D.data); % demodulated 465 raw data, sometimes the data format may appear as singles;
r405r = double(data.streams.x405D.data); % demodulated 405 raw data;
figure
hold on
plot(tim2,r560) % plot blue signals, gCamP
hold on;
plot(tim2,r405r) % plot uv signals.


%% check raw demodulated signals, if it is GCaMP6, essentially using r465 for places of r560
% also plot with TTL pulses timing, for simple direct plotting
start1 = data.streams.x465C.startTime; % session start time
tim1=start1:1/(data.streams.x465C.fs):(start1 + (length(data.streams.x465C.data)-1)/(data.streams.x465C.fs)); % time for streamed data
r465 = double(data.streams.x465C.data); % demodulated 465 raw data, sometimes the data format may appear as singles;
r405 = double(data.streams.x405C.data); % demodulated 405 raw data;
figure
plot(tim1,r465) % plot blue signals, gCamp
hold on;
plot(tim1,r405) % plot uv signals

%% first step doing some smooth filter for the reference signal, 3rd order butterknive filter
F = 3; %corner frequency
Fs = data.streams.x405D.fs; %sampling rate
[y, x] = butter(3, F/(Fs/2)); % create how pass filter, 3rd order, if other filter time, define Ftype
inputSignal = r405; % create input signal
r405 = filter(y, x, inputSignal); % created filtered signals
% [Passing the input signal as an input to the butterworth filter created]
figure
plot(r405) % plot filtered signals to check

%% first step of post processing, calculating dF/F using raw 10th percenptile with 5s sliding window, -1.28 based on z score table for 10th percentile
% z-score table based on percentile can be found here: https://www.statology.org/z-table/
i_win = 5; % use 5 second window for moving average computation
i_er = -1.28; % F0 percentile, if using 10th percentile, using value -1.28, based on z score table
kk = 1000; % typically exclude first 1000 samples, to get rid of 1st second large ringing artifacts when turning on the recording
F0_1 = movmean(r560(kk:end),i_win*data.streams.x560D.fs) - 1.28*movstd(r560(kk:end),i_win*data.streams.x560D.fs); %smoothed 560 signal, with 10th percentile adjustment
F0_2 = movmean(r405r(kk:end),i_win*data.streams.x560D.fs) - 1.28*movstd(r405r(kk:end),i_win*data.streams.x560D.fs);% smoothed UV signal

[r,p] = corrcoef(F0_1,F0_2); % calculate cross correlation between uv and gcamp signals, only r<0.6 will be used


%% %% Signal normalization, Method#0, fit 560 with 405 data, using RANSAC least square fit
sampleSize = 2; % number of points to sample per trial
maxDistance = 2; % max allowable distance for inliers
points = []; 
% points(:,1) = dF_r405_1_f; points(:,2) = dF_r560_1;

points(:,1) = F0_2; % the uv
points(:,2) = F0_1; % the 560 signal


fitLineFcn = @(points) polyfit(points(:,1),points(:,2),1); % fit function using polyfit
evalLineFcn = ...   % distance evaluation function
  @(model, points) sum((points(:, 2) - polyval(model, points(:,1))).^2,2);

[modelRANSAC, inlierIdx] = ransac(points,fitLineFcn,evalLineFcn, ...
  sampleSize,maxDistance);

% refit the model
modelInliers = polyfit(points(inlierIdx,1),points(inlierIdx,2),1);

inlierPts = points(inlierIdx,:);
x = [min(inlierPts(:,1)) max(inlierPts(:,1))];
y = modelInliers(1)*points(:,1) + modelInliers(2); %UV fitted signal,'F' in normalization 
figure
plot(tim1, y, 'g-')
% legend('Noisy points','Least squares fit','Robust fit');
hold off


%% compute normalized dF/F for the 560 signal after post-processing
figure;
dF_560_f = (F0_1-y')./y'; % normalized 560 signal, dF/F
plot(tim2(kk1:end), dF_560_f, 'g-')
set(gcf,'color','w')
%% parsing into optic sti

i_dur =15; % duration of stimulation, change that in each session if duration is different
TDT_TrialHead = data.epocs.PC0_.onset; % rearing start 
TDT_TrialTail = data.epocs.PC0_.offset + i_dur; % duration of stimulation is 5 sec
DAF.r560 = dF_560_f-mean(dF_560_f); %subtract mean to offset individual variabilities among subjects
tim1 = tim2(kk:end);
% then parsing sorted spikes into trials, typically add 10s after stim off,if can be longer if needed 
    DAF.time_Trials =  arrayfun(@(x,y)...
                   tim1(tim1 >= x-10 & tim1 <= y+10),...
                   TDT_TrialHead, TDT_TrialTail,...
                   'UniformOutput',false)
    
    DAF.r560_Trials =  arrayfun(@(x,y)...
                   DAF.r560(tim1 >= x-10 & tim1 <= y+10),...
                   TDT_TrialHead, TDT_TrialTail,...
                   'UniformOutput',false);

    % plotting trial by trial data and generate cell arrrays for average
    figure
    trs1_1 = {};
for i = 1:size(DAF.time_Trials,1)
    trim1_t = DAF.time_Trials{i} - TDT_TrialHead(i); % aligned on trial start
    f560_t = DAF.r560_Trials{i};
    trs1_1{i} = f560_t(trim1_t >=-10 & trim1_t <=10+i_dur);% each individual trials with stimulation
    plot(trim1_t,f560_t)
    hold on
    xlabel('Time from trial onset (s)')
    ylabel('Normalized fluorescence signal (dF/F)')

end


%% generate average response across trials and figures with error bar (SEM)
figure
set(gcf,'color','w')
len0=min(cell2mat(cellfun(@(x) length(x), trs1_1, 'UniformOutput', false))); % find minmum length of a trial
tem0=trim1_t(trim1_t >=-10 & trim1_t <=10+i_dur);
array1=cell2mat(cellfun(@(x) x(1:len0)', trs1_1, 'UniformOutput', false));


% % calculating dF/F with % and erros
x00=tem0(1:len0);
y00=mean(array1,2)'*100;
sd00  = std(array1,0,2)'*100; n00 = size(array1,2); er00  = sd00/(n00^0.5);
errorSchmear(x00,y00(1:end),er00(1:end),'color',[0.4 0.4 0.4], 'LineWidth', 1.5)

xlabel('Time from stimulation onset (sec)')
ylabel('F560 Mean dF/F (%)')
title('3mW 15s SNr stimulation') % change stim parameters in figure making
box off
set(gca,'tickdir','out')



