%% find locomotion, acceleration onset or offset based on video tracking, if use DLC output, use resampled locomotion data
% velocity is always positive, all behavioral data (i.e.velocity) stored in the array called spdR1 

spdR1 = [tim2;data.streams.Wav1.data(kk1:end)];% first column timestamp, 2nd column analog velocity output from treadmill, kk1=1000;
spdR1(:,2) = spdR1(:,2)*13; % multiply scaling factor according to manufactuer instruction to convert unit to cm/s
spdR1(:,3) = [nan;diff(spdR1(:,2))./diff(spdR1(:,1))];% acceleration
spdR1 = double(spdR1); % somehow the original data format was singles  

sFs = data.streams.x465C.fs; % photometry data sampling rate

% with treadmill, same rample rate in adc and photometry
sigv = spdR1(2:end,2);% velocity
siga = spdR1(2:end,3);% acceleration
Ty = tim2(2:end)';% timestamps
r465 = dF_465_f(2:end); % use normalized GCaMP signals

x1=zeros(size(sigv,1),1); % array for movements
x1(sigv>=0.2)=1;% using 0.2cm/s as movement threshold for head-fixed treadmill, can adjust
x1_i0 = diff([0;x1]);
x1_i = find(x1_i0==1); % find index of movement bouts
x1_j = find(x1_i0==-1); % find index of movement bouts end

if length(x1_i)>length(x1_j) % in case the last data point has positive velocity
    x1_i = x1_i(1:end-1);
else
    x1_j = x1_j;
end

Tr1=Ty(x1_i((x1_j-x1_i)>=round(sFs/2.0)),1); % timestamp for movement start, use 500ms minimum locomotion duration
Tr2=Ty(x1_j((x1_j-x1_i)>=round(sFs/2.0)),1); % timestamp for movement end, use 500ms minimum locomotion duration

x1_i0=x1_i(2:end);
Tr3=Ty(x1_i0((x1_i(2:end)-x1_j(1:end-1))>=round(sFs/1.0)),1); % movement bout start, using 1 sec no movement period to separate movement bouts, as some movements were jittering
Tr4=Ty(x1_j((x1_i(2:end)-x1_j(1:end-1))>=round(sFs/1.0)),1); % movement bout ends
Tr4=[Tr4(2:end);Ty(x1_j(end))];
% Parsing into trial cell arrays
% find trial start and end

TDT_TrialHead = Tr3; % event start
TDT_TrialTail = Tr4; % event end



 

% 
    CAF.time_Trials =  arrayfun(@(x,y)...
                   Ty(Ty >= x & Ty <= y),...
                   TDT_TrialHead, TDT_TrialTail,...
                   'UniformOutput',false);
   
    CAF.r465_Trials =  arrayfun(@(x,y)...
                   r465(Ty >= x & Ty <= y),...
                   TDT_TrialHead, TDT_TrialTail,...
                   'UniformOutput',false);
% velocity
    CAF.vel_Trials =  arrayfun(@(x,y)...
                   sigv(Ty >= x & Ty <= y),...
                   TDT_TrialHead, TDT_TrialTail,...
                   'UniformOutput',false);




%% parsing into trial segments around event onset for individual dataset
% then parsing sorted signals into trials, for correlation analysis

% % 
r465 = (r465-mean(r465))./std(r465); % z-score the array
r465 = r465'; % transpose the array

r405 = (r405-mean(r405))./std(r405); % z-score the array
r405 = r405'; % transpose the array

    XAF.time_Trials =  arrayfun(@(x,y)...
                   Ty(Ty >= x-2.0 & Ty <= y+2.0),...
                   TDT_TrialHead, TDT_TrialTail,...
                   'UniformOutput',false);

% photometry signals   
    XAF.r465_Trials =  arrayfun(@(x,y)...
                   r465(Ty >= x-2.0 & Ty <= y+2.0),...
                   TDT_TrialHead, TDT_TrialTail,...
                   'UniformOutput',false);
    XAF.r405_Trials =  arrayfun(@(x,y)...
                   r405(Ty >= x-2.0 & Ty <= y+2.0),...
                   TDT_TrialHead, TDT_TrialTail,...
                   'UniformOutput',false);
% velocity
    XAF.vel_Trials =  arrayfun(@(x,y)...
                   sigv(Ty >= x-2.0 & Ty <= y+2.0),...
                   TDT_TrialHead, TDT_TrialTail,...
                   'UniformOutput',false);
% acceleration
    XAF.acc_Trials =  arrayfun(@(x,y)...
                   siga(Ty >= x-2 & Ty <= y+2),...
                   TDT_TrialHead, TDT_TrialTail,...
                   'UniformOutput',false);



    trs3_1 = {};trs3_2 = {};trs3_3 = {};trs3_0 = {};
    kkk = 0; % counter
for i = 3:size(XAF.time_Trials,1)
    if mean(XAF.vel_Trials{i})>0.20 % only includes locomotion bouts with average velocity above threshold
        kkk = kkk + 1;
    trim1_t = XAF.time_Trials{i} - TDT_TrialHead(i); % aligned on locomotion onset, use TrialHead, if aligned with locomotion offset, use TrialTail
    f465_t = XAF.r465_Trials{i}; % photometry signals
    trs3_1{kkk} = f465_t(trim1_t >=-2.0 & trim1_t <=2.0);
    f465_t0 = XAF.r405_Trials{i}; % photometry signals
    trs3_0{kkk} = f465_t0(trim1_t >=-2 & trim1_t <=2);
    f465_t1 = double(XAF.vel_Trials{i}); % velocity
    trs3_2{kkk} = f465_t1(trim1_t >=-2 & trim1_t <=2);
    end

end

a_cog_0 = double(cell2mat(trs3_1)'); % gCaMP or other main fluorophor
a_cou_0 = double(cell2mat(trs3_0)'); % uv
a_cov_0 = double(cell2mat(trs3_2)'); % velocity
tem0=trim1_t(trim1_t >=-2 & trim1_t <=2); % time

%% make plot of the mean trace with errorbars
figure;

x00 = tem0; %timestamps
y00=nanmean(a_cog_0,1); %mean across trials
sd00  = nanstd(a_cog_0,0,1); n00 = size(a_cog_0,1); er00  = sd00/(n00^0.5);% calculating errorbars with SEM

%  %for locomotion offset data
% y01=nanmean(af_cog_0,1);
% sd01  = nanstd(af_cog_0,0,1); n00 = size(af_cog_0,1); er01  = sd01/(n00^0.5);

hold on
errorSchmear(x00,y00(1:end),er00(1:end),'color',[0.4 0.4 0.8], 'LineWidth', 2)
xlim([-1 1]) % unit of x-axis here is samples
% 
% hold on
% errorSchmear(x00,y01(1:end),er01(1:end),'color',[0.4 0.4 0.4], 'LineWidth', 0.75)
% xlim([-1 1]) % unit of x-axis here is samples


set(gcf,'color','w');set(gca,'tickdir','out');box off

%% create variable for generate heatmap of the gCamp activity for each mouse
% create cell for gcamp signal aligned on locomotion onset, later rename as
% a_SPC
SPC=[mean(a_cog_1,1);mean(a_cog_2,1);mean(a_cog_3,1);mean(a_cog_4,1);mean(a_cog_5,1);mean(a_cog_6,1);mean(a_cog_7,1);mean(a_cog_8,1);mean(a_cog_9,1);mean(a_cog_10,1);mean(a_cog_11,1);mean(a_cog_12,1)];
% for data aligned on locomotion offset, use af_cog_x instead and rename as af_SPC;

figure;
imagesc(sortrows(bsxfun(@rdivide,SPC1,max(abs(SPC),[],2)),1018:1:2035,'ascend')) %normalize the data with maximum value, and sorting data based on values within the window   
set(gcf,'color','w');set(gca,'tickdir','out');box off

%% get the rolling window stats, to calculate activity rise time, maximum activity within a window

GR={}; %create a cell variable for loading data from individual animals 
GR{1}=a_cog_1;GR{2}=a_cog_2;GR{3}=a_cog_3;GR{4}=a_cog_4;GR{5}=a_cog_5;GR{6}=a_cog_6;GR{7}=a_cog_7;GR{8}=a_cog_8;GR{9}=a_cog_9;GR{10}=a_cog_10;GR{11}=a_cog_11;

TM1 = NaN(size(GR,2),1);
TM_o=zeros(size(GR,2),1);

kk = data.streams.x465C.fs; %TDT system sampling rate, or hard coded as 1017

for i=1:size(GR,2)

a_cog = GR{i};

% get rolling window analysis for onset time

t1 = a_cog(:,kk+1:3*kk);% get trunk of 2s data

t1_1= mean(t1(:,1:round(kk/2)),2);% baseline 500ms across trials, 509 samples
z=NaN(501,1);

w1 = 20; %20ms sliding window
q=0;
for j = 0:2:1000
    q=q+1;
  t2 = mean(t1(:,510+w1+j:520+w1+j),2);
 [h,p]=ttest(t2-t1_1,0,'tail','right'); %do a one-tailed t-test to determine if activity exceed baseline
  z(q)=h;
end

x0=find(z==1,1,'first');% find index of activity increase first
ot0 = (- 500 + (x0-1)*2 + w1/2)/1000; % activity increase onset time in ms

TM1(i) = ot0; %tally it from individual animals

% find timing of max activity

AAt = tem0(kk+1:3*kk); %time vector within the range
AA1 = mean(t1,1);

TM_o(i)=AAt(AA1==max(AA1)); %finding the timestamp for maximum activity


end


%% get the peak and slope, mostly for offset

% TM1 = NaN(size(GR,2),1);
TM_f=zeros(size(GF,2),1);
TM_s=zeros(size(GF,2),1);% slope

%sampling rate is 1017, hard-code numbers here

for i=1:size(GF,2)

af_cog = GF{i};

% get rolling window analysis for onset time

t1 = af_cog(:,1018:3051);% get trunk of 2s data


% get the slope

AAt2 = tem0(1527:2034); %time vector within the range (500ms interval)
AA2 = mean(af_cog(:,1527:2034),1);


coef1=polyfit(AAt2', AA2, 1);
TM_s(i)=coef1(1);% coefficient 1 is slope


% find timing of max activity

AAt = tem0(1018:3051); %time vector within the range
AA1 = mean(t1,1);


TM_f(i)=AAt(AA1==max(AA1));


end



%% get the mean movement velocity, mostly for offset and onset

GF1={}; %create a cell variable for loading data from individual animals, GF is for offset, GR is for onset 
GF1{1}=af_cov_1;GF1{2}=af_cov_2;GF1{3}=af_cov_3;GF1{4}=af_cov_4;GF1{5}=af_cov_5;GF1{6}=af_cov_6;GF1{7}=af_cov_7;GF1{8}=af_cov_8;GF1{9}=af_cov_9;GF1{10}=af_cov_10;GF1{11}=af_cov_11;


VM_f=zeros(size(GF1,2),1);%offset
VM_o=zeros(size(GR1,2),1);% onset

kk = data.streams.x465C.fs; %TDT system sampling rate, or hard coded as 1017

for i=1:size(GF1,2)

af_cov = GF1{i};
ao_cov = GR1{i};

% get rolling window analysis for onset time

t0 = mean(mean(ao_cov(:,2*kk+1:3*kk)));% get trunk of 1s data
t1 = mean(mean(af_cov(:,kk+1:2*kk)));% get trunk of 1s data


VM_o(i)=t0;%onset velocity
VM_f(i)=t1;%offset velocity


end





