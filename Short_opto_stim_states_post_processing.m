%%
%% Parsing into individual trial, cell arrays, only use the normalized data
% find trial start and end, read from Excel files of post video analysis
% from EthoVision
% Kremen1-ChR2 stim data, later renamed as pS
S.sub{1}=X1_R657;S.sub{2}=X1_R658;S.sub{3}=X1_R659;S.sub{4}=X1_R805;
S.sub{5}=X1_R807;S.sub{6}=X1_R816;S.sub{7}=X1_R883;S.sub{8}=X1_R944;S.sub{9}=X1_R945;

% Kremen1-YFP stim data, later renamed as pG
S1.sub{1}=X2_R543;S1.sub{2}=X2_R646;S1.sub{3}=X2_R660;S1.sub{4}=X2_R661;
S1.sub{5}=X2_R662;S1.sub{6}=X2_R529;S1.sub{7}=X2_R542;S1.sub{8}=X2_R586;S1.sub{9}=X2_R588;

%%
pS.sub{1}=cell2mat(Da(:,2));pS.sub{2}=cell2mat(Da(:,3));pS.sub{3}=cell2mat(Da(:,4));% Kremen1-ChR2 2nd cohort
cS.sub{1}=cell2mat(Da(:,5));cS.sub{2}=cell2mat(Da(:,6));cS.sub{3}=cell2mat(Da(:,7));cS.sub{4}=cell2mat(Da(:,8));% Calb1-ChR2 2nd cohort
tim_1 = cell2mat(Da(:,1));% time stamps
%%

%% Separate trials based on locomotion state prior or during optogenetic stimulation

% critical variables: data structure for locomotion velocity: pS
% (kremen1-ChR2-stim);pG (kremen1-YFP-stim); cS(Calb1-ChR2-stim);
% cG(Calb1-YFP-stim); within the structure, each cell is data from a single
% subject
% Timestamps for stimulation onset time,doubles array, each row is a trial,
% each column is a subject. TP(kremen1-ChR2-stim); TPg(kremen1-YFP-stim);
% TC(calb1-ChR2-stim);TCg(calb1-YFP-stim);
 
  p0 = 0; p1 = 0; p2 = 0; p3 = 0;
     trs1_1 = {};trs1_2 = {};trs2_1 = {};trs2_2 = {};pp={};
for i = 1:size(pS.sub,2)

x1=zeros(size(pS.sub{i},1),1); % array for movements
x2=zeros(size(pS.sub{i},1),1); % array for stationary
x1(pS.sub{i}(:,1)>=2)=1;% using 2.0cm/s as movement threshold
x2(pS.sub{i}(:,1)<0.5)=1;% using 0.5cm/s as stationary threshold

TDT_TrialHead = cell2mat(TP(:,i)); % trial starts 
TDT_TrialTail = TDT_TrialHead + 10; % trial ends

tim1 = pS.sub{i}(:,2); % time
da1 = pS.sub{i}(:,1); % velocity
% then parsing sorted spikes into trials, taking 2min before and 1min after 

    DAF.time_Trials =  arrayfun(@(x,y)...
                   tim1(tim1 >= x-5 & tim1 <= y+5),...
                   TDT_TrialHead, TDT_TrialTail,...
                   'UniformOutput',false);
    
    DAF.ve_Trials =  arrayfun(@(x,y)...
                   da1(tim1 >= x-5 & tim1 <= y+5),...
                   TDT_TrialHead, TDT_TrialTail,...
                   'UniformOutput',false);

    DAF.mv_Trials =  arrayfun(@(x,y)...
                   x1(tim1 >= x-5 & tim1 <= y+5),...
                   TDT_TrialHead, TDT_TrialTail,...
                   'UniformOutput',false);       % movement bouts
    DAF.im_Trials =  arrayfun(@(x,y)...
                   x2(tim1 >= x-5 & tim1 <= y+5),...
                   TDT_TrialHead, TDT_TrialTail,...
                   'UniformOutput',false);       % immobility

   
      p4 = 0;pp1=[];
for j = 1:size(DAF.time_Trials,1)
    trim1_t = DAF.time_Trials{j} - TDT_TrialHead(j); % aligned on trial start
    f0_t = DAF.ve_Trials{j};
    f0_t1 = DAF.mv_Trials{j};
    f0_t2 = DAF.im_Trials{j};
    if sum(f0_t1(trim1_t >=-1 & trim1_t <=0))>=15 % movement before stim start
        p0 = p0 + 1;
        trs1_1{i}(:,p0) = f0_t(trim1_t >=-2 & trim1_t <=5);
    elseif sum(f0_t2(trim1_t >=-1 & trim1_t <=0))>=15 % immob before stim start
        p1 = p1 + 1;
        trs1_2{i}(:,p1) = f0_t(trim1_t >=-2 & trim1_t <=5);
    end

    trim2_t = DAF.time_Trials{j} - TDT_TrialTail(j); % aligned on trial end

    if sum(f0_t1(trim2_t >=-1 & trim2_t <=0))>=15 % movement bef stim end
        p2 = p2 + 1; p4 = p4+1;
        trs2_1{i}(:,p2) = f0_t(trim2_t >=-2 & trim2_t <=5);
    elseif sum(f0_t2(trim2_t >=-0.6 & trim2_t <=0))>=8 % immob before stim end
        p3 = p3 + 1;
        trs2_2{i}(:,p3) = f0_t(trim2_t >=-2 & trim2_t <=5);
    end

pp1=[pp1;p4];

end


pp{i}=pp1;

end

% later added a prefix to these trs cell arrays, p-
% (kremen1-ChR2);pg- (kremen1-YFP);c- (calb1-ChR2);cg- (calb1-YFP);


%% making plots with errorbars 
figure
hold on
len0=min(cell2mat(cellfun(@(x) length(x), trs1_2(1:end), 'UniformOutput', false))); % find minmum length of a trial
tem0=trim1_t(trim1_t >=-2 & trim1_t <=5);
array1=cell2mat(cellfun(@(x) x(1:len0), trs1_2(1:end), 'UniformOutput', false));
x00=tem0(1:len0);
y00=mean(array1,2)'*1;
sd00  = std(array1,0,2)'*1; n00 = size(array1,2); er00  = sd00/(n00^0.5);
errorSchmear(x00,y00(1:end),er00(1:end),'color',[0.4 0.4 0.4], 'LineWidth', 2)
% plot(downsample(x00,10),downsample(y00,10));

xlabel('Time from stimulation onset (sec)')
ylabel('Velocity cm/s')
% title('Immobility prior to stimulation')
% ylim([-4.5 4.5])
set(gcf,'color','w')
set(gca,'tickdir','out');
box off



%% calculating mean velocity across trials within the defined interval, pre or post stimulation
% the generated data were used in Prism for plotting 

tem1 = tem0(tem0>=-2&tem0<=5);
am_1 = NaN(size(p_trs1_1,2),1);
am_2 = NaN(size(p_trs1_1,2),1);

for i=1:size(p_trs1_1,2)
 
    a1=mean(p_trs1_2{1,i},2);
    am_1(i) = mean(a1(tem1>=-1.0 & tem1<=-0.0));
    am_2(i) = mean(a1(tem1>0.0 & tem1<=1.0));

end

