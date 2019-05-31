%% script to calculate delay between haptic touch and DCS onset
% Compares between DCS onset and haptic touch for a desired condition
% following the stimulation response timing experiment
% David.J.Caldwell 7-19-2018

%% load in subject
clear all; clc

% ui box for input
% have to have a value assigned to the file to have it wait to finish
% loading...mathworks bug

[filename,filepath] = uigetfile('.mat');
load([filepath filename],'ECO1','Sing','Stim','Tact','Tra1','Tra2','Valu')

prompt = {'Which block on TDT front panel was this? e.g 1',...
    'Plot intermediate voltages? 0 for no, 1 for yes',...
    'input the condition of interest e.g. 2,3,4,5,6 (-1,0,1 are haptic, null, off target)',...
    'new tactor ? 0 for no, 1 for yes'};
dlgTitle = 'block';
numLines = 1;
defaultans = {'4','0','2','1'};
answer = inputdlg(prompt,dlgTitle,numLines,defaultans);
s = str2num(answer{1});
plotIt = str2num(answer{2});
condInt = str2num(answer{3});
newTactor = str2num(answer{4});

%% load in data of interest
stim = Stim.data;
fsStim = Stim.info.SamplingRateHz;
clear Stim

fsData = ECO1.info.SamplingRateHz;
sing = Sing.data;
fsSing = Sing.info.SamplingRateHz;
clear Sing

tact = Tact.data;
fsTact = Tact.info.SamplingRateHz;
clear Tact

clear Valu

%% figure out stim times
% vector of condition type - for first subject, looks like condition type
% is what was used , rather than test_condition
if s == 4
    condType = dlmread('rxnTime_condition_1_priming_5_cond.txt');
    train = dlmread('rxnTime_stimTrainDelivery_1_priming_5_cond.txt');
    priming =  dlmread('rxnTime_primingOption_1_priming_5_cond.txt');
elseif s == 3
    
    condType = dlmread('C:\TDT\OpenEx\MyProjects\StimResponseTiming\UserFiles\unblocked\rxnTime_condition_1_shortened.txt');
    train = dlmread('C:\TDT\OpenEx\MyProjects\StimResponseTiming\UserFiles\unblocked\rxnTime_stimTrainDelivery_1_shortened.txt');
    % stim cue from file
    
elseif s==2
    %     condType= dlmread('rxnTime_condition_2_modified_9_23_2016.txt');
    condType= dlmread('rxnTime_condition_2.txt');
    train = dlmread('rxnTime_stimTrainDelivery_2.txt');
    
elseif s==1
    %     condType= dlmread('rxnTime_condition_1_modified_9_23_2016.txt');
    condType= dlmread('rxnTime_condition_1.txt');
    train = dlmread('rxnTime_stimTrainDelivery_1.txt');
    % shrink condition type to be 140
end
%%
% stim cue from file
stimFromFile = tact(:,3);

% button press, start being onset of stimulation marker
trainTimes = find(stimFromFile~=0);

% pick condition type where stimulation was delivered
uniqueCond = unique(condType);
trainTimesTotal = {};

for i = 1:length(uniqueCond)
    trainTimesTotal{i} = trainTimes(condType==uniqueCond(i));
end


%% plot stim
if plotIt
    figure
    hold on
    for i = 1:size(stim,2)
        
        t = (0:length(stim)-1)/fsStim;
        subplot(3,2,i)
        plot(t*1e3,stim(:,i))
        title(sprintf('Channel %d',i))
        
    end
    
    xlabel('Time (ms)')
    ylabel('Amplitude (V)')
    subtitle('Stimulation Channels')
end

%% Sing looks like the wave to be delivered, with amplitude in uA
% 1st stim channel
Sing1 = sing(:,1);
% 2nd stim channel
Sing2 = sing(:,4);

samplesOfPulse = round(2*fsStim/1e3);

% build a burst table with the timing of stimuli
bursts = [];
bursts(1,:) = condType;
bursts(2,:) = trainTimes;
bursts(3,:) = trainTimes + samplesOfPulse;

stims1 = squeeze(getEpochSignal(Sing1,(bursts(2,condType>=2)-1),(bursts(3,condType>=2))+1));
t = (0:size(stims1,1)-1)/fsSing;
t = t*1e3;
if plotIt
    figure
    plot(t,stims1,'b','linewidth',2)
    xlabel('Time (ms)');
    ylabel('Current to be delivered (mA)')
    ylim([(min(stims1(:))-100) (max(stims1(:))+100)])
    title('Current to be delivered for all trials  on 1st channel')
end

stims2 = squeeze(getEpochSignal(Sing2,(bursts(2,condType==1)-1),(bursts(3,condType==1))+1));
t = (0:size(stims1,1)-1)/fsSing;
t = t*1e3;

if plotIt
    figure
    plot(t,stims2,'b','linewidth',2)
    xlabel('Time (ms)');
    ylabel('Current to be delivered (mA)')
    ylim([(min(stims2(:))-100) (max(stims2(:))+100)])
    title('Current to be delivered for all trials  on 2nd channel')
end

% delay loks to be 0.2867 ms from below.

%% Plot stims with info from above

% 1st stimulation channel
stim1 = stim(:,1);
stim1Epoched = squeeze(getEpochSignal(stim1,(bursts(2,condType>=2)-1),(bursts(3,condType>=2))+1));
t = (0:size(stim1Epoched,1)-1)/fsStim;
t = t*1e3;

if plotIt
    figure
    plot(t,stim1Epoched)
    xlabel('Time (ms)');
    ylabel('Voltage (V)');
    title('Finding the delay between current output and stim delivery - 1st stim channel')
    
    % hold on
    plot(t,stims1)
end
% get the delay in stim times
delay = round(0.2867*fsStim/1e3);

% plot the appropriately delayed signal
stimTimesBegin = bursts(2,condType>=2)-1+delay;
stimTimesEnd = bursts(3,condType>=2)-1+delay;
stim1Epoched = squeeze(getEpochSignal(stim1,stimTimesBegin,stimTimesEnd));
t = (0:size(stim1Epoched,1)-1)/fsStim;
t = t*1e3;
if plotIt
    figure
    plot(t,stim1Epoched)
    xlabel('Time (ms)');
    ylabel('Voltage (V)');
    title('Stim voltage monitoring with delay added in - 1st stim channel ')
end
% 2nd stimulation channel
stim2= stim(:,4);
stim2Epoched = squeeze(getEpochSignal(stim2,(bursts(2,condType==1)-1),(bursts(3,condType==1))+1));
t = (0:size(stim2Epoched,1)-1)/fsStim;
t = t*1e3;
if plotIt
    figure
    plot(t,stim2Epoched)
    xlabel('Time (ms)');
    ylabel('Voltage (V)');
    title('Finding the delay between current output and stim delivery - 2nd stim channel')
    
    % hold on
    plot(t,stims2)
end
% get the delay in stim times
delay = round(0.2867*fsStim/1e3);

% plot the appropriately delayed signal
stimTimesBegin = bursts(2,condType==1)-1+delay;
stimTimesEnd = bursts(3,condType==1)-1+delay;
stim2Epoched = squeeze(getEpochSignal(stim2,stimTimesBegin,stimTimesEnd+5));
t = (0:size(stim2Epoched,1)-1)/fsStim;
t = t*1e3;
if plotIt
    figure
    plot(t,stim2Epoched)
    xlabel('Time (ms)');
    ylabel('Voltage (V)');
    title('Stim voltage monitoring with delay added in - 2nd stim channel')
end
%% extract data
% try and account for delay for the stim times
stimTimes = bursts(2,:)+delay;
trainTimes=stimTimes;

% DJC 7-7-2016, changed presamps and post samps to 1 second
presamps = round(1 * fsData); % pre time in sec
postsamps = round(1 * fsData); % post time in sec, % modified DJC to look at up to 300 ms after

% sampling rate conversion between stim and data
fac = fsStim/fsData;

% find times where stims start in terms of data sampling rate
sts = round(stimTimes / fac);

%% look at tactor
tactorData = tact(:,1);

% threshold for new tactor for rough approximation
if newTactor
    tactorData(tactorData > 0.2) = 0.2;
end

tTact = (0:length(tactorData)-1)/fsTact;
figure
plot(tTact,tactorData);

title('tactor data')

% look at button press
buttonData = tact(:,2);
tButton = (0:length(buttonData)-1)/fsTact;
figure
plot(tButton,buttonData);

title('button data')

% look at stim from file saved
tStimFile = (0:length(stim)-1)/fsTact;
figure
plot(tStimFile,stimFromFile);
title('stim from file')

% look at all 3 + stim waveform
figure
ax1 = subplot(5,1,1);
plot(tTact,tactorData)
title('tactor data')

ax2 = subplot(5,1,2);
plot(tButton,buttonData);
title('button data')

ax3 = subplot(5,1,3);
plot(tStimFile,stimFromFile);
title('stim from file')

% assuming stim1 here is the channel where stim was being delivered
ax4 = subplot(5,1,4);
plot(tStimFile,stim1);
title('S1 Stim Channel')

%
ax5 = subplot(5,1,5);
plot(tStimFile,stim2);
title('Off Target Stim Channel')

%link axis
linkaxes([ax1,ax2,ax3,ax4,ax5],'x')

%% button press
% this is a coarse approximation, better is to use something like
% findchangepts
% find button data peaks
% set above certain threshold to 0.009

buttonDataClip = buttonData;
buttonDataClip(buttonData >= 0.009) = 0.009;
[buttonPks,buttonLocs] = findpeaks(buttonDataClip,tButton,'MinpeakDistance',2,'Minpeakheight',0.008);

figure
findpeaks(buttonDataClip,tButton,'MinpeakDistance',2,'Minpeakheight',0.008);
hold on
plot(tStimFile,stimFromFile,'g');
plot(tStimFile,stim1,'r');
%plot(t_stimFile,stim2,'m');

legend({'Button Data','Button Press Onset Peaks','Stimulation Times From File','S1 stim output','Off Target Stim Output'})
%% QUANTIFY RXN TIME TO CORTICAL STIM
sampsEnd = round(3.5*fsStim);

% epoched button press
epochedButton = {};

for i = 1:length(uniqueCond)
    epochedButton{i} = squeeze(getEpochSignal(buttonDataClip,trainTimesTotal{i},(trainTimesTotal{i} + sampsEnd)));
end

epochedTactor = squeeze(getEpochSignal(tactorData,trainTimesTotal{1},trainTimesTotal{1} + sampsEnd));
tEpoch = [0:size(epochedTactor,1)-1]/fsStim;

buttonPks = {};
buttonLocs = {};

buttonPksTempVec = [];
buttonLocsTempVec = [];

%%
for i = 1:length(uniqueCond)
    
    % for stimulation conditions
    if uniqueCond(i)~=-1
        for j = 1:length(trainTimesTotal{i})
            [buttonPksTemp,buttonLocsTemp] = findpeaks(epochedButton{i}(:,j),tEpoch,'NPeaks',1,'Minpeakheight',0.008);
            if isempty(buttonPksTemp)
                buttonPksTemp = NaN;
                buttonLocsTemp = NaN;
            end
            buttonPksTempVec(j) = buttonPksTemp;
            buttonLocsTempVec(j) = buttonLocsTemp;
        end
        buttonPks{i} = buttonPksTempVec;
        buttonLocs{i} = buttonLocsTempVec;
        
        % for tactor target condition
    elseif uniqueCond(i)==-1
        for j = 1:length(trainTimesTotal{i})
            
            [buttonPksTemp,buttonLocsTemp] = findpeaks((epochedButton{i}(tEpoch>1,j)),tEpoch(tEpoch>1),'NPeaks',1,'Minpeakheight',0.008);
            sprintf(['button ' num2str(buttonLocsTemp)])
            
            if newTactor
                [tactorPksTemp,tactorLocsTemp] = findpeaks((epochedTactor(:,j)),tEpoch,'NPeaks',1,'Minpeakheight',0.2);
            else
                [tactorPksTemp,tactorLocsTemp] = findpeaks((epochedTactor(:,j)),tEpoch,'NPeaks',1,'Minpeakheight',1);
                
            end
            
            sprintf(['tactor ' num2str(tactorLocsTemp)])
            
            if isempty(buttonPksTemp)
                buttonPksTemp = NaN;
                buttonLocsTemp = NaN;
            end
            
            if isempty(tactorPksTemp)
                tactorPksTemp = NaN;
                tactorLocsTemp = NaN;
            end
            
            buttonPksTempVec(j) = buttonPksTemp;
            buttonLocsTempVec(j) = buttonLocsTemp;
            
            tactorPksVec(j) = tactorPksTemp;
            tactorLocsVec(j) = tactorLocsTemp;
        end
        
        buttonPks{i} = buttonPksTempVec;
        buttonLocs{i} = buttonLocsTempVec;
        
    end
end

%%
% calculate differences - MAKE SURE YOU ONLY DO THIS ONCE
figure
buttonTactDiff = (buttonLocs{uniqueCond==-1} - tactorLocsVec);
buttonStimDiff = buttonLocs{uniqueCond==condInt};

respLo = 0.150;
respHi = 1;
subplot(1,2,1)
plot(buttonTactDiff,'b.','MarkerSize',50)
keep = find(buttonTactDiff >= respLo & buttonTactDiff <= respHi);
buttonTactDiffTrimmed = buttonTactDiff(keep);
hold on; plot(keep,buttonTactDiffTrimmed,'r.','MarkerSize',30)
hold on; plot([1 length(buttonTactDiff)],[respLo respLo],'k--')
hold on; plot([1 length(buttonTactDiff)],[respHi respHi],'k--')
avgRT = mean(buttonTactDiffTrimmed);
hold on; plot([1 length(buttonTactDiff)],[avgRT avgRT],'k','LineWidth',3)
ylabel('Response Time (s)')
xlabel('Trial')
title(['Tactor Average Response Time: ' num2str(mean(buttonTactDiffTrimmed)) ' s'])
ylim([-1.5 4])

respLo = 0.15;
respHi = 1;
subplot(1,2,2)
plot(buttonStimDiff,'b.','MarkerSize',50)
keep = find(buttonStimDiff >= respLo & buttonStimDiff <= respHi);
buttonStimDiffTrimmed = buttonStimDiff(keep);
hold on; plot(keep,buttonStimDiffTrimmed,'r.','MarkerSize',30)
hold on; plot([1 length(buttonStimDiff)],[respLo respLo],'k--')
hold on; plot([1 length(buttonStimDiff)],[respHi respHi],'k--')
avgRT = mean(buttonStimDiffTrimmed);
hold on; plot([1 length(buttonStimDiff)],[avgRT avgRT],'k','LineWidth',3)
ylabel('Response Time (s)')
xlabel('Trial')
title(['DCS Stim Average Response Time: ' num2str(mean(buttonStimDiffTrimmed)) ' s for condition = ' num2str(condInt)])
ylim([-1.5 4])

subtitle(['Time to Perception (TTP) = ' num2str(abs(mean(buttonStimDiffTrimmed-mean(buttonTactDiffTrimmed))*1000)) ' ms'])

