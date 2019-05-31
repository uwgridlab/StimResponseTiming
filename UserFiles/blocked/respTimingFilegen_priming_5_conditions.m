%% - 5-11-2018 - DJC - response timing file generator for 5 conditions of pilot
% stim train delivery file - this has the sample of delivery
%
% condition file - this has the conditions of interest 
% -1 = tactor
%  0 = null
%  1 = off-target
%  2 = train condition 1
%  3 = train condition 2
%  4 = train condition 3 
%  5 = train condition 4 (priming is on)
%  6 = train condition 5 (priming is on) 
% priming condition file - this tells the TDT whether or not priming is on
%  0 - priming off
%  1 - priming on
%
% all should be the same length 


prompt = {'Enter subject name','What is the range of ITI?', 'What is the sample rate of the TDT?','Which file number is this?'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'rxnTime','[1.5,2.5]','24414','1'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
sid = answer{1};
ITI = str2num(answer{2});
fs = str2num(answer{3});
fileNum = answer{4};


%% make the timing file
% add 1 to ITI - changed from 0.8 which we originally thought. This is to
% account for 200 ms of pulse, then the 800 ms of downtime, then the ITI

% so the range of the ITI is actually 2.5-3.5 s

ITIlo = ITI(1)+1;
ITIhi = ITI(2)+1;

% number of trials

numTrials = 40;
randTimes = unifrnd(ITIlo,ITIhi,numTrials,1);


% here the vector is converted to the sample number where the stimulus
% train should start to be delivered 
sample = 1; % start with sample 1 
pts = [];
for i = 1:length(randTimes)
    sample = round(sample + randTimes(i)*fs);
    pts = [pts; sample];
end

pts1 = pts;
clear pts


%% make the conditions file

% tactor 
tact = repmat(-1,20,1);

% no stim
noStim = repmat(0,10,1);

% stim conditions

%numbers of stimuli to delivery 
numEachStim = 20;

offTarg = repmat(1,numEachStim,1);
stim1 = repmat(2,numEachStim,1);
stim2 = repmat(3,numEachStim,1);
stim3 = repmat(4,numEachStim,1);
stim4 = repmat(5,numEachStim,1);
stim5 = repmat(6,numEachStim,1);

vectorCond = [stim4;stim5;noStim];
%vectorCond = [stim1;stim2;stim3;stim4;stim5;noStim;offTarg];
%vectorTact = [tact;noStim];
vectorTact = [];

vectorCondRand1 = vectorCond(randperm(length(vectorCond)));
vectorTactRand1 = vectorTact(randperm(length(vectorTact)));

ptsTotal = [pts1];
vectorCondRandTotal = [vectorCondRand1; vectorTactRand1];

%% account for latency in beeps for tactor touch of experimenter
% assuming 200 ms hi, 300 ms lo for each beep, and a train of 3 beeps where
% the experimenter is delivering a touch on the 3rd beep, subtract 1000 ms
% from each response time pt where the condition is for tactor touch 

timeToDelay = 1; % 1000 ms = 1 s
logicalValues = (vectorCondRandTotal == -1);
ptsTotal(logicalValues) = (ptsTotal(logicalValues)-(round(fs*timeToDelay)));


%% write these times to file for stim train delivery

filename = sprintf('%s_stimTrainDelivery_%s_priming_5_cond.txt',sid,fileNum);
fileID = fopen(filename,'w+');
fprintf(fileID,'%d\r\n',ptsTotal);
fclose(fileID);

%% write these times to file for condition 

filename = sprintf('%s_condition_%s_priming_5_cond.txt',sid,fileNum);
fileID = fopen(filename,'w+');
fprintf(fileID,'%d\r\n',vectorCondRandTotal);
fclose(fileID);
%% deal with priming 
% want to have stim3 and stim 4 have the priming

primingTotal = zeros(size(vectorCondRandTotal));
primingTotal(vectorCondRandTotal == 5 | vectorCondRandTotal == 6) = 1;

%% write these times for priming

filename = sprintf('%s_primingOption_%s_priming_5_cond.txt',sid,fileNum);
fileID = fopen(filename,'w+');
fprintf(fileID,'%d\r\n',primingTotal);
fclose(fileID);
