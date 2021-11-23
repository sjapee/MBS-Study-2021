%This script reads in behavioral data (Foils excluded) for the thresholding
%tasks (emotion and control) and plots and models the data.

%% initializations and setup
% cd to directory with preprocessed images you want to use
clear
inDir = '/misc/data12/sjapee/MoebiusSyndrome/Behavioral_Data/';
outDir = inDir;
cd(inDir);

task = 1; %1 for Dynamic Faces and 2 for Dynamic Bodies
if task == 1
    taskName = 'DynamicFaceMorphs'; %DynamicFaceMorphs or %Threshold_AFHN_BodyMovements
    numDurs = 8;
else
    taskName = 'Threshold_AFHN_BodyMovements'; %DynamicFaceMorphs or %Threshold_AFHN_BodyMovements
    numDurs = 12;
end

inCtrlFname = [taskName '/ControlTask/FinalGroupJan2020_1.txt']; %%This means we are only using 1 run for everyone%%
inEmotFname = [taskName '/FinalGroupJan2020.txt'];

inDataTask = 'C'; %'C' or 'E';
if inDataTask == 'C'
    inFname = inCtrlFname
else
    inFname = inEmotFname;
end
if task == 1
    tasktag = 'DFaces';
else
    tasktag = 'DBodies';
end
if inDataTask == 'E'
    datatag = 'Emotion';
else
    datatag = 'Control';
end
analysisTag = 'add0Level';
outxlfname = ['Jan2020_' tasktag '_' datatag 'Task-FitDataClean-' analysisTag '.xlsx'];
outmatname = ['Jan2020_' tasktag '_' datatag 'Task-AllDataClean-' analysisTag '.mat'];
outctrlsfigname = ['Jan2020_Controls_' tasktag '_' datatag 'Task-AllDataClean-' analysisTag '.fig'];
outmosfigname = ['Jan2020_MoS_' tasktag '_' datatag 'Task-AllDataClean-' analysisTag '.fig'];



%% read in data file
opts = detectImportOptions(inFname);
%opts = setvartype(opts,{'Age'},'double');
%inRawTable = readtable(inCtrlFname);
inRawTable = readtable(inFname,opts);
summary(inRawTable);
SScodes = unique(inRawTable.Subject);
tempSScodes = unique(inRawTable.Subject);
temp = ~contains(tempSScodes,'AnCh');
SScodes = SScodes(temp);
Ecodes = unique(inRawTable.EmotionType);
Ecodes = {'H';'F';'A'};
Gcodes = unique(inRawTable.Group)
Gcodes = {'Control'; 'Patient'};

AllLevels = [1/(numDurs):1/(numDurs):1]';
nLevels = size(AllLevels,1);
totalSS = size(SScodes,1);
numG = size(Gcodes,1);


%% now separate Control subject's data from Patients
if inDataTask == 'C'
    findTaskType = 'Control';
else
    findTaskType = 'Emotion';
end
CtrlsRawTable = inRawTable(find(strcmp(cellstr(inRawTable.Group),'Control')),:);
PtsRawTable = inRawTable(find(strcmp(cellstr(inRawTable.Group), 'Patient')),:);
GoodCtrlsRawTable = inRawTable(find(strcmp(cellstr(inRawTable.Group), 'Control')),:);
GoodPtsRawTable = inRawTable(find(strcmp(cellstr(inRawTable.Group), 'Patient')),:);


%% do the splitting/grouping and plotting and fitting
for g = 1:size(Gcodes,1)
    g;
    if g == 1
        currGrpTable = CtrlsRawTable;
    else
        currGrpTable = PtsRawTable;
    end
    tempSScodes = unique(currGrpTable.Subject);
    temp = ~contains(tempSScodes,'AnCh');
    SScodes = tempSScodes(temp);
    numSS = size(SScodes,1)
    if g == 1
        numCtrls = numSS;
        CtrlSScodes = SScodes;
    else
        numPts = numSS;
        PtSScodes = SScodes;
    end
    numplotsx = 4;
    numplotsy = ceil(numSS/numplotsx);
    numploti = 1;
    currG = Gcodes{g}
    for s = 1:size(SScodes,1)
        currSS = SScodes{s};
        for e = 1:size(Ecodes,1)
            currE = Ecodes{e};
            if strcmp(currE , 'A')
                currColor = 'r';
            elseif strcmp(currE , 'F')
                currColor = 'b';
            elseif strcmp(currE , 'H')
                currColor = 'g';
            end
            currSS;
            %find all the data for the given group, given SS and given emotion
            currRawTable = currGrpTable(strcmp(cellstr(currGrpTable.Group), currG) & ...
                strcmp(cellstr(currGrpTable.Subject), currSS) & ...
                strcmp(cellstr(currGrpTable.EmotionType),currE),:);
            size(currRawTable);
            if ~isempty(currRawTable)
                currUseTable = currRawTable(find(currRawTable.AllCorrectRespType >=-1),:);
                %%%%%%%%%%%%%%%%
                %fit the model and get param estimates
                currMLs = currUseTable.DurationCode;
                currG = currUseTable.Group{1};
                %%%%%%%%%%%%%%%%
                currUseTable = currRawTable(find(currRawTable.AllCorrectRespType >=0),:);
                [MorphLoc,MorphLevels] = findgroups(currUseTable.DurationCode);
                MorphLevels = MorphLevels/(size(AllLevels,1));
                [currPCData] = splitapply(@mean,currUseTable.AllCorrectRespType,MorphLoc);
                [currNData] = splitapply(@numel,currUseTable.AllCorrectRespType,MorphLoc);
                Xmeasured = [MorphLevels];
                Fmeasured = [currPCData];
                cleanCurrNData = currNData>5;
                cleanCurrNData = currNData(currNData>5);
                cleanCurrPCData = currPCData(currNData>5);
                cleanMorphLevels = MorphLevels(currNData>5);
                cleanXmeasured = [cleanMorphLevels];
                cleanFmeasured = [cleanCurrPCData];
                if inDataTask == 'E'
                    currPCData = cleanCurrPCData;
                    MorphLevels = cleanMorphLevels;
                    Xmeasured = cleanXmeasured;
                    Fmeasured = cleanFmeasured;
                end
                figure(g+2); subplot(numplotsx,numplotsy,s); plot(Xmeasured,Fmeasured,[currColor '*-']); hold on;
                [ParamEstimates, SSE, Fit, Thresh50, Thresh79] = FitSModel_alphaDynamic([0;Xmeasured], [0;Fmeasured], [0;AllLevels]);
                outSS{g,s,:} = currSS;
                outParams(g,s,e,:) = ParamEstimates;
                outSSE(g,s,e) = SSE;
                outFit(g,s,e,:) = Fit;
                outThresh(g,s,e,:) = [Thresh50 , Thresh79];
                nValues = size(MorphLevels,1);
                %need to pad -1s here so size is the same if there are missing MLs due to all trials being misses
                nPad = size(AllLevels,1) - nValues;
                tempML = [MorphLevels];
                outML(g,s,e,:) = padarray(tempML, nPad, -1, 'post');
                tempPC = [currPCData];
                outPC(g,s,e,:) = padarray(tempPC, nPad, -1, 'post');
                tempX = Xmeasured;
                outXmeasured(g,s,e,:) = padarray(tempX, nPad, -1, 'post');
                tempF = Fmeasured;
                outFmeasured(g,s,e,:) = padarray(tempF, nPad, -1, 'post');
                %%%%%%%%
                figure(g+2); subplot(numplotsx,numplotsy,s);plot([0;AllLevels],Fit,[currColor '--']);
                currG = currUseTable.Group{1};
                title([currG ' - ' currSS]);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %add a part here to compute the mean across the three emotions and
        %then do the fitting for these data
        e = 4;
        currE = 'M';
        currColor = 'k';
        %find all the data for the given group, given SS and given emotion
        currRawTable = currGrpTable(strcmp(cellstr(currGrpTable.Group), currG) & ...
            strcmp(cellstr(currGrpTable.Subject), currSS),:);
        size(currRawTable);
        if ~isempty(currRawTable)
            currUseTable = currRawTable(find(currRawTable.AllCorrectRespType >=-1),:);
            %%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%
            %fit the model and get param estimates
            currMLs = currUseTable.DurationCode;
            currG = currUseTable.Group{1};
            currUseTable = currRawTable(find(currRawTable.AllCorrectRespType >=0),:);
            [MorphLoc,MorphLevels] = findgroups(currUseTable.DurationCode);
            MorphLevels = MorphLevels/(size(AllLevels,1));
            [currPCData] = splitapply(@mean,currUseTable.AllCorrectRespType,MorphLoc);
            [currNData] = splitapply(@numel,currUseTable.AllCorrectRespType,MorphLoc);
            cleanCurrNData = currNData>5;
            cleanCurrNData = currNData(currNData>5);
            cleanCurrPCData = currPCData(currNData>5);
            cleanMorphLevels = MorphLevels(currNData>5);
            cleanXmeasured = [cleanMorphLevels];
            cleanFmeasured = [cleanCurrPCData];
            currPCData = cleanCurrPCData;
            MorphLevels = cleanMorphLevels;
            Xmeasured = [cleanXmeasured];
            Fmeasured = [cleanFmeasured];
            figure(g+2); subplot(numplotsx,numplotsy,s); plot(Xmeasured,Fmeasured,[currColor '*-']); hold on;
            [ParamEstimates, SSE, Fit, Thresh50, Thresh79] = FitSModel_alphaDynamic([0;Xmeasured], [0;Fmeasured], [0;AllLevels]);
            outSS{g,s,:} = currSS;
            outParams(g,s,e,:) = ParamEstimates;
            outSSE(g,s,e) = SSE;
            outFit(g,s,e,:) = Fit;
            outThresh(g,s,e,:) = [Thresh50 , Thresh79];
            nValues = size(MorphLevels,1);
            %need to pad -1s here so size is the same if there are missing MLs due to all trials being misses
            nPad = size(AllLevels,1) - nValues;
            tempML = [MorphLevels];
            outML(g,s,e,:) = padarray(tempML, nPad, -1, 'post');
            tempPC = [currPCData];
            outPC(g,s,e,:) = padarray(tempPC, nPad, -1, 'post');
            tempX = Xmeasured;
            outXmeasured(g,s,e,:) = padarray(tempX, nPad, -1, 'post');
            tempF = Fmeasured;
            outFmeasured(g,s,e,:) = padarray(tempF, nPad, -1, 'post');
            %%%%%%%%
            figure(g+2); subplot(numplotsx,numplotsy,s);plot([0;AllLevels],Fit,[currColor '--']);
            currG = currUseTable.Group{1};
            title([currG ' - ' currSS]);
        end
    end
    if g == 1
        numCtrls = numSS;
    else
        numPts = numSS;
    end
end

%%
%write out the data into tables
T = table;
T.Group = [repmat(['Control'],numCtrls*4,1);repmat(['Patient'],numPts*4,1)];
T.Name =  [repmat([CtrlSScodes],4,1); repmat([PtSScodes],4,1)];
T.Emotion = [repmat(cellstr('1Happy'),numCtrls,1); repmat(cellstr('2Fearful'),numCtrls,1);repmat(cellstr('3Angry'),numCtrls,1);repmat(cellstr('4Average'),numCtrls,1); ...
    repmat(cellstr('1Happy'),numPts,1); repmat(cellstr('2Fearful'),numPts,1);repmat(cellstr('3Angry'),numPts,1);repmat(cellstr('4Average'),numPts,1)];
T.SSEs = [squeeze(outSSE(1,1:numCtrls,1))'; squeeze(outSSE(1,1:numCtrls,2))'; squeeze(outSSE(1,1:numCtrls,3))'; squeeze(outSSE(1,1:numCtrls,4))'; ...
    squeeze(outSSE(2,1:numPts,1))'; squeeze(outSSE(2,1:numPts,2))'; squeeze(outSSE(2,1:numPts,3))'; squeeze(outSSE(2,1:numPts,4))'];
T.Thresh50 = [squeeze(outThresh(1,1:numCtrls,1,1))'; squeeze(outThresh(1,1:numCtrls,2,1))'; squeeze(outThresh(1,1:numCtrls,3,1))'; squeeze(outThresh(1,1:numCtrls,4,1))'; ...
    squeeze(outThresh(2,1:numPts,1,1))'; squeeze(outThresh(2,1:numPts,2,1))'; squeeze(outThresh(2,1:numPts,3,1))'; squeeze(outThresh(2,1:numPts,4,1))'];
T.Slope = [squeeze(outParams(1,1:numCtrls,1,2))'; squeeze(outParams(1,1:numCtrls,2,2))'; squeeze(outParams(1,1:numCtrls,3,2))'; squeeze(outParams(1,1:numCtrls,4,2))'; ...
    squeeze(outParams(2,1:numPts,1,2))'; squeeze(outParams(2,1:numPts,2,2))'; squeeze(outParams(2,1:numPts,3,2))'; squeeze(outParams(2,1:numPts,4,2))'];
T.Fits = [squeeze(outFit(1,1:numCtrls,1,:)); squeeze(outFit(1,1:numCtrls,2,:)); squeeze(outFit(1,1:numCtrls,3,:)); squeeze(outFit(1,1:numCtrls,4,:)); ...
    squeeze(outFit(2,1:numPts,1,:)); squeeze(outFit(2,1:numPts,2,:)); squeeze(outFit(2,1:numPts,3,:)); squeeze(outFit(2,1:numPts,4,:))];

%write out xls file, mat file and figures with appropriate names
writetable(T,outxlfname);
savefig(3,outctrlsfigname);
savefig(4,outmosfigname);
save(outmatname);

