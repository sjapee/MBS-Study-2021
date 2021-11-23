%This script reads in behavioral data (Foils excluded) for the thresholding
%tasks (emotion and control) and plots and models the data.

%% initializations and setup
% cd to directory with preprocessed images you want to use
clear
inDir = '/misc/data12/sjapee/MoebiusSyndrome/Behavioral_Data/';
outDir = inDir;
cd(inDir);

taskName = 'StaticFace';
tasktag = 'SFaces';


inCtrlFname = 'ControlTask_Data/FinalGroupJan2020-withfoils.txt';
inEmotFname = 'Threshold_Durs_Catchboth_Data/FinalGroupJan2020-withfoils.txt';



inDataTask = 'C'; %'C' or 'E';
if inDataTask == 'C'
    inFname = inCtrlFname;
    datatag = 'Control';
else
    inFname = inEmotFname;
    datatag = 'Emotion';
end
analysisTag = '';
outxlfname = ['Jan2020_' tasktag '_' datatag 'Task-FitDataClean' analysisTag '.xlsx'];
outmatname = ['Jan2020_' tasktag '_' datatag 'Task-AllDataClean' analysisTag '.mat'];
outctrlsfigname = ['Jan2020_Controls_' tasktag '_' datatag 'Task-AllDataClean' analysisTag '.fig'];
outmosfigname = ['Jan2020_MoS_' tasktag '_' datatag 'Task-AllDataClean' analysisTag '.fig'];


%% read in data file
opts = detectImportOptions(inFname);
opts = setvartype(opts,{'Age'},'double');
inRawTable = readtable(inFname,opts);
summary(inRawTable);
SScodes = unique(inRawTable.Subject);
Ecodes = unique(inRawTable.MorphType);
Ecodes = {'Happy';'Fear';'Angry'};
Gcodes = unique(inRawTable.Group);
Gcodes = {'Control'; 'Patient'};

AllLevels = [0:0.05:1]';
nLevels = size(AllLevels,1);
totalSS = size(SScodes,1);
numG = size(Gcodes,1);


%% now separate Control subject's data from Patients
if inDataTask == 'C'
    findTaskType = 3
else
    findTaskType = 2;
end
CtrlsRawTable = inRawTable(find(strcmp(cellstr(inRawTable.Group),'Control')),:);
PtsRawTable = inRawTable(find(strcmp(cellstr(inRawTable.Group), 'Patient')),:);
GoodCtrlsRawTable = inRawTable(find(strcmp(cellstr(inRawTable.Group), 'Control') & ...
    inRawTable.TaskType == findTaskType & inRawTable.Age >= 13),:);
GoodPtsRawTable = inRawTable(find(strcmp(cellstr(inRawTable.Group), 'Patient') & ...
    inRawTable.TaskType ==findTaskType & inRawTable.Age >= 13),:);


%% do the splitting/grouping and plotting and fitting

for g = 1:size(Gcodes,1)
    g;
    if g == 1
        currGrpTable = CtrlsRawTable;
    else
        currGrpTable = PtsRawTable;
    end
    SScodes = unique(currGrpTable.Subject);
    numSS = size(SScodes,1)
    if g == 1
        numCtrls = numSS;
        CtrlSScodes = SScodes;
    else
        numPts = numSS;
        PtSScodes = SScodes;
    end
    numplotsx = 5;
    numplotsy = ceil(numSS/numplotsx);
    numploti = 1;
    currG = Gcodes{g};
    for s = 1:size(SScodes,1)
        currSS = SScodes{s};
        for e = 1:size(Ecodes,1)
            currE = Ecodes{e};
            if strcmp(currE , 'Angry')
                currColor = 'r';
            elseif strcmp(currE , 'Fear')
                currColor = 'b';
            elseif strcmp(currE , 'Happy')
                currColor = 'g';
            end
            currSS;
            %find all the data for the given group, given SS and given emotion
            currRawTable = currGrpTable(find(strcmp(cellstr(currGrpTable.Group), currG) & ...
                strcmp(cellstr(currGrpTable.Subject), currSS) & ...
                strcmp(cellstr(currGrpTable.MorphType),currE) & ...
                currGrpTable.TaskType == findTaskType),:);
            size(currRawTable);
            if ~isempty(currRawTable)
                currUseTable = currRawTable(find(currRawTable.CorrectRespType >=-1),:);
                %%%%%%%%%%%%%%%%
                % find the stable thresholds
                currMLs = currUseTable.Morph_;
                currG = currUseTable.Group{1};
                %%%%%%%%%%%%%%%%
                %fit the model and get param estimates
                currUseTable = currRawTable(find(currRawTable.CorrectRespType >=0),:);
                [MorphLoc,MorphLevels] = findgroups(currUseTable.Morph_);
                [currPCData] = splitapply(@mean,currUseTable.CorrectRespType,MorphLoc);
                [currNData] = splitapply(@numel,currUseTable.CorrectRespType,MorphLoc);
                Xmeasured = [MorphLevels];
                Fmeasured = [currPCData];
                Fmeasured(1,1) = 1 - Fmeasured(1,1); % invert the PC at level 0 since that is the neutral foil
                cleanCurrNData = currNData>5;
                cleanCurrNData = currNData(currNData>5);
                cleanCurrPCData = currPCData(currNData>5);
                cleanMorphLevels = MorphLevels(currNData>5);
                cleanXmeasured = [cleanMorphLevels];
                cleanFmeasured = [cleanCurrPCData];
                cleanFmeasured(1,1) = 1 - cleanFmeasured(1,1);
                if cleanFmeasured(1,1) > .3 %if the False alarm rate for the neutral foil is more than 30%, then do not use this point as the 0
                    cleanFmeasured(1) = [];
                    cleanXmeasured(1) = [];
                end
                currPCData = cleanCurrPCData;
                MorphLevels = cleanMorphLevels;
                Xmeasured = cleanXmeasured;
                Fmeasured = cleanFmeasured;
                %Xmeasured = [0; MorphLevels; 1];
                %Fmeasured = [0; currPCData; 1];
                figure(g+2); subplot(numplotsx,numplotsy,s); plot(Xmeasured,Fmeasured,[currColor '*-']); hold on;
                [ParamEstimates, SSE, Fit, Thresh50, Thresh79] = FitSModel_alpha(Xmeasured, Fmeasured, AllLevels);
                outSS{g,s,:} = currSS;
                outParams(g,s,e,:) = ParamEstimates;
                outSSE(g,s,e) = SSE;
                outFit(g,s,e,:) = Fit;
                outThresh(g,s,e,:) = [Thresh50 , Thresh79];
                nValues = size(MorphLevels,1);
                outML(g,s,e,:) = padarray(MorphLevels,[nLevels-nValues,0],-1,'post');
                outPC(g,s,e,:) = padarray(currPCData,[nLevels-nValues,0],-1,'post');
                newXnValues = size(Xmeasured,1);
                outXmeasured(g,s,e,:) = padarray(Xmeasured,[nLevels-newXnValues,0],-1,'post');
                outFmeasured(g,s,e,:) = padarray(Fmeasured,[nLevels-newXnValues,0],-1,'post');
                %%%%%%%%
                figure(g+2); subplot(numplotsx,numplotsy,s);plot(AllLevels,Fit,[currColor '--']);
                currG = currUseTable.Group{1};
                title([currG ' - ' currSS]);
            end;
        end;
    end;
    if g == 1
        numCtrls = numSS;
    else
        numPts = numSS;
    end
end


%%
%% write out the parameter estimates into a table file
T = table;
T.Group = [repmat(['Control'],numCtrls*3,1);repmat(['Patient'],numPts*3,1)];
T.Name =  [repmat([CtrlSScodes],3,1); repmat([PtSScodes],3,1)];
T.Emotion = [repmat(cellstr('1Happy'),numCtrls,1); repmat(cellstr('2Fearful'),numCtrls,1);repmat(cellstr('3Angry'),numCtrls,1); ...
    repmat(cellstr('1Happy'),numPts,1); repmat(cellstr('2Fearful'),numPts,1);repmat(cellstr('3Angry'),numPts,1)];
T.SSEs = [squeeze(outSSE(1,1:numCtrls,1))'; squeeze(outSSE(1,1:numCtrls,2))'; squeeze(outSSE(1,1:numCtrls,3))';  ...
    squeeze(outSSE(2,1:numPts,1))'; squeeze(outSSE(2,1:numPts,2))'; squeeze(outSSE(2,1:numPts,3))'];
T.Thresh50 = [squeeze(outThresh(1,1:numCtrls,1,1))'; squeeze(outThresh(1,1:numCtrls,2,1))'; squeeze(outThresh(1,1:numCtrls,3,1))';  ...
    squeeze(outThresh(2,1:numPts,1,1))'; squeeze(outThresh(2,1:numPts,2,1))'; squeeze(outThresh(2,1:numPts,3,1))'];
T.Slope = [squeeze(outParams(1,1:numCtrls,1,2))'; squeeze(outParams(1,1:numCtrls,2,2))'; squeeze(outParams(1,1:numCtrls,3,2))'; ...
    squeeze(outParams(2,1:numPts,1,2))'; squeeze(outParams(2,1:numPts,2,2))'; squeeze(outParams(2,1:numPts,3,2))'];
T.Fits = [squeeze(outFit(1,1:numCtrls,1,:)); squeeze(outFit(1,1:numCtrls,2,:)); squeeze(outFit(1,1:numCtrls,3,:));  ...
    squeeze(outFit(2,1:numPts,1,:)); squeeze(outFit(2,1:numPts,2,:)); squeeze(outFit(2,1:numPts,3,:))];

%write out xls file, mat file and figures with appropriate names
writetable(T,outxlfname);
savefig(3,outctrlsfigname);
savefig(4,outmosfigname);
save(outmatname);








