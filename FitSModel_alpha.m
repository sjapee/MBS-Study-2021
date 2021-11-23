%% function [ParaEstimates, SSE, Fit] = FitSModel_alpha(Xmeasured, Fmeasured)
% reads in Morph Levels and percentage yes responses and fits a sigmoidal
% function and returns param estimates, SSE and Fit

function [ParamEstimates, SSE, Fit, Thresh50 Thresh79] =  FitSModel_alpha(Xmeasured, Fmeasured, inLevels)
%% Setup the function to fit the data for each group, for each subject, for each task type (control and emotion detection) and for each morph type
% sigmoid:
% F(x) = a + (b-a)/(1+exp((alpha-x)/beta))

% Initializing [lb ub alpha beta]
% lb and ub, lower and upper boundaries
% alpha, inflexion point
% beta, steepness of the curve (beta close to 0: increased steepness)

alphaInit = 0.5;
betaInit = 0.03;
ParamInit = [alphaInit,betaInit];


% lb = min(inLevels);
% ub = max(inLevels);
lb = 0;
ub = 1;


SigFun = @(params,X) (lb + (ub - lb))./(1+exp((params(1) - X)./params(2)));
SigFunSum = @(params) sum((SigFun(params,Xmeasured)-Fmeasured).^2);
options = optimset('MaxFunEvals',100000);
[ParamEstimates , SSE] = fminsearch(SigFunSum,ParamInit,options);
Fit = SigFun(ParamEstimates,inLevels);
RevSigFun = @(params,Y) params(1) - (params(2).*log(((lb + (ub - lb))/Y)-1));
Thresh50 = RevSigFun(ParamEstimates,0.5);
Thresh79 = RevSigFun(ParamEstimates,0.794);
end
