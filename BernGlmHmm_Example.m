clear
close all force
try
 rmdir('Test_Fit_Name', 's');
catch

end
% Fit Name (Where Data is Saved to)
settingsBernGlmHmm.FitName = "Test_Fit_Name";

% Predictors Of Interest
settingsBernGlmHmm.RelevantInputs = ["Stim","Bias","PrevChoice","WSLS"]; 

% Model Settings
settingsBernGlmHmm.stateMin = 2;
settingsBernGlmHmm.stateMax = 5;
settingsBernGlmHmm.prior_alphas = [1,2]; % Prior on Transition Matrix. High alpha =  results in more uniform probability distribution for state transition
settingsBernGlmHmm.prior_sigmas = [0.5 0.75 1 2 3]; % Prior on Weights in Glm. Higher sigma allows greater standard deviation of weights
settingsBernGlmHmm.AinitType = "Sticky"; % Controls Transition Matrix Initialisation 

% Fit Settings
settingsBernGlmHmm.folds = 10;
settingsBernGlmHmm.MaxIter = 150;
settingsBernGlmHmm.Tolerance = 10^-4;
settingsBernGlmHmm.DataParallelism = false;
settingsBernGlmHmm.ModelParallelism = false;

% Create Subject Data 
settingsBernGlmHmm.SubjectFlag = true;

% Fake Data
settingsBernGlmHmm.FakeData= true;
settingsBernGlmHmm.FakeTrialsPerSess= 97;
settingsBernGlmHmm.FakeProbL= 0.5;
settingsBernGlmHmm.FakeSessions = 11;
settingsBernGlmHmm.FakeNumberOfSubjects = 4;


[DataInput,settingsBernGlmHmm] = PrepareBernGlmHmm(settingsBernGlmHmm);
% If you have real data, Comment out above and uncomment below. 
% assertions will check to see if data is acceptable.
% Check FakeData Generated for an Example

%[DataOutput,settingsBernGlmHmm] = PrepareBernGlmHmm(settingsBernGlmHmm,data); 

FitType = "fullfit";
FullFitData = FitBernGlmHmm( DataInput, settingsBernGlmHmm, FitType);

PlotData = PlotBernGlmHmm(DataInput, settingsBernGlmHmm, FullFitData);