function [DataOutput,settingsBernGlmHmm] = PrepareBernGlmHmm(settingsBernGlmHmm,data)


    assert(isstruct(settingsBernGlmHmm),"Argument 1 must be a Struct") 

    FakeDataLogic = isfield(settingsBernGlmHmm,"FakeData");
    % Fake Data
    if FakeDataLogic == false
        settingsBernGlmHmm.FakeData = false;
    end

    if settingsBernGlmHmm.FakeData == false && isempty(data)
        error("No Data Table Provided")
    end
    
    % --------- Checking Fake Data--------------
    if settingsBernGlmHmm.FakeData == true
        
        FakeTrialsPerSessLogic = isfield(settingsBernGlmHmm,"FakeTrialsPerSess") &&...
            all((mod(settingsBernGlmHmm.FakeTrialsPerSess, 1) == 0) == 1) &&...
            all(settingsBernGlmHmm.FakeTrialsPerSess >= 1);
        assert(FakeTrialsPerSessLogic,"settingsBernGlmHmm.FakeTrialsPerSess = Positive Integer ");

        FakeProbLLogic = isfield(settingsBernGlmHmm,"FakeProbL") &&...
            settingsBernGlmHmm.FakeProbL > 0 &&...
            settingsBernGlmHmm.FakeProbL < 1;
        assert(FakeProbLLogic,"settingsBernGlmHmm.FakeProbL = 0 < Float < 1 ");

        FakeSessionsLogic = isfield(settingsBernGlmHmm,"FakeSessions") &&...
            all((mod(settingsBernGlmHmm.FakeSessions, 1) == 0) == 1) &&...
            all(settingsBernGlmHmm.FakeSessions >= 1);
        assert(FakeSessionsLogic,"settingsBernGlmHmm.FakeSessions = Positive Integer ");

        FakeNumberOfSubjectsLogic =  isfield(settingsBernGlmHmm,"FakeSessions") &&...
            all((mod(settingsBernGlmHmm.FakeNumberOfSubjects, 1) == 0) == 1) &&...
            all(settingsBernGlmHmm.FakeNumberOfSubjects >= 1);
        assert(FakeNumberOfSubjectsLogic,"settingsBernGlmHmm.FakeNumberOfSubjects = Positive Integer ");

    end
    


    if settingsBernGlmHmm.FakeData == true
        data =  CreateFakeMouseData(settingsBernGlmHmm);
    end

    assert(istable(data),"Argument 2 must be a Table")
    
    
    settingsBernGlmHmm = CreateParamStruct_ver2(data,settingsBernGlmHmm);
    
    DataOutput = PrepareHMMdata_ver2(data,settingsBernGlmHmm);

    BinaryCheck = all(cellfun(@(x) all(x ==1 | x == 0) ,DataOutput.TrainData.Global.OutputFoldCell) == 1);
    assert(BinaryCheck == true,"Only two catgeories allowed in the output Vector (0 and 1)") 
    obs_dim = all(cellfun(@(x) size(x,2) ,DataOutput.TrainData.Global.OutputFoldCell) == 1);
    assert(obs_dim == true,"Observation diemnsion should be 1")
    
    CovariateRangeCheck = all(cellfun(@(x)all(all(x >= -1 & x <= 1)),DataOutput.TrainData.Global.InputFoldCell) == 1);
    if CovariateRangeCheck == 0
        warning("Covariates May be difficult to interpret")
    end

end




function settingsBernGlmHmm = CreateParamStruct_ver2(data,settingsBernGlmHmm)
        
    % CHECK INPUTS BEGIN
    
    %---------------------------------------------------------------------
    %-------------------------- Check Data Table -------------------------
    %---------------------------------------------------------------------
    TableVariables = data.Properties.VariableNames;
    TableLogic = all(ismember(["SubjectName","Trial","TrialTag","Session","Intervention"],TableVariables) == 1);
    assert(TableLogic,"Table requires columns 'SubjectName','Trial','TrialTag','Session','Intervention'")
    
    TableSubjectLogic = all(cellfun(@(x)ischar(x) ,cellstr(data.SubjectName)) == 1 );
    assert(TableSubjectLogic,"In Data Table Provided, SubjectNames Not All Strings ");
    
    TableTrialLogic = all((mod(data.Trial, 1) == 0) == 1) && all(data.Trial >= 1);
    assert(TableTrialLogic,"In Data Table Provided, All Trial Values Should Be Positive Integers");
    
    TableTrialTagLogic = all((mod(data.TrialTag, 1) == 0) == 1)  && ( length(data.TrialTag) == length(unique(data.TrialTag)) ) && all(data.Trial >= 1);
    assert(TableTrialTagLogic,"In Data Table Provided, All TrialTag Values Should Be Unique Positive Integers");
    
    TableSessionLogic = all((mod(data.Session, 1) == 0) == 1) && all(data.Trial >= 1);
    assert(TableSessionLogic,"In Data Table Provided, All Sessions Values Should Be Positive Integers");
    
    TableInterventionLogic = all((mod(data.Intervention, 1) == 0) == 1) && all(data.Intervention >= 0); 
    assert(TableInterventionLogic,"In Data Table Provided, All Intervention Values Should Be Greater Than Zero And Integers");

    %---------------------------------------------------------------------
    %--------------------- Check settingsBernGlmHmm ----------------------
    %---------------------------------------------------------------------
    % Check Fit Name (Where Data is Saved to)
    assert(isfield(settingsBernGlmHmm,"FitName") && ~isempty(settingsBernGlmHmm.FitName),"settingsBernGlmHmm.FitName = [Insert Fit Name]");
    
    % Check Predictors Of Interest
    assert(isfield(settingsBernGlmHmm,"RelevantInputs"),"settingsBernGlmHmm.RelevantInputs = [String Array]"); 
    
    PossibleNames = data.Properties.VariableNames;
    RelevantInputs = ["Stim","Bias","PrevChoice","WSLS"];
    assert(all(ismember(RelevantInputs,PossibleNames)==1),"Relevant Inputs not present in Data Table Provided");
    
    % Check Model Settings
    assert(settingsBernGlmHmm.stateMin >= 2,"Minimum States must be greater than 2. settingsBernGlmHmm.stateMin = 2");
    assert(settingsBernGlmHmm.stateMax <= 6,"Maximum States must be less than 6. settingsBernGlmHmm.stateMin = 6 ");
    
    PriorAlphaLogic = isfield(settingsBernGlmHmm,"prior_alphas") &&  all((mod(settingsBernGlmHmm.prior_alphas, 1) == 0) == 1) && all(settingsBernGlmHmm.prior_alphas >= 1);
    assert(PriorAlphaLogic,"settingsBernGlmHmm.prior_alphas = [Positive Integer Array] >= 1")
    
    PriorSigmaLogic = isfield(settingsBernGlmHmm,"prior_sigmas")  && all(settingsBernGlmHmm.prior_sigmas > 0);
    assert(PriorSigmaLogic,"settingsBernGlmHmm.prior_sigmas = [Positive Float Array] >= 0")
    
    AinitTypeLogic = isfield(settingsBernGlmHmm,"AinitType") && ismember(settingsBernGlmHmm.AinitType,["Dir","Uni","Sticky"]);
    assert(AinitTypeLogic,"settingsBernGlmHmm.AinitType = 'Dir' OR 'Uni' OR 'Sticky' ")
    
    % Check Fit Settings
    NumberOfSessions = numel(unique(data.Session));
    FoldLogic = isfield(settingsBernGlmHmm,"folds") &&...
        all(settingsBernGlmHmm.folds <= NumberOfSessions) &&...
        all((mod(settingsBernGlmHmm.folds, 1) == 0) == 1) &&...
        all(settingsBernGlmHmm.folds >= 1);
    assert(FoldLogic,"settingsBernGlmHmm.folds = [ Postive Integer <= Number Of Sessions ] ")
    
    MaxIterLogic = isfield(settingsBernGlmHmm,"MaxIter") && all((mod(settingsBernGlmHmm.MaxIter, 1) == 0) == 1) && settingsBernGlmHmm.MaxIter >= 10;
    assert(MaxIterLogic,"settingsBernGlmHmm.MaxIter = [Positive Integer >= 10]")
    
    ToleranceLogic = isfield(settingsBernGlmHmm,"Tolerance") && settingsBernGlmHmm.Tolerance > 0 && settingsBernGlmHmm.Tolerance <= 0.01;
    assert(ToleranceLogic,"settingsBernGlmHmm.Tolerance = [Positive Float <= 0.01] ")
    
    DataParallelismLogic = isfield(settingsBernGlmHmm,"DataParallelism") && islogical(settingsBernGlmHmm.DataParallelism);
    assert(DataParallelismLogic,"settingsBernGlmHmm.DataParallelism = 'true' OR 'false' ")
    
    ModelParallelismLogic = isfield(settingsBernGlmHmm,"ModelParallelism") && islogical(settingsBernGlmHmm.ModelParallelism);
    assert(ModelParallelismLogic,"settingsBernGlmHmm.ModelParallelism = 'true' OR 'false' ")

    ExclusiveParallelismLogic = not(settingsBernGlmHmm.DataParallelism == true &&  settingsBernGlmHmm.ModelParallelism == true);
    assert(ExclusiveParallelismLogic,"Model Parallelism or Data Parallelism, but not both")
    
    % Create Subject Data 
    SubjectFlagLogic = isfield(settingsBernGlmHmm,"SubjectFlag") && islogical(settingsBernGlmHmm.SubjectFlag);
    assert(SubjectFlagLogic,"settingsBernGlmHmm.SubjectFlag = 'true' OR 'false' ");

    if settingsBernGlmHmm.SubjectFlag == true
        
        Subjects = unique(data.SubjectName);
        for sub = 1:1:length(Subjects)  
            MaxFolds = max(unique(data{data.SubjectName == Subjects(sub),"Session"}));
            SubjectFoldsLogic = settingsBernGlmHmm.folds <= MaxFolds;
            assert(SubjectFoldsLogic, Subjects(sub) + ": Number Of Folds Must be less than Maximum Number Of Sesions" )
        end

        
    end
    

   
    % CHECK INPUTS OVER

    settingsBernGlmHmm.SubjectNames = unique(data.SubjectName);
    settingsBernGlmHmm.NumberOfSubjects = numel(data.SubjectName);
    
    stateMin = settingsBernGlmHmm.stateMin;
    stateMax = settingsBernGlmHmm.stateMax;
    settingsBernGlmHmm.States = stateMin:1:stateMax;
    
    settingsBernGlmHmm.NumberOfFolds = settingsBernGlmHmm.folds;
    settingsBernGlmHmm.NumberOfModels = numel(settingsBernGlmHmm.States);
    settingsBernGlmHmm.NumberOfAlphas = numel(settingsBernGlmHmm.prior_alphas);
    settingsBernGlmHmm.NumberOfSigmas = numel(settingsBernGlmHmm.prior_sigmas);

    settingsBernGlmHmm.InputDim = length(settingsBernGlmHmm.RelevantInputs);



end




function [DataOut] = PrepareHMMdata_ver2(data,settingsBernGlmHmm)

    Folds = settingsBernGlmHmm.folds;
    FitName = settingsBernGlmHmm.FitName;
    RelevantInputs = settingsBernGlmHmm.RelevantInputs;
    SubjectFlag = settingsBernGlmHmm.SubjectFlag;
    


    [TrainSessionFold, TestSessionFold] = CreateFolds(data,Folds,FitName,"Global");
    FoldData.Global.TrainSessionFold = TrainSessionFold;
    FoldData.Global.TestSessionFold = TestSessionFold;
    [InputFoldCell,OutputFoldCell,SessionFoldCell,LaserCell,TrialTagFoldCell] = DataForHmm(data,TrainSessionFold,RelevantInputs,Folds,FitName,"Global","Train");
    TrainData.Global.InputFoldCell = InputFoldCell;
    TrainData.Global.OutputFoldCell = OutputFoldCell;
    TrainData.Global.SessionFoldCell = SessionFoldCell;
    TrainData.Global.LaserCell = LaserCell;
    TrainData.Global.TrialTagFoldCell = TrialTagFoldCell;
    
    [InputFoldCellTest,OutputFoldCellTest,SessionFoldCellTest,LaserCellTest,TrialTagFoldCellTest] = DataForHmm(data,TestSessionFold,RelevantInputs,Folds,FitName,"Global","Test");
    TestData.Global.InputFoldCellTest = InputFoldCellTest;
    TestData.Global.OutputFoldCellTest = OutputFoldCellTest;
    TestData.Global.SessionFoldCellTest = SessionFoldCellTest;
    TestData.Global.LaserCellTest = LaserCellTest;
    TestData.Global.TrialTagFoldCellTest = TrialTagFoldCellTest;


    
    if SubjectFlag == true

        SubjectNames = settingsBernGlmHmm.SubjectNames;
        for sub = 1:1:length(SubjectNames)
            SubjectName = SubjectNames{sub};
            dataSubject = data(ismember(data.SubjectName,SubjectName),: );
            [TrainSessionFold, TestSessionFold] = CreateFolds(dataSubject,Folds,FitName,SubjectName);
            FoldData.(SubjectName).TrainSessionFold = TrainSessionFold;
            FoldData.(SubjectName).TestSessionFold = TestSessionFold;
            
            [InputFoldCell,OutputFoldCell,SessionFoldCell,LaserCell,TrialTagFoldCell] = DataForHmm(dataSubject,TrainSessionFold,RelevantInputs,Folds,FitName,SubjectName,"Train");
            TrainData.(SubjectName).InputFoldCell = InputFoldCell;
            TrainData.(SubjectName).OutputFoldCell = OutputFoldCell;
            TrainData.(SubjectName).SessionFoldCell = SessionFoldCell;
            TrainData.(SubjectName).LaserCell = LaserCell;
            TrainData.(SubjectName).TrialTagFoldCell = TrialTagFoldCell;
        
            [InputFoldCellTest,OutputFoldCellTest,SessionFoldCellTest,LaserCellTest,TrialTagFoldCellTest] = DataForHmm(dataSubject,TestSessionFold,RelevantInputs,Folds,FitName,SubjectName,"Test");
            TestData.(SubjectName).InputFoldCellTest = InputFoldCellTest;
            TestData.(SubjectName).OutputFoldCellTest = OutputFoldCellTest;
            TestData.(SubjectName).SessionFoldCellTest = SessionFoldCellTest;
            TestData.(SubjectName).LaserCellTest = LaserCellTest;
            TestData.(SubjectName).TrialTagFoldCellTest = TrialTagFoldCellTest;      
        end

    end
    
    % Return this data
    DataOut.FoldData = FoldData;
    DataOut.TrainData = TrainData;
    DataOut.TestData = TestData;

end



function [InputFoldCell,OutputFoldCell,SessionFoldCell,InterventionCell,TrialTagFoldCell] = DataForHmm(data,SessionFold,RelevantInputs,Folds,FitName,SubjectName,Type)

    InputFoldCell  = cell(1,Folds);
    OutputFoldCell = cell(1,Folds);
    SessionFoldCell = cell(1,Folds);
    TrialTagFoldCell = cell(1,Folds);
    for f = 1:1:Folds
        Input = data{ismember(data.Session,SessionFold{f}),RelevantInputs};
        Output = data{ismember(data.Session,SessionFold{f}),"Output"};
        Session = data{ismember(data.Session,SessionFold{f}),"Session"};
        TrialTag =  data{ismember(data.Session,SessionFold{f}),"TrialTag"};
        
        InputFoldCell{f}  = Input;
        OutputFoldCell{f} = Output;
        SessionFoldCell{f} = Session;
        TrialTagFoldCell{f} = TrialTag;

    end
    Sessions = unique(data.Session);
    NumberOfSessions = numel(Sessions);
    InterventionCell = cell(1,NumberOfSessions);
    for sess = 1:1:NumberOfSessions   
        InterventionIdx = logical( data{data.Session == sess,"Intervention"}' );
        SessionLength = numel(InterventionIdx);
        x = 1:1:SessionLength;
        InterventionCell{sess} = x(InterventionIdx);
    end

    DataFolder = fullfile(pwd,FitName,"Prepared_Data",SubjectName,Type);    
    if ~exist("Prepared_DataFolder","dir")
        mkdir(DataFolder)
    end
    
    DataHmm = fullfile(DataFolder, "InputFoldCell" + ".mat");
    save(DataHmm,'InputFoldCell');

    DataHmm = fullfile(DataFolder, "OutputFoldCell" + ".mat");
    save(DataHmm,'OutputFoldCell');

    DataHmm = fullfile(DataFolder, "SessionFoldCell" + ".mat");
    save(DataHmm,'SessionFoldCell');

    DataHmm = fullfile(DataFolder, "InterventionCell" + ".mat");
    save(DataHmm,'InterventionCell');

    DataHmm = fullfile(DataFolder, "TrialTagFoldCell" + ".mat");
    save(DataHmm,'TrialTagFoldCell');

end




function [TrainSessionFold, TestSessionFold] = CreateFolds(data,Folds,FitName,SubjectName)
    disp("Creating Folds For: " + SubjectName)
    Session = unique(data.Session)'; % NB Transpose to Row Vector
    MaxSession = numel(Session);
        
    % indices_or_sections is a scalar, not an array.
    Nsections = int64(Folds);
    assert(Nsections > 0,"number sections must be larger than 0.")
    
    assert(Nsections > 0,"number sections must be larger than 0.")
    assert(MaxSession>=Nsections,"Folds must be less than the number of Sessions")

    Neach_section = floor(MaxSession/Folds);
    extras = mod(MaxSession,Folds);
    
    ExtraDiv = repelem(Neach_section+1,extras);
    NormalDiv = repelem(Neach_section,Nsections-extras);
    
    section_sizes = [0 ExtraDiv NormalDiv];
    div_points = cumsum(section_sizes) + 1;
    
    
    FoldArray = cell(1,Nsections);
    for i = 1:1:Nsections
        st = div_points(i);
        et = div_points(i+1);
        FoldArray{i} = Session(st:et-1); 
    end
    
    folds = 1:1:Folds;
    
    TrainFolds = nchoosek(folds,Folds-1);
    TrainFoldsNan = [TrainFolds,nan(Folds,1)];
    TestFolds = zeros(1,Folds);
    for f = 1:1:Folds
        TestFolds(f) = setdiff(folds,TrainFoldsNan(f,:));
    end
    
    TrainSessionFold = cell(1,Folds);
    TestSessionFold = cell(1,Folds);
    for f = 1:1:Folds
        TrainSessionFold{f} = cell2mat( FoldArray(TrainFolds(f,:)) );
        TestSessionFold{f} = cell2mat( FoldArray(TestFolds(1,f)) );
    end

    FoldFolder = fullfile(pwd,FitName,"Prepared_FoldData",SubjectName);    
    if ~exist("Prepared_FoldFolder","dir")
        mkdir(FoldFolder)
    end

    FoldsHmm = fullfile(FoldFolder, "TrainSessionFold" + ".mat");
    save(FoldsHmm,'TrainSessionFold');

    FoldsHmm = fullfile(FoldFolder, "TestSessionFold" + ".mat");
    save(FoldsHmm,'TestSessionFold');

end

























function data =  CreateFakeMouseData(settingsBernGlmHmm)


TrialsPerSess = settingsBernGlmHmm.FakeTrialsPerSess;
ProbL = settingsBernGlmHmm.FakeProbL;
Sessions = settingsBernGlmHmm.FakeSessions;
NumberOfSubjects = settingsBernGlmHmm.FakeNumberOfSubjects;

NT = repmat(TrialsPerSess,[1 Sessions]);
ProbLeftT = repmat(ProbL,[1 Sessions]);

output = [];
stimOutput = [];
bias = [];
prevChoice = [];
reward = [];
wSLS = [];
state = [];
Target = [];
trialNum = [];
session = [];
LaserOni = [];
SubjectName = [];

assert(length(NT) == length(ProbLeftT),"N ~= ProbLeft")
NumberOfSessions = length(NT);
    
Nstart = 1;
Nend = 0;

for subi = 1:1:NumberOfSubjects

    suboutput = [];
    substimOutput = [];
    subbias = [];
    subprevChoice = [];
    subreward = [];
    subwSLS = [];
    substate = [];

    subTarget = [];
    subtrialNum = [];
    
    subsession = [];
    subLaserOni = [];



    for sess = 1:1:Sessions
    
        N = NT(sess);
        ProbLeft = ProbLeftT(sess);
    
        Targeti = zeros(1,N);
        for t = 1:1:N
            SideProb = rand(1,1);
            if SideProb < ProbLeft
                Targeti(t) = 3;
            else
                Targeti(t) = 4;
            end
        end
        %histcounts(Targeti)
        
        Cuei = zeros(1,length(Targeti));
        for i = 1:1:length(Targeti)   
            if Targeti(i) == 3
                Cuei(i) = randi([0,49])/100;
            else
                Cuei(i) = randi([50,100])/100;
            end
        end
       
        
        % Create Choice Data 
        
        NumTrials = length(Targeti);
        N1 = (NumTrials/4)*1;
        N2 = (NumTrials/4)*2;
        N3 = (NumTrials/4)*3;
        N4 = (NumTrials/4)*4;
        
        trialNumi = zeros(1,NumTrials);
        statei = zeros(1,NumTrials);
        Choicei = zeros(1,NumTrials);
        Laseri = randsample([0 1],NumTrials,true,[0.8 0.2]);
        
        
        for t = 1:1:NumTrials 
            trialNumi(t) = t;
            if t <= N1
                % ENGAGED STATE
                %disp("Enagaged")
                statei(t) = 1;
                
                if Targeti(t) == 3 
                    if Cuei(t) < 0.35 % Correct Trial
                        Choicei(t) = 3;
                    else
                        Choicei(t) = 4;
                    end
                end
                if Targeti(t) == 4
        
                    if Cuei(t) > 0.65 % Correct Trial
                        Choicei(t) = 4;
                    else
                        Choicei(t) = 3;
                    end
        
                end
            elseif t <= N2
                
                %disp("Bias")
                Choicei(t) = 3; % Left Bias
                if rand > 0.95
                    Choicei(t) = 4;
                end
    
                %{
                if t < N1 + (N2-N1)/2 
                    statei(t) = 2;
                    if Cuei(t) < 1.0
                        Choicei(t) = 3;
                    else
                        Choicei(t) = 4;
                    end
                end
        
                if  t <= N2 && t >= N1 + (N2-N1)/2 
                    statei(t) = 3;
                    if Cuei(t) > 0.00
                        Choicei(t) = 4; 
                    else
                        Choicei(t) = 3;
                    end
        
                end
                %}
            elseif t <= N3
                statei(t) = 4;
                %disp("PrevChoice: " + num2str(Choicei(t-1)))
                Choicei(t) = Choicei(t-1);
                if Cuei(t) < 0.2
                    %disp("set Left Prev Choice")
                    Choicei(t) = 3;
                end
        
                if Cuei(t) > 0.8
                    %disp("set Right Prev Choice")
                    Choicei(t) = 4;
                end
        
            else
                statei(t) = 5;
                PreviousCorrect = Targeti(t-1) == Choicei(t-1);
                 if PreviousCorrect
                    % Win Stay
                    if Choicei(t-1) == 3
                        %disp("WSLS Win Left")
                        Choicei(t) = 3;
                    else
                        %disp("WSLS Win Right")
                        Choicei(t) = 4;
                    end
                else
                    % Lose Switch
                    if Choicei(t-1) == 3
                        %disp("WSLS Lose Right")
                        Choicei(t) = 4;
                    else
                        %disp("WSLS Lose Left")
                        Choicei(t) = 3;
                    end
                
                end
            end
        end
        
        outputi = Choicei;
        outputi(outputi == 4) = 1;
        outputi(outputi == 3) = 0;
        
        stimOutputi = (Cuei  - 0.5).*2;
    
        biasi = repelem(1,length(Targeti)); % Final Value
    
        prevChoicei = 2*[0.5, outputi(1:end-1)] - 1;
        
        rewardi = double(Targeti(1:end) == Choicei(1:end));
        rewardi(rewardi==0) = -1;
        prevRewardi = [0, rewardi(1:1:end-1)];
        wSLSi = prevRewardi.*prevChoicei;
    
        suboutput = [suboutput,outputi];
        substimOutput = [substimOutput,stimOutputi];
        subbias = [subbias,biasi];
        subprevChoice = [subprevChoice,prevChoicei];
        subreward = [subreward,rewardi];
        subwSLS = [subwSLS,wSLSi];
        substate = [substate,statei];
    
        subTarget = [subTarget,Targeti];
        subtrialNum = [subtrialNum,trialNumi];
        
        sessioni = repmat(sess,[1 N]);
        subsession = [subsession,sessioni];
        subLaserOni = [subLaserOni,Laseri];
    
    end
    
   

    output = [output,suboutput];
    stimOutput = [stimOutput,substimOutput];
    bias = [bias,subbias];
    prevChoice = [prevChoice,subprevChoice];
    reward = [reward,subreward];
    wSLS = [wSLS,subwSLS];
    state = [state,substate];

    Target = [Target,subTarget];
    trialNum = [trialNum,subtrialNum];
    
    session = [session,subsession];
    LaserOni = [LaserOni,subLaserOni];

    SubjectNamei = repmat("FakeMouse" + subi,length(subTarget),1);
    SubjectName = [SubjectName; SubjectNamei];
end
%size(SubjectName)
TargetSide = Target';
%size(TargetSide)
Output = output';
Stim = stimOutput';
Bias = bias';
PrevChoice = prevChoice';
WSLS = wSLS';
Reward = reward';
States = state';

Trial = trialNum';
Session = session';
TrialTag = ( 1:1:length(TargetSide) )' ;
Intervention = LaserOni';

data = table(Output,TargetSide,Stim,Bias,PrevChoice,WSLS,States,Reward,SubjectName,Trial,TrialTag,Session,Intervention);
end