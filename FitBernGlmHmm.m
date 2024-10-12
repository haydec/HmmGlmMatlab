function FullFitData = FitBernGlmHmm(DataInput, settingsBernGlmHmm, FitType, FullFitData)
    
    FitName = settingsBernGlmHmm.FitName;
    assert(isstring(FitName),"First Input Argument (FitName) Should Be A String")
    assert(isstruct(DataInput),"Second Input Argument (DataInput) Should Be A Struct. Please use 'PrepareHMMdata' Function")
    assert(isstruct(settingsBernGlmHmm),"Third Input Argument (ParamStruct) Should Be A Struct")
        
    if isempty(FitType)
        FitType = "fullfit"; 
    end

    if exist("FullFitData","var")
        assert( isstruct(FullFitData),"Fifth Input Argument (FullfitData) Should Be A Struct" ) 
    end
  
    FitType = convertCharsToStrings(FitType);
    FitType = lower(FitType);
    assert(any(strcmpi(FitType,["one","global","subject","fullfit"])),"choose 1 of the 4 options: one, global, subject, fullfit")
    

    %----------------------------------------------------------------------------------------------------------------
    % -------------------------------- One State MLE Fit 
    %----------------------------------------------------------------------------------------------------------------
    if FitType == "one" || (FitType == "global" && ~exist("FullFitData","var") ) || ( FitType == "subject" && ~exist("FullFitData","var") ) ||FitType == "fullfit" 
        disp("---------------------------------")
        disp("-Maximum Likelihood (One State)-")
        disp("---------------------------------")
        FullFitData = fitOneStateMLE(DataInput,settingsBernGlmHmm);
        
    end
 
    %----------------------------------------------------------------------------------------------------------------
    % -------------------------------- Multi State MAP Fit 
    %----------------------------------------------------------------------------------------------------------------
    if FitType == "global" || (FitType == "subject" && ~isfield(FullFitData,"GlobalMap")) || FitType == "fullfit"
        
        disp("---------------------------------")
        disp("-Maximum Aposterior (Multi State)-")
        disp("---------------------------------")
        FullFitData = fitHmmGlmFitMultiState_MAP(DataInput,settingsBernGlmHmm,FullFitData);
        
    end

    %----------------------------------------------------------------------------------------------------------------
    % -------------------------------- Multi State MAP Fit Per Subject
    %----------------------------------------------------------------------------------------------------------------
    if FitType == "subject"  || FitType == "fullfit"
        
        disp("------------------------------------------")
        disp("-Maximum Aposterior (Multi State Subject)-")
        disp("------------------------------------------")
        FullFitData = fitHmmGlmFitMultiState_MAP_Subject(DataInput,settingsBernGlmHmm, FullFitData);

    end
    
    FullFitData.DataInput = DataInput;
    FullFitData.settingsBernGlmHmm = settingsBernGlmHmm;

    FolderHmm = fullfile(pwd,FitName);
    if ~exist(FolderHmm ,"dir")
        mkdir(FolderHmm)
    end

    save(FolderHmm,"FullFitData")
end


function FullFitData = fitOneStateMLE(DataInput,settingsBernGlmHmm)
    
    SubjectName = "Global";

    mle_model =  HmmGlmFitOneState_MLE(DataInput,settingsBernGlmHmm,SubjectName);
    FullFitData.OneStateMle.rawdata = mle_model;

    FitName = settingsBernGlmHmm.FitName;
    Input = DataInput.TrainData.Global.InputFoldCell;
    Output = DataInput.TrainData.Global.OutputFoldCell;
    [WkSort,TransMatSort, EstStateViterbiSort, GammaSort, EpsilonSort] = ...
        SortModels(mle_model,Input,Output,SubjectName,FitName);
    DataSort.WkSort = WkSort;
    DataSort.TransMatSort = TransMatSort;
    DataSort.EstStateViterbiSort = EstStateViterbiSort;
    DataSort.GammaSort = GammaSort;
    DataSort.EpsilonSort = EpsilonSort;

    FullFitData.OneStateMle.sortdata = DataSort;
end


function FullFitData = fitHmmGlmFitMultiState_MAP(DataInput,settingsBernGlmHmm,FullFitData)
    
    SubjectName = "Global";
    FitName = settingsBernGlmHmm.FitName;

    map_model = HmmGlmFitMultiState_MAP(DataInput,settingsBernGlmHmm,SubjectName,FullFitData);
    FullFitData.GlobalMap.rawdata = map_model;

    [BestState, BestPriorAlpha, BestPriorSigma,I1,I2,I3] =...
        BestSinglSubjectModel(map_model,settingsBernGlmHmm);

    MultiDataSort.BestModel.BestState = BestState;
    MultiDataSort.BestModel.BestPriorAlpha = BestPriorAlpha;
    MultiDataSort.BestModel.BestPriorSigma = BestPriorSigma;
    MultiDataSort.BestModel.I1 = I1;
    MultiDataSort.BestModel.I2 = I2;
    MultiDataSort.BestModel.I3 = I3;

    
    
    Input = DataInput.TrainData.Global.InputFoldCell;
    Output = DataInput.TrainData.Global.OutputFoldCell;
    [WkSort, TransMatSort, EstStateViterbiSort, GammaSort, EpsilonSort] = ...
        SortModelsMultiState(map_model,Input,Output,I2,I3,settingsBernGlmHmm,SubjectName,FitName);
    MultiDataSort.WkSort = WkSort;
    MultiDataSort.TransMatSort = TransMatSort;
    MultiDataSort.EstStateViterbiSort = EstStateViterbiSort;
    MultiDataSort.GammaSort = GammaSort;
    MultiDataSort.EpsilonSort = EpsilonSort;        
    FullFitData.GlobalMap.sortdata = MultiDataSort;
end


function FullFitData = fitHmmGlmFitMultiState_MAP_Subject(DataInput,settingsBernGlmHmm, FullFitData)
    
    FitName = settingsBernGlmHmm.FitName;

    subject_map_model =  HmmGlmFitMultiState_MAP_Subject(DataInput,settingsBernGlmHmm, FullFitData);
    FullFitData.SubjectMap.rawdata =  subject_map_model;

    BestSubject = BestSinglSubjectModel_Subject(subject_map_model,settingsBernGlmHmm);
    FullFitData.SubjectMap.BestSubject = BestSubject;  
    
    SubjectDataSort = SortModelsMultiState_Subject(subject_map_model,BestSubject,DataInput,settingsBernGlmHmm,FitName);      
    FullFitData.SubjectMap.sortdata = SubjectDataSort;

end


function Estimate_glmhmmCell =  HmmGlmFitOneState_MLE(DataInput,settingsBernGlmHmm,SubjectName)

    
    %------------Load Global Data-----------------  
    try
        InputFoldCell = DataInput.TrainData.Global.InputFoldCell;
        OutputFoldCell = DataInput.TrainData.Global.OutputFoldCell;
    catch
        disp("Global Train Data does not exist")
        return
    end
      
    NumberOfFolds = settingsBernGlmHmm.NumberOfFolds;
    
    Estimate_glmhmmCell = cell(1,NumberOfFolds);

    obs_dim = size(OutputFoldCell{1},2);
    num_categories = numel(unique(OutputFoldCell{1}));
    BinaryCheck = all(unique(OutputFoldCell{1}) == [0 ;1]);
    
    assert(obs_dim == 1,"Observation Dimension (Columns in Output) must be equal to one")
    assert(num_categories == 2,"Only two catgeories allowed in the output Vector (0 and 1)")
    assert(BinaryCheck == true,"Category numbers must be 0 and 1")
    
    %------Settings for MLE One State Fit---------
    ParamInput.num_states = 1;
    ParamInput.prior_sigma = 100;
    ParamInput.prior_alpha = 1;
    ParamInput.prior_mean = 0;

   
    %--------------Load Paramters-----------------
    input_dim = settingsBernGlmHmm.InputDim;
    categories = 2; % Hard Coded
    AinitType = settingsBernGlmHmm.AinitType;
    num_states = ParamInput.num_states;
    Wmle = [];

    %----------Fit Settings-----------------------
    ParamInput.MaxIter = settingsBernGlmHmm.MaxIter;
    ParamInput.Tolerance = settingsBernGlmHmm.Tolerance;

    %----------Parallel Use-----------------------
    ParamInput.ParallelUse = settingsBernGlmHmm.DataParallelism;

    %-------Setup Progress GUI-------------------
    TotalFits = NumberOfFolds;
    counter = 0;
    wh = waitbar(counter*(1/TotalFits),SubjectName+": Please wait...");

    for f = 1:1:NumberOfFolds

        %-------------Initialize Variables------------
        [cinit,Ainit,Winit] = InitializeHmmGlmParameters(num_states,input_dim,categories,AinitType,Wmle);

        ParamInput.Winit = Winit;
        ParamInput.cinit = cinit;
        ParamInput.Ainit = Ainit;
    
        input = InputFoldCell{f};
        output = OutputFoldCell{f};    
        DataOutput = GlmHMM(output,input,ParamInput);
        Estimate_glmhmmCell{f} = DataOutput;
        counter = counter + 1;
        waitbar(counter*(1/TotalFits),wh,SubjectName+": Fits Left: " + num2str(TotalFits-counter) );
    end
    waitbar(1,wh,SubjectName+": All Fits Complete");


end


function Estimate_glmhmmCell =  HmmGlmFitMultiState_MAP(DataInput,settingsBernGlmHmm,SubjectName,FullFitData)


    %------------Load Global Data-----------------  
    try
        InputFoldCell = DataInput.TrainData.Global.InputFoldCell;
        OutputFoldCell = DataInput.TrainData.Global.OutputFoldCell;
    catch
        disp("Global Train Data does not exist")
        return
    end

    %-----------MLE Weights for Initialization---------------
    WinitMean = mean(FullFitData.OneStateMle.sortdata.WkSort,3);
    
    obs_dim = size(OutputFoldCell{1},2);
    num_categories = numel(unique(OutputFoldCell{1}));
    BinaryCheck = all(unique(OutputFoldCell{1}) == [0 ;1]);
    
    assert(obs_dim == 1,"Observation Dimension (Columns in Output) must be equal to one")
    assert(num_categories == 2,"Only two catgeories allowed in the output Vector (0 and 1)")
    assert(BinaryCheck == true,"Category numbers must be 0 and 1")    

    %---------------Prepare GlmHmm Fits-----------------------
    States = settingsBernGlmHmm.States;
    NumberOfModels = settingsBernGlmHmm.NumberOfModels;

    Prior_Alphas = settingsBernGlmHmm.prior_alphas;
    NumberOfAlphas = settingsBernGlmHmm.NumberOfAlphas;

    Prior_Sigmas = settingsBernGlmHmm.prior_sigmas;
    NumberOfSigmas = settingsBernGlmHmm.NumberOfSigmas;

    NumberOfFolds = settingsBernGlmHmm.NumberOfFolds;

    Estimate_glmhmmCell = cell(NumberOfModels,NumberOfAlphas,NumberOfSigmas,NumberOfFolds);
    TotalSubjectFits = NumberOfModels*NumberOfAlphas*NumberOfSigmas*NumberOfFolds;
    
    %-------Setup Progress GUI-------------------
    TSCi = 1/TotalSubjectFits;
    counter = 0;
    wh = waitbar(counter*TSCi,SubjectName+": Please wait...");

    for k = 1:1:NumberOfModels
        for a = 1:1:NumberOfAlphas
            for s = 1:1:NumberOfSigmas

                %------Settings for MAP K State Fit---------
                ParamInput.num_states = States(k);
                ParamInput.prior_alpha = Prior_Alphas(a);
                ParamInput.prior_sigma = Prior_Sigmas(s);
                ParamInput.prior_mean = 0;

                %----------Fit Settings-----------------------
                ParamInput.MaxIter = settingsBernGlmHmm.MaxIter;
                ParamInput.Tolerance = settingsBernGlmHmm.Tolerance;

                %----------Parallel Use-----------------------
                ParamInput.ParallelUse = settingsBernGlmHmm.DataParallelism;

                for f = 1:1:NumberOfFolds

                    %--------------Load Paramters-----------------
                    input_dim = settingsBernGlmHmm.InputDim;
                    categories = 2;
                    AinitType = settingsBernGlmHmm.AinitType;
                    num_states = ParamInput.num_states;


                    %-------------Initialize Variables------------


                    [cinit,Ainit,Winit] = InitializeHmmGlmParameters(num_states,input_dim,categories,AinitType,WinitMean);

                    ParamInput.Winit = Winit;
                    ParamInput.cinit = cinit;
                    ParamInput.Ainit = Ainit;
                    
                     
                    input = InputFoldCell{f};
                    output = OutputFoldCell{f};    
                    
                    DataOutput = GlmHMM(output,input,ParamInput);
                    Estimate_glmhmmCell{k,a,s,f} = DataOutput;
                    
                    counter = counter + 1;
                    waitbar(counter*TSCi,wh,SubjectName+": Fits Left: " + num2str(TotalSubjectFits-counter) );

                end
            end
        end
    end
    waitbar(1,wh,SubjectName+": All Fits Complete");

end


function Estimate_glmhmmSubMapStruct =  HmmGlmFitMultiState_MAP_Subject( DataInput, settingsBernGlmHmm, FullFitData)
    
    SubjectName = settingsBernGlmHmm.SubjectNames;
    
    if length(SubjectName) > 1
        Estimate_glmhmmSubMapCell = cell(1,length(SubjectName));
        % Parfor does not allow stucutures to be updated

        for sub = 1:length(SubjectName)
            
            

            % - Get Data for each subject to fit to            
            try
                InputFoldCell = DataInput.TrainData.(SubjectName(sub)).InputFoldCell;
                OutputFoldCell = DataInput.TrainData.(SubjectName(sub)).OutputFoldCell;
            catch
                disp("Subject Train Data does not exist: " + SubjectName(sub))
                continue
            end 
            
            %-----------MLE Weights for Initialization---------------
            WinitMean = cellfun(@(x)mean(x,3),FullFitData.GlobalMap.sortdata.WkSort,"UniformOutput",false);
            
            obs_dim = size(OutputFoldCell{1},2);
            num_categories = numel(unique(OutputFoldCell{1}));
            BinaryCheck = all(unique(OutputFoldCell{1}) == [0 ;1]);
            
            assert(obs_dim == 1,"Observation Dimension (Columns in Output) must be equal to one")
            assert(num_categories == 2,"Only two catgeories allowed in the output Vector (0 and 1)")
            assert(BinaryCheck == true,"Category numbers must be 0 and 1")    
        
            %---------------Prepare GlmHmm Fits-----------------------
            States = settingsBernGlmHmm.States;
            NumberOfModels = settingsBernGlmHmm.NumberOfModels;
        
            Prior_Alphas = settingsBernGlmHmm.prior_alphas;
            NumberOfAlphas = settingsBernGlmHmm.NumberOfAlphas;
        
            Prior_Sigmas = settingsBernGlmHmm.prior_sigmas;
            NumberOfSigmas = settingsBernGlmHmm.NumberOfSigmas;
        
            NumberOfFolds = settingsBernGlmHmm.NumberOfFolds;
        
            Estimate_glmhmmCell = cell(NumberOfModels,NumberOfAlphas,NumberOfSigmas,NumberOfFolds);
            TotalSubjectFits = NumberOfModels*NumberOfAlphas*NumberOfSigmas*NumberOfFolds;
            
            %-------Setup Progress GUI-------------------
            TSCi = 1/TotalSubjectFits;
            counter = 0;
            wh = waitbar(counter*TSCi,SubjectName(sub) + ": Please wait...");
        
            for k = 1:1:NumberOfModels
                for a = 1:1:NumberOfAlphas
                    for s = 1:1:NumberOfSigmas
        
                        %------Settings for MAP K State Fit---------
                        ParamInput.num_states = States(k);
                        ParamInput.prior_alpha = Prior_Alphas(a);
                        ParamInput.prior_sigma = Prior_Sigmas(s);
                        ParamInput.prior_mean = 0;
        
                        %----------Fit Settings-----------------------
                        ParamInput.MaxIter = settingsBernGlmHmm.MaxIter;
                        ParamInput.Tolerance = settingsBernGlmHmm.Tolerance;

                        %----------Parallel Use-----------------------
                        ParamInput.ParallelUse = settingsBernGlmHmm.DataParallelism;
        
                        for f = 1:1:NumberOfFolds
        
                            %--------------Load Paramters-----------------
                            input_dim = settingsBernGlmHmm.InputDim;
                            categories = 2;
                            AinitType = settingsBernGlmHmm.AinitType;
                            num_states = ParamInput.num_states;
        
        
                            %-------------Initialize Variables------------
        
        
                            [cinit,Ainit,Winit] = InitializeHmmGlmParameters(num_states,input_dim,categories,AinitType,WinitMean{k});
        
                            ParamInput.Winit = Winit;
                            ParamInput.cinit = cinit;
                            ParamInput.Ainit = Ainit;
                            
                             
                            input = InputFoldCell{f};
                            output = OutputFoldCell{f};    
                            
                            DataOutput = GlmHMM(output,input,ParamInput);
                            Estimate_glmhmmCell{k,a,s,f} = DataOutput;
                            
                            counter = counter + 1;
                            waitbar(counter*TSCi,wh,SubjectName(sub) +": Fits Left: " + num2str(TotalSubjectFits-counter) );
        
                        end
                    end
                end
            end
                      
            Estimate_glmhmmSubMapCell{sub} = Estimate_glmhmmCell; % Cell because parfor of implementation issues

        end
        waitbar(1,wh,SubjectName(sub)+": All Fits Complete")
        
        %
        for sub = 1:1:length(SubjectName)
            Estimate_glmhmmSubMapStruct.(SubjectName(sub)) = Estimate_glmhmmSubMapCell{sub};% Cell because parfor of implementation issues
        end

    else

        disp("Only One Subject. No Subject Fit Required")

    end

    

end


function DataOutput = GlmHMM(output,input,DataInput)
    
    ParallelUse = DataInput.ParallelUse;
    
    if isstruct(DataInput)
        % -----------Initialize Parameters------------
        prior_mean = DataInput.prior_mean;
        prior_sigma = DataInput.prior_sigma;
        prior_alpha = DataInput.prior_alpha; 
        num_states = DataInput.num_states;
        MaxIter = DataInput.MaxIter;
        Tolerance = DataInput.Tolerance;
        
        %-------------Initialize Variables------------
        Atrans = DataInput.Ainit;
        c = DataInput.cinit;
        weights = DataInput.Winit;
        GlmProb = BernGlm(weights,input,output);
    elseif iscell(DataInput)


    else
        error("Incorrect Data Type")
        
    end
       
    % -------Begin Fitting-----------
    for iter = 1:1:MaxIter
        
        % -----------Expectation Step------------
        [Alpha, negLogLike(iter),cs] = forward_algorithm_Vec(output,Atrans,c,GlmProb); % Forward Pass
        
        [Beta] = backward_algorithm_Vec(output,Atrans,GlmProb,cs); % Back Pass
        
        % -----------Expectated State Estimation------------
        [Gama,StateMax] = ComputeGama(Alpha,Beta); % Calculate State Probability 
    
        [Xis,XisN,XisD] = ComputeEpsilonVec(output,Atrans,Alpha,Beta,GlmProb,cs); % Calculate joint state probability
        
        % -----------Maximization Step--------------------
        XisNa = (prior_alpha-1) + XisN;
        XisDa = num_states*(prior_alpha-1) + XisD;    
        AtransNew = XisNa./XisDa; % Calculate Transition Matrix
    
        cNew = Gama(1,:)/sum(Gama(1,:)); % Calculate Initial Distribution
    
        wMaxF = @(w) -(sum(Gama.*logBernGlm(w,input,output),"all") + sum(priorLogMVN(w,prior_mean,prior_sigma),"all") ) ;
        options = optimoptions("fminunc",'Algorithm',"quasi-newton","Display","none","UseParallel",ParallelUse);
        [wMax,f(iter)] = fminunc(wMaxF,weights,options); % Calculate Weights
        
        % ------- Update Values-----------
        Atrans = AtransNew;
        c = cNew;
        weights = wMax;
        GlmProb = BernGlm(weights,input,output);
        
        % -------Check Convergence-----------
        if iter > 5
            tol = all(abs(negLogLike(end-4:end) - negLogLike(end)) < Tolerance);
            if tol == true
                break
            end
        end
        
    end

    ValidDist = CheckStochasity(Atrans,Gama,Xis);
    if ValidDist == false
        
        disp("Distirbutions may not be valid")
    
    end
    
    
    
    DataOutput.Atrans = Atrans; % Transition Matrix
    DataOutput.weights = weights; % Weights from Bernoulli Glm
    DataOutput.initState =  c; % Initial State Distribution
    DataOutput.negLogLike  = negLogLike; % Negative Logliklihood from Forward Pass
    DataOutput.ExpectedState = Gama; %  State Probability at each t-step
    DataOutput.ExpectedState_Max = StateMax;% Max State Probability at each t-step
    DataOutput.JointExpectedState = Xis; % Joint State Probability at each t-step
    DataOutput.GlmProb = GlmProb; % GLM Probabilities
    DataOutput.K = size(Atrans,1); % Number Of States
    
end



function GlmProb = BernGlm(Winit,inpts,y)
    
    K = size(Winit,1);
    N = size(y,1);
    
    b = Winit';
    Logit0 = 1./(1 + exp(-inpts*b));
    Logit1 = 1 - Logit0;
    
    GlmProb = zeros(N,K,2);

    GlmProb(:,:,1) = Logit0;
    GlmProb(:,:,2) = Logit1;

    %p = repmat((y==0),1,size(b,2)).*Logit0 + repmat( (y==1),1,size(b,2)).*Logit1;
end


function [ll] = logBernGlm(Winit,inpts,y)
    
    b = Winit';
    Logit0 = -log1p(exp(-inpts*b));
    Logit1 = -inpts*b -log1p(exp(-inpts*b));
    ll = repmat((y==0),1,size(b,2)).*Logit0 + repmat( (y==1),1,size(b,2)).*Logit1;
end


function [alpha, negLogLike, cs] = forward_algorithm_Vec(y,Atrans,c,GlmProb)
    
    % Compute Forward Probability
    % Reference Bishop Chapter 13 page 628, Equation 13.59
    % https://www.microsoft.com/en-us/research/publication/pattern-recognition-machine-learning/

    K = size(Atrans,1);
    N = size(y,1);
    
    alpha = zeros(N,K);
    alpha_prior = zeros(N,K);
    cs = zeros(N,1);
    
    % Initialization
    p = GlmProb(1,:,y(1)+1).*c;
    cs(1) = sum(p);
    alpha(1,:) = p/cs(1);
    alpha_prior(1,:) = 1/K;
    
    % Recursion
    for t = 2:1:size(alpha,1)
    
        alpha_prior(t,:) = alpha(t-1,:)*Atrans; 
        p = GlmProb(t,:,y(t)+1).*alpha_prior(t,:);
        cs(t) = sum(p);
        alpha(t,:) = p/cs(t);
    
    end
    
    negLogLike = sum(log(cs));

end


function [beta] = backward_algorithm_Vec(y,Atrans,GlmProb,cs)

    % Compute Backward Probability
    % Reference Bishop Chapter 13 page 628, Equation 13.62
    % https://www.microsoft.com/en-us/research/publication/pattern-recognition-machine-learning/
    
    K = size(Atrans,1);
    N = size(y,1);
    
    beta = zeros(N,K);
    
    % Initialization
    beta(end,:) = 1;

 
    % Recursion
    for t = N-1:-1:1
    
        beta_prior = beta(t+1,:).*GlmProb(t+1, :, y(t+1)+1); 
        p = Atrans*beta_prior';
        beta(t,:) = p/cs(t+1);
    
    end
        
end

function [Gama,StateMax] = ComputeGama(alpha,beta)
        
    % Compute State Probability 
    % Reference Bishop Chapter 13 page 628, Equation 13.64
    % https://www.microsoft.com/en-us/research/publication/pattern-recognition-machine-learning/
    Gama = alpha.*beta;
    [~ ,StateMax] = max(Gama,[],2); 

end

function [Xis,XisN,XisD] = ComputeEpsilonVec(y,Atrans,Alpha,Beta,GlmProb,cs)

    % Compute Joint State Probability 
    % Reference Bishop Chapter 13 page 628, Equation 13.65
    % https://www.microsoft.com/en-us/research/publication/pattern-recognition-machine-learning/
    N = length(y);
    K = size(Atrans, 1);
    Xis = zeros(N - 1, K, K);
    
    for t = 1:1:N-1
        BetaPhi = GlmProb(t+1,:,y(t+1)+1).* Beta(t+1,:); % Current State
        AlphaReshaped = reshape(Alpha(t,:),K,1); % Previous State
        Xis(t,:,:) = (BetaPhi.*AlphaReshaped.*Atrans)/cs(t+1);
    end
    
    XisN = reshape(sum(Xis, 1), K, K);
    XisD = reshape(sum(Xis, [1, 3]), K, 1);
end



function ValidDist = CheckStochasity(Atrans,Gama,Xis)
    
    K = size(Atrans,1);
     
    stochasityTol = ismembertol(sum(Atrans,2)', ones(1,size(Atrans,1)), 0.05, 'ByRows', true);% Check that Stochasity cosntraint was maintained
    Gama1 = Gama(1,:);
    XN = reshape(Xis(1,:,:),K,K);
    Gama2 = sum(XN,2)'; % Sum of Previous to states to currents sates should equal probability of current states
    iswithintol = ismembertol(Gama1, Gama2, 0.05, 'ByRows', true);

    ValidDist = all([stochasityTol,iswithintol] == 1);

end



% 

function p = priorLogMVN(w,mu,sig)
    % Prior for multiple states
    States = size(w,1);
    p = zeros(1,States);
    for s = 1:1:States 
    
        p(s) = log_multivariate_normal_pdf(w(s,:), mu, sig);
    
    end
end

function log_pdf = log_multivariate_normal_pdf(x, mu, sig)
    % sig: Standard Deviation
    % mu: mean
    k = size(x,2);  % Dimensionality of the distribution
    Sigma = eye(k)*sig^2; % Covariance Matrix
    Mu = repmat(mu,1,k); % Mean Vector

    c = -0.5 * (k * log(2*pi) + log(det(Sigma)));  % Constant term
    
    % Computing the exponent part
    exponent = -0.5 * (x - Mu) * pinv(Sigma) * (x - Mu)';
    
    log_pdf = c + exponent;  % Logarithm of the PDF
end   






function [cNew,ATrans,Weights] = InitializeHmmGlmParameters(num_states,input_dim,categories,AinitType,Winit)
    
    initOptions = ["Dir","uni","Sticky"];
    assert(contains(AinitType,initOptions),"Options for Transition Matrix Initiation are 'Dir', 'Uni' and 'Sticky'")

    if AinitType == "Dir"
        %disp("Transition Matrix initialized with Dirichlet Distribution") 
        a = 1; % Concentration Parameter fixed to 1
        alphaVec = repmat(a,1,num_states);
        ATrans = drchrnd(alphaVec,num_states);
    end
    if AinitType == "Uni"
        %disp("Transition Matrix initialized with Uniform Distribution")  
        ATrans = ones(num_states,num_states)/num_states; % Uniform Init
    end
    if AinitType == "Sticky"
        %disp("Transition Matrix initialized with 'Sticky' Distribution")
        A = eye(num_states)*0.95 + rand(num_states,num_states)* 0.05;
        Anorm = sum(A,2);
        ATrans = A./Anorm;
    end
    
    %disp("Prior State Vector initialized with Uniform Distribution")
    cNew = ones(1,num_states)/num_states;
    
    %disp("Weights Initialized for Bernoulli GLM")
    if isempty(Winit)
        Weights = randn(1,input_dim,categories - 1 );
    elseif size(Winit,1) == 1
        sigma = 0.2; % std
        mu = 0; % Mean
        GaussNoise = sigma.*randn(num_states,input_dim) + mu;
        WkInit =  repmat(Winit,num_states,1);
        Weights = WkInit + GaussNoise;
    else
        sigma = 0.2; % std
        mu = 0; % Mean
        GaussNoise = sigma.*randn(num_states,input_dim) + mu;
        WkInit =  Winit;
        Weights = WkInit + GaussNoise;

    end
    
    function r = drchrnd(a,n)
        % take a sample from a dirichlet distribution
        % https://cxwangyi.wordpress.com/2009/03/18/to-generate-random-numbers-from-a-dirichlet-distribution/
        p = length(a);
        r = gamrnd(repmat(a,n,1),1,n,p);
        r = r ./ repmat(sum(r,2),1,p);
    end

end




function [WkSort,TransMatSort,EstStateViterbiSort,EstExpectatedStateSort,EstJointExpectatedStateSort] = SortModels(Estimate_glmhmmCell,InputFoldCell,OutputFoldCell,SubjectName,FitName)

    
    

    Folds = length(Estimate_glmhmmCell);
    num_states = size(Estimate_glmhmmCell{1}.weights,1);
    input_dim = size(Estimate_glmhmmCell{1}.weights,2);

    Wk = zeros(num_states,input_dim,Folds);
    distMat = zeros(num_states,num_states,Folds);
    M = zeros(num_states,2,Folds);
    WkSort = zeros(num_states,input_dim,Folds);
    
    for f = 1:1:Folds
        Wk(:,:,f) = Estimate_glmhmmCell{f}.weights;
    end
    
    for f = 1:1:Folds
        distMat(:,:,f) = pdist2(Wk(:,:,1),Wk(:,:,f),"cityblock");
        M(:,:,f)  = sortrows(matchpairs(distMat(:,:,f),10000));
        WkSort(:,:,f) = Wk(M(:,2,f),:,f);
    end

    WeightsFolder = fullfile(pwd,FitName,"WeightsData",SubjectName);    
    if ~exist("WeightsFolder","dir")
        mkdir(WeightsFolder)
    end

    WeightsHmm = fullfile(WeightsFolder,"BestHyperWeights_State" + num2str(num_states) + ".mat");
    save(WeightsHmm,'WkSort');

    
    TransMatSort = zeros(num_states,num_states,Folds);
    for f = 1:1:Folds
        TransMat = double(Estimate_glmhmmCell{f}.Atrans);
        TransMatSort(:,:,f) = TransMat(M(:,2,f),M(:,2,f));
    end


    TransMatFolder = fullfile(pwd,FitName,"TransMatData",SubjectName);    
    if ~exist("TransMatFolder","dir")
        mkdir(TransMatFolder)
    end

    TransMatHmm = fullfile(TransMatFolder,"BestHyperTransMat_State" + num2str(num_states) + ".mat");
    save(TransMatHmm,'TransMatSort');
    
    EstStateViterbi = cell(1,Folds);
    for f = 1:1:Folds
    
        Input = InputFoldCell{f};
        Output = OutputFoldCell{f};
        A = Estimate_glmhmmCell{f}.Atrans;
        W = Estimate_glmhmmCell{f}.weights;
        GlmProb = BernGlm(W,Input,Output);
        c = Estimate_glmhmmCell{f}.initState;
        EstStateViterbi{f} = viterbi_algorithm(A,GlmProb,Output,c) ;
    end
    
    IDXSwap = cell(Folds,num_states);
    EstStateViterbiSort = EstStateViterbi;
    for f = 1:1:Folds
    
        for s = 1:1:num_states
            IDXSwap{f,s} = ismember( EstStateViterbi{f}, M(s,2,f) );
        end
    
        for s = 1:1:num_states
            EstStateViterbiSort{f}( IDXSwap{f,s}  ) = M(s,1,f);
        end
    
    end

    
    ViterbiFolder = fullfile(pwd,FitName,"ViterbiData",SubjectName);    
    if ~exist("ViterbiFolder","dir")
        mkdir(ViterbiFolder)
    end

    ViterbiHmm = fullfile(ViterbiFolder, "BestHyperViterbi_State" + num2str(num_states) + ".mat");
    save(ViterbiHmm,'EstStateViterbiSort');

    EstExpectatedStateSort = cell(1,Folds);
    K = double(Estimate_glmhmmCell{f}.K);
    for f = 1:1:Folds
    
        EstExpectatedStateSort{f} = Estimate_glmhmmCell{f}.ExpectedState;
        EstExpectatedStateSort{f}(:,1:1:K) = EstExpectatedStateSort{f}(:,M(:,2,f));

    end

    GammaFolder = fullfile(pwd,FitName,"GammaData",SubjectName);    
    if ~exist("GammaFolder","dir")
        mkdir(GammaFolder)
    end

    GammaHmm = fullfile(GammaFolder, "BestHyperGamma_State" + num2str(num_states) + ".mat");
    save(GammaHmm,'EstExpectatedStateSort');


    EstJointExpectatedStateSort = cell(1,Folds);
    for f = 1:1:Folds
    
        EstJointExpectatedStateSort{f} = Estimate_glmhmmCell{f}.JointExpectedState;      
        EstJointExpectatedStateSort{f}(:,:,:) = EstJointExpectatedStateSort{f}(:,M(:,2,f),M(:,2,f));

    end

    EpsilonFolder = fullfile(pwd,FitName,"EpsilonData",SubjectName);    
    if ~exist("EpsilonFolder","dir")
        mkdir(EpsilonFolder)
    end

    EpsilonHmm = fullfile(EpsilonFolder, "BestHyperEpsilon_State" + num2str(num_states) + ".mat");
    save(EpsilonHmm,'EstJointExpectatedStateSort');

end






function[ WkSortCell,TransMatSortCell,EstStateViterbiSortCell,EstGammaSortCell,EstEpsilonSortCell] = SortModelsMultiState(Estimate_glmhmmMapCell,InputFoldCell,OutputFoldCell,I2,I3,settingsBernGlmHmm,SubjectName,FitName)
    
    NumberOfModels = settingsBernGlmHmm.NumberOfModels;
    WkSortCell = cell(1,NumberOfModels);
    TransMatSortCell = cell(1,NumberOfModels);
    EstStateViterbiSortCell = cell(1,NumberOfModels);
    EstGammaSortCell = cell(1,NumberOfModels);
    EstEpsilonSortCell = cell(1,NumberOfModels);
    
    for k = 1:1:NumberOfModels
        
        Estimate_glmhmmCell = permute(Estimate_glmhmmMapCell(k,I2,I3,:),[4 2 3 1]);
        [WkSortCell{k},TransMatSortCell{k},EstStateViterbiSortCell{k},EstGammaSortCell{k},EstEpsilonSortCell{k}] = SortModels(Estimate_glmhmmCell,InputFoldCell,OutputFoldCell,SubjectName,FitName);
    
    end

end

function SubjectSort = SortModelsMultiState_Subject(Estimate_glmhmmSubMapStruct,BestSubject,DataInput,settingsBernGlmHmm,FitName)
    
    SubjectNames = settingsBernGlmHmm.SubjectNames;
    for sub = 1:1:length(SubjectNames)
        Estimate_glmhmmMapCell = Estimate_glmhmmSubMapStruct.(SubjectNames(sub));
        InputFoldCell = DataInput.TrainData.(SubjectNames(sub)).InputFoldCell;
        OutputFoldCell = DataInput.TrainData.(SubjectNames(sub)).OutputFoldCell;
        I2 = BestSubject.(SubjectNames(sub)).I2;
        I3 = BestSubject.(SubjectNames(sub)).I3;
        SubjectName = SubjectNames(sub);
        [ WkSortCellGmap,TransMatSortCell,EstStateViterbiSortCell,EstGammaSortCell,EstEpsilonSortCell] = SortModelsMultiState(Estimate_glmhmmMapCell,InputFoldCell,OutputFoldCell,I2,I3,settingsBernGlmHmm,SubjectName,FitName);
        SubjectSort.(SubjectNames(sub)).WkSortCellGmap = WkSortCellGmap;
        SubjectSort.(SubjectNames(sub)).TransMatSortCell = TransMatSortCell;
        SubjectSort.(SubjectNames(sub)).EstStateViterbiSortCell= EstStateViterbiSortCell;
        SubjectSort.(SubjectNames(sub)).EstGammaSortCell = EstGammaSortCell;
        SubjectSort.(SubjectNames(sub)).EstEpsilonSortCell = EstEpsilonSortCell;

    end
end




function [BestState,BestPriorAlpha,BestPriorSigma,I1,I2,I3] = BestSinglSubjectModel(Estimate_glmhmmMapCell,settingsBernGlmHmm)
    
    % Determine Which Hyperparameters (alpha and sigma) are the best and
    % use those models
                      
    LogProb = LogProbability(Estimate_glmhmmMapCell,settingsBernGlmHmm);

    meanLogProb = mean(LogProb,4); % Mean LogProb across Folds
    [~,idx] = max(meanLogProb,[],[1,2,3],"linear");
    [I1,I2,I3] = ind2sub(size(meanLogProb),idx);
    BestState = settingsBernGlmHmm.States(I1);
    BestPriorAlpha = settingsBernGlmHmm.prior_alphas(I2);
    BestPriorSigma = settingsBernGlmHmm.prior_sigmas(I3);
end



function BestSubject = BestSinglSubjectModel_Subject(Estimate_glmhmmSubMapStruct,settingsBernGlmHmm)
    
    SubjectName = settingsBernGlmHmm.SubjectNames;
    for sub = 1:1:length(SubjectName)
        Estimate_glmhmmMapCell = Estimate_glmhmmSubMapStruct.(SubjectName(sub));
        [BestState,BestPriorAlpha,BestPriorSigma,I1,I2,I3] = BestSinglSubjectModel(Estimate_glmhmmMapCell,settingsBernGlmHmm);
        BestSubject.(SubjectName(sub)).BestState = BestState;
        BestSubject.(SubjectName(sub)).BestPriorAlpha = BestPriorAlpha;
        BestSubject.(SubjectName(sub)).BestPriorSigma = BestPriorSigma;
        BestSubject.(SubjectName(sub)).I1 = I1;
        BestSubject.(SubjectName(sub)).I2 = I2;
        BestSubject.(SubjectName(sub)).I3 = I3;
    end
end






function LogProb = LogProbability(Estimate_glmhmmMapCell,settingsBernGlmHmm)
    
    NumberOfModels = settingsBernGlmHmm.NumberOfModels;
    NumberOfAlphas = settingsBernGlmHmm.NumberOfAlphas;
    NumberOfSigmas = settingsBernGlmHmm.NumberOfSigmas;
    NumberOfFolds = settingsBernGlmHmm.NumberOfFolds;

    LogProb = zeros(NumberOfModels,NumberOfAlphas,NumberOfSigmas,NumberOfFolds);
    
    for k = 1:1:NumberOfModels
    
        for a = 1:1:NumberOfAlphas
    
            for s = 1:1:NumberOfSigmas
    
                for f = 1:1:NumberOfFolds   
                    alpha = settingsBernGlmHmm.prior_alphas(a);
                    sig = settingsBernGlmHmm.prior_sigmas(s);
                    weights = Estimate_glmhmmMapCell{k,a,s,f}.weights;
                    mean = 0;
                    Atrans = Estimate_glmhmmMapCell{k,a,s,f}.Atrans;
                    c = Estimate_glmhmmMapCell{k,a,s,f}.initState;

                    logPrior = LogPriorGlmHmm(weights,mean,sig,Atrans,c,alpha);
                    negLogLike= Estimate_glmhmmMapCell{k,a,s,f}.negLogLike(end);
                    LogProb(k,a,s,f) = negLogLike + logPrior;
                end
            end
        end
    end
   
end

function q=viterbi_algorithm(A,GlmProb,O,c) 
    % Viterbi Algorithm for discrete hidden Markov Models with ’m’ hidden 
    % states ,  n  observable states and  N  observations . 
    % A - mxm (state transition matrix) 
    % B - mxn (confusion matrix) 
    % O - 1xN (observations vector)
    % c - 1xm (initial probabilities vector) 
   
    [N,m,~]=size(GlmProb);
    %N=length(O); 
    
    %% Initialization 
    delta=zeros(N,m);phi=zeros(N,m); 
    t=1; 
    for k=1:m     
        delta(t,k)=c(k)*GlmProb(t,k,O(t)+1); 
        phi(t,k)=0;
    end 
    
    %% Recursion 

    for t=2:N
        for k=1:m
        
            for l=1:m 
                tmp(l)=delta(t-1,l)*A(l,k)*GlmProb(t,k,O(t)+1); 
            end 
            [delta(t,k),phi(t,k)]=max(tmp); 
        end
    end 
    %% Path finding
    q=zeros(N,1);
    [~,Inx]=max(delta(N,:));
    q(N)=Inx;
    for k=N-1:-1:1
    
        q(k)=phi(k+1,q(k+1)); 
    end
    
end

function logPrior = LogPriorGlmHmm(weights,mean,sig,Atrans,c,alpha)
    pW = priorLogMVN(weights,mean,sig);
    pA = priorLogDirA(Atrans,alpha);
    pC = priorLogDirC(c,alpha);
    logPrior = sum(pW)+sum(pA)+sum(pC);
end

  
function p = dirichletpdf(x, alpha)
    
    if length(alpha) ~= length(x)
        error('The alpha and x vectors must have the same length.');
    end

    % Add a small positive value to x to avoid taking the logarithm of zero
    %epsilon = 1e-10;
    %x = x + epsilon;
    
    % Calculate the normalization constant
    C =  prod(gamma(alpha))/gamma(sum(alpha));
    
    % Calculate the probability density function
    p = C * prod(x.^(alpha-1));

end

function p = priorLogDirA(A,alpha)
    
    States = size(A,1);
    p = zeros(1,States);
    for s = 1:1:States
        x = A(s,:);
        a = repmat(alpha,1,States);
        p(s) = log(dirichletpdf(x, a));
    end

end

function p = priorLogDirC(c,alpha)
    
    States = size(c,2);
    p = zeros(1,States);
    for s = 1:1:States        
        a = repmat(alpha,1,States);
        p(s) = log1p(dirichletpdf(c, a));
    end

end