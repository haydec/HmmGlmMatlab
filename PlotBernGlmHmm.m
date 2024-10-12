function PlotData = PlotBernGlmHmm(DataInput, settingsBernGlmHmm, FullFitData)

    if isfield(FullFitData,"OneStateMle")
        
        Alpha = []; Sigma = []; I1 = [];        
        RelevantInputs = settingsBernGlmHmm.RelevantInputs;
        FitName = settingsBernGlmHmm.FitName;
        SubjectName = "Global";
        
        WkSort = FullFitData.OneStateMle.sortdata.WkSort;
        TransMatSort = FullFitData.OneStateMle.sortdata.TransMatSort;

        [Fw,Wkmean] = PlotWeights(WkSort,RelevantInputs,Alpha,Sigma,FitName,SubjectName,settingsBernGlmHmm,I1);
        [Ft, MeanTransMat] = PlotMatrix(TransMatSort,settingsBernGlmHmm,FitName,SubjectName,I1);
        
        PlotData.OneStateMle.Fw = Fw;
        PlotData.OneStateMle.Wkmean = Wkmean;
        PlotData.OneStateMle.Ft = Ft;
        PlotData.OneStateMle.MeanTransMat = MeanTransMat;
    end

    if isfield(FullFitData,"GlobalMap") 
        
        SessionFoldCell = DataInput.TrainData.Global.SessionFoldCell;
        LaserCell = DataInput.TrainData.Global.LaserCell;
        FitName = settingsBernGlmHmm.FitName;
        SubjectName = "Global";


        BestPriorAlpha = FullFitData.GlobalMap.sortdata.BestModel.BestPriorAlpha;
        BestPriorSigma = FullFitData.GlobalMap.sortdata.BestModel.BestPriorSigma;
        I1 = FullFitData.GlobalMap.sortdata.BestModel.I1;

        WkSortCellGmap = FullFitData.GlobalMap.sortdata.WkSort;
        TransMatSortCell = FullFitData.GlobalMap.sortdata.TransMatSort;
        EstGammaSortCell = FullFitData.GlobalMap.sortdata.GammaSort;
        EstStateViterbiSortCell = FullFitData.GlobalMap.sortdata.EstStateViterbiSort;

        [FkCell,WkmeanCell] = PlotWeightsMultiState(WkSortCellGmap,settingsBernGlmHmm,BestPriorAlpha,BestPriorSigma,FitName,SubjectName,I1);
        [FtCell, MeanTransMatCell] = PlotMatrixMultiState(TransMatSortCell,settingsBernGlmHmm,FitName,SubjectName,I1);
              
        [GammaSessionCell,ViterbiSessionCell] = GammaVertCollect(EstGammaSortCell,EstStateViterbiSortCell,SessionFoldCell,settingsBernGlmHmm,FitName,SubjectName);
        [Fmode,Fmode_T,normRCell,HistCell] = ViterbiPlot(ViterbiSessionCell,settingsBernGlmHmm,FitName,SubjectName,I1);
        FmuE = GammaPlot(GammaSessionCell,LaserCell,settingsBernGlmHmm,FitName,SubjectName,I1);
        
        PlotData.GlobalMap.FkCell = FkCell;
        PlotData.GlobalMap.WkmeanCell = WkmeanCell;
        PlotData.GlobalMap.FtCell = FtCell;
        PlotData.GlobalMap.MeanTransMatCell = MeanTransMatCell;

        PlotData.GlobalMap.GammaSessionCell = GammaSessionCell;
        PlotData.GlobalMap.ViterbiSessionCell = ViterbiSessionCell;

        PlotData.GlobalMap.Fmode = Fmode;
        PlotData.GlobalMap.Fmode_T = Fmode_T;
        PlotData.GlobalMap.normRCell = normRCell;
        PlotData.GlobalMap.HistCell = HistCell;

        PlotData.GlobalMap.FmuE = FmuE;

    end

    if isfield(FullFitData,"SubjectMap")

        SubjectSort = FullFitData.SubjectMap.sortdata;
        BestSubject = FullFitData.SubjectMap.BestSubject;
        FitName = settingsBernGlmHmm.FitName;

        TrainData = DataInput.TrainData;

        WeightSubject = PlotWeightsMultiState_Subject(SubjectSort,settingsBernGlmHmm,BestSubject,FitName);
        SubjectTrans =  PlotMatrixMultiState_Subject(SubjectSort,settingsBernGlmHmm,BestSubject,FitName);

        GamaVertSubject = GammaVertCollect_Subject(SubjectSort,settingsBernGlmHmm,TrainData,FitName);
        ViterbiSubject = ViterbiPlot_Subject(GamaVertSubject,settingsBernGlmHmm,BestSubject,FitName);
        GammaSubject = GamaPlot_Subject(GamaVertSubject,settingsBernGlmHmm,BestSubject,TrainData,FitName);
       
        PlotData.SubjectMap.Weights = WeightSubject;
        PlotData.SubjectMap.TransMat = SubjectTrans;
        PlotData.SubjectMap.GammaVit = GamaVertSubject;
        PlotData.SubjectMap.Viterbi = ViterbiSubject;
        PlotData.SubjectMap.Gamma = GammaSubject;
    end

    if(~isfield(FullFitData,"OneStateMle") || ~isfield(FullFitData,"GlobalMap") || ~isfield(FullFitData,"SubjectMap"))
        warning("No Fit Data Provided")
    end

    disp("Plotting: Done")

end


function [F,Wkmean] = PlotWeights(WkSort,RelevantInputs,Alpha,Sigma,FitName,SubjectName,ParamStruct,I1)
    
    
    num_states = size(WkSort,1);
    input_dim = size(WkSort,2);
    Folds = size(WkSort,3);


    Wkmean = -mean(WkSort,3); % Flip for Interpretation ( Consequence of Minimising the NEGATIVE Log-Likehood in the Function GlmHMM.m )
    Wkstd = std(WkSort,1,3);    
    StatesC = ParamStruct.States;
    if isempty(I1)

        FigureName = "GLM Weights - States: " + num2str(num_states);
        F = figure("Name",FigureName,"Visible","on","WindowState","minimized");

    else
        if num_states == StatesC(I1)
        
            FigureName = "GLM Weights - States: " + num2str(num_states);
            F = figure("Name",FigureName,"Visible","on","WindowState","minimized");
        else
            
            FigureName = "GLM Weights - States: " + num2str(num_states);
            F = figure("Name",FigureName,"Visible","off","WindowState","minimized");
        end
    end
    A = axes(F);
    
    x = 1:1:input_dim;
    xconf = repmat([x x(end:-1:1)],num_states,1) ;         
    WkUp = Wkmean(:,1:1:end) + Wkstd;
    WkDown = Wkmean(:,1:1:end) - Wkstd;
    yconf = [WkUp,WkDown(:,end:-1:1)];
    
    Colours = [0.8500 0.3250 0.0980; % Orange
               0.4660 0.6740 0.1880; % Green
               0.0000 0.4470 0.7410; % Blue
               0.4940 0.1840 0.5560; % Purple
               0.9290 0.6940 0.1250; % Yellow
               0.6350 0.0780 0.1840  % Red
               ];
    
    States = size(Wkmean,1);
    
    StateArray = strings(1,size(Wkmean,1));
    for s = 1:1:size(Wkmean,1)
        StateArray(s) = "State: " + num2str(s);
    end
    
    p = gobjects(1, num_states);
    hold on
    for s = 1:1:size(yconf,1)
        p(s) = plot(A,x,Wkmean(s,:),"Color",Colours(s,:));
        fill(A,xconf(s,:),yconf(s,:),Colours(s,:),"FaceAlpha",0.2)
    end
    
    if isempty(Alpha)
        Alpha = "NA";
    else
        Alpha = num2str(Alpha);
    end

    if isempty(Sigma)
        Sigma = "NA";
    else
        Sigma = num2str(Sigma);
    end
    yline(A,0,"LineStyle","--","Color",[.7 .7 .7])
    title(A,"States: "+ num2str(States) + "|| Prior Alpha: "+ Alpha +"|| Prior Sigma: "+Sigma)  
    legend(A,p,StateArray)
    xticks(A,1:1:input_dim)
    xticklabels(A,RelevantInputs)
    ylabel(A,"GLM Weights")
   
    hold off

    FolderHmm = fullfile(pwd,FitName,"WeightsPlot",SubjectName);
    if ~exist(FolderHmm ,"dir")
        mkdir(FolderHmm)
    end
    Name = "GLM Weights-States-" + num2str(num_states)+"-"+SubjectName;
    FigName = fullfile(FolderHmm,Name + ".fig");
    % Before saving, set the figure to be visible
    
    StatesC = ParamStruct.States;
    num_states = size(WkSort,1);
    set(F, 'Visible', 'on');
    savefig(F,FigName); 
    if num_states ~= StatesC(I1)
        close(F)
        pause(0.01)
        disp("Close Weights Plot: " + num_states)
    end

end

function [F, MeanTransMat] = PlotMatrix(TransMatSort,ParamStruct,FitName,SubjectName,I1)
    
    num_states = size(TransMatSort,1);
    MeanTransMat = mean(TransMatSort,3);

    StatesC = ParamStruct.States;
    if num_states == StatesC(I1)

        FigureName = "Transition Matrix - States: " + num2str(num_states);
        F = figure("Name",FigureName,"Visible","on","WindowState","minimized");
    else
        FigureName = "Transition Matrix - States: " + num2str(num_states);
        F = figure("Name",FigureName,"Visible","off","WindowState","minimized");
    end
    A = axes(F);
    
    axis off
    
    cOffSet = 0.05*size(MeanTransMat,1); % use to centre values in the horizontal direction 
    imagesc(A,MeanTransMat,"AlphaData",0.65)
    colormap bone
    
    StateArray = strings(1,size(MeanTransMat,1));
    for s = 1:1:size(MeanTransMat,1)
        StateArray(s) = "States " + num2str(s);
    end
    
    xticks(1:1:num_states)
    xticklabels(StateArray)
    yticks(1:1:num_states)
    yticklabels(StateArray)
    for K1 = 1:1:size(MeanTransMat,1)
        for K2 = 1:1:size(MeanTransMat,2)
            text(A,K1-cOffSet,K2, sprintfc('%0.3f',MeanTransMat(K2,K1) ) ); 
        end
    end
    title(A, SubjectName)
    xlabel(A,"State K")
    ylabel(A,"State K - 1")

    FolderHmm = fullfile(pwd,FitName,"TransMatPlot",SubjectName);
    if ~exist(FolderHmm ,"dir")
        mkdir(FolderHmm)
    end
    Name = "TransitionMat-" + num2str(num_states)+"-"+SubjectName;
    FigName = fullfile(FolderHmm,Name + ".fig");
    
    StatesC = ParamStruct.States;
    num_states = size(TransMatSort,1);
    set(F, 'Visible', 'on');
    savefig(F,FigName);
    if num_states ~= StatesC(I1)
        close(F)
        pause(0.01)
        disp("Close Transition Mat: " + num_states)
    end

end

function [FkCell,WkmeanCell] = PlotWeightsMultiState(WkSortCell,ParamStruct,BestPriorAlpha,BestPriorSigma,FitName,SubjectName,I1)

    %input_dim = ParamStruct.InputDim; % Not Used In Calculation Scope
    NumberOfModels = ParamStruct.NumberOfModels;
    RelevantInputs = ParamStruct.RelevantInputs;
    FkCell = cell(1,NumberOfModels);
    WkmeanCell = cell(1,NumberOfModels);
    
    for k = 1:1:NumberOfModels
        %num_states = ParamStruct.States(k); % Not Used In Calculation Scope
        [FkCell{k},WkmeanCell{k}] = PlotWeights(WkSortCell{k},RelevantInputs,BestPriorAlpha,BestPriorSigma,FitName,SubjectName,ParamStruct,I1);
    end

end

function [FtCell, MeanTransMatCell] = PlotMatrixMultiState(TransMatSortCell,ParamStruct,FitName,SubjectName,I1)

    NumberOfModels = ParamStruct.NumberOfModels;
    FtCell = cell(1,NumberOfModels);
    MeanTransMatCell = cell(1,NumberOfModels);
    
    for k = 1:1:NumberOfModels
        [FtCell{k}, MeanTransMatCell{k}] = PlotMatrix(TransMatSortCell{k},ParamStruct,FitName,SubjectName,I1);
    end
end


function [GammaSessionCell,ViterbiSessionCell] = GammaVertCollect(EstGammaSortCell,EstStateViterbiSortCell,SessionFoldCell,ParamStruct,FitName,SubjectName)
    
    NumberOfModels = ParamStruct.NumberOfModels;
    NumberOfSessions = numel( unique( cat(1,SessionFoldCell{:}) ) );
    NumberOfFolds = ParamStruct.NumberOfFolds;
    
    GammaSessionCell = cell(NumberOfModels,NumberOfSessions);
    ViterbiSessionCell = cell(NumberOfModels,NumberOfSessions);

    UniqueSessions =  unique( cat(1,SessionFoldCell{:}) );
    
    for sessi = 1:1:NumberOfSessions

        sess = UniqueSessions(sessi);    
        SessionCell(sessi,:) = cellfun(@(x)x == sess,SessionFoldCell,"UniformOutput",false);
    
    end



    for st = 1:1:NumberOfModels

    
        for sess = 1:1:NumberOfSessions

            if exist("GammaSession","var")
                clear GammaSession
            end
            if exist("ViterbiSession","var")
                clear ViterbiSession
            end
        
            g = 0;
            for f = 1:1:NumberOfFolds
                
                Idx = SessionCell{sess,f};
                
                if all(Idx == 0)
                    continue        
                end
            
                g = g + 1;    
                GammaSession(:,:,g) = EstGammaSortCell{st}{f}(Idx,:);
                ViterbiSession(g,:) = EstStateViterbiSortCell{st}{f}(Idx);
            end
        
            GammaSessionCell{st,sess} = GammaSession;
            ViterbiSessionCell{st,sess} = ViterbiSession;
        end
    end
    
    DataFolder = fullfile(pwd,FitName,"VertCollect",SubjectName);    
    if ~exist("DataFolder","dir")
        mkdir(DataFolder)
    end
    
    DataHmm = fullfile(DataFolder,"GammaSessionCell" + ".mat");
    save(DataHmm,'GammaSessionCell');
    
    DataHmm = fullfile(DataFolder,"ViterbiSessionCell" + ".mat");
    save(DataHmm,'ViterbiSessionCell');

end



function [Fmode,Fmode_T,normRCell,HistCell] = ViterbiPlot(ViterbiSessionCell,ParamStruct,FitName,SubjectName,I1)

    FolderHmm = fullfile(pwd,FitName,"ViterbiMode",SubjectName);
    if ~exist(FolderHmm ,"dir")
        mkdir(FolderHmm)
    end
    
    [NumberOfModels,NumberOfSessions] = size(ViterbiSessionCell);
    
    Fmode   = cell(NumberOfModels,NumberOfSessions);
    Fmode_T = cell(NumberOfModels);
    normRCell = cell(NumberOfModels,NumberOfSessions);
    HistCell = cell(NumberOfModels,NumberOfSessions);
    for model = 1:1:NumberOfModels
            
        NumberOfStates = double(ParamStruct.States(model));
            
        for state = 1:1:NumberOfStates
            StateArray(state) = "State: "+num2str(state);
        end
        
        Name_T = "vMode_Total"+ "-States-" + num2str(NumberOfStates) + "AllSessions" ;
        if model == I1
            fmode_T = figure("name",Name_T,"Visible","on","WindowState","minimized");
        else
            fmode_T = figure("name",Name_T,"Visible","off","WindowState","minimized");
        end
        amode_T = axes(fmode_T);



      
        for sess = 1:1:NumberOfSessions
           
            Name = "vMode"+ "-States-" + num2str(NumberOfStates) + "-Session-" + num2str(sess) ;
            if model == I1 && NumberOfSessions < 20
                fmode = figure("name",Name,"Visible","on","WindowState","minimized");
            else
                fmode = figure("name",Name,"Visible","off","WindowState","minimized");
            end
            amode = axes(fmode);
           
    
            ViterbiSession = ViterbiSessionCell{model,sess};
        
            Edges = 0.5:1:NumberOfStates+0.5;
            SesionLength = size(ViterbiSession,2);
            
            if exist("HistX","var")
                clear HistX
            end
    
            for i = 1:1:SesionLength 
                HistX(i,:) = histcounts(ViterbiSession(:,i),Edges);
            end
            

            HistMax = max(HistX,[],2);
            normR = HistX./HistMax;

            %normR = normr(HistX);
            normRCell{model,sess} = normR; 
            HistCell{model,sess} = HistX;
                
            NumberOfTrials = size(ViterbiSession,2);

            hold(amode,"on")
            for i = 0:1:NumberOfStates-1
                
                lengthx = 1:1:NumberOfTrials;
                x = [lengthx lengthx(end:-1:1)];
                lengthy1 = repelem(0.5+i,NumberOfTrials);
                lengthy2 = repelem(1.5+i,NumberOfTrials);
                y = [lengthy1 lengthy2];
                c = [normR(:,i+1)' normR(end:-1:1,i+1)']';
                fill(amode,x,y,c);
                ylim([0.5 NumberOfStates+0.05]);
                xlim([1 NumberOfTrials]);
            end
            title("States: " + num2str(NumberOfStates) + " Session: " + num2str(sess));
            subtitle("Modal Values of Viterbi Algorithim");
            xlabel("Trials");
            yticks([1:1:NumberOfStates]);
            yticklabels(StateArray);
            cbh = colorbar(amode);
            NumberOfLabels = length(cbh.TickLabels);
            cbh.TickLabels = num2cell(linspace(0,10,NumberOfLabels));
            ylabel(cbh,'Folds','FontSize',8,'Rotation',270);
            cbh.Label.Position(1) = 3;
            hold(amode,"off")
            
            Fmode{model,sess} = fmode;
            FigName = fullfile(FolderHmm,Name + ".fig");
            
            num_states = double(ParamStruct.States(model));
            StatesC = ParamStruct.States;
            set(fmode, 'Visible', 'on');
            savefig(fmode,FigName);  
            if num_states ~= StatesC(I1)
                close(fmode)
                pause(0.01)
                disp("Close Viterbi: " + num_states)
            end
          
        end
       

        normC_T = vertcat(normRCell{model,:}); % Color Intensity of Viterbi Plot
        NumberOfTrials_T = size(normC_T,1);
        lengthx = 1:1:NumberOfTrials_T;
        xT = [lengthx lengthx(end:-1:1)];
        
        hold(amode_T,"on") 
        for j = 0:1:NumberOfStates-1
            lengthy1 = repelem(0.5+j,NumberOfTrials_T);
            lengthy2 = repelem(1.5+j,NumberOfTrials_T);
            yT = [lengthy1 lengthy2];
            cT = [normC_T(:,j+1)' normC_T(end:-1:1,j+1)']';    
            fill(amode_T,xT,yT,cT);
            ylim(amode_T,[0.5 NumberOfStates+0.05]);
            xlim(amode_T,[1 NumberOfTrials_T]);
        end
        title(amode_T,"States: " + num2str(NumberOfStates) + " Total");
        subtitle(amode_T,"Modal Values of Viterbi Algorithim");
        xlabel(amode_T,"Trials");
        yticks(amode_T,[1:1:NumberOfStates]);
        yticklabels(amode_T,StateArray);
        cbh_T = colorbar(amode_T);
        NumberOfLabels = length(cbh_T.TickLabels);
        cbh_T.TickLabels = num2cell(linspace(0,10,NumberOfLabels));
        ylabel(cbh_T,'Folds','FontSize',8,'Rotation',270);
        cbh_T.Label.Position(1) = 3;
        hold(amode_T,"off") 

        Fmode_T{model} = fmode_T;
        FigName_T = fullfile(FolderHmm,Name_T + ".fig");
        
        StatesC = ParamStruct.States;
        num_states = double(ParamStruct.States(model));        
        set(fmode_T, 'Visible', 'on');
        savefig(fmode_T,FigName_T);
        if num_states ~= StatesC(I1)
            close(fmode_T)
            pause(0.01)
            disp("Close Viterbi T: " + num_states)
        end



    end
end


function FmuE = GammaPlot(GammaSessionCell,LaserCell,ParamStruct,FitName,SubjectName,I1)

    [NumberOfModels,NumberOfSessions] = size(GammaSessionCell);
    NumberOfFolds = ParamStruct.NumberOfFolds;
    
    FolderHmm = fullfile(pwd,FitName,"ExpectedStates",SubjectName);
    if ~exist(FolderHmm ,"dir")
        mkdir(FolderHmm)
    end
    FmuE = cell(NumberOfModels,NumberOfSessions);
    for model = 1:1:NumberOfModels
        for sess = 1:1:NumberOfSessions
    
            GammaSession = GammaSessionCell{model,sess};
            muState = mean(GammaSession,3);
            sdState = std(GammaSession,0,3);      
            NumberOfStates = size(muState,2);
            Name = "Gamma"+ "-States-" + num2str(NumberOfStates) + "-Session-" + num2str(sess);
            
            if model == I1 && NumberOfSessions < 20
                fmuE = figure("name",Name,"Visible","on","WindowState","minimized");
            else
                fmuE = figure("name",Name,"Visible","off","WindowState","minimized");
            end

            amuE = axes(fmuE);
    
            Colours = [0.8500 0.3250 0.0980; % Orange
                       0.4660 0.6740 0.1880; % Green
                       0.0000 0.4470 0.7410; % Blue
                       0.4940 0.1840 0.5560; % Purple
                       0.9290 0.6940 0.1250; % Yellow
                       0.6350 0.0780 0.1840  % Red
                       ];
    
            if exist("StateArray","var")
                clear StateArray
            end
            
            for state = 1:1:NumberOfStates
                StateArray(state) = "State: "+num2str(state);
            end
            
            Trials = 1:1:size(muState,1);
            z = 1.96;
            if exist("p","var")
                clear p
            end
    
    
            hold on
            for state = 1:1:NumberOfStates
                p(state) = plot(amuE,muState(:,state),"Color",Colours(state,:));
                xconf = [Trials Trials(end:-1:1)] ;  
            
                yconf1 = [muState(:,state) + z*sdState(:,state)/sqrt(NumberOfFolds - 1)]';
                yconf2 = [muState(end:-1:1,state) - z*sdState(end:-1:1,state)/sqrt(NumberOfFolds -1)]';
                yconf = [ yconf1,yconf2];
            
                fill(amuE,xconf,yconf,Colours(state,:),"FaceAlpha",0.1,"EdgeColor","none");
                try
                    xline(amuE,LaserCell{sess},"Color","y","LineWidth",0.3,"Alpha",0.9);
                    %plot(amuE,LaserCell{sess},1.10,'*y')
                catch
                    % No Laser Trial
                end
            end
            title(amuE,"Expected States: "+"States: " +  num2str(NumberOfStates) + " Session: " + num2str(sess) );
            subtitle(amuE,"\gamma - Normality assumed for CI");
            hold off
            legend(p,StateArray);
    
            FmuE{model,sess} = fmuE;
            FigName = fullfile(FolderHmm,Name + ".fig");
            
            num_states = double(ParamStruct.States(model));
            StatesC = ParamStruct.States;
            set(fmuE, 'Visible', 'on');
            savefig(fmuE,FigName)
            if num_states ~= StatesC(I1)
                close(fmuE)
                pause(0.01)
                disp("Close Expected States: " + num_states)
            end
            
        end
    end
end




function WeightSubject = PlotWeightsMultiState_Subject(SubjectSort,ParamStruct,BestSubject,FitName)

    SubjectNames = ParamStruct.SubjectNames;
    for sub = 1:1:length(SubjectNames)

        SubjectName = SubjectNames(sub);

        BestPriorAlpha = BestSubject.(SubjectName).BestPriorAlpha;
        BestPriorSigma = BestSubject.(SubjectName).BestPriorSigma;
        I1 = BestSubject.(SubjectName).I1;

        WkSortCell = SubjectSort.(SubjectName).WkSortCellGmap;

        [FkCell,WkmeanCell] = PlotWeightsMultiState(WkSortCell,ParamStruct,BestPriorAlpha,BestPriorSigma,FitName,SubjectName,I1);

        WeightSubject.(SubjectName).FkCell = FkCell;
        WeightSubject.(SubjectName).WkmeanCell = WkmeanCell;

    end

end


function SubjectTrans =  PlotMatrixMultiState_Subject(SubjectSort,ParamStruct,BestSubject,FitName)

    SubjectNames = ParamStruct.SubjectNames;
    for sub = 1:1:length(SubjectNames)
        SubjectName = SubjectNames(sub);      
        TransMatSortCell = SubjectSort.(SubjectName).TransMatSortCell;
        I1 = BestSubject.(SubjectName).I1;

        [FtCell, MeanTransMatCell] = PlotMatrixMultiState(TransMatSortCell,ParamStruct,FitName,SubjectName,I1);
        SubjectTrans.(SubjectName).FtCell = FtCell;
        SubjectTrans.(SubjectName).MeanTransMatCell = MeanTransMatCell;

    end

end


function GamaVertSubject = GammaVertCollect_Subject(SubjectSort,ParamStruct,TrainData,FitName)

    SubjectNames = ParamStruct.SubjectNames;
    for sub = 1:1:length(SubjectNames)
        SubjectName = SubjectNames(sub);

        EstGammaSortCell = SubjectSort.(SubjectName).EstGammaSortCell;
        EstStateViterbiSortCell = SubjectSort.(SubjectName).EstStateViterbiSortCell;

        SessionFoldCell = TrainData.(SubjectName).SessionFoldCell;
        
        [GammaSessionCell,ViterbiSessionCell] = GammaVertCollect(EstGammaSortCell,EstStateViterbiSortCell,SessionFoldCell,ParamStruct,FitName,SubjectName);

        GamaVertSubject.(SubjectName).GammaSessionCell = GammaSessionCell;
        GamaVertSubject.(SubjectName).ViterbiSessionCell = ViterbiSessionCell;
    end
end

function ViterbiSubject = ViterbiPlot_Subject(GamaVertSubject,ParamStruct,BestSubject,FitName)

    SubjectNames = ParamStruct.SubjectNames;
    for sub = 1:1:length(SubjectNames)
        SubjectName = SubjectNames(sub);
        
        ViterbiSessionCell = GamaVertSubject.(SubjectName).ViterbiSessionCell;
        I1 = BestSubject.(SubjectName).I1;

        [Fmode,Fmode_T,normRCell,HistCell] = ViterbiPlot(ViterbiSessionCell,ParamStruct,FitName,SubjectName,I1);
        ViterbiSubject.(SubjectName).Fmode = Fmode;
        ViterbiSubject.(SubjectName).Fmode_T = Fmode_T;
        ViterbiSubject.(SubjectName).normRCell = normRCell;
        ViterbiSubject.(SubjectName).HistCell = HistCell;
    end
    
end

function GamaSubject = GamaPlot_Subject(GamaVertSubject,ParamStruct,BestSubject,TrainData,FitName)
    
    SubjectNames = ParamStruct.SubjectNames;
    for sub = 1:1:length(SubjectNames)
        SubjectName = SubjectNames(sub);

        GammaSessionCell = GamaVertSubject.(SubjectName).GammaSessionCell;
        I1 = BestSubject.(SubjectName).I1;
        LaserCell = TrainData.(SubjectName).LaserCell;

        FmuE = GammaPlot(GammaSessionCell,LaserCell,ParamStruct,FitName,SubjectName,I1);
        GamaSubject.(SubjectName).FmuE = FmuE;
    end

end