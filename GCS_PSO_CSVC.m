classdef GCS_PSO_CSVC< ALGORITHM
    % <multi> <real> <expensive>
    % wmax --- 20 --- The maximum number of internal evluation
    
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [educational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------
    
    % This function is written by Qi-Te Yang
    
    methods
        function main(Algorithm,Problem)
            %% Parameter setting 
            [wmax] = Algorithm.ParameterSet(20);  
            %% Generate initial population based on Latin hypercube sampling
            InitN          = 11*Problem.D-1;
            P          = UniformPoint(InitN,Problem.D,'Latin');
            Population = Problem.Evaluation(repmat(Problem.upper-Problem.lower,InitN,1).*P+repmat(Problem.lower,InitN,1));            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                %% database
                [~,index]  = unique(Population.decs,'rows');
                Population = Population(index);
                randIndex = randperm(length(Population));
                Population = Population(randIndex);
                PopDec = Population.decs;
                PopObj = Population.objs;
                [N,~] = size(PopObj);
                
                div = 30;
                %% Calculate the grid location of each solution
                fmax = max(PopObj,[],1);
                fmin = min(PopObj,[],1);
                lb   = fmin-(fmax-fmin)/2/div;
                ub   = fmax+(fmax-fmin)/2/div;
                d    = (ub-lb)/div;
                lb   = repmat(lb,N,1);
                d    = repmat(d,N,1);
                GLoc = floor((PopObj-lb)./d);
                GLoc(isnan(GLoc)) = 0;
                
                %% train a SVM for each objective
                Model = cell(1,Problem.M);
                for i = 1 : Problem.M
                        Model{i} = svmtrain(GLoc(:,i),PopDec,'-s 0 -t 1 -c 1 -g 2.8');
                end           
                [front,~] = NDSort(PopObj,inf);
                FrontValue = front';
                
                %% Grid PF
                GPF = GLoc(FrontValue==1,:);    
                OptDec = PopDec(FrontValue==1,:);               
                %% optimization
                choose = EnvironmentalSelectionNSGA2(Population,Problem.N);
                ArcDec = PopDec(choose,:);
                ArcGLoc = GLoc(choose,:);
                w = 1;
                while w < wmax
                    %% update GR and NCI
                    [Problem.N,~] = size(ArcDec);
                    ArcGR = sum(ArcGLoc,2);   % ArcGR
                    ArcNCI = zeros(Problem.N,1);
                    for i = 1 : Problem.N
                        ArcNCI(i) = NCI_Cal(ArcNCI,i,Problem);
                    end
                    
                    %% SLPSO
                    Xmean = mean(ArcDec);
                    OffDec = SLPSOmyself(ArcDec,OptDec,Xmean);
                    OffGLoc = zeros(Problem.N,Problem.M);
                    for i = 1 :Problem.M
                        [OffGLoc(:,i),~,~] = svmpredict(OffGLoc(:,i),OffDec,Model{i},'-q 0');
                    end
                    %% GR and NCI calculation of offsprings
                    OffGR = sum(OffGLoc,2);
                    OffNCI = zeros(Problem.N,1);
                    for i = 1 : Problem.N
                        TempGLoc = ArcGLoc;
                        TempGLoc(i,:) = OffGLoc(i,:);
                        OffNCI(i) = NCI_Cal(TempGLoc,i,Problem);
                    end
                    
                    %% update
                    % update swarm
                    replace = (OffGR < ArcGR) | (OffGR == ArcGR & OffNCI < ArcNCI);
                    ArcDec(replace,:) = OffDec(replace,:);
                    ArcGLoc(replace,:) = OffGLoc(replace,:);
                    % update GPF
                    for i = 1 : Problem.N
                        [isDom,DomNo] = GDominated(ArcGLoc(i,:),GPF);
                        if  isDom
                            OptDec(DomNo,:) = [];
                            GPF(DomNo,:) = [];
                            OptDec = [OptDec;ArcDec(i,:)];
                            GPF = [GPF;ArcGLoc(i,:)];
                        end
                    end
                    [~,index]  = unique(OptDec,'rows');
                    OptDec = OptDec(index,:);
                    GPF = GPF(index,:);
                    w = w + 1;
                end
                
                %% infill solutions criterion
                % remove solutions that have been real evaluated
                OldDec = Population.decs;
                Loc = [];
                for i = 1 : size(ArcDec,1)
                    for j = 1 : N
                        if isequal(ArcDec(i,:),OldDec(j,:))
                            Loc = [Loc,i];
                            break;
                        end
                    end
                end
                ArcDec(Loc,:) = [];
                ArcGLoc(Loc,:) = [];
                
                % infill criterion
                [n,~] = size(ArcDec);
                ArcGR = sum(ArcGLoc,2);   % ArcGR
                if n <= 5
                    NewArc = ArcDec;
                else
                    choose = EnvironmentalSelectionNSGA2(Population,5); % choose reference points
                    Ref = PopDec(choose,:);
                    Location = zeros(1,n);
                    for i = 1 : n
                        X = (Ref - repmat(ArcDec(i,:),5,1));
                        dis = sum(abs(X).^2,2).^(1/2);
                        [~,Location(i)] = min(dis);
                    end
                    for i = 1 : 5
                        if isempty(find(Location==i))
                            X = (ArcDec - repmat(Ref(i,:),n,1));
                            [~,loc] = min(sum(abs(X).^2,2).^(1/2));
                            Location(loc) = i;
                        end
                    end
                    NewArc = [];
                    for i = 1 : 5
                        Ci = find(Location==i);
                        CA = ArcDec(Ci,:);
                        [~,loc] = min(ArcGR(Ci));
                        NewArc = [NewArc;CA(loc,:)];
                    end 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
                
                if ~isempty(NewArc)
                    PopNew = Problem.Evaluation(NewArc);
                    Population = [Population,PopNew];
                end
                
            end
        end
    end
end