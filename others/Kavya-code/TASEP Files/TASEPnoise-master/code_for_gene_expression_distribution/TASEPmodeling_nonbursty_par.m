function ans = TASEPmodeling_nonbursty_par(pauseProfile,promoter,parRun,kLoading,avgSpeed,pauseSite,pauseDuration,pauseProb, mRNALL, mRNALL2, proteinLL)
%{
-About- (edited by Kavya Vaidya 2020-21)
This function runs TASEP modeling of transcription, translation and mRNA degradation based on parfor function
Nonursty promoter case: same as TASEPmodeling_bursty_par, except
transcription initiation

-Version changes-
Editor: Kavya Vaidya. 5/19/2021
Edited to include two-phase degradation (separating co-transcriptional and post-transcriptional rates. mRNALL and mRNALL2 respectively) so it needs file "exp2rnd.m" to
run. It also has no Par Runs.
Usual variable changes to explore - geneLength, avgSpeed, N_totalloci,
stoptime (Stop loading of RNAP, around 90 seconds for induction-repression
and ~720 seconds for steady-state studies(depends on mRNALL values))

The output file name in this version is - 
('m',string(mRNALL),string(mRNALL2),'.', string(stopTime),'.',string(geneLength), '.',string(avgSpeed), 'bps','.mat')


-Inputs-
pauseProfile: name of the pause profile (this will be included in the output file name)
promoter: name of the promoter (this will be included in the output file name)
parRun: integer, indicating # of files in a series (only to name output
file name)
kLoading: RNAP loading rate (transcription initiation rate) in 1/sec
avgSpeed: average speed of RNAP (and ribosome)
pauseSite: location of the pause (nt)
pauseDuration: duration of the pause (sec)
pauseProb: probability of pausing at the site (%)
mRNALL: mRNA lifetime 1; (co-transcriptional)
mRNALL 2: mRNA lifetime 2; (post-transcriptional
proteinLL: protein lifetime. If 0, do not consider protein degradation.

-varargin-

-Outputs-
output: matlab data file containing RNAP loading number per template,
headway between RNAPs at initiation, headway between RNAPs at termination,
mRNA number, protein number (steady state and accumulation), temporal profile of protein production
(rEndStampA), RNAP number density on the DNA template, number of ribosomes
per transcript (translation efficiency), mRNA lifetime   

-Example-
TASEPmodeling_nonbursty_par('flat','P1nb',1,1/15,30,0,0,0,90,1200)
   
-Supplementary-

-Keywords-
TASEP modeling, gene expression, parfor, nonbursty

-Dependencies-
masterscript: see masterscript for how to use this function

-References-

-Author-
Sangjin Kim, 2017 September 30
This version: Kavya Vaidya. 5/19/2021
%}

%--------------------------------------------------------------------------

%% Parameters
% All changeable parameters are listed here. Feel free to modify depending
% on the purpose

dx = 1; % lattice site = 1 nt
geneLength = 1000; %3075; % total gene legnth in nt 
polWidth = 35;  % RNAP size in nt
N_totalloci = 1000; % total number of iterations when this function is called
simTime = 0:1:40*60; % simulation time in seconds (only the first and last time point is important)
stopTime = 720; %simTime(end) %edit this to stop loading
riboWidth = 30; % ribosome size in nt

% Input parameters for ribosome initiation and elongation
kRiboLoading = 0.2; % initiation rate of ribosome (1/sec); 
riboSpeed = avgSpeed; % ribosome speed in nt/sec. Same as RNAP

% Nonbursty initiation parameters: kLoading = input value
 
% Parameters for analysis
sampleWindow = [15*60,30*60]; % simulation time window when the analysis happened
loadingSample = zeros(1,N_totalloci); % # of loading per loci during sampling

sampleT = [20*60,1201];
presentonDNA = zeros(1, N_totalloci); %for # of RNAPs present on a DNA

% mRNA life-time
binLT = 0:5:1000;
lifeTimeHist1 = zeros(length(binLT),N_totalloci); %for mRNA lifetime distribution
lifeTimeHist2 = zeros(length(binLT),N_totalloci); %for mRNA lifetime distribution 
lifeTimeAvg = zeros(2,N_totalloci); %for mRNA lifetime average
criticalTime = zeros(N_totalloci); %in case of two-phase degradation

% Number of ribosomes per mRNA (translation efficiency)
binFine = 0:1:400;
riboNumHistFine = zeros(length(binFine),N_totalloci);

% bin for time points for translation completion (rEndStampA)
% temporal profile of protein production
riboTbin = 600:10:2400;
rEndStampA = zeros(length(riboTbin),N_totalloci);

%Distribution of mRNA numbers (FISH)
probe1 = 1;%40:40:960; %1 % probing 5'-end mRNA
%probe12 = 40:40:960;
probe2 = geneLength; %2000:40:2940;%geneLength; % probing 3'-end mRNA
%probe22 = 2000:40:2940;%760:40:960;
fishTime = 0:10:1800;%0:10:1800;%max(simTime); %I think making this 0:30:max should be enough?
fishSignal1 = zeros(length(fishTime),N_totalloci); % #mRNA per loci
fishSignal12 = zeros(length(fishTime),N_totalloci); % #mRNA per loci
fishSignal2 = zeros(length(fishTime),N_totalloci); % #mRNA per loci
fishSignal22 = zeros(length(fishTime),N_totalloci); % #mRNA per loci

% for steady-state (SS) protein amount
proteinSS =  zeros(length(fishTime),N_totalloci); % for #protein per loci

%for distribution of RNAP interval (t-headway)
binT = 0:5:200;
tDiffStartHist = zeros(length(binT),N_totalloci);
tDiffEndHist = tDiffStartHist;
bindT = -50:5:50;
tDiffDHist = zeros(length(bindT),N_totalloci);

%% Define RNAP dwelltime profile 
avgDwelltime = dx/avgSpeed; %sec per nucleotide
if strcmp(pauseProfile, 'flat')
    specificDwelltime = avgDwelltime * ones(geneLength/dx,1); % no pause, flat profile
    runCondition = strcat(pauseProfile,'-NO-');
elseif strcmp(pauseProfile, 'OnepauseAbs')
    specificDwelltime = avgDwelltime * ones(geneLength/dx,1); % no pause, flat profile
    specificDwelltime(pauseSite) = pauseDuration;
    runCondition = strcat(pauseProfile,sprintf('%01.0f',pauseSite),'x',sprintf('%01.0f',pauseDuration),'xp',sprintf('%01.0f',pauseProb),'-NO-');
elseif strcmp(pauseProfile, 'MultipauseAbs')
    specificDwelltime = avgDwelltime * ones(geneLength/dx,1); % no pause, flat profile
    specificDwelltime(pauseSite) = pauseDuration;
    runCondition = strcat(pauseProfile,sprintf('%01.0f',length(pauseSite)),'x',sprintf('%01.0f',pauseDuration(1)),'x',sprintf('%01.0f',pauseDuration(2)),'xp',sprintf('%01.0f',pauseProb(1)),'-NO-');
end;


%% par run
%parallelobject = parpool(12);
tic 
for lociID = 1:N_totalloci %each loci is independent
    
    %---RNAP loading attempts
    % For each RNAP, there are...
    % [1] loading time = time point of loading "attempt"
    % [2] exitTime = stepping time at each nt of an RNAP
    % [3] polStatus = status of RNAP at each nt
    % 0 not yet loaded or elongating (0.5 for collision); 3 loading inhibited;
    % This is the section that only differs from bursty mode of run
    loadTime = [];
    pauseRNAP = [];
    exitTime = -inf(geneLength/dx+polWidth/dx,1); 
    polStatus = zeros(geneLength/dx+polWidth/dx,1);
    
    j = 0;
    t = simTime(1) + exprnd(1/kLoading)*rand;
    % because previous loading (exprnd(1/kLoading) before can be before
    % simTime(1)
    while t < stopTime 
        j = j+ 1; %loaded RNAP on a locus
        loadTime(j) = t;
        pauseRNAP(j) = rand<= (pauseProb/100); % 1= yes pause 0 = no pause
        if j == 1
            % for the first RNAP, determine stepping through entire length 
            % because there is none ahead.
            if pauseRNAP(j) == 0
                specificDwelltime1 = avgDwelltime * ones(geneLength/dx,1); % no pause, flat profile
            else
                specificDwelltime1 = specificDwelltime;
            end;
            exitTime(1:geneLength/dx,j) = t + cumsum(exprnd(specificDwelltime1));
        else
            exitTime(1:geneLength/dx,j) = zeros(geneLength/dx,1);
        end;
        polStatus(1:geneLength/dx,j) = zeros(geneLength/dx,1);
        t = t + exprnd(1/kLoading);
    end;
    
    %---Elongation/ RNAP translocation
    % Iterate over RNAP (p = polID) -> iterate over x = 1:geneLength
    for p = 2:length(loadTime)
        % check which pol is ahead on the template
        ahead = 0;
        for i = p:-1:2
            if polStatus(1,i-1)<1
                ahead = i-1; 
                break;
            end;
        end;
        
        % check if the loading attempt (loadTime) is hindered by RNAP ahead
        % loadTime = arrivalTime at x = 1 = (exitTime at x = 0)
        if loadTime(p) <= exitTime(polWidth/dx,ahead) 
            %loading inhibited
            polStatus(:,p) = 3*ones(size(polStatus,1),1); % Loading failed
            exitTime(1:geneLength/dx,p) = -inf(geneLength/dx,1); % reset (no loading)
        else
            % translocate the RNAP 
           if pauseRNAP(p) == 0 % no pause
                specificDwelltime1 = avgDwelltime * ones(geneLength/dx,1); % no pause, flat profile
            else
                specificDwelltime1 = specificDwelltime; %original input, with pause
            end;
            for x = 1:geneLength
                % arrival time at x = exitTime at x-1
                if x == 1
                    arrival = loadTime(p);
                else
                    arrival = exitTime(x-1,p);
                end;

                % exitTime at x = arrival time at x + dwell time at x (=shift)
                shift = exprnd(specificDwelltime1(x));
                exitTime(x,p) = arrival + shift;

                % check for collision with RNAP ahead
                % collision -> trailing RNAP polStatus = 0.5; leading RNAP
                % polStatus = 0.2
                % collision -> wait until leading RNAP steps
                if arrival + shift <=  exitTime(x+polWidth/dx,ahead) 
                    polStatus(x,p) = polStatus(x,p)+0.5; 
                    polStatus(x+polWidth/dx,ahead) = polStatus(x+polWidth/dx,ahead)+ 0.2;
                    exitTime(x,p) = exitTime(x+polWidth/dx,ahead) + shift;
                end;
            end;
        end;
    end;
    
    %---Translation
    % Initialize parameters
    decayTime=zeros(size(exitTime));
    ribo = 0;
    ribotEnd = zeros(1,p); 
    lifeTimePol = [];
    
    % Iterate over elongated RNAP = mRNAs
    for polID = min(find(exitTime(1,:)>0)):max(find(exitTime(1,:)>0)) %checks if RNAP loaded on first position(5') and how many.
        if exitTime(1,polID)>0
            % Determine life-time of 5'-end mRNA (lifeTime)
            %if exitTime(geneLength, polID) > 0
            %    tcrit = exitTime(geneLength, polID);
            %    lifeTime = exp2rnd(mRNALL,mRNALL2,tcrit,1,1);
            %else
            %    lifeTime = exprnd(mRNALL,1,1);
            %end
            %lifeTime = exprnd(mRNALL,1,1);
            tcrit = exitTime(geneLength, polID)-loadTime(polID);
            lifeTime = exp2rnd(mRNALL,mRNALL2,tcrit,1,1);
            lifeTimePol(polID) = lifeTime;
            
            % For all ribo's for this mRNA:
            % Each ribo requires [1] ribo loadTime = Time point of load attempt 
            % [2] ribo exitTime = ribo translocation on mRNA
            % [3] ribo status = ribo status at each nt position
            RiboloadTime = []; 
            RiboexitTime = zeros(size(exitTime,1),floor(2*lifeTime*kRiboLoading)); 
            RiboStatus = zeros(size(RiboexitTime));
            
            % First ribo j = 1 show couping to RNAP for loading and
            % translocation
            coupleWidth = round((polWidth + riboWidth)/2);
            j = 1;
            t = exitTime(coupleWidth,polID);
            RiboloadTime(j) = t;
            % if the first t was after lifeTime, no ribo. 
            % keep riboexitTime for degradation though (delayed by 40 nt from RNAP)
            noriboFlag = 0;
            if t >= exitTime(1,polID)+lifeTime
                noriboFlag = 1;
                RiboexitTime(1:geneLength,j) = exitTime(1:geneLength,polID) + lifeTime;
            else
                RiboDwelltime = exprnd(1/riboSpeed,coupleWidth,1);
                RiboexitTime(1:(geneLength-coupleWidth),j) = exitTime(coupleWidth+1:geneLength,polID);
                RiboexitTime((geneLength-coupleWidth+1):geneLength,j) = RiboexitTime(geneLength-coupleWidth,j)+cumsum(RiboDwelltime);
            end;
            
            % loadTime for j>1 ribo's
            % load continues until 5'-end mRNA is gone (lifeTime)
            t = t + exprnd(1/kRiboLoading);
            while t < exitTime(1,polID) + lifeTime; 
                j = j+ 1; 
                RiboloadTime(j) = t;
                RiboexitTime(1:geneLength/dx,j) = zeros(geneLength/dx,1);
                RiboStatus(1:geneLength/dx,j) = zeros(geneLength/dx,1);
                t = t + exprnd(1/kRiboLoading);
            end;            
            
            % exitTime for j>1 ribo's
            % Iterate over ribos (r = riboID)-> iterate over x = 1:geneLength
            for r = 2:length(RiboloadTime)
                % check which ribo is ahead on the mRNA template
                ahead = 0;
                for i = r:-1:2
                    if RiboStatus(1,i-1)<1
                        ahead = i-1; 
                        break;
                    end;
                end;
                % check if the loading attempt (loadTime) is hindered by ribo ahead
                % loadTime = arrivalTime at x = 1 = (exitTime at x = 0)
                if RiboloadTime(r) <= RiboexitTime(riboWidth/dx,ahead) 
                    %loading inhibited
                    RiboStatus(:,r) = 3*ones(size(RiboStatus,1),1); % No loading "3"
                    RiboexitTime(1:geneLength/dx,r) = -inf(geneLength/dx,1); % reset (no loading)
                else
                    % translocate ribo
                    for x = 1:geneLength
                        if x == 1
                            arrival = RiboloadTime(r);
                        else
                            arrival = RiboexitTime(x-1,r);
                        end;
                        % exitTime at x = arrival time at x + dwell time at x (=shift)
                        shift = exprnd(1/riboSpeed,1,1);
                        RiboexitTime(x,r) = arrival + shift;
                        % check for collision with ribo ahead
                        if arrival + shift <=  RiboexitTime(x+riboWidth/dx,ahead) 
                            RiboStatus(x,r) = RiboStatus(x,r)+0.5; 
                            RiboStatus(x+riboWidth/dx,ahead) = RiboStatus(x+riboWidth/dx,ahead)+ 0.2;
                            RiboexitTime(x,r) = RiboexitTime(x+riboWidth/dx,ahead) + shift;
                        end;
                    end;
                end;
            end; 
            
            % Check last ribo's trajectory and take it as a degradation of
            % mRNA
            lastribo = max(find(RiboexitTime(1,:)>0));
            decayTime(:,polID) = RiboexitTime(:,lastribo);
            
            % Number of ribosomes loaded per mRNA
            ribo(polID) = length(find(RiboexitTime(1,:)>0));
            % Record time point of ribosome arrival at the end of gene
            % = time of protein production/release
            for h = 1:size(RiboexitTime,2)
                if RiboexitTime(1,h)>0
                    ribotEnd(h,polID) = RiboexitTime(geneLength,h);
                end;
            end;
            
            if noriboFlag == 1, 
                ribo(polID) = 0; 
                ribotEnd(1,polID) = 0;
            end;
        end; 
    end;
    %filling the critical time
    criticalTime(1:size(exitTime,2),lociID) = (exitTime(geneLength,:));;
    % Check # of protein produced at certain times = trajectory of protein
    % production number
    for polID = 1: size(ribotEnd,2);
        rEndStampA(:,lociID) = rEndStampA(:,lociID) + hist(ribotEnd(find(ribotEnd(:,polID)>0),polID),riboTbin)';   
    end;

    
    %---Ribo Analysis
    % Check proteins accumulated per DNA during t = 1200-1800 (10 min)
    % Check proteins at each time point
    if isempty(find(ribotEnd>0)) 
        proteinLoci2(lociID) = 0;
        proteinSS(:,lociID) = 0;
    else
        proteinLoci2(lociID) = length(find(ribotEnd>1200 & ribotEnd<=1800));
        proteinmadeTime = sort(ribotEnd(:)'); 
        proteinmadeTime(find(proteinmadeTime==0)) = [];
        % one protein is made at the end of translation elongation
        % protein decays proteinLL after one protein is made.
        proteindecayTime = proteinmadeTime + exprnd(proteinLL,1,length(proteinmadeTime));
        % for protein SS (FISH)
        for tt = 1:length(fishTime)
            made1 = length(find(proteinmadeTime < fishTime(tt)));
            decayed1 = length(find(proteindecayTime <fishTime(tt)));
            loc = zeros(length(fishTime),1); loc(tt) = 1;
            proteinSS(:,lociID) = proteinSS(:,lociID) + (made1-decayed1)*loc;
        end
    end;
    
    
    %---Analysis during sampleWindow time
    rNumber = []; 
    i = 0; 
    lastPol = max(find(exitTime(1,:)>0));
    lifeTimeSite = [0,0];
    tStart = []; tEnd = []; dheadway = [];
    if isempty(lastPol)
        loadingSample(lociID) = 0;
    else
        for polID = 1: lastPol
            if loadTime(polID)>=sampleWindow(1) && loadTime(polID)<sampleWindow(2)     
                if exitTime(1,polID) >=0 % successfully loaded RNAP
                    i = i+1; % traj index in a lociID
                    % number of successful RNAP initiation per loci
                    loadingSample(lociID) = loadingSample(lociID) + 1;
                    % 5'-end mRNA life time
                    lifeTimeSite(i,1) = lifeTimePol(polID);
                    % 3'-end mRNA life time
                    lifeTimeSite(i,2) = decayTime(geneLength,polID)-exitTime(geneLength,polID); %based on exit at the last base
                    % Number of RNAPs on DNA at sample T (or t = 20 min)
                    presentonDNA(1,lociID) = presentonDNA(1,lociID) + ~isempty(find(exitTime(:,polID)>=sampleT(1) & exitTime(:,polID)<sampleT(2)));
                    % Number of ribosomes loaded per mRNA
                    rNumber(i) = ribo(polID);
                    % For the interval between RNAPs upon loading
                    tStart(i) = loadTime(polID); % or exitTime(1,polID)
                    % For the interval between RNAPs upon completion
                    tEnd(i) = exitTime(geneLength,polID);

                end;
            end;
        end;
        lifeTimeHist1(:,lociID) = hist(lifeTimeSite(:,1),binLT)*100/size(lifeTimeSite,1)';
        lifeTimeHist2(:,lociID) = hist(lifeTimeSite(:,2),binLT)*100/size(lifeTimeSite,1)';%gives #RNAP x iterations(loops,loci)
        lifeTimeAvg(:,lociID) = [mean(lifeTimeSite(:,1)),mean(lifeTimeSite(:,2))];
        riboNumHistFine(:,lociID) = hist(rNumber,binFine)*100/length(rNumber)';
        tDiffStartHist(:,lociID) = hist(diff(tStart),binT)*100/max((length(tStart)-1),1);
        tDiffEndHist(:,lociID) = hist(diff(tEnd),binT)*100/max((length(tEnd)-1),1);
        tDiffDHist(:,lociID) = hist(diff(tEnd)-diff(tStart),bindT)*100/max((length(tEnd)-1),1);
    end;
       
    % quantify mRNA number distributions at different simulation time points
    % similar to fluorescence in situ hybridization (FISH) for quantifying mRNA 
    for tt = 1:length(fishTime)
        for polID = 1: size(exitTime,2)
            if exitTime(1,polID) >0 %realTraj
               made1 =  length(find(exitTime(probe1,polID)<fishTime(tt)))*1;
               decayed1 = length(find(decayTime(probe1,polID)<fishTime(tt)))*1;
               made2 =  length(find(exitTime(probe2,polID)<fishTime(tt)))*1;
               decayed2 = length(find(decayTime(probe2,polID)<fishTime(tt)))*1;
               %made12 =  length(find(exitTime(probe12,polID)<fishTime(tt)))*1;
               %decayed12 = length(find(decayTime(probe12,polID)<fishTime(tt)))*1;
               %made22 =  length(find(exitTime(probe22,polID)<fishTime(tt)))*1;
               %decayed22 = length(find(decayTime(probe22,polID)<fishTime(tt)))*1;
                             
               loc = zeros(length(fishTime),1); loc(tt) = 1;
               fishSignal1(:,lociID) = fishSignal1(:,lociID) + (made1-decayed1)/length(probe1)*loc;
               fishSignal2(:,lociID) = fishSignal2(:,lociID) + (made2-decayed2)/length(probe2)*loc;
               %fishSignal12(:,lociID) = fishSignal12(:,lociID) + (made12-decayed12)/length(probe12)*loc;
               %fishSignal22(:,lociID) = fishSignal22(:,lociID) + (made22-decayed22)/length(probe22)*loc;
               
            end;
        end;
    end;
end;
%delete(parallelobject);
toc

fileName = strcat('m',string(mRNALL),string(mRNALL2),'.', string(stopTime),'.',string(geneLength), '.',string(avgSpeed), 'bps','.mat');
save(fileName,'specificDwelltime','loadingSample','presentonDNA','fishSignal1','fishSignal2','geneLength',...
    'criticalTime', 'fishTime', 'mRNALL', 'mRNALL2','avgSpeed','stopTime','riboNumHistFine','rEndStampA','proteinLoci2','proteinSS',...
    'lifeTimeHist1','lifeTimeHist2','lifeTimeAvg',...
    'tDiffStartHist','tDiffEndHist','tDiffDHist',...
    'runCondition','promoter',...
    '-v7.3');


