function ansAll = TASEPmodeling_par_analysis(pauseProfile,promoter,totalparRun,pauseSite,pauseDuration,pauseProb)

%{
-About-
This function analyzes results of TASEPmodeling_(non)bursty_par in
case of the cells in the population contain only 1 gene copy per cell
This code combine results from separate parfor runs and analyze 
This generates data that is plotted in the paper

-Inputs-
pauseProfile: name of the pause profile (this will be included in the output file name)
promoter: name of the promoter (this will be included in the output file name)
totalparRun: total par file number (the output files should be numbered as
1, 2, 3, ... totalparRun. If only out output file, use 1)
pauseSpecific = pauseSite: location of the pause (nt)
pauseDuration: duration of the pause (sec)
pauseProb: probability of pausing at the site (%)

-varargin-

-Outputs-
ansAll: Transcription/translation initiation rate (or intervals), headway
distribtuion, RNAP density on the DNA template, Distribution of mRNAs and
proteins per cell (mean, fano factor, CV2), duration with no new protein
production

-Example-
   
-Supplementary-

-Keywords-
TASEP modeling, gene expression

-Dependencies-
see masterscript

-References-

-Author-
Sangjin Kim, 2017 September 30
%}

%--------------------------------------------------------------------------

if strcmp(pauseProfile, 'flat')
    runCondition = strcat(pauseProfile,'-NO-');
elseif strcmp(pauseProfile, 'OnepauseAbs')
    runCondition = strcat(pauseProfile,sprintf('%01.0f',pauseSite),'x',sprintf('%01.0f',pauseDuration),'xp',sprintf('%01.0f',pauseProb),'-NO-');
elseif strcmp(pauseProfile, 'MultipauseAbs')
    runCondition = strcat(pauseProfile,sprintf('%01.0f',length(pauseSite)),'x',sprintf('%01.0f',pauseDuration(1)),'x',sprintf('%01.0f',pauseDuration(2)),'xp',sprintf('%01.0f',pauseProb(1)),'-NO-');
end;

loadingSample = []; presentonDNA = [];
lifeTimeHist1 = []; lifeTimeHist2= []; lifeTimeAvg = []; 
fishSignal1 = []; fishSignal2 = []; fishSignal3 = []; 
rEndStampA = []; riboNumHistFine = []; 
proteinLoci2 = []; proteinSS = [];
tDiffStartHist = []; tDiffEndHist = []; tDiffDHist = [];
% Load all par runs and combine together
for i = 1:totalparRun
    fileName = strcat(runCondition,promoter,'-',sprintf('%01.0f',i),'par.mat'); 
    tmp = load(fileName);
    loadingSample = [loadingSample, tmp.loadingSample];
    presentonDNA = [presentonDNA, tmp.presentonDNA];
    
    tDiffStartHist = [tDiffStartHist, tmp.tDiffStartHist];
    tDiffEndHist = [tDiffEndHist, tmp.tDiffEndHist];
    tDiffDHist = [tDiffDHist, tmp.tDiffDHist];
    
    lifeTimeHist1 = [lifeTimeHist1 tmp.lifeTimeHist1]; 
    lifeTimeHist2 = [lifeTimeHist2 tmp.lifeTimeHist2]; 
    lifeTimeAvg = [lifeTimeAvg tmp.lifeTimeAvg];
    
    fishSignal1 = [fishSignal1 tmp.fishSignal1];
    fishSignal2 = [fishSignal2 tmp.fishSignal2];
    fishSignal3 = [fishSignal3 tmp.fishSignal3];
    
    rEndStampA = [rEndStampA tmp.rEndStampA];
    riboNumHistFine = [riboNumHistFine tmp.riboNumHistFine];
    
    proteinLoci2 = [proteinLoci2 tmp.proteinLoci2];
    proteinSS = [proteinSS tmp.proteinSS];
    
    specificDwelltime = tmp.specificDwelltime;
end;

sampleWindow = [15*60,30*60];

%% ans1 loadingSample = Number of successful initiation per DNA during sampleWindow
ans1 = BootstrapMeanNoise(loadingSample/(15*60),3000); %This is for effective loading rate (sec-1)
ans1(2,:) = BootstrapMeanNoise((15*60)./loadingSample,3000); %This is for effective loading interval (sec)
ansAll.tsxInitiationrate = ans1(:,1:2);

%% ans2 presentonDNA = Number of RNAP on a DNA at an instant (t = 1200 s)
ans2 = presentonDNA';
ansAll.RNAPnumDNA = ans2;

%% ans3 tDiffStart(End)Hist = temporal separation between RNAPs at the start and emd
binT = 0:5:200;
ans3(:,1) = binT';
ans3(:,2) = mean(tDiffStartHist,2);
ans3(:,3) = mean(tDiffEndHist,2);
ans3(end+1,:) = [205,0,0];
ansAll.RNAPheadway = ans3;

%% ans30 tDiffDHist = headway change from initiation to elongation (headway END - headway START)
bindT = -50:5:50;
ans30(:,1) = bindT';
ans30(:,2) = mean(tDiffDHist,2);
ansAll.RNAPheadwaychange = ans30;

%% ans4 mRNA life-time
%% DNA with no RNAP initiation during sampleWindow will show mRNA life-time = 0
%% and hence should be eliminated from the analysis
tmp2 = find(loadingSample>0);
binLT = 0:5:1000;
ans4 = [];
ans4(:,1) = binLT';
ans4(:,2) = mean(lifeTimeHist1(:,tmp2),2); %life-time of 5'-end mRNA
ans4(:,3) = mean(lifeTimeHist2(:,tmp2),2); %life-time of 3'-end mRNA
ans4(:,4) = hist(lifeTimeAvg(1,tmp2),binLT)*100/length(tmp2); %life-time of 5'-end mRNA
ans4(:,5) = hist(lifeTimeAvg(2,tmp2),binLT)*100/length(tmp2); %life-time of 3'-end mRNA
ans4stat = BootstrapMeanNoise(lifeTimeAvg(1,tmp2),3000);
ans4stat(10:17) = BootstrapMeanNoise(lifeTimeAvg(2,tmp2),3000);
ansAll.mRNAlifetime = ans4;
ansAll.mRNAlifetimestat = ans4stat(:,1:2);

%% ans5/6 mRNA FISH signal
%% probe1 = 5'-end; probe2 = 3'-end; probe3 = 72 tiling probes
ssTime = 30; %steady-state time index of the fish time
ssfish1 = []; ssfish2 = []; ssfish3 = []; 
for i = ssTime
    ssfish1 = [ssfish1, fishSignal1(i,:)]; 
    ssfish2 = [ssfish2, fishSignal2(i,:)];
    ssfish3 = [ssfish3, fishSignal3(i,:)]; 
end;

% For mRNA distribution at steady-state
ans5 = [];
ans5(:,1) = (0:1:50)'; %binN';
ans5(:,2) = hist(ssfish1', ans5(:,1))*100/length(ssfish1);
ans5(:,3) = hist(ssfish2', ans5(:,1))*100/length(ssfish2);
ans5(:,4) = hist(ssfish3', ans5(:,1))*100/length(ssfish3);
ansAll.mRNAdist = ans5;

% For mean, fano, CV, CV^2
ans6(1,:) = BootstrapMeanNoise(ssfish1',3000);
ans6(2,:) = BootstrapMeanNoise(ssfish2',3000);
ans6(3,:) = BootstrapMeanNoise(ssfish3',3000);
ansAll.mRNAstat = ans6;


%% ans7/8 = Total protein number per DNA during 10 min
% % For distribution
% ans7(:,1) = (0:100:2000)'; %binP
ans7(:,1) = (0:200:4000)'; %binP
ans7(:,2) = hist(proteinLoci2',ans7(:,1))*100/length(proteinLoci2);
ansAll.proteindist = ans7;
% For mean and fano (mean+/-ste; fano +/- ste; CV +/-ste, CV^2 +/- ste)
ans8 = BootstrapMeanNoise(proteinLoci2,3000);
ansAll.proteinstat = ans8;


%% ans9 = # of ribosomes per mRNA
binFine = 0:1:400;
binFinemat = repmat(binFine',1,size(riboNumHistFine,2));
riboNum = sum(binFinemat.*riboNumHistFine,1)/100;
riboNum(find(isnan(riboNum))) = []; %Delete NaN which is due to 0 mRNA on DNA
ans9 = BootstrapMeanNoise(riboNum,3000); % only for mean +/-ste
% Translation initiation rate is ribosome number divided by mRNA lifetime
ansAll.tslefficiency = ans9(1:2);

%%  New protein production rate per DNA
% Example trajectories of protein production (Fig. 3)
% are from rEndStampA(t,DNA), which contains
% # of new protein produced per DNA, every 10sec during t = (600:10:2400).

% Time duration with zero protein production
X = rEndStampA(2:end-1,:);
% X = cellRibostamp(2:end-1,:);
k = 0; zeroDall = [];
for i = 1:size(X,2)
    % time point with 0 will get 1 time points >0 will get 0
    zeroTs = X(:,i) == 0;
    if isempty(find(X(:,i) ==0)) % no zeros 
        zeroDall = [zeroDall 0];
    else
        h = 0; k = 0; zeroD = 0;
        for j = 1:length(zeroTs)
            if zeroTs(j)==1
                h = h+1; %count 1's 
                if (j<length(zeroTs) && zeroTs(j+1) == 0)
                    k = k+1;
                    zeroD(k) = h;
                end;
            else
                h = 0;
            end;
        end;
        if zeroD == 0
            zeroDall = [zeroDall 0];
        else
            zeroDall = [zeroDall zeroD(2:end)];
        end;
    end;
end;
ans10 = BootstrapMeanNoise(zeroDall,3000);
ansAll.durationzeroprotein = ans10(1:2);

%% examine steady-state protein levels
fishTime = 0:10:(size(fishSignal1,1)-1)*10;
figure();
%yyaxis left
xlabel('Time (sec)'); 
ylabel('mRNA per loci');
plot(fishTime, mean(fishSignal1,2), 'r-o'); hold on; plot(fishTime, mean(fishSignal2,2), 'g-o');% plot(fishTime, mean(fishSignal3,2), 'b-');
legend("5' signal", "3' signal", "Location", "northeast");
%yyaxis right
%ylabel('Protein per loci'); % plot(fishTime, mean(proteinSS,2), 'k-'); 
%hold off;
yyaxis left 
ylabel('mRNA per loci')

ssTime = size(fishSignal1,1); %steady-state time index of the fish time
ssProtein = [];
for i = ssTime
    ssProtein = [ssProtein, proteinSS(i,:)]; 
end;
ans11(1,:) = BootstrapMeanNoise(ssProtein',3000);
ansAll.steadystateprotein = ans11;
