%{
-About-
This script contains command lines to run
TASEPmodeling_bursty_par or TASEPmodeling_nonbursty_par and to analyze the
run results

-Inputs-
conditions for promoter activity (bursty, nonbursty, kinetic parameters), 
conditions for RNAP elongation (RNAP speed, pause site, duration, probability)

-varargin-

-Outputs-W
output1: results of the TASEPmodeling_par runs will be saved as independent
%matlab data files. See the last lines of TASEPmodeling_(non)bursty_par.m
for what is contained in those data files
output2: results of TASEPmodeling_par_analysis will be in the Workspace. See
TASEPmodeling_par_analysis for the details

-Example-
masterscript
   
-Supplementary-

-Keywords-
TASEP modeling, gene expression

-Dependencies-

-References-

-Author-
Sangjin Kim, 2017 September 30
%}

%--------------------------------------------------------------------------
% This "masterscript" is to run the TASEP model in a parallel fashion 
% This parfor-based code is advantageous for running 10000 simulations in a timely manner
% This parfor-based code only records statistics from each simulation and does not save raw trajectories
% For raw trajectories of individual motors (RNA polymerases and ribosomes), see code_for_RNAP_trafficking
% -----------------------------------------------------------------------------------------------------
%% Example scenario: 
% A bursty promoter, average RNAP loading interval ~15sec, no pause
close all; clear;
pauseProfile = 'flat'; 
promoter='P1'; %arbitrary name, but will be included the output file name
fOn = 0.25; %fOn constant for promoter
avgSpeed =30; pauseSite = 0; pauseDuration = 0; pauseProb = 0; %for RNAP elongation profile
mRNALL = 90; proteinLL = 0;
for i = 1:10 %each run does 100 iterations
     TASEPmodeling_bursty_par(pauseProfile,promoter,i,fOn,avgSpeed,pauseSite, pauseDuration, pauseProb, mRNALL, proteinLL)
end;
%% Result files are saved in the folder.
%% For analysis of the par results:
%% per DNA template or per simulations
TASEPmodeling_par_analysis('flat','P1',10,pauseSite,pauseDuration,pauseProb)


%%-------------
%% More example scenarios: result files are not save in the folder
%% Example scenario: 
%% A bursty promoter, average RNAP loading interval ~15sec, pause at xp = 1500 for tp = 10sec
close all; clear;
pauseProfile = 'OnepauseAbs'; 
promoter='P1'; 
fOn = 0.25;
avgSpeed =30;  pauseSite = 1500; pauseDuration = 10; pauseProb = 80;
mRNALL = 90; proteinLL = 0;
for i = 1:10
     TASEPmodeling_bursty_par(pauseProfile,promoter,i,fOn,avgSpeed,pauseSite, pauseDuration, pauseProb, mRNALL, proteinLL)
end;

%%
TASEPmodeling_par_analysis('OnepauseAbs','P1',10,pauseSite,pauseDuration,pauseProb)

%% Example scenario: 
%% A bursty promoter, average RNAP loading interval ~15sec, pause at xp = 1500 and 2500 for tp = 10, 15 sec respectively
close all; clear;
pauseProfile = 'MultipauseAbs'; 
promoter='P1'; 
fOn = 0.25;
avgSpeed =30; pauseSite = [1500,2500]; pauseDuration = [10,15]; pauseProb = [100,100];
mRNALL = 90; proteinLL = 0;
for i = 1:1
     TASEPmodeling_bursty_par(pauseProfile,promoter,i,fOn,avgSpeed,pauseSite, pauseDuration, pauseProb, mRNALL, proteinLL)
end;

%%
TASEPmodeling_par_analysis('MultipauseAbs','P1',10,pauseSite,pauseDuration,pauseProb)

%% Example scenario: 
%% A nonbursty promoter, average RNAP loading interval ~15sec, no pause
close all; clear;
pauseProfile = 'flat'; 
promoter='P1nb';
kLoading = 1/15;
avgSpeed =30;  pauseSite = 0; pauseDuration = 0; pauseProb = 0;
mRNALL = 90; mRNALL2 = 40; proteinLL = 0;
for i = 1:1
     TASEPmodeling_nonbursty_par(pauseProfile,promoter,i,kLoading,avgSpeed,pauseSite, pauseDuration, pauseProb, mRNALL, mRNALL2, proteinLL)
end;

close all; clear;
pauseProfile = 'OnepauseAbs'; 
promoter='P1nb';
kLoading = 1/15;
avgSpeed =30;  pauseSite = 1500; pauseDuration = 10; pauseProb = 80;
mRNALL = 90; proteinLL = 180;
for i = 1:10
     TASEPmodeling_nonbursty_par(pauseProfile,promoter,i,kLoading,avgSpeed,pauseSite, pauseDuration, pauseProb, mRNALL, proteinLL)
end;

TASEPmodeling_par_analysis('OnepauseAbs','P1nb',10,pauseSite,pauseDuration,pauseProb)
