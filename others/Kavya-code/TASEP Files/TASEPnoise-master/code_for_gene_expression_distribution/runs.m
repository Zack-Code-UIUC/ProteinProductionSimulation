%%
%Kavya Vaidya: This to run and make files based on various parameters
pauseProfile = 'flat';
promoter='P1nb';
kLoading = 1/10;
avgSpeed =10;  pauseSite = 0; pauseDuration = 0; pauseProb = 0;
mRNALL = 40; mRNALL2 = 50; proteinLL = 0;
for i = 1:1
     TASEPmodeling_nonbursty_par(pauseProfile,promoter,i,kLoading,avgSpeed,pauseSite, pauseDuration, pauseProb, mRNALL, mRNALL2, proteinLL)
end;
%%
pauseProfile = 'flat'; 
promoter='P1nb';
kLoading = 1/10;
avgSpeed =10;  pauseSite = 0; pauseDuration = 0; pauseProb = 0;
mRNALL = 300; mRNALL2 = 200; proteinLL = 0;
for i = 1:1
     TASEPmodeling_nonbursty_par(pauseProfile,promoter,i,kLoading,avgSpeed,pauseSite, pauseDuration, pauseProb, mRNALL, mRNALL2, proteinLL)
end;
%%
pauseProfile = 'flat'; 
promoter='P1nb';
kLoading = 1/10;
avgSpeed =10;  pauseSite = 0; pauseDuration = 0; pauseProb = 0;
mRNALL = 140; mRNALL2 = 120; proteinLL = 0;
for i = 1:1
     TASEPmodeling_nonbursty_par(pauseProfile,promoter,i,kLoading,avgSpeed,pauseSite, pauseDuration, pauseProb, mRNALL, mRNALL2, proteinLL)
end;
