%% A nonbursty promoter, average RNAP loading interval ~10(was15)sec, no pause
close all; clear;
pauseProfile = 'OnepauseAbs'; 
promoter='P1nb';
kLoading = 1/10;
avgSpeed = 10; %30;    
pauseSite = 1500; pauseDuration = 10; pauseProb = 80; 
mRNALL = 120; mRNALL2 = 90; proteinLL = 0;
for i = 1:1
     TASEPmodeling_nonbursty_par(pauseProfile,promoter,i,kLoading,avgSpeed,pauseSite, pauseDuration, pauseProb, mRNALL, mRNALL2, proteinLL)
     %avgSpeed = avgSpeed + 10;
end;

%%
pauseProfile = 'flat'; 
promoter='P1nb';
kLoading = 1/10;
avgSpeed =10;  pauseSite = 0; pauseDuration = 0; pauseProb = 0;
mRNALL = 70; mRNALL2 = 90; proteinLL = 0;
for i = 1:1
     TASEPmodeling_nonbursty_par(pauseProfile,promoter,i,kLoading,avgSpeed,pauseSite, pauseDuration, pauseProb, mRNALL, mRNALL2, proteinLL)
end;
%TASEPmodeling_par_analysis('flat','P1nb',1,0,0,0);
%speedcrit = zeros(6,4);
%%
% Kavya Vaidya: Analysis begins here
load('m4590.1000.10.90.mat');
clearvars -except fishSignal1 fishSignal2 mRNALL mRNALL2 criticalTime speedcrit
%%
fishTime = (0:10:(size(fishSignal2,1)-1)*10)'; %fishTime (fishSig at pos 1 <=> time = 0
%%
fishSignal = mean(fishSignal2,2); % Change what signal you want to look at
%% to print log plot
figure();
plot(fishTime,log(mean(fishSignal1,2)), 'r-o');
hold on;
%%
% To print normal plot
figure();
plot(fishTime,fishSignal, 'ro');
xlabel("FISH time course (seconds)");
ylabel("5' signal");
title("Simulation of 5' signal for mu = 90 seconds")
%% 
% get more information about data
Last_RNAP_exit = max(criticalTime, [], 'All');
criticalTime(criticalTime < 1) = NaN; %omitting all none time related values
avg_First_RNAP_exit = mean(criticalTime(1,:), 'omitnan');
First_RNAP_exit = min(criticalTime, [], 'All', 'omitnan');
%i = 2;
%speedcrit(i,1) = i*10;
%speedcrit(i,2) = tcrit;
%speedcrit(i,3) = tcrit_max;
%speedcrit(i,4) = tcrit_min;
%% 
% manually select regions to fit two-phase degradation in log plot (so
% linear equation)
% fit first half
a1 = 14;
b1 = 18;
x1 = fishTime(a1:b1);
y1 = log(mean(fishSignal1,2));
f1 = fit(x1, y1(a1:b1), 'poly1');
m1 = -(1/(f1.p1));
xf1 = linspace(fishTime(a1), fishTime(b1), 100);
yf1 = f1.p1*xf1 + f1.p2;
plot(xf1, yf1, 'k', 'linewidth', 1);
hold on;
%fit second half
a2 = 21;
b2 = 35;
x2 = fishTime(a2:b2);
y2 = log(mean(fishSignal1,2));
f2 = fit(x2, y2(a2:b2), 'poly1');
m2 = -(1/(f2.p1));
xf2 = linspace(fishTime(a2), fishTime(b2), 100);
yf2 = f2.p1*xf2 + f2.p2;
plot(xf2, yf2, 'b', 'linewidth', 1);
xlabel("FISH time course (seconds)");
ylabel("log(fishSignal1)");
title(sprintf("5' signal fitting, mu1 = %g, mu2 = %g", mRNALL, mRNALL2));
legend({"5' data",sprintf('m1 = %g (%g - %g)', m1, x1(1), x1(end) ),sprintf('m2 = %g (%g - %g)', m2, x2(1), x2(end))}, 'Location', 'northeast');
axis([0 800 -6 2]);
hold off;
