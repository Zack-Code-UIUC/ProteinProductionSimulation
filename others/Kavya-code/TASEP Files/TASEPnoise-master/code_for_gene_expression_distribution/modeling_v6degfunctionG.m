%%%%%%%
% use this to check Z5 degradation
%%%%%%%%
function ans = modeling_v6degfunctionG(s2,s2after,s4,processingDiff,speed, DF, pulse)

% Y = 0;
timeDuration = 14*60; % total time period 12 min or 720 seconds.
dt = 10; % time step in 10 second
N_timestep = timeDuration/dt; % total number of time steps
geneLength = 5000; % bases
enzymeWidth = 40; % bases (physical size + random min distance that two pols to not come close together)
%% input parameters for elongation models
% speed = 30;
AF = 1; %DF = 0.4;
processingRate = 0; % Z/Y processing 0 = no processing; 40 = 40% of chance of processing
%% input parameters from data fitting
% pulse = 'late';
s2 = s2 *60/dt;%(in min) 
s2 = 1/s2; %(in min-1)
s2after = s2after*60/dt; %(in min) 
s2after = 1/s2after;
s1 = 2; s3 = 0.25; % per min
% s1 = 4.1209; s3 = 1.0656; %high IPTG
% s1 = 3.18; s3 = 1.12; %high IPTG -short
% s1 = 1.2; s3 = 0.24; %B8 short
% s1 = 1.5; s3 = 0.24; %B8 long
%% RE: promoter-proximal peak, initially slow.
slowPeriod = 20; % in seconds
slowPeriod = slowPeriod /dt; % during 30 seconds after loading, speed is low.
initiationSpeedReduction = 1; %during initiation, speed is 33% of regular speed;

if strcmp(pulse,'early')
%     speed = 6.417*3; % 14.7; %nt/sec %6.417*3;
    Toff = 1.5*60; % time when glucose was added to shut off (in second)
elseif strcmp(pulse,'late')
%     speed = 9.91*3;%22.8;
    Toff = 3.5*60;
end;
speed = speed*dt; %nt/dt
s3 = s3/60*dt;
loadingInterval = 60/s1/dt; %timestep interval between RNAP loading (initiation 1/min => 1/6 timestep)
s4 = s4/60*dt; %LacI rebinding
% activation_rate = 0.0275; % per timestep

%% Populate cells with 2+gene loci
N_totalcell = 1000; % total number of cells in simulation
doublegeneloci = 57.8; % % of cells with 2 gene loci
cell_geneloci = ones(N_totalcell,1);
cell_geneloci(1:floor(doublegeneloci*N_totalcell/100)) = 2;
% cell_geneloci = cell_geneloci(randperm(length(cell_geneloci)),:);

gl = 1;
for i = 1:N_totalcell
    genelociTable(gl:(gl+cell_geneloci(i)-1),1) = i; %cell number
    genelociTable(gl:(gl+cell_geneloci(i)-1),2) = 1:cell_geneloci(i); %loci number 1st or 2nd loci in a cell
    gl = gl + cell_geneloci(i);
end;
% genelociTable = genelociTable(randperm(length(genelociTable)),:);

N_totalloci = gl-1;
for i = 1:N_totalloci, lociList(i).ontime = 0; end; %initialize as 0.
for i = 1:N_totalloci, lociList(i).ontime = exprnd(1/s3,1,1); end; %initialize as 0.

%% Elongation - for simplicity, assume very processive + no pause or backtracking during elongation
initiationblockcount = 0; collisioncount = 0;
for lociID = 1:N_totalloci % loci are sorted such that all on loci are in 1:N_cumOnloci
    lociOntime = lociList(lociID).ontime; %in timestep
%     lociOntime = exprnd(1/s3,1,1);
%     lociList(lociID).ontime = lociOntime;
    lociOfftime = Toff/dt + exprnd(1/s4,1,1);
    lociList(lociID).offtime = lociOfftime;
    % % determine # and timing of RNAP loading
    polID = 0; initiationPol = 0; 
    t = lociOntime + exprnd(loadingInterval); %%%EXPONENTIAL DIST!!
    while t < lociOfftime
        polID = polID + 1;
        initiationPol(polID) = t; % Pol initiation time (in timestep)
        t = t +  exprnd(loadingInterval); % in each loading, get random loading interval.
    end; % % if lociOntime was close to Toff, no loading can happen.
    
    initiationPol = sort(initiationPol); %XXdue to randn, sorting is necessary to have initiation time in ascending order.
    lociList(lociID).numPols = polID;
    lociList(lociID).initiationTime = initiationPol;
    
    positionPol = zeros(N_timestep,polID); polStatus = zeros(N_timestep,polID); finishedRNA = 0;
    % % move RNAPs over time
    for t = 1 : (N_timestep-1)
        for p = 1 : polID
            if t+1 > initiationPol(p) && polStatus(t,p)==0 % initiation always at the timestep edge
                thisPos = positionPol(t,p);
                nextPos = positionPol(t+1,p);
                dts = min(t+1-initiationPol(p),1);
                if t < initiationPol(p) + slowPeriod
                    initiationSpeed = initiationSpeedReduction;
                else
                    initiationSpeed = 1;
                end;
 
                % % % % count the number of trailing RNAPs
                count = 0;
                if p<polID && positionPol(t,p+1)>speed*initiationSpeedReduction*slowPeriod % p is not the last RNAP, 1+ behind pols is elongating.
                    for i = p:1:(polID-1)
                        if positionPol(t,i+1)>speed*initiationSpeedReduction*slowPeriod
                            count = count+1;
                        end;
                    end;
                end;
                 % % % %
                if t>=lociOfftime %&& t<=Toff/dt+90/dt
%                     if t>=Toff/dt && count<10% 
%                 if t>=Toff/dt %&& thisPos>1500 %t>=Toff/dt && count <10 % thisPos>1500 && count<3% 
                    deccelerationFactor = DF;
                    accelerationFactor = 1;
                    processingRate = processingDiff;
                else
                    deccelerationFactor = 1;
                     accelerationFactor = AF^(count);
                     processingRate = 0;
%                         accelerationFactor = 1 + (AF-1)*(count);
                end;
                
                speedFactor = 1 + randn * 0.2; 
                shift = speedFactor * speed * initiationSpeed * accelerationFactor * deccelerationFactor * dts;
                nextPos = thisPos + shift;
%                 plot(t*dt/60,nextPos,'o'); hold on;
                if t <= initiationPol(p) % detect collision before initiation
                    for i = p : -1 : 2
                        if polStatus(t,i-1) == 0
                            x1 = positionPol(t,i-1) + (positionPol((t+1),(i-1))-positionPol(t,i-1)) * (initiationPol(p)-t); %interpolate position at firing moment
                            if thisPos + enzymeWidth > x1
                                nextPos = 0;
                                polStatus(:,p) = 3; % 3 = RNAP loading was inhibited
%                                 display(['initiation inhibited for', num2str(p)]);
                                inititationblockcount = initiationblockcount + 1;
                            end;
                            break;
                        end;
                    end;
                end;
                
                if nextPos > 0  %avoid collisions by checking first existing polymerase ahead of me
                    for i = p : -1 : 2
                        if polStatus(t+1,i-1) == 0 %% instead of if (polExistence[i-1][t+1]=='m') 
                            if nextPos + enzymeWidth > positionPol(t+1,i-1)
                                nextPos = positionPol(t+1,i-1) - enzymeWidth;
%                                 display(['collision for', num2str(p)]);
                                collisioncount = collisioncount + 1;
                            end;
%                             break; %needed?
                        end;
                    end;
                end;
                positionPol(t,p) = thisPos;
                positionPol(t+1,p) = nextPos;
                
                %% early release at the end of Z. ProcessingRate = 40 means, 40% of time, it is early released. Keep this 0 to bypass.
                if nextPos >= 3000 && thisPos <3000 && randi(100,1) < processingRate
                    polStatus(t+1:N_timestep,p) = 2; % 2 means released
                    finishedRNA = finishedRNA + 1;
                    positionPol(t+1:N_timestep,p) = 3000;
                end;
                
                %% at the end of gene (ZYA), RNAP is released, and does not move this polid any more by scoring status as "2".
                if nextPos >= geneLength
                    polStatus(t+1:N_timestep,p) = 2; % 2 means released
                    finishedRNA = finishedRNA + 1;
                    positionPol(t+1:N_timestep,p) = geneLength;
                end;
            end;
        end;
    end;
    lociList(lociID).positionPol = positionPol;
    lociList(lociID).polStatus = polStatus;
%     Y(lociID,1) = finishedRNA;
%     if finishedRNA <1;
%         display('..');
%     end;
end;

%% 1. lociList --> cellList and calculate signal per cell, etc
% % assign loci into cells. Shuffle/random
% % 2. Degradation = inactivation (s2) & ribosome protection
%%%%%%%%%%%%%%%%%%%%%%%%
probe1Region5 = 0; probe1Region3 = 1000;  
probe2Region5 = 2000; probe2Region3 = 3000;  
probe1Signal = zeros(N_timestep,N_totalloci);
probe2Signal = zeros(N_timestep,N_totalloci); %nascent only
probe1SignalAll = zeros(N_timestep,N_totalloci); % released + nascent
probe2SignalAll = zeros(N_timestep,N_totalloci);

x = (1:N_timestep)'*dt/60;
totalpol = 0;
Z5spot = zeros(N_timestep,1); Z3spot = zeros(N_timestep,1);

% figure,
for lociID = 1 : length(lociList)
    if ~isempty(find(lociList(lociID).polStatus == 2)) % loci with .ontime>0 might have no loading if ontime is close to Toff
        polTraj = lociList(lociID).positionPol;
        polStat = lociList(lociID).polStatus;
        fiveTraj = zeros(size(polTraj)); % 5'end of mRNA
        initiationTime = lociList(lociID).initiationTime;
        lociOfftime = lociList(lociID).offtime;
%         Z5inactivationTime = initiationTime + exp2rnd(1/s2,1/s2after,lociOfftime-initiationTime);

%         % % search the time when the leading RNAP enters Y for
%         % 'translocation'
%         enterYtime = [];
%         for polID = 1:size(polTraj,2)
%             if polStat(1,polID) == 0
%                 if ~isempty(find(polTraj(:,polID)>3200,1))
%                     enterYtime = find(polTraj(:,polID)>3200,1); %time when the loci moves due to the leader
%                     break;
%                 end;
%             end;
%         end;
%         Z5inactivationTime = initiationTime + exp2rnd(1/s2,1/s2after,enterYtime-initiationTime);
        % % 2regime by lacY translocation
        releasedTime = zeros(size(polTraj,2),1);     
        for polID = 1:size(polTraj,2)
            if polStat(1,polID) ==0 && ~isempty(find(polStat(:,polID)==2,1))%initiation not inhibited
                releasedTime(polID) = find(polStat(:,polID)==2,1);
                Z5inactivationTime(polID) = initiationTime(polID) + exp2rnd(1/s2,1/s2after,releasedTime(polID)-initiationTime(polID));
                    % % 2regime by cotsx/posttsx
                if floor(Z5inactivationTime(polID))+1 < releasedTime(polID) 
                    if floor(Z5inactivationTime(polID))+1 <= 1, Z5inactivationTime(polID) = 1; end;
                    %co-transcriptional
                    %follow ribo
                    for t = floor(Z5inactivationTime(polID))+1 : releasedTime(polID)
                        fiveTraj(t,polID) = fiveTraj(t-1,polID) + polTraj(t,polID) - polTraj(t-1,polID);
                    end;
                    for t = releasedTime(polID)+1 :N_timestep
                        fiveTraj(t,polID) = fiveTraj(t-1,polID) + 40*dt;
                    end;
                    polStat(floor(Z5inactivationTime(polID))+1:releasedTime(polID),polID) = 4; %assign co-tsx degradation
                    polStat(releasedTime(polID)+1: end, polID) = 5; %post-tsx degradation
                else %post-transcriptional, follow ribo
                    for t = floor(Z5inactivationTime(polID))+1 : N_timestep
                        fiveTraj(t,polID) = fiveTraj(t-1,polID) + 40*dt;
                    end;
                    polStat(floor(Z5inactivationTime(polID))+1:end, polID) = 5;
                end;
            end;
                
            for t = 1:N_timestep
                if polStat(t,polID) == 0 || polStat(t,polID) == 4 %elongating
%                         plot(t*dt/60,polTraj(t,polid),'o'); hold on;
                    if polTraj(t,polID)>probe1Region5 && fiveTraj(t,polID)<probe1Region3
                        covered = min(probe1Region3,polTraj(t,polID))-fiveTraj(t,polID);
                        probe1Signal(t,lociID) = probe1Signal(t,lociID) + covered/1000;
                        probe1SignalAll(t,lociID) = probe1SignalAll(t,lociID) + covered/1000;
                    end;
                    if polTraj(t,polID)>probe2Region5 && fiveTraj(t,polID)<probe2Region3
                        covered = min(probe2Region3,polTraj(t,polID))-max(fiveTraj(t,polID),probe2Region5);
                        probe2Signal(t,lociID) = probe2Signal(t,lociID) + covered/1000;
                        probe2SignalAll(t,lociID) = probe2SignalAll(t,lociID) + covered/1000;
                    end;
                elseif polStat(t,polID) == 2 || polStat(t,polID) == 5 %released
                    if fiveTraj(t,polID)<probe1Region3
                       covered = min(probe1Region3,polTraj(t,polID))-fiveTraj(t,polID); 
                        probe1SignalAll(t,lociID) = probe1SignalAll(t,lociID) + covered/1000;
                    end;
                    if fiveTraj(t,polID)<probe2Region3
                        covered = min(probe2Region3,polTraj(t,polID))-max(fiveTraj(t,polID),probe2Region5);
                        probe2SignalAll(t,lociID) = probe2SignalAll(t,lociID) + covered/1000;
                    end;
                end;
            end;
        end;
        lociList(lociID).fiveTraj = fiveTraj;
        lociList(lociID).polStatusdecay = polStat;
    end;
    
%     %loci that has Z5 or Z3 signal are counted.
%     Z5spot(find(probe1SignalAll(:,lociID)>0)) = Z5spot(find(probe1SignalAll(:,lociID)>0)) + 1;
%     Z3spot(find(probe2SignalAll(:,lociID)>0)) = Z3spot(find(probe2SignalAll(:,lociID)>0)) + 1;
end;

SMthrsh = 0;
% for i = 1:N_timestep
%     probe1DetachedSignal(i,:) = probe1SignalAll(i,:)-probe1Signal(i,:);
%     probe2DetachedSignal(i,:) = probe2SignalAll(i,:)-probe2Signal(i,:);
%     noNascent1(i,1) = length(find(probe1Signal(i,:)>SMthrsh)); %nascent spot #
%     noReleased1(i,1) = sum(probe1DetachedSignal(i,:)); % released spot # assume individual molecule separates into individual spots
%     noNascent2(i,1) = length(find(probe2Signal(i,:)>SMthrsh));
%     noReleased2(i,1) = sum(probe2DetachedSignal(i,:));
%     probe1Signalperspot(i,1) = sum(probe1SignalAll(i,:))/(noNascent1(i,1)+noReleased1(i,1));
%     probe2Signalperspot(i,1) = sum(probe2SignalAll(i,:))/(noNascent2(i,1)+noReleased2(i,1));
%     probe1Signalperloci(i,1) = mean(probe1Signal(i,:)); %nascent signal
%     probe2Signalperloci(i,1) = mean(probe2Signal(i,:));
%     probe1SignalAllperloci(i,1) = mean(probe1SignalAll(i,:)); % all signal
%     probe2SignalAllperloci(i,1) = mean(probe2SignalAll(i,:));
% end;

% x = (1:N_timestep)'*dt/60;
% Xans(:,1) = x;
% Xans(:,2) = probe1Signalperspot; %spot intensity
% Xans(:,3) = probe2Signalperspot;
% 
% Z5spot = Z5spot/N_totalloci*100;
% Z3spot = Z3spot/N_totalloci*100;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for shuffle, combine loci into a cell

for shufflei = 1:100
    
    genelociTable = genelociTable(randperm(length(genelociTable)),:); %(cell number, 1st/2nd loci)
    genelociTable(:,3) =(1:N_totalloci)'; %give locus index
    
    probe1SignalIncell = zeros(N_timestep,N_totalcell); probe2SignalIncell = zeros(N_timestep,N_totalcell); 
    probe1SpotnumIncell = zeros(N_timestep,N_totalcell); probe2SpotnumIncell = zeros(N_timestep,N_totalcell);
    for cellID = 1:N_totalcell
        clear lociNum_tmp;
        lociNum_tmp = genelociTable(find(genelociTable(:,1) == cellID),3); % e.g., [#3,#1000]
        for lociID = 1:length(lociNum_tmp)
            probe1SignalIncell(:,cellID) = probe1SignalIncell(:,cellID)+ probe1SignalAll(:,lociNum_tmp(lociID)); %all signal per cell
            probe2SignalIncell(:,cellID) = probe2SignalIncell(:,cellID) + probe2SignalAll(:,lociNum_tmp(lociID));
%             probe1Nascent = probe1Signal(:,lociNum_tmp(lociID))>SMthrsh; probe2Nascent = probe2Signal(:,lociNum_tmp(lociID))>SMthrsh;
%             % # of loci containing signal = 1 or 0
%             probe1SpotnumIncell(:,cellID) = probe1SpotnumIncell(:,cellID) + probe1SignalAll(:,lociNum_tmp(lociID)) - probe1Signal(:,lociNum_tmp(lociID)) +probe1Nascent;
%             probe2SpotnumIncell(:,cellID) = probe2SpotnumIncell(:,cellID) + probe2SignalAll(:,lociNum_tmp(lociID)) - probe2Signal(:,lociNum_tmp(lociID)) +probe2Nascent;
        end;
    end;
    
    for i = 1:N_timestep
        probe1Signalpercell(i,shufflei) = mean(probe1SignalIncell(i,:)); % all signal per cell
        probe2Signalpercell(i,shufflei) = mean(probe2SignalIncell(i,:));
%         probe1Spotnumpercell(i,shufflei) = mean(probe1SpotnumIncell(i,:)); % spot number per cell
%         probe2Spotnumpercell(i,shufflei) = mean(probe2SpotnumIncell(i,:));
%         probe1Percentcell(i,shufflei) = length(find(probe1SignalIncell(i,:)>0))/N_totalcell*100;
%         probe2Percentcell(i,shufflei) = length(find(probe2SignalIncell(i,:)>0))/N_totalcell*100;
%         probe1Spotnumpercell2(i,shufflei) = mean(probe1SpotnumIncell(i,find(probe1SpotnumIncell(i,:)>0))); % spot number per cell (>0)
%         probe2Spotnumpercell2(i,shufflei) = mean(probe2SpotnumIncell(i,find(probe2SpotnumIncell(i,:)>0)));
    end;
end;

for i = 1:N_timestep
    XansShuffle(i,1) = mean(probe1Signalpercell(i,:));
    XansShuffle(i,2) = mean(probe2Signalpercell(i,:));
%     XansShuffle(i,3) = Xans(i,2);
%     XansShuffle(i,4) = Xans(i,3);
%     XansShuffle(i,5) = mean(probe1Percentcell(i,:));
%     XansShuffle(i,6) = mean(probe2Percentcell(i,:));
%     XansShuffle(i,7) = mean(probe1Spotnumpercell(i,:));
%     XansShuffle(i,8) = mean(probe2Spotnumpercell(i,:));
%     XansShuffle(i,9) = mean(probe1Spotnumpercell2(i,:));
%     XansShuffle(i,10) = mean(probe2Spotnumpercell2(i,:));
end;


Xans(:,1) = x;
Xans(:,2) = XansShuffle(:,1);
Xans(:,3) = XansShuffle(:,2);

MGearlyspotDATA = [0,0.0410681000000000,0.0428664000000000;1,0.202405000000000,0.0113691000000000;2,0.614529000000000,0.0381357000000000;3,0.525201000000000,0.162129000000000;4,0.439032000000000,0.461142000000000;5,0.284526000000000,0.490169000000000;6,0.143265000000000,0.227940000000000;7,0.117927000000000,0.170424000000000;8,0.0689753000000000,0.101503000000000;9,0.0221215000000000,0.0409212000000000;10,0.0686922000000000,0.0930030000000000;12,0.0307892000000000,0.0231767000000000;];
MGlatespotDATA = [0,0.0148500000000000,0.0307189000000000;1,0.173199000000000,0.0199088000000000;2,0.920316000000000,0.157376000000000;3,1.44341000000000,0.850463000000000;4,1.88367000000000,1.29346000000000;5,1.56962000000000,1.44383000000000;6,0.626673000000000,1.10983000000000;7,0.407285000000000,0.872549000000000;8,0.227583000000000,0.548508000000000;9,0.153310000000000,0.277760000000000;10,0.0814297000000000,0.131522000000000;11,0.0533351000000000,0.128123000000000;12,0.0369998000000000,0.0728880000000000;14,0.0480840000000000,0.133206000000000;];
MGearlyintDATA = [0,-2.05117000000000e-05,-1.85171000000000e-05;1,0.212811000000000,-0.0606560000000000;2,0.822809000000000,0.0259615000000000;3,0.732416000000000,0.213399000000000;4,0.589036000000000,0.667880000000000;5,0.382955000000000,0.638346000000000;6,0.177199000000000,0.286242000000000;7,0.139322000000000,0.224066000000000;8,0.127647000000000,0.153899000000000;9,-0.0242904000000000,0.100817000000000;10,0.0615288000000000,0.0820574000000000;12,0.0535092000000000,0.00259389000000000;];
MGlateintDATA = [0,-0.000173119000000000,-3.25718000000000e-05;1,0.385416000000000,0.214730000000000;2,1.57367000000000,0.661704000000000;3,2.25131000000000,1.59036000000000;4,2.74765000000000,2.33183000000000;5,2.39896000000000,2.58055000000000;6,1.32193000000000,2.35691000000000;7,1.00142000000000,1.96133000000000;8,0.710183000000000,1.35341000000000;9,0.497579000000000,0.753505000000000;10,0.221720000000000,0.463888000000000;11,-0.000197473000000000,0.223018000000000;12,0.114521000000000,0.176315000000000;14,0.00369274000000000,0.226152000000000;];

if strcmp(pulse,'early')
    XansDATA1 = MGearlyspotDATA;
    XansDATA2 = MGearlyintDATA;
elseif strcmp(pulse,'late')
    XansDATA1 = MGlatespotDATA;
    XansDATA2 = MGlateintDATA;
end;
figure, subplot(1,2,1), plot(Xans(:,1),Xans(:,2),'r-'); hold on; plot(Xans(:,1),Xans(:,3),'g-');
plot(XansDATA1(:,1),XansDATA1(:,2),'ro'); plot(XansDATA1(:,1),XansDATA1(:,3),'go'); hold off;
title(['s2 ',num2str(dt/60/s2), ' s2after ', num2str(dt/60/s2after)]); 
subplot(1,2,2), plot(Xans(:,1),Xans(:,2),'r-'); hold on; plot(Xans(:,1),Xans(:,3),'g-');
plot(XansDATA2(:,1),XansDATA2(:,2),'ro'); plot(XansDATA2(:,1),XansDATA2(:,3),'go'); hold off;

ans = [0,0,0,0];
for i = 2: size(XansDATA1,1);
    ti = find(x == XansDATA1(i,1));
    ans(1) = ans(1) + (Xans(ti,2) - XansDATA1(i,2))^2; %compare Z5
    ans(2) = ans(2) + (Xans(ti,3) - XansDATA1(i,3))^2; %compare Z3
    ans(3) = ans(3) + (Xans(ti,2) - XansDATA2(i,2))^2; %compare Z5
    ans(4) = ans(4) + (Xans(ti,3) - XansDATA2(i,3))^2; %compare Z5
end;

