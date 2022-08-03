function ans = BootstrapMeanNoise(dataX,bsN)
%{
-About-
This function computes statistics of dataX (e.g., mean and noise) 
with standard error by the bootstrapping method.
This was used to analyze TASEP simulation data by Sangjin Kim.

-Inputs-
dataX: One dimensional array of raw data
bsN: number of iterations (e.g. 3000)

-Outputs-
ans: mean, ste(standard error) of mean, fano, ste of fano, CV, ste of CV,
CV^2, ste of CV2

-Example-
BootstrapMeanNoise(X,3000)   

-Supplementary-

-Keywords-
bootstrap, error estimation, 

-Dependencies-
This function is called from TASEPmodeling_par_analysis

-References-

-Author-
Sangjin Kim, 2017 September 16
%}

% For bootstrapping: randomly pick elements of raw data in each iteration
for n = 1:bsN
    tmp_i = unidrnd(length(dataX),length(dataX),1);
    dataStat(n,1) = mean(dataX(tmp_i)); %mean
    dataStat(n,2) = var(dataX(tmp_i))/mean(dataX(tmp_i)); %fano
    dataStat(n,3) = std(dataX(tmp_i))/mean(dataX(tmp_i)); %CV
    dataStat(n,4) = (dataStat(n,3))^2; %CV^2
end

%ans = [mean_mean, mean_ste, fano_mean, fano_ste, CV_mean, CV_ste, CV2_mean, CV2_ste];
ans = [mean(dataStat(:,1)),std(dataStat(:,1)),...
    mean(dataStat(:,2)),std(dataStat(:,2)),...
    mean(dataStat(:,3)),std(dataStat(:,3)),...
    mean(dataStat(:,4)),std(dataStat(:,4))];