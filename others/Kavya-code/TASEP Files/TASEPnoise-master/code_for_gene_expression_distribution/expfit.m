function [fitresult, gof] = createFit(fishTime, fishSignal, mRNALL, mRNALL2)
%CREATEFIT(FISHTIME,FISHSIGNAL)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : fishTime
%      Y Output: fishSignal
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 21-Jul-2020 05:07:31


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( fishTime, fishSignal );

% Set up fittype and options. - m
ft = fittype( 'exp1' );
excludedPoints = (xData < 90) | (xData > 120);
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [1 -0.1];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Set up fittype and options. - m1
ft = fittype( 'exp1' );
excludedPoints1 = (xData < 80) | (xData > 230);
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [1 -0.1];
opts.Exclude = excludedPoints1;

% Fit model to data.
fitresult1 = fit( xData, yData, ft, opts );

% Set up fittype and options. - m2
ft = fittype( 'exp1' );
excludedPoints2 = (xData < 90) | (xData > 110);
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [1 -0.1];
opts.Exclude = excludedPoints2;

% Fit model to data.
fitresult2 = fit( xData, yData, ft, opts );

% Set up fittype and options. - m3
ft = fittype( 'exp1' );
excludedPoints3 = (xData < 170) | (xData > 500);
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [1 -0.1];
opts.Exclude = excludedPoints3;

% Fit model to data.
fitresult3 = fit( xData, yData, ft, opts );

% Plot fit with data.
m = -(1/(fitresult.b));
m1 = -(1/(fitresult1.b));
m2 = -(1/(fitresult2.b));
m3 = -(1/(fitresult3.b));
figure( 'Name', 'Exponential fit for 5 data' );
h1 = plot( fitresult1,'c', xData, yData,'oc', excludedPoints1, '.b' );
hold on;
h2 = plot( fitresult2,'m', xData, yData,'dm', excludedPoints2, '.b' );
h3 = plot( fitresult3,'g', xData, yData,'xg', excludedPoints3, '.b' );
h = plot( fitresult, xData, yData,'.r', excludedPoints, '.b' );
legend([h(1) h(3) h1(3) h2(3) h3(3) h(2)], 'Fit FISHsignal(m)',sprintf('Fit m_1 = %g', m), sprintf('Fit m_1.1 = %g', m1), sprintf('Fit m_1.2 = %g', m2), sprintf('Fit m_2 = %g', m3), 'Excluded FISHsignal', 'Location', 'SouthWest', 'Interpreter', 'none' );
%legend('Fit FISHsignal(m)', 'Excluded FISHsignal', sprintf('Fit m = %g', m), 'Fit FISHsignal(m_1)', 'Excluded FISHsignal', sprintf('Fit m_1 = %g', m1), 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'fishTime', 'Interpreter', 'none' );
ylabel( 'fishSignal', 'Interpreter', 'none' );
title(sprintf("5' signal fitting, mu1 = %g, mu2 = %g", mRNALL, mRNALL2));
axis([0 1200 0 5]);
hold off;
grid on

