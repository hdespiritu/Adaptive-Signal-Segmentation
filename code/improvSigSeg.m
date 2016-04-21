%% My implementation of 
%
%An Improved Signal Segmentation Using Moving Average and Savitzky-Golay Filter
%Azami et al. 2012
%
%Inputs:
%   eegSignal: vector of EEG signal values
%   time: vector of time values
%   fs: sampling frequency
%Outputs:
%   savGolSegments: segment boundaries for Savitzky-Golay filtered signal
%   movAvgSegments: segment boundaries for Moving Average filtered signal
%
%% 
function [savGolSegments, movAvgSegments] = improvSigSeg(eegSignal, time, fs)
%Moving Avg Filter
% x(n) = 0.2[x(n-4) + x(n-3) + x(n-2) + x(n-1) + x(n)]
mafInput = tsmovavg(eegSignal.', 's', 5);%did transpose b/c tsmovavg only takes row vectors

    
%Savitzky-Golay Filter
sgfInput = sgolayfilt(eegSignal, 3, 51);
sgfInput = sgfInput.'; %did transpose to match dimensions of tsmovavg 

%Calculate the G variations for the filtered signal
G_movAvg = modifiedVarri(mafInput, fs);
G_savGol = modifiedVarri(sgfInput, fs);

%G threshold: mean(G)
threshold_movAvg = mean(G_movAvg, 'omitnan');
threshold_savGol = mean(G_savGol);

%Find local maxima
[maxima_movAvg, locMax_movAvg]= findpeaks(G_movAvg);
[maxima_savGol, locMax_savGol]= findpeaks(G_savGol);

%Scale window time frame back to signal time frame
gScaledTime = time(1:256:end);

%Local maxima above the threshold are segment boundaries
movAvgSegments = [];
savGolSegments = [];
for idx = 1:size(maxima_movAvg,2)
    if(maxima_movAvg(idx) > threshold_movAvg)
        movAvgSegments = [movAvgSegments gScaledTime(locMax_movAvg(idx))];
        
    end
end
for idx = 1:size(maxima_savGol,2)
    if(maxima_savGol(idx) > threshold_savGol)
        savGolSegments = [savGolSegments gScaledTime(locMax_savGol(idx))];
    end
end

end

%Calculates A_dif for a given window
function[A_dif] = adif(window)
A_dif = 0;
for idx = 1:size(window,2)
    A_dif = A_dif + abs(window(idx));
end
end

%Calculates F_dif for a given window
function [F_dif] = fdif(window)
F_dif = 0;
for idx = 2:size(window, 2)
    F_dif = F_dif + abs(window(idx) - window(idx - 1));
end
end

%Modified Varri Method
function [G_m] = modifiedVarri(signal, fs)
winLen = 2; %window length (seconds)
signalLen = size(signal,2)/fs; %total signal length (seconds)
numWin = ceil(signalLen/winLen); %Num of windows for signal
winNumElem = winLen*fs; %Num elements in the first numWin-1 windows
lastWinNumElem = size(signal,2) - ((numWin-1)*winNumElem); %Num elements in the last window

A1 = 7; F1 = 1; %Modified Varri constants
G_m = zeros(1,numWin);    %initialize measure difference function

%Calculate G_m
for win = 1:numWin-1

    %handles the last window case(not necessarily same len as other windows)
    if(win == numWin-1)
        winNumElemNext = size(signal,2)/(win+1);
    else
        winNumElemNext = winNumElem;
    end
    
    currWin = signal((win-1)*winNumElem + 1 : win*winNumElem);
    nextWin = signal(win*winNumElem + 1: (win+1)*winNumElemNext);
    G_m(win) = A1 * abs(adif(nextWin) - adif(currWin)) + F1 * abs(fdif(nextWin) - fdif(currWin));
end
end