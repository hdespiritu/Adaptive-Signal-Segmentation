%% This program plots the signal segments computed for the sample EEG data


%Import patient data
patientData = importdata('eegData.dat');
n = length(patientData);
Time = patientData(1:n,1);
EEG1 = patientData(1:n,2);
EEGSignal1_42sec = EEG1(1:5400);


%Perform adpative signal segmentation
[savGolSegments, movAvgSegments] = improvSigSeg(EEGSignal1_42sec, Time, 128);


%Plot original EEG signal
figure; 
hold;
plot(Time(1:5400),EEGSignal1_42sec);
xlabel('Time(seconds)');
ylabel('Amplitude');
title('Original EEG Signal');


%Plot moving average filtered EEG signal
figure; 
hold; 
plot(Time(1:5400), tsmovavg(EEGSignal1_42sec.','s',5)); 
xlabel('Time(seconds)');
ylabel('Amplitude');
title('Moving Average Filtered Signal');
for idx = 1:size(movAvgSegments,2)
plot([movAvgSegments(idx) movAvgSegments(idx)], [-0.2 0.2],'r');
end


%Plot Savitzky-Golay filtered EEG signal
figure; 
hold;
plot(Time(1:5400), sgolayfilt(EEGSignal1_42sec,3,51));
title('Savitzky-Golay Filtered Signal');
xlabel('Time(seconds)');
ylabel('Amplitude');
for idx = 1:size(savGolSegments,2)
plot([savGolSegments(idx) savGolSegments(idx)], [-0.2 0.2],'r');
end
