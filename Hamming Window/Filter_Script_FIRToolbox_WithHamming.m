% %% Obtain Data
 [fname1, pname1] = uigetfile('.xlsx'); %Prompt user for filename %change to txt/xlxs later
 filename = fullfile(pname1, fname1); %Create fully-formed filname as a string
 assert(exist(filename,'file')==2, '%s does not exist', filename); %Check that file exists
% 
 data = readtable(filename); %Read in the Data
% 
 channel = input('Enter the number of the channel: '); %Input 1st channel for which data was measured
 raw_data = data(:,channel); %change to data_wo_time/data later
 raw_data = table2array(raw_data); %omit later

%raw_data = xlsread('megboth.xlsx');
%raw_data = raw_data(:,1);

%% Plot Raw Data and FFT 
ws_5 = 256;
N_5 = length(raw_data);
dw_5 = 0:ws_5/N_5:ws_5;
t=0:1/256:(length(raw_data)-1)/256;
figure(1);
subplot(2,1,1);
plot(t, raw_data); title('Raw Data'); xlabel('Time(s)'); ylabel('Microvolts(mV)')

y = fft(raw_data);
figure(1);
subplot(2,1,2); suptitle('Index Finger; Channel 6')
dw_5(15346)=[];
plot(dw_5, y); title('FFT'); xlabel('Frequency(Hz)'); ylabel('Amplitude');

%% Filter Data: Bandstop

Hd = FIR_Bandstop_Toolbox_HammingWin
channel = step(Hd,raw_data);

t=0:1/256:(length(raw_data)-1)/256;
figure(2);
subplot(2,1,1);
plot(t,channel); title('Bandstop Filter Applied; Hamming Window'); xlabel('Time'); ylabel('Microvolts');

%% Filter Data: Bandpass

Hd = FIR_Bandpass_Toolbox_HammingWin
new = step(Hd,channel);

mean = sum(new)/length(new);
new1 = new - mean;
figure(2);
subplot(2,1,2);
plot(t,new1); title('Bandstop and Bandpass Filter Applied; Hamming Window'); xlabel('Time'); ylabel('Microvolts');
suptitle('Index Finger Data; Channel 6')





% subplot(2,1,1);
% plot(t,channel3);
% xlabel('time(sec)');
% ylabel('microvolts');
% title('EMG At Rest raw data (0) - Channel 3 ');
% subplot(2,1,2);
% plot(t,channel3y);
