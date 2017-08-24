%% Clear Workspace and Command Window 
clc;
clearvars;

%% Obtain Raw Data File
%[fname1, pname1] = uigetfile('.txt'); %Prompt user for filename %change to txt/xlxs later
%filename = fullfile(pname1, fname1); %Create fully-formed filname as a string
%assert(exist(filename,'file')==2, '%s does not exist', filename); %Check that file exists
filename = 'C:/Users/Aradhana/Downloads/Ga_label_1.txt';

data = readtable(filename); %Read in the Data
%temporarily omit

channel = input('Enter the number of the channel: '); %Input 1st channel for which data was measured
raw_data = data(:,channel); %change to data_wo_time/data later
raw_data.Properties.VariableNames = {'Temp'}; % Rename single column
raw_data = raw_data.Temp;
%raw_data = table2array(raw_data); %omit later

%% Implement DFT Routine on Raw Data to convert it from time domain --> frequency domain
t_5 = 0:1/250:(length(raw_data)-1)/250; %250-->Sampling Frequency
ws_5 = 250;%Hz
N_5 = length(raw_data);
dw_5 = 0:ws_5/N_5:ws_5;%divide by the #of points in data in order to get an evenly-spaced set of frequencies up to the sampling frequency
dw_5(19289)=[];

% for k=1:N_5-1
%     raw_data_freq(k+1)=0;
%         for n=1:N_5-1
%             raw_data_freq(k+1)=raw_data_freq(k+1)+raw_data(n+1)*(exp(-2*pi*1i/N_5))^(n*k);
%         end 
% end

    %OR can use fft() function to acheive the same result
    raw_data_freq = fft(raw_data);
    
%% Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
    P2 = abs(raw_data_freq/N_5);
    P1 = P2(1:N_5/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
%% Define the frequency domain f
    f = ws_5*(0:(N_5/2))/N_5;
    
%% Plot the single-sided amplitude spectrum P1, and other stuff 
figure(1);      
%subplot(3,1,1); suptitle(['Channel ' num2str(channel_2) ' General Information']);
plot(t_5,raw_data); 
xlabel('Time'); ylabel('MicroVolts'); title('EMG Data: Time Domain');
 
figure(2);
%subplot(3,1,2); 
plot(dw_5,abs(raw_data_freq)); xlabel('Frequency'); ylabel('|H(w)|'); title('Frequency Response');

figure(3);
%subplot(3,1,3);
plot(f,P1); xlabel('f(Hz)'); ylabel('|P1(f)|'); title('Single-Sided Spectrum X(t)');

%% %% Define Filter Order and Implement Zero-Padding. 
M=35; %order

data_processed = zeros(length(raw_data),1);

zero_pad = zeros(length(raw_data) + 2*M, 1);
zero_pad(M+1:end-M) = raw_data; 

%% Remove Bias: Make sure the average of data is zero

%b = isnan(raw_data);
%mean_data1 = mean(raw_data); 
meanValue = mean(raw_data,'omitnan');
raw_data = raw_data - meanValue;

%% Implement Bandpass Filter ; without using toolbox ; without hamming window
wc1 = 50*2*pi; %units: rad/s NOT Hz
wc2 = 125*2*pi;
sf = 250; 
delta_t = 1/sf; 
time = length(raw_data);

[data_processed] = FIR_Bandpass(wc1, wc2, sf, delta_t, M, time, raw_data, zero_pad, data_processed);

figure(4);
plot(t_5,data_processed); xlabel('Time'); ylabel('MicroVolts'); title('Bandpass filter;wc_l=0;wc_u=100');

%% Implement Bandstop FIR filter ; without using toolbox; without hamming window
M=35; 
zero_pad(M+1:end-M) = data_processed; 

wc1=55*2*pi;
wc2 = 65*2*pi;
sf = 250;
delta_t = 1/sf;
time = length(raw_data);

[data_processed] = FIR_Bandstop(wc1, wc2, sf, delta_t, M, time, raw_data, zero_pad, data_processed);

figure(5);
plot(t_5, data_processed); title('Bandstop filter;wc_l=55;wc_u=65'); xlabel('Time'); ylabel('Microvolts'); 
%suptitle(['Filtered Channel ' num2str(channel_1) ' Data']); hold off;

%% Feature Extraction 

avg = mean(data_processed) %average value
Y = rms(data_processed) % eqn: sqrt[(x-xavg)/n] --> done for each point
energy = sum(data_processed)^2 %energy carries by the signal 
%power_spec = periodogram(abs(raw_data)); %power spectral density
% OR abs(raw_data_freq).^2; power spectral density proportional to the absolute value squared of the DFT (calculated by fft())
%figure(6); dw_6 = 0:ws_5/length(power_spec):ws_5; dw_6(8193)= []; plot(dw_6,power_spec);
title('Power Spectral Density'); xlabel('Frequency(Hz)'); ylabel('Power (db/Hz)');

zero_cross = 0; %number of times your filtered data crosses zero-axis (How many time its amplitude is zero)
for i = 1:N_5-1
    if data_processed(i) > 0 && data_processed(i+1) < 0 
        zero_cross = zero_cross + 1;
    elseif data_processed(i) < 0 && data_processed(i+1) > 0 
        zero_cross = zero_cross + 1;
    end
end

variance = var(data_processed) %eqn: (1/N-1)*summation i = 0 to N-1 of [(x-xavg)^2] 

data_table = [avg Y energy zero_cross variance]'; 
%feature_matrix = horzcat(data_table, feature_matrix);
fprintf('Average Root Mean Square Zero-Crossing Variance\n')
fprintf('%5.0f\t%10.3f\t%13.3f\t%7.3f\t%7.3f\n' , data_table')

