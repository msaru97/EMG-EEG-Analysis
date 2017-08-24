function [data_processed] = FIR_Highpass(wc, ~, delta_t, M, time, ~, zero_pad, data_processed)
%wc:cutoff frquency (in rad/s)
%sf:sampling frequency
%delta_t:period
%M:filter order
%time:number of discrete points in raw data
%raw_data:data obtained from signal 

for n = 1:time
    for k = -M:M
     h_k = ((((-delta_t*wc)/pi)* sinc(wc*k*delta_t))); %filter coefficient for low-pass filter
     y_k = zero_pad(n+M-k)*h_k; %output
     data_processed(n) = data_processed(n) + y_k;
    end
end
end