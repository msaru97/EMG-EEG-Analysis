function [data_processed] = FIR_Bandstop(wc1, wc2, ~, delta_t, M, time, ~, zero_pad, data_processed)
%wc1:lower cutoff frquency (in rad/s)
%wc2:upper cutoff frquency (in rad/s)
%sf:sampling frequency
%delta_t:period
%M:filter order
%time:number of discrete points in raw data
%raw_data:data obtained from signal 

for n = 1:time 
    for k = -M:M
     h_k = ((((delta_t*wc1)/pi)* sinc(wc1*k*delta_t)) - (((delta_t*wc2)/pi)* sinc(wc2*k*delta_t))); %filter coefficient for low-pass filter
     y_k = zero_pad(n+M-k)*h_k; %output
     data_processed(n) = data_processed(n) + y_k;
    end
end
end