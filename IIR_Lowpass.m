function [data_processed] = IIR_Lowpass(wc, sf, ~, N, ~, raw_data, zero_pad, data_processed)
%% Background 
    %IIR filters are causal filters: filter output depends only on past and present inputs 
    
%% Variables used in function
%wc:cutoff frquency
%sf:sampling frequency
%delta_t:period
%N:filter order
%raw_data:data obtained from signal
%
    
%% First generate the lowpass analogue Butterworth prototypes
[zs,ps,ks]=buttap(N);
%buttap() --> Creates an Mth order normalized Butterworth analog lowpass filter prototype 
%M --> Filter Order
%zs --> zeroes of filter prototype
%ps --> poles of filter prototype 
%ks --> gain of filter prototype 

%% Put the filter prototype in terms of a transfer function in the s-domain
[nums, dens] = zp2tf(zs,ps,ks);
%zp2tf() --> forms the transfer function H(s) = NUM(s)/DEN(s) given a set of zero locations in vector Z
%, a set of pole locations in vector P, and a gain in scalar K
%nums, dens --> returned with numerator and denominator coefficients in descending powers of s

%% Transform the transfer function to an analog lowpass filter with cut-off frequency wc
[nums, dens] = lp2lp(nums,dens,wc*2*pi);
%Lowpass to lowpass analog filter transformation. Transforms the lowpass filter prototype NUM(s)/DEN(s)
%with unity cutoff frequency of 1 rad/sec to a lowpass filter with cutoff frequency Wo (rad/sec)
%nums, dens --> returned with numerator and denominator coefficients in descending powers of s

%% Use Bilinear transformation to transform the analog lowpass filter from the s-domain to the z-domain
[numzb, denzb] = bilinear(nums, dens, sf);
%bilinear() --> where nums and dens are row vectors containing numerator and denominator transfer function 
%coefficients, NUM(s)/DEN(s), indescending powers of s, transforms to z-transform coefficients NUMd(z)/DENd(z) 
%sf --> sampling frequency 

%% Divide the coefficients of the numerator and denominator by first denominator coefficient in order to yeild a canonical form for the filter in accordance with Equation (11) in Barr Chan paper

numzb = numzb/denzb(1);
denzb = denzb/denzb(1);

%% Design the recursion formula with vectors numzb and denzb in the time-domain and calculate the output; Equation (30) in Barr Chan Paper
k = length(numzb);
k1 = length(denzb);

m=0; %both for loop below implement Equation (30) in Barr Chan paper 
for l = 1:1:k  
    for n=1:length(raw_data)
    
    sub1=n-m;
    
    if sub1>0
       x = zero_pad(sub1);
    end
   
    if sub1<=0 
        x=0;   
    end
   
    data_processed(n) = data_processed(n) + (numzb(l)*x);

    end
    m=m-1;
end

m=1;

for l = 1:1:k1
    for n=1:length(raw_data)
        sub2=n-m;
         
    if sub2>0
      x1 = data_processed(sub2);
      
    if sub2<=0
      x1=0;
    end
    
    data_processed(n) = data_processed(n) - (denzb(l)*x1);
    
    end
    m=m+1;  
    end
    
end

