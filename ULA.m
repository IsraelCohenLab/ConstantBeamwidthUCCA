clear all; close all; clc
%% Define Parameters
N_bin = 75;         %Number of frequency bins in the upper band
L = 37;
M = 7;         %Number of microphones on one side of the ULA (there are 2*M+1)
d = 3.5*0.01;   %Microphone spacing [m]
c = 340;        %Speed of sound [m/s]
theta_BW = 30*pi/180;  %Desired beamwidth [rad]
Fs = 16*1e3;    %Sampling frequency [Hz]
Ts = 1/Fs;      %Delay between adjacent elements
desired_bw_point = db2mag(-3);

%% find minimal frequency of the relevant configuration
mic_locations = linspace(-d*M,d*M,2*M+1);       %microphones locations

Df=@(f) (abs(sum(cos(2*pi/c*sin(theta_BW/2)*f.*mic_locations))/(2*M+1)-desired_bw_point));        %Directivity function
f0 = fminsearch(Df,100);        %find minimal frequency at BW

%plot directivity function% 
figure;
fplot(Df,[0 8e3])
grid
xlabel('Frequency [Hz]')
ylabel('Directivity pattern for BW/2')
hold on
plot(f0,Df(f0),'ro')
hold off
%% find the limits for each geometry set
x = mic_locations(M+2:end);       %right side of mics (without middle)
f_low = 0;
frequency_limits = [];
%%%%% low frequency band - start from the outer most ring %%%%%
flipped_x_vec = fliplr(x);
for curr_allPass_mics = 1:length(x)-1
    % define directivity pattern depending on the number of rings
    directivity_function_handle = @(f) (abs((2*sum(cos(2*pi/c*sin(theta_BW/2)*f.*flipped_x_vec(1:curr_allPass_mics))))/(2*curr_allPass_mics))-desired_bw_point) ;
    % find the frequency for which the directivity handle function = 0 
     reference_frequency_vec = f_low:10000;
    directivity_value = cell2mat(arrayfun(directivity_function_handle,reference_frequency_vec ,'UniformOutput', false));
    plot(reference_frequency_vec, directivity_value)
    [~, lower_idx] = find(directivity_value < 0, 1 , 'first');
    [~, higher_idx] = find(directivity_value > 0, 1 , 'first');
    if isempty(lower_idx) || isempty(higher_idx)
        continue;
    end
    f_band = [reference_frequency_vec(lower_idx), reference_frequency_vec(higher_idx)];
    intersection_freq = bisection(directivity_function_handle,f_band(1),f_band(2));
% % %     intersection_freq = bisection(directivity_function_handle,f_low,f_high);
    frequency_limits = [frequency_limits, intersection_freq];
end

%%%%% higher frequency band - start from all the rings and reduce %%%%%
for curr_allPass_mics = length(x)-1:-1:1
    % define directivity pattern depending on the number of rings
    directivity_function_handle = @(f) (abs((1+2*sum(cos(2*pi/c*sin(theta_BW/2)*f.*x(1:curr_allPass_mics))))/(2*curr_allPass_mics+1))-desired_bw_point) ;
    % find the frequency for which the directivity handle function = 0 
     reference_frequency_vec = f_low:10000;
    directivity_value = cell2mat(arrayfun(directivity_function_handle,reference_frequency_vec ,'UniformOutput', false));
    plot(reference_frequency_vec, directivity_value)
    [~, lower_idx] = find(directivity_value < 0, 1 , 'first');
    [~, higher_idx] = find(directivity_value > 0, 1 , 'first');
    if isempty(lower_idx) || isempty(higher_idx)
        continue;
    end
    f_band = [reference_frequency_vec(lower_idx), reference_frequency_vec(higher_idx)];
    intersection_freq = bisection(directivity_function_handle,f_band(1),f_band(2));
% % %     intersection_freq = bisection(directivity_function_handle,0,8000);
    frequency_limits = [frequency_limits, intersection_freq];
end
%% Design filters for microphones
W = ones(N_bin,M);    % magnitudes of the frequency responses of the filters

%design filter for the middle microphone (All-Pass)% 
Hd_g0 = design_g0(L,Fs);        %descrete time filter object
[h_g0,w] = freqz(Hd_g0.Numerator,1,N_bin);  % w is frequency on [0,pi]  ;   h_g0 is the impulse responce
freq_vec = w/pi*Fs/2;       %from [0,pi] to [0,Fs/2]

%design filter for the rest of microphones (lower effective number of mics)%
start_ind = find(freq_vec>f0, 1, 'first');
x = mic_locations(M+2:end);       %right side of mics (without middle)
curr_allPass_mics = length(x)-1;      %at first iteration we look for the wheight of the right mic

for freq_bin = start_ind:N_bin
    disp(freq_bin)
    curr_freq = freq_vec(freq_bin);
    Df = @(a)  (abs((1+2*sum(cos(2*pi/c*sin(a)*curr_freq.*x(1:curr_allPass_mics))))/(2*curr_allPass_mics+1)-desired_bw_point)); % beampattern as a function of the angle phi
    a_min = fminsearch(Df,theta_BW/8);
    while a_min<theta_BW/2 && curr_allPass_mics>=1
        curr_allPass_mics = curr_allPass_mics-1;   % remove one mic to increas the beamwidth
        Df = @(a) (abs((1+2*sum(cos(2*pi/c*sin(a)*curr_freq.*x(1:curr_allPass_mics))))/(2*curr_allPass_mics+1) -desired_bw_point)); % beampattern as a function of the angle phi
        a_min = fminsearch(Df,theta_BW/8);
    end
    % now find the weight g for the (curr_allPass_mics+1)th mic
    Df=@(a,g) (abs((1+2*sum(cos(2*pi/c*sin(a)*curr_freq.*x(1:curr_allPass_mics)))+2*g*cos(2*pi/c*sin(a)*curr_freq.*x(1+curr_allPass_mics)))./(2*curr_allPass_mics+1+2*g) -desired_bw_point)); % beampattern as a function of the angle phi
    
    g=W(freq_bin-1,curr_allPass_mics+1);
    a_min = fminsearch(@(a) Df(a,g),theta_BW/8);
    while a_min<theta_BW/2 && g>=0
        g=g-0.01;
        a_min = fminsearch(@(a) Df(a,g),a_min);
    end
    W(freq_bin,curr_allPass_mics+1)=g;
    W(freq_bin,curr_allPass_mics+2:M)=0;
end
%%
% % % W(1:start_ind,:)=0;
% % % % W(1:start_ind,0.5*(M-1):M)=1;
% % % W(1:start_ind,M)=1;
% % % n_mics = 1;
% % % for i = start_ind-5:start_ind
% % %     W(i-2:i,end-n_mics+1:end) = 1;
% % %     if n_mics<M
% % %     W(i-2:i-1,end-n_mics) = 0.5;
% % %     end
% % %     n_mics = n_mics+1;
% % % end
% count = 1;
% for idx = 5:M
%     W(1:start_ind,idx)=0.5*count;
%     count = count+1;
% end
figure;
imagesc(1:M,freq_vec,W)
colorbar
xlabel('Filter index')
ylabel('Frequency [Hz]')
title('Frequency response magnitude')

total_W = [fliplr(W) ,ones(N_bin,1), W];
%% Calculate the impulse response of all the filters and normalization
%%%%% without weights %%%%%
% W = ones(N_bin,M); 
%%%%% without weights %%%%%

total_filter_mat = zeros(N_bin,2*M+1);      %each column is the impulse response of the corresponding mic
total_filter_mat(:,M+1) = h_g0;     %the middle microphone

for mic_id = 1:M
    Hd_g1 =design_g1(freq_vec,W(:,mic_id),Fs,L);
    h_g= freqz(Hd_g1.Numerator,1,w);
    total_filter_mat(:,mic_id+M+1) = h_g;    
    total_filter_mat(:,M-mic_id+1) = h_g;      %the symmetric mic 
end

%calculate normalization filter%
total_norm_coeff = 1./( sum(abs(total_filter_mat),2 ) );
Hd_normalization = design_g1(freq_vec,total_norm_coeff,Fs,2*L);
h_normalization= freqz(Hd_normalization.Numerator,1,w);
%
% total_filter_mat = total_W;
% h_normalization = 1./( sum(abs(total_filter_mat),2 ) );
figure;
imagesc(1:M,freq_vec,abs(total_filter_mat(: , 8:end)))
colorbar
xlabel('Filter index')
ylabel('Frequency [Hz]')
title('Frequency response magnitude')

%% Compute the actual beampattern of the broadside beam
phi = (-90:1:90)';
phi_rad = phi/180*pi;

phi_mat = repmat(phi_rad,[1,2*M+1]);
location_mat = repmat(mic_locations,[length(phi_rad),1]);


Da = zeros(length(phi),length(freq_vec));     % Actual beampattern

for freq_bin = 1:length(freq_vec)
    f = freq_vec(freq_bin);
    
    curr_freq_filter = total_filter_mat(freq_bin,:);
    curr_freq_filter_mat = repmat(curr_freq_filter,[length(phi_rad),1]);
    curr_phase_mat = curr_freq_filter_mat.*exp(-1j*2*pi*f/c*location_mat.*sin(phi_mat));
    
    Da(:,freq_bin) = abs(h_normalization(freq_bin)*sum(curr_phase_mat,2));
    
end
figure;
imagesc(phi,freq_vec,20*log10(Da'));
colormap jet
contourf(phi,freq_vec,20*log10(Da'),-18:3:-3);
imagesc(phi,freq_vec,20*log10(Da'));

grid
axis 'xy'
xlabel('\theta   [deg]')
ylabel('Frequency   [Hz]')
title('Actual Beampattern')
% % % % colorbar