% This code implements the constant beamwidth beamformers presented in the
% paper "Constant-Beamwidth Beamforming with Concentric Ring Arrays", published
% in Sensors. Please note the following before running the script:
%       1. To calculate the ULA beampattern - please run ULA.m
%       2. To calculate the original constant BW UCCA (LP) - set is_modified =
%       false , is_improved = false
%       3. To calculate the modified UCCA (BP) beampattern - set is_modified = true,
%       is_improved = false
%       4. To calculate the improved modified beampattern - set is_modified
%       = true, is_improved = true.
%       5. If you'd like to save the filters - set should_save_filters =
%       true
% !!!!! IMPORTANT NOTE !!!!!
%    The improved beamformer requires an external library called CVX. Please
%    download and install fron http://cvxr.com/cvx/download/ before running
%    the improved configuration
clear all; close all; clc;
%% Choose Desig specifications
is_modified = true;
is_improved = false;
should_save_filters = false;
%% Define constants
R_vec = [ 2.5,(5:5:25)]*1e-2;        % possible radii
theta_BW = 30*pi/180;               % Desired beamwidth [rad]
f_low = 0;                          % lower frequency limit [Hz]
f_high = 8000;                      % high frequency limit [HZ]
desired_bw_point = db2mag(-3);      % the -3dB point defines the beam-width
c = 340;
num_of_rings = length(R_vec);


BW = theta_BW;     % mainlobe width in rad
L = 37;      % length of the filters (the length of the normalization filter is 2L)
Fs = 16000;  % Sampling Frequency

G=ones(L,num_of_rings);    % magnitudes of the frequency responses of the filters

Hd_g0 = design_g0(L,Fs);
[h_g0,w] = freqz(Hd_g0.Numerator,1,75);  % w is frequency on [0,pi]
coefs=Hd_g0.Numerator;
fw=w/pi*Fs/2;   % frequencies for designing the filters
frequency_vec = fw;

%% check discrete setting
N_mic_vec = [16,16,16*ones(1,num_of_rings-2)];
N_mic_vec = [16*ones(1,num_of_rings)];
M = sum(N_mic_vec);

mic_pos_mat = CalculateArrayGeometry(R_vec,N_mic_vec);
mic_locations = mic_pos_mat(:,2)';

%%%%% low frequency band - start from the outer most ring %%%%%
low_frequencies_limit = [];
flipped_num_mics = fliplr(N_mic_vec);
flipped_mic_locations = fliplr(mic_locations);
for ring_num = 1:num_of_rings
    num_of_mics = sum(flipped_num_mics(1:ring_num));
    active_mic_locations = flipped_mic_locations(1:num_of_mics);
    % define directivity pattern depending on the number of rings
    directivity_function_handle = @(f) ( (1/num_of_mics)*sum(  exp( -1j*2*pi/c*f*sin(theta_BW/2).*active_mic_locations  ) )  - desired_bw_point);
    reference_frequency_vec = f_low:10000;
    directivity_value = cell2mat(arrayfun(directivity_function_handle,reference_frequency_vec ,'UniformOutput', false));
    plot(reference_frequency_vec, real(directivity_value))
    [~, lower_idx] = find(directivity_value < 0, 1 , 'first');
    [~, higher_idx] = find(directivity_value > 0, 1 , 'first');
    if isempty(lower_idx) || isempty(higher_idx)
        continue;
    end
    f_band = [reference_frequency_vec(lower_idx), reference_frequency_vec(higher_idx)];
    intersection_freq = bisection(directivity_function_handle,f_band(1),f_band(2));

    low_frequencies_limit = [low_frequencies_limit, intersection_freq];
end
%%%%% higher frequency band - start from all the rings and reduce %%%%%
high_frequencies_limit = [];

for ring_num = num_of_rings:-1:1
    num_of_mics = sum(N_mic_vec(1:ring_num));
    active_mic_locations = mic_locations(:,1:num_of_mics);

    % define directivity pattern depending on the number of rings
    directivity_function_handle = @(f) ( (1/num_of_mics)*sum(  exp( -1j*2*pi/c*f*sin(theta_BW/2).*active_mic_locations  ) )  - desired_bw_point);
    
    % find the frequency for which the directivity handle function = 0 
    reference_frequency_vec = f_low:10000;
    directivity_value = cell2mat(arrayfun(directivity_function_handle,reference_frequency_vec ,'UniformOutput', false));
    plot(reference_frequency_vec, real(directivity_value))
    [~, lower_idx] = find(directivity_value < 0, 1 , 'first');
    [~, higher_idx] = find(directivity_value > 0, 1 , 'first');
    if isempty(lower_idx) || isempty(higher_idx)
        continue;
    end
    f_band = [reference_frequency_vec(lower_idx), reference_frequency_vec(higher_idx)];
    intersection_freq = bisection(directivity_function_handle,f_band(1),f_band(2));
    high_frequencies_limit = [high_frequencies_limit, intersection_freq];
end
%%
frequency_vec = fw;
%%%%% lower frequency band %%%%%
low_frequency_vec = frequency_vec(frequency_vec < high_frequencies_limit(1));
low_frequencies_weights = zeros(length(low_frequency_vec),sum(N_mic_vec));
for f_bin = 1:length(low_frequency_vec)
    f = low_frequency_vec(f_bin);
    current_freq_band = find(f>low_frequencies_limit,1,'last');
    if isempty(current_freq_band)
        % only the outer ring is active %
        low_frequencies_weights(f_bin,1:flipped_num_mics(1)) = 1;
        continue;
    end
    % define mic number and locations %
    num_of_active_mics = sum(flipped_num_mics(1:current_freq_band));
    num_of_weighted_mics = flipped_num_mics(current_freq_band+1);
    active_mic_locations = flipped_mic_locations(1:num_of_active_mics);
    weighted_mic_locations = flipped_mic_locations(num_of_active_mics+1: num_of_active_mics+ num_of_weighted_mics);
    % define directivity function and find the corresponding weight %
    active_directivity = sum(  exp( -1j*2*pi/c*f*sin(theta_BW/2).*active_mic_locations  ) );
    weighted_deirectivity = sum(  exp( -1j*2*pi/c*f*sin(theta_BW/2).*weighted_mic_locations  ) );
    current_alpha = ( desired_bw_point*num_of_active_mics -  active_directivity) / ...
                    ( weighted_deirectivity -  desired_bw_point*num_of_weighted_mics );
    
    low_frequencies_weights(f_bin,1:num_of_active_mics) = 1;
    low_frequencies_weights(f_bin,num_of_active_mics+1: num_of_active_mics+ num_of_weighted_mics) = abs(current_alpha);
end
% flip the weight matrix %
low_frequencies_weights = fliplr(low_frequencies_weights);

%%%%% higher frequency band %%%%%
high_frequency_vec = frequency_vec(frequency_vec >= low_frequencies_limit(end));
high_frequencies_weights = zeros(length(high_frequency_vec),sum(N_mic_vec));
for f_bin = 1:length(high_frequency_vec)
    f = high_frequency_vec(f_bin);
    current_freq_band = find(f>high_frequencies_limit,1,'last');
    if isempty(current_freq_band)
        % all ring are active %
        high_frequencies_weights(f_bin,:) = 1;
        continue;
    end
    total_active_rings = num_of_rings - current_freq_band;
    % define mic number and locations %
    num_of_active_mics = sum(N_mic_vec(1:total_active_rings));
    num_of_weighted_mics = N_mic_vec(total_active_rings+1);
    active_mic_locations = mic_locations(1:num_of_active_mics);
    weighted_mic_locations = mic_locations(num_of_active_mics+1: num_of_active_mics+ num_of_weighted_mics);
    % define directivity function and find the corresponding weight %
    active_directivity = sum(  exp( -1j*2*pi/c*f*sin(theta_BW/2).*active_mic_locations  ) );
    weighted_deirectivity = sum(  exp( -1j*2*pi/c*f*sin(theta_BW/2).*weighted_mic_locations  ) );
    current_alpha = ( desired_bw_point*num_of_active_mics -  active_directivity) / ...
                    ( weighted_deirectivity -  desired_bw_point*num_of_weighted_mics );
    
    high_frequencies_weights(f_bin,1:num_of_active_mics) = 1;
    high_frequencies_weights(f_bin,num_of_active_mics+1: num_of_active_mics+ num_of_weighted_mics) = abs(current_alpha);
end

if is_modified
    total_weights_mat = [low_frequencies_weights ; high_frequencies_weights];
else
    total_weights_mat = [ones(size(low_frequencies_weights)) ; high_frequencies_weights];
end
%% Add lower frequency optimization
if is_improved
    optimize_DI;
end
clc;
%% design time domain filters
G = total_weights_mat(:,1:16:end);
% find the FIR impulse responses of the filters
Hd_g1 =design_g1(fw,G(:,1),Fs,L);
Hd_g2 =design_g1(fw,G(:,2),Fs,L);
Hd_g3 =design_g1(fw,G(:,3),Fs,L);
Hd_g4 =design_g1(fw,G(:,4),Fs,L);
Hd_g5 =design_g1(fw,G(:,5),Fs,L);
Hd_g6 =design_g1(fw,G(:,6),Fs,L);

h_g1= freqz(Hd_g1.Numerator,1,w);
h_g2= freqz(Hd_g2.Numerator,1,w);
h_g3= freqz(Hd_g3.Numerator,1,w);
h_g4= freqz(Hd_g4.Numerator,1,w);
h_g5= freqz(Hd_g5.Numerator,1,w);
h_g6= freqz(Hd_g6.Numerator,1,w);

A1= 1./ (16*(abs(h_g1)+ abs(h_g2)+ abs(h_g3)+ abs(h_g4)+ abs(h_g5)+ abs(h_g6)));
Hd_normalization = design_g1(fw,A1,Fs,16*L);
h_normalization= freqz(Hd_normalization.Numerator,1,w);
%% save filters
if should_save_filters
    if is_improved
        main_path = 'FilterValues\Improved_Modified\';
    elseif is_modified
        main_path = 'FilterValues\Modified\';
    else
        main_path = 'FilterValues\Low_Pass\';
    end
    
    coefs=Hd_g1.Numerator;
    fid=fopen([main_path,'filter_g1.txt'],'wt'); fprintf(fid,'%f\n',coefs); fclose(fid);

    coefs=Hd_g2.Numerator;
    fid=fopen([main_path,'filter_g2.txt'],'wt'); fprintf(fid,'%f\n',coefs); fclose(fid);

    coefs=Hd_g3.Numerator;
    fid=fopen([main_path,'filter_g3.txt'],'wt'); fprintf(fid,'%f\n',coefs); fclose(fid);

    coefs=Hd_g4.Numerator; % take into acoount the division by 2
    fid=fopen([main_path,'filter_g4.txt'],'wt'); fprintf(fid,'%f\n',coefs); fclose(fid);

    coefs=Hd_g5.Numerator;
    fid=fopen([main_path,'filter_g5.txt'],'wt'); fprintf(fid,'%f\n',coefs); fclose(fid);

    coefs=Hd_g6.Numerator; % take into acoount the division by 2
    fid=fopen([main_path,'filter_g6.txt'],'wt'); fprintf(fid,'%f\n',coefs); fclose(fid);

    coefs=Hd_normalization.Numerator;
    fid=fopen([main_path,'filter_gt.txt'],'wt'); fprintf(fid,'%f\n',coefs); fclose(fid);
end

%% Combine filters into a matrix
ideal_weights = total_weights_mat;
total_weights_mat = [ repmat(h_g1 , 1 , 16) , ...
    repmat(h_g2 , 1 , 16) , ...
    repmat(h_g3 , 1 , 16) , ...
    repmat(h_g4 , 1 , 16) , ...
    repmat(h_g5 , 1 , 16) , ...
    repmat(h_g6 , 1 , 16) ];
%% plot weights mat
figure;
imagesc(1:16,frequency_vec,abs(total_weights_mat(:,1:16:end)))
axis xy
colorbar
ax = gca;
new_label = arrayfun( @(x) num2str(x/1000) , ax.YTick, 'UniformOutput', false );
ax.YTickLabel = new_label;
ax.XTick = 1:3:16;
new_label = arrayfun( @(x) num2str(x) , 1:6, 'UniformOutput', false );
ax.XTickLabel = new_label;
caxis([0 , 1])

%% Calculate preformance measures
Calculate_DI_and_WNG_discrete_CCA;
%% Plot actual beampattern
phi=(-90:1:90)';
phi_rad=phi/180*pi;

Da=zeros(length(phi),length(fw));% Actual beampattern
for k=1:length(fw)
    f=fw(k);
    Da(:,k) = h_g6(k)  * sum(  exp( -1j*2*pi/c*f*sin(phi_rad).*mic_locations(16*5 + 1 :end)  ),2 ) + ...
        h_g5(k) *sum(  exp( -1j*2*pi/c*f*sin(phi_rad).*mic_locations(16*4 + 1 :16*5)  ),2 ) +...
        h_g4(k) * sum(  exp( -1j*2*pi/c*f*sin(phi_rad).*mic_locations(16*3 + 1 :16*4)  ) ,2)  +...
        h_g3(k)  * sum(  exp( -1j*2*pi/c*f*sin(phi_rad).*mic_locations(16*2 + 1 :16*3)  ),2 )  +...
        h_g2(k)  * sum(  exp( -1j*2*pi/c*f*sin(phi_rad).*mic_locations(16 + 1 :16*2)  ) ,2)  +...
        h_g1(k) * sum(  exp( -1j*2*pi/c*f*sin(phi_rad).*mic_locations(1 :16)  ),2 ) ;
    Da(:,k)=abs(Da(:,k)*h_normalization(k));
end

ideal_Da=zeros(length(phi),length(fw));% Actual beampattern
ideal_h_normalization = 1./( sum(abs(ideal_weights),2 ) );
phi_mat = repmat(phi_rad,[1,sum(N_mic_vec)]);
location_mat = repmat(mic_locations,[length(phi_rad),1]);
for freq_bin = 1:length(frequency_vec)
    f = fw(freq_bin);
    
    curr_freq_filter = ideal_weights(freq_bin,:);
    curr_freq_filter_mat = repmat(curr_freq_filter,[length(phi_rad),1]);
    curr_phase_mat = curr_freq_filter_mat.*exp(-1j*2*pi*f/c*location_mat.*sin(phi_mat));
    
    ideal_Da(:,freq_bin) = abs(ideal_h_normalization(freq_bin)*sum(curr_phase_mat,2));
    
end

figure;
imagesc(phi,fw,20*log10(Da'));
axis 'xy'
ax = gca;
new_label = arrayfun( @(x) num2str(x/1000) , ax.YTick, 'UniformOutput', false );
ax.YTickLabel = new_label;
colorbar;
colormap jet
caxis([-40 0])
%%%%%%
new_label = arrayfun( @(x) num2str(x/1000) , ax.YTick, 'UniformOutput', false );
ax.YTickLabel = new_label;
v = [-3,-3];
hold on;
contour(phi,fw,20*log10(Da'),v,'LineColor','white','LineWidth',4,'LineStyle',':')

%% get BW at each frequency
line_width = 2;
line_style = '--';
marker = '*';
line_color = 'b';

difference_mat = abs(Da-desired_bw_point);
[min_points_values,min_points_id] = min(difference_mat,[],1);
figure;
plot(frequency_vec,2*abs(phi(min_points_id)),'LineStyle',line_style,'Color',line_color,'LineWidth',line_width,'Marker',marker)