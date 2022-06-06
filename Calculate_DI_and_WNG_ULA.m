phi = linspace(0,pi,100);
phi_rad = phi;
phi_mat = repmat(phi',[1,2*M+1]);
location_mat = repmat(mic_locations,[length(phi_rad),1]);

theta_vec = linspace(0,pi,100);
phi_vec = linspace(0,2*pi,100);
d_theta = theta_vec(2)-theta_vec(1);
d_phi = phi_vec(2)-phi_vec(1);
[T,P] = meshgrid(theta_vec,phi_vec);
alpha_x_mat = sin(T).*cos(P);
%% claculate broadsight
broad_theta = 0;
broad_phi = pi/2;
alpha_broad = sin(broad_theta).*cos(broad_phi);
d_f_broad = zeros(2*M+1,length(freq_vec)); 
for freq_bin = 1:length(freq_vec)
    f = freq_vec(freq_bin);
    d_f_broad(:,freq_bin) = exp(-1j*2*pi*f/c.*(mic_locations').*alpha_broad);
end
%%
total_gamma_mat = zeros(2*M+1,2*M+1,length(freq_bin));
for freq_bin = 1:length(freq_vec)
    f = freq_vec(freq_bin);
    total_sum = 0;
    for t_id = 1:length(theta_vec)
        curr_theta = theta_vec(t_id);
        for p_id = 1:length(phi_vec)
            curr_phi = phi_vec(p_id);
            alpha = sin(curr_theta)*cos(curr_phi);
            curr_d = exp(-1j*2*pi*f/c.*(mic_locations').*alpha);
            curr_mat = curr_d*curr_d';
            total_sum = total_sum + curr_mat.*sin(curr_theta);
        end
    end
    total_gamma_mat(:,:,freq_bin) = total_sum.*d_theta.*d_phi;

end


%%
per_freq_denum = zeros(1,length(freq_vec));
per_freq_num = zeros(1,length(freq_vec));

per_freq_denum_WNG = zeros(1,length(freq_vec));
per_freq_num_WNG = zeros(1,length(freq_vec));
% h_normalization = h_normalization./(max_vec');
for freq_bin = 1:length(freq_vec)
    w = h_normalization(freq_bin).*total_filter_mat(freq_bin,:);
    n = sum(w>0);
    curr_gamma = total_gamma_mat(:,:,freq_bin);
    
    curr_gamma_WNG = eye(2*M+1);

    per_freq_denum(freq_bin) = w*curr_gamma*w';
    per_freq_num(freq_bin) = w*d_f_broad(:,freq_bin);
    
        per_freq_denum_WNG(freq_bin) = w*curr_gamma_WNG*w';

    per_freq_num_WNG(freq_bin) = w*d_f_broad(:,freq_bin);
end
% per_freq_denum = abs(per_freq_denum);
per_freq_denum = (1/(4*pi)).*per_freq_denum;

%%
total_DF = (abs(per_freq_num).^2)./abs(per_freq_denum);

total_DF_WNG = (abs(per_freq_num_WNG).^2)./abs(per_freq_denum_WNG);

figure;
plot(freq_vec,10*log10(total_DF))

figure;
plot(freq_vec,10*log10(total_DF_WNG))