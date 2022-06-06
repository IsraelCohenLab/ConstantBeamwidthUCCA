%% parameters definition
c = 340;
bw_mag = db2mag(-3);
theta_BW = 30*pi/180;
phi_broadsight = pi/2;
angles_to_BW = linspace(0,theta_BW/2,20);
R_vec = [ 2.5,(5:5:25)]*1e-2;        % possible radii
R_vec = [0.2, 0.25];
num_of_points = 100;
phi = linspace(0,pi,num_of_points);
d_theta = phi(2)-phi(1);
%%
freq_vec = frequency_vec;
location_mat = repmat(R_vec,[length(phi),1]);
phi_mat = repmat(phi',[1,length(R_vec)]);

d_f = zeros(length(phi),length(R_vec),length(freq_vec));     % Actual beampattern
for freq_bin = 1:length(freq_vec)
    f = freq_vec(freq_bin);
        d_f(:,:,freq_bin) =(besselj(0,2*pi.*location_mat*f/c.*sin(phi_mat)));
%     d_f(:,:,freq_bin) = location_mat.*(besselj(0,2*pi.*location_mat*f/c.*sin(phi_mat)));
end

total_gamma_mat = zeros(length(R_vec),length(R_vec),length(freq_bin));

for freq_bin = 1:length(freq_vec)
    curr_df = d_f(:,:,freq_bin);
    total_sum = 0;
    
    for t_id = 1:length(phi)
        curr_phi = phi(t_id);
        d = curr_df(t_id,:)';
        d_mult = d*d';
        d_mult_theta = d_mult.*sin(curr_phi);
        total_sum = total_sum+d_mult_theta;
    end
    total_gamma_mat(:,:,freq_bin) = total_sum.*d_theta;
end
%%
very_low_f = frequency_vec(frequency_vec<low_frequencies_limit(2));

%%
max_g = 1;
total_W = [];
DI_vec = [];
for f_id = length(very_low_f):-1:2
    f = very_low_f(f_id);
    cvx_begin
        variable weight_vector(1,2)
        
        curr_gamma = total_gamma_mat(:,:,f_id);
        DI_denum = weight_vector*curr_gamma*weight_vector';
        beampattern_BW_value =  besselj(0,2*pi*f/c*R_vec*sin(theta_BW/2))*weight_vector';
        rings_Da_vec = [];
        for i = 1:length(angles_to_BW)
           rings_Da_vec = [rings_Da_vec , besselj(0, 2*pi*f/c*R_vec*sin(angles_to_BW(i)))*weight_vector' ] ;
        end
        
        minimize(DI_denum);
        
        subject to
        rings_Da_vec(2:end) <= rings_Da_vec(1:end-1);
        rings_Da_vec >= bw_mag*ones(1,length(angles_to_BW));
%         0 <= weight_vector
%         sum(weight_vector) == 1;
        0<= weight_vector(1) <= weight_vector(2);
%         weight_vector(1) <= max_g
    cvx_end
    total_W = [total_W ; weight_vector];
    max_g = weight_vector(1);
    DI_vec = [DI_vec, 2/DI_denum];
end
%%
norm_total_W = total_W./(repmat(total_W(:,2),[1,2]));
optimized_w = [0 , 1; flipud(norm_total_W)];
for i = 1:length(very_low_f)
    current_filter = zeros(1,M);
    current_filter(M-N_mic_vec(num_of_rings)+1:M) = optimized_w(i,2);
    current_filter(M-N_mic_vec(num_of_rings)-N_mic_vec(num_of_rings-1)+1:M-N_mic_vec(num_of_rings)) = optimized_w(i,1);
    total_weights_mat(i,:) = current_filter;
end