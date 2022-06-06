function [total_pos_mat] = CalculateArrayGeometry(Radii_vec,N_mic_vec)
%Recieve the microphone positions acoording to the total array
%configuration
%Input:
%   N_mic_vec --> Number of microphones at each ring
%   Radii_vec --> Radius of each ring (usually uniform differences)
%Output:
%   mic_pos_mat --> a NumOfMic x 2 matrix where each row is the xy position
%                   of the microphone

total_num_of_mic = sum(N_mic_vec);
total_pos_mat = zeros(total_num_of_mic,2);      %Assuming the array lays on the xy plane

curr_ring_id = 1;
curr_ring_mic_count = 0;
for mic_id = 1:total_num_of_mic
    if curr_ring_mic_count >= N_mic_vec(curr_ring_id)
        curr_ring_id = curr_ring_id+1;
        curr_ring_mic_count = 0;

    end
    curr_ring_mic_count = curr_ring_mic_count+1;
    curr_ring_radi = Radii_vec(curr_ring_id);       %Rh the radius of the ring
    if mod(N_mic_vec(curr_ring_id),2) == 0 
        curr_mic_phase = 2*pi*(curr_ring_mic_count)/N_mic_vec(curr_ring_id);      %Phi_hk = 2*pi*(m-1)/M
    else
        curr_mic_phase = 2*pi*(curr_ring_mic_count-1)/N_mic_vec(curr_ring_id);      %Phi_hk = 2*pi*(m-1)/M
    end
    total_pos_mat(mic_id,1) = curr_ring_radi*cos(curr_mic_phase);
    total_pos_mat(mic_id,2) = curr_ring_radi*sin(curr_mic_phase);
end
   
    
       


end

