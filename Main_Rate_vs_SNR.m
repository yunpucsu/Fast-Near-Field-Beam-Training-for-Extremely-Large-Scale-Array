clc; clear all; close all;

%% System Parameters
M = 256; % number of antennas
f = 100e9; % 100Ghz
c = 3e8;
lambda = c/f;
d = lambda/2;
Rayleigh_dist = 2*(M*d)^2/lambda;
Fresnel_dist = 0.62*sqrt((M*d)^3/lambda);
noise_pow = 1e-10;
beta_ref = (lambda/(4*pi))^2;

num_samples = 10;
max_rounds = 100;

SNR_dB = [21:2:39];  
SNR_linear = 10.^(SNR_dB/10.);
Num_SNR = length(SNR_linear);



theta = -1 + 2/M : 2/M : 1;

Min_distance = 5; Max_distance = 0.1*Rayleigh_dist;
Lef_angle = -0.3; right_angle = 0.3;

%% exhaustive search based codebook generation
rho_min = 3; rho_max = 64; beta = 1.2; 
[Un, label] = PolarCodeBook(M, d, lambda, beta, rho_min, rho_max);

%% main loop
         
mean_perfect_rate = zeros(1,Num_SNR);
mean_proposed_rate_K3 = zeros(1,Num_SNR);
mean_proposed_rate_K1 = zeros(1,Num_SNR);
mean_ES_polar_rate = zeros(1,Num_SNR);
mean_far_field_rate = zeros(1,Num_SNR);
mean_LS_rate = zeros(1,Num_SNR);
for iter=1:max_rounds
 iter
rate_achievable = zeros(num_samples,Num_SNR);  
proposed_rate_K3= zeros(num_samples,Num_SNR);
proposed_rate_K1= zeros(num_samples,Num_SNR);
uni_rate_achievable= zeros(num_samples,Num_SNR);
rate_real = zeros(num_samples,Num_SNR);
LS_rate = zeros(num_samples,Num_SNR);
far_field_rate = zeros(num_samples,Num_SNR);
ES_polar_rate = zeros(num_samples,Num_SNR);

% randomly generate 10 user angle-distance pair
user_angle_pro = zeros(1,num_samples);
user_distance = zeros(1,num_samples);
for i =1:num_samples
    user_angle_pro(1,i) = Lef_angle + (right_angle-Lef_angle)*rand;
    user_distance(1,i) = Max_distance-rand*(Max_distance-Min_distance);
end
for t = 1:num_samples
    user_angle_sca = user_angle_pro(1,t);
    user_angle = asin(user_angle_sca);
    user_loc = user_distance(1,t);
    alpha_n=sqrt(beta_ref)/user_loc;
  
    hn_real_ini = near_field_manifold(M,d,f,user_loc,user_angle);
    hn_real = exp(-1j*2*pi*f*user_loc/c)*alpha_n*hn_real_ini*sqrt(M);
    an_real = near_field_manifold(M,d,f,user_loc,user_angle)';
    
    for snr = 1:Num_SNR
        trans_power = (SNR_linear(snr)*user_loc^2*noise_pow)/(beta_ref*M);
        
        % First_phase: angle estimation
        three_choose_index = zeros(1,3);
        angle_snr = zeros(1,M);
        noise = sqrt(noise_pow)*(randn(1)+1i*randn(1))/sqrt(2);
        for m = 1:M
            af = far_field_mainfold(M,d,f,asin(theta(m))); 
            yn = sqrt(trans_power)*hn_real*af + noise;
            angle_snr(1,m) = norm(yn)^2;
    
        end
        [ang_max_gain, ang_max_index] = max(angle_snr);   
        indices = find(angle_snr > (sqrt(2)/2)^2 * (ang_max_gain)); 
        % find the middle location
        med_idx = median(indices);
        if med_idx == fix(med_idx)   
            est_angle = theta(1,med_idx);
            three_choose_index(1,1) = theta(1,med_idx);
            three_choose_index(1,2) = theta(1,med_idx-1);
            three_choose_index(1,3) = theta(1,med_idx+1);
        else
            est_angle = theta(1,floor(med_idx));
            three_choose_index(1,1) = theta(1,floor(med_idx));
            three_choose_index(1,2) = theta(1,ceil(med_idx));
            three_choose_index(1,3) = theta(1,floor(med_idx)-1);
        end
    
        % Second_phase: distance estimation 
        % distance estimation for one candidate angle
        [num_dist_sample,dist_sample_set] = generate_dist_samples(est_angle,label);
        dist_snr_profile = zeros(1,num_dist_sample);
        for k = 1:num_dist_sample
            an = near_field_manifold(M, d, f, dist_sample_set(k), asin(est_angle))';
            y_eff = sqrt(trans_power)*hn_real*an + noise;
            dist_snr_profile(1,k) = norm(y_eff)^2;
        end
        [~,idx_max] = max(dist_snr_profile);
        est_distance = dist_sample_set(idx_max);
        
        % distance estimation for three candidate angles
        dist_snr_profile_K3 = zeros(1,length(three_choose_index));
        est_distance_set = zeros(1,length(three_choose_index));

        for jj = 1:length(three_choose_index)
            cur_angle = three_choose_index(jj);
            [num_dist_sample_K3,dist_sample_set_K3] = generate_dist_samples(cur_angle,label);
            dist_snr_profile_K3_tmp =zeros(1,num_dist_sample_K3);
            
            for mm = 1:num_dist_sample_K3
                an = near_field_manifold(M, d, f, dist_sample_set_K3(mm), asin(cur_angle))';
                yf = sqrt(trans_power)*hn_real*an+noise;
                dist_snr_profile_K3_tmp(1,mm)=norm(yf)^2;
            end
            [max_value,max_idx] = max(dist_snr_profile_K3_tmp); 
            dist_snr_profile_K3(1,jj) = max_value; 
            est_distance_set(1,jj) = dist_sample_set_K3(max_idx); 
        end
        [~,final_max_idx] = max(dist_snr_profile_K3);
        est_distance_K3 = est_distance_set(final_max_idx);
        est_angle_K3 = three_choose_index(final_max_idx);
   
    
        %% achievable rate
        % Proposed two-phase beam training (K=3)
        estimated_bf_K3 = near_field_manifold(M,d,f,est_distance_K3,asin(est_angle_K3))';
        y_est_eff_K3 = sqrt(trans_power)*hn_real*estimated_bf_K3;
        estimated_bf_gain_K3 = norm(y_est_eff_K3)^2;
        proposed_rate_K3(t,snr) = log2(1+estimated_bf_gain_K3/noise_pow);


        % Proposed two-phase beam training (K=1)
        estimated_bf = near_field_manifold(M,d,f,est_distance,asin(est_angle))';
        y_est_eff = sqrt(trans_power)*hn_real*estimated_bf;
        estimated_bf_gain = norm(y_est_eff)^2;
        proposed_rate_K1(t,snr) = log2(1+estimated_bf_gain/noise_pow);


        % Perfect CSI based beamforming
        real_gain = norm(sqrt(trans_power)*hn_real*an_real)^2;
        rate_real(t,snr) = log2(1+real_gain/noise_pow);
    

        % LS channel estimation 
        DFT = (1/sqrt(M))*exp(-1i*pi*[0:M-1]'*[-(M-1)/2:1:(M/2)]*(2/M));
        y_LS = sqrt(trans_power)*hn_real*DFT+sqrt(noise_pow)*(randn(1,M)+1i*randn(1,M))/sqrt(2);
        h_LS1 = y_LS*inv(DFT);
        an_LS = exp(1j*phase(h_LS1'))/sqrt(M);
        y_LS_eff = sqrt(trans_power)*hn_real*an_LS;
        LS_gain = norm(y_LS_eff)^2;
        LS_rate(t,snr) = log2(1+LS_gain/noise_pow);
    

        % far-field
        far_beamf = far_field_mainfold(M,d,f,asin(est_angle));
        y_far = sqrt(trans_power)*hn_real *far_beamf;
        far_field_rate_gain = norm(y_far)^2;
        far_field_rate(t,snr) = log2(1+far_field_rate_gain/noise_pow);
          

        % Exhaustive-search over polar-domain codebook
        S = size(Un, 2);
        ES_polar_gain = zeros(1,S);
        for s=1:S
        cur_beamform_vec = Un(:,s);
        cur_beamform_vec = conj(cur_beamform_vec);
        y_polar = sqrt(trans_power)*hn_real*cur_beamform_vec;
        ES_polar_gain(1,s) = norm(y_polar)^2;
        end
        [ES_polar_max_gain,~] = max(ES_polar_gain);
        ES_polar_rate(t,snr) = log2(1+ES_polar_max_gain/noise_pow);
    end   
end
mean_perfect_rate = mean_perfect_rate + mean(rate_real);
mean_proposed_rate_K3 = mean_proposed_rate_K3 + mean(proposed_rate_K3);
mean_proposed_rate_K1 = mean_proposed_rate_K1 + mean(proposed_rate_K1);
mean_ES_polar_rate = mean_ES_polar_rate + mean(ES_polar_rate);
mean_far_field_rate = mean_far_field_rate + mean(far_field_rate);
mean_LS_rate = mean_LS_rate + mean(LS_rate);
end



%% Plot achievable rate
figure;

plot(SNR_dB,mean_perfect_rate/max_rounds,'k--','linewidth',1.6); hold on;
plot(SNR_dB,mean_proposed_rate_K3/max_rounds,'ro-','linewidth',1.6); hold on;
plot(SNR_dB,mean_proposed_rate_K1/max_rounds,'+-.','color',[0.39 0.83 0.07],'linewidth',1.6); hold on;
plot(SNR_dB,mean_ES_polar_rate/max_rounds,'b--','linewidth',2); hold on;
plot(SNR_dB,mean_far_field_rate/max_rounds,'m^-','linewidth',1.6); hold on;
plot(SNR_dB,mean_LS_rate/max_rounds,'s-','color',[0.8500 0.3250 0.0980],'linewidth',1.6); hold on;
grid on;

legend("Perfect CSI based beamforming","Proposed two-phase beam training (K=3)","Proposed two-phase beam training (K=1)","Exhaustive-search based near-field beam trainin","Far-field beam training","LS channel estimation",'location','SouthEast')
xlabel('SNR (dB)');
ylabel('Achievable Rate (bit/s/Hz)');
xlim([SNR_dB(1) SNR_dB(end)])