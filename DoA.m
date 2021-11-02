clc;
clear all;
close all;

%spcom_10 for dataset with 10 degree seperation.
load('spcom_50sep');

Rx = cov(X.');
theta_deg = (0:0.01:180).';
theta_rad = 2*pi/360*theta_deg;
res = length(theta_deg);
dist = (0:M-1)*Delta;
cand_CB = zeros(1,res);
cand_MVDR = zeros(1,res);
cand_MUSIC = zeros(1,res);
cand_corr = zeros(1,res);

% for MUSIC: eigenvalue decomposition and ordering in decreasing order
% two nonzero eigenvalues and three zero eigenvalues as expected
[U, D] = eig(Rx);
[B, I] = sort(diag(D), "descend");
U = U(:,I(num_sources + 2:end));
Cl = zeros(2*M-1,1);

for i=1:res
    a_theta = exp(-2*1i*pi*cos(theta_rad(i))*dist).';

    cand_corr(end-i+1) = mean(abs(a_theta' * X));
    cand_CB(i) = abs(a_theta'*Rx*a_theta/(a_theta'*a_theta));
    cand_MVDR(i) = abs((a_theta'/Rx*a_theta)^-1);
    cand_MUSIC(i) = 1./(a_theta'*(U*U')*a_theta);
end

C = (U*U');
for l = 0:M-1
    Cl(l+M) = sum(diag(C,l));
end

corr_delta_w = half_freq(cand_corr, res);
CB_delta_w = half_freq(cand_CB, res);
MVDR_delta_w = half_freq(cand_MVDR, res);
MUSIC_delta_w = half_freq(cand_MUSIC, res);

disp("Half Frequency Ranges (\Delta w) of Correlation Method are ")
disp(corr_delta_w)
disp("Half Frequency Ranges (\Delta w) of Classical Beamforming are ")
disp(CB_delta_w)

Cl(1:M-1) = flip(Cl(M+1:end));
rts = roots(Cl);
rts(abs(rts)>1) = 0;

[~,ind] = maxk(abs(rts),2);
true_rts = rts(ind);
true_angles = acos(Delta^-1*(2*pi)^-1*imag(log(true_rts)))/pi*180;

figure;
subplot(4,1,1);
plot(theta_deg,10*log10(abs(cand_CB)),'LineWidth',2);
title('Output Power vs Theta (Correlation Method)');
xlabel('Angle (dg)');ylabel('Power');

subplot(4,1,2);
plot(theta_deg,10*log10(cand_CB),'LineWidth',2);
title('Output Power vs Theta (Classical Beamformer Approach)');
xlabel('Angle (dg)');ylabel('Power');

subplot(4,1,3);
plot(theta_deg,10*log10(cand_MVDR),'LineWidth',2);
title('Output Power vs Theta (MVDR)');
xlabel('Angle (dg)');ylabel('Power');

subplot(4,1,4);
plot(theta_deg,10*log10(abs(cand_MUSIC)),'LineWidth',2);
title('Output Power vs Theta (MUSIC)');
xlabel('Angle (dg)');ylabel('Power');

function [delta_w] = half_freq(powers, res)
    max_power = max(powers);
    freq_range = find(powers>max_power/2);
    indices = [freq_range(1), 0; 0, freq_range(end)];
    for i = 1:length(freq_range)-1
        if freq_range(i)+1 < freq_range(i+1)
            indices(1,2) = freq_range(i);
            indices(2,1) = freq_range(i+1);
        end
    end
    
    if indices(1,2) == 0 && indices(2,1) == 0
        indices = sum(indices);
    end
    
    freqs = indices*180/res;
    delta_w = freqs(:,2) - freqs(:,1);
end
