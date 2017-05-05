clc, clear, close all;
% M_sr = load('sridhar_adbadd.txt'); 
% M_ht = load('sridhar_adbadd_ours.txt');
% M_ht = load('sridhar_adbadd_ours_.txt'); 

M_sr = load('E:/OneDrive/EPFL/Code/honline-mat/honline-figures/htrack_results_on_dexter/sridhar/flexex1.txt'); 
M_ht = load('E:/OneDrive/EPFL/Code/honline-mat/honline-figures/htrack_results_on_dexter/htrack_reinit/flexex1.txt');

% M_ht = load('sridhar_flexex1_ours.txt');


% M_sr = load('sridhar_pinch.txt'); 
% M_ht = load('sridhar_pinch_ours.txt'); 
% M_ht = load('sridhar_pinch_ours_.txt'); 

% M_sr = load('sridhar_fingercount.txt'); 
% M_ht = load('sridhar_fingercount_ours.txt'); 
% M_ht = load('sridhar_fingercount_ours_.txt'); 

% M_sr = load('sridhar_tigergrasp.txt'); 
% M_ht = load('sridhar_tigergrasp_ours.txt'); 
% M_ht = load('sridhar_tigergrasp_ours_.txt'); 

% M_sr = load('sridhar_fingerwave.txt'); 
% M_ht = load('sridhar_fingerwave_ours.txt');
% M_ht = load('sridhar_fingerwave_ours_.txt');

% M_sr = load('sridhar_random.txt'); 
% M_ht = load('sridhar_random_ours_.txt');
% M_ht = load('sridhar_random_ours__.txt');


num_frames = size(M_ht,1);
num_joints = size(M_ht,2)/3;
M_sr_B = M_sr;

% M_sr_B = M_sr(1:12,:)
% M_sr = reshape(M_sr(:)',2,6*3);

M_sr = nan(num_frames, 6*3);
for i=0:num_frames-1
    B = M_sr_B((i*6)+1:((i+1)*6),1:3 )';
    M_sr(i+1,:) = B(:);
end

figure(1); hold on;
e_ht = zeros(num_frames,1);
for i=1:num_frames
    cla
    sr = reshape(M_sr(i,:),3,6); sr(1,:) = -sr(1,:);
    ht = reshape(M_ht(i,:),3,num_joints);

    % subsample indexes (8 thumb)
    sr = sr(:,[1  2  3  4  5]);
    ht = ht(:,[8 24 20 16 12]);

    res = sqrt( sum( (ht-sr).^2, 1 ) );
    res( res>40 ) = []; 
    res( res>40 ) = 40;

	e_ht(i) = mean(res);
        
    cla
 	scatter3(sr(1, :), sr(2, :), sr(3, :),  20, '.k');
    scatter3(ht(1, :), ht(2, :), ht(3, :), 20, '.r' );
    view(-180,-90);
    drawnow;
end

figure(2)
e_ht = e_ht( ~isnan(e_ht) );
plot(e_ht);
mean(e_ht)
% median(e_ht)