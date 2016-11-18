function [history] = get_last_frame_covariance(H, history, settings, N)

B = 3; T = 3;

%% for displaying covariance

if settings.display_covariance
    beta_only_covariance = true;
    beta_theta_covariance = false;
    
    invert_block = false;
    invert_all = true;
    
    if (settings.quadratic_one || settings.kalman_two)
        invert_block = true;
        invert_all = false;
    end
    
    %% beta covariance
    if (beta_only_covariance)
        beta_indices = repmat([ones(B, 1); zeros(B, 1)], min(N, settings.batch_size), 1);
        H_beta = H(beta_indices == 1, beta_indices == 1);
        Sigma_beta = inv(H_beta);
        %Sigma_beta = Sigma_beta / norm(Sigma_beta);
        
        if (invert_block)
            sigma_last_frame = inv(H_beta(end - B + 1:end, end - B + 1:end));
        end
        
        if (invert_all)
            sigma_last_frame = Sigma_beta(end - B + 1:end, end - B + 1:end);
        end
        
    end
    
    %% beta theta covariance
    if (beta_theta_covariance)
        Sigma = inv(H);
        sigma_last_frame = Sigma(end - B - T + 1:end - T, end - B  - T + 1:end - T);
    end
    
    history.covariance(N, :, :) = sigma_last_frame;
end


if settings.display_full_covariance
    history.full_covariance(N, :, :) = H(end - B - T + 1:end, end - B - T + 1:end);
end

