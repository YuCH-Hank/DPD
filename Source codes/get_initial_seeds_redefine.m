%% GO through (run) this code step by step and see how each parameter change.20131231

%% >>>> Get Initial Seeds
function  seeds_data = get_initial_seeds_redefine(seeds_data_in, im_base_name)
% OUTPUT ===================================================================
%
% seeds_data.seeds_
% seeds_data.initial_seeds_shape
% seeds_data.min_allowed_points_for_initial_seeding
% seeds_data.im_d
% seeds_data.places_to_check_for_initial_seeding
% seeds_data.divergance_from_plan_mat
%
% =========================================================================

% >>>> Image
im_d = seeds_data_in.im_d;
[r, c] = size(im_d); % the size of the input image

% >>>> seeds_data_in >>>>>>>>>>>>>>>>>>>>>
places_to_check_for_initial_seeding     = seeds_data_in.places_to_check_for_initial_seeding;
initial_seeds_shape                     = seeds_data_in.initial_seeds_shape;
min_allowed_points_for_initial_seeding  = seeds_data_in.min_allowed_points_for_initial_seeding; %(jz)it's decided by the seeds shape
seed_shape                              = seeds_data_in.seed_shape;
[r_seed, c_seed]                        = size(initial_seeds_shape); % row of seed and col of seed
divergance_from_plan_mat = NaN * ones(r, c); %(initialise the matrix of divergance plan with value NaN)

% >>>> Initial Parameter >>>>>>>>>>>>>>>>>>>>>
idx = 0;
[ri_idx_base, ci_idx_base] = find( initial_seeds_shape ~= 0 ); % change into index;
truth = ones(r,c);
cand  = ones(r,c,r_seed * c_seed ,3); 
equation = ones(r,c,5);

% >>>> Check Each point in patch  >>>>>>>>>>>>>>>>>>>>>
for ci = 1 : (c - c_seed - 1) % don't change the order of scaning the matrix %(jz) the scaning order is from top to bottom and then from left to right
                              % column index = column of image - column of patch
    clc;
    fprintf('Finished %.2f %% >>>>>>> \n', ci/(c - c_seed - 1)*100);
    fprintf('Processing column %d \n', ci);
    
    parfor ri = 1 : (r - r_seed - 1)    % candidate_points can be treated as a butter
                                        % row index = row of image - row of patch
        candidate_points = zeros(r_seed * c_seed,3);
        candidate_points( : , 1 ) = ri_idx_base + ri - 1;
        candidate_points( : , 2 ) = ci_idx_base + ci - 1;
        
        % >>>> using the function 'sub2ind' to find the exact point in the im_d
        % >>>> if the input image is depth map then it's Z value
        candidate_points( : , 3 ) = im_d( sub2ind(  size(im_d), ...
                                                    candidate_points( : , 1 ) , ...
                                                    candidate_points( : , 2 )   )   ); 
        
        % >>>> Check Validation
        [valid_points, candidate_points] = are_points_valid_candidate_for_initial_seeding(...
                                                candidate_points, ...
                                                places_to_check_for_initial_seeding, ...
                                                min_allowed_points_for_initial_seeding);
        % >>>> Fit plane
        if valid_points
            truth(ri,ci) = true;
            
            try 
                cand(ri,ci,:,:) = candidate_points;
                
                ri_idx = candidate_points(:, 1);    % >>>> row
                ci_idx = candidate_points(:, 2);    % >>>> cols
                z_data = candidate_points(:, 3);    % >>>> depth
                warning('off')
                
                sf = fit( [ri_idx, ci_idx], z_data, 'poly11', 'Robust', 'Bisquare'); %Matlab built-in function
                equation(ri,ci,:) = [sf.p10, sf.p01, -1, sf.p00, 0];    %(jz) the "0" here is left for future recording the update times, since it could have have more fitting surface functions.

                divergance_from_plan_mat(ri,ci) = mean((z_data - feval(sf, ri_idx, ci_idx)).^2); %(jz)find the average divergance of each point in the seeds shape to the fitting surface.
            
            catch
                print('%d,%d',ri,ci);
                
            end
        end
    end
end

%% seperate
for ci = 1 : (c - c_seed - 1)   % don't change the order of scaning the matrix %(jz) the scaning order is from top to bottom and then from left to right
                                % column index = column of image - column of patch
    for ri = 1 : (r - r_seed + 1) % candidate_points can be treated as a butter
                                    % row index = row of image - row of patch
        clc;
        fprintf('Finished %.2f %% >>>>>>> \n', ci/(c - c_seed - 1)*100);
        fprintf('Processing column %d \n', ci);
        
        if truth(ri,ci)
            idx = idx + 1;
            data_point(:,:) = cand(ri,ci,:,:);
            equations_{ 1 } = [equation(ri,ci,1), equation(ri,ci,2), equation(ri,ci,3), equation(ri,ci,4), 0];
            seeds_data.seeds_{idx}.equations_of_plan_ = equations_;
            seeds_data.seeds_{idx}.points = data_point(:,:);
            seeds_data.seeds_{idx}.divergance_from_plan = divergance_from_plan_mat(ri,ci);
        end
        
    end
end

%%
% >>>> Record File
% >>>> copy data from seeds_data_in
seeds_data.initial_seeds_shape                      = seeds_data_in.initial_seeds_shape;
seeds_data.min_allowed_points_for_initial_seeding   = seeds_data_in.min_allowed_points_for_initial_seeding;
seeds_data.im_d                                     = im_d;
seeds_data.places_to_check_for_initial_seeding      = places_to_check_for_initial_seeding;
seeds_data.divergance_from_plan_mat                 = divergance_from_plan_mat;
seeds_data.order_of_generating_seeds                = seeds_data_in.order_of_generating_seeds;
seeds_data.seed_shape                               = seed_shape;

show_im = 1; fig_idx = 1;

if show_im
    figure(fig_idx); image(double(divergance_from_plan_mat)); colormap(gray(255)); title(sprintf('the divergance from plan for seeds'));
end

% >>>> Record File
results_file = sprintf('%s_SeedName_%s_SeedSize_%03d.mat', im_base_name, upper(seeds_data.seed_shape), sum(sum(initial_seeds_shape)));
save(results_file);

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

%% >>>> Check Validation
function [valid_points, candidate_points_2nd_filtered] = are_points_valid_candidate_for_initial_seeding(candidate_points, places_to_check_for_initial_seeding, min_allowed_points_for_initial_seeding)

% >>>> check if depth is not zero
% >>>> return 1 or 0
idx_valid_points_1st_filtered = (candidate_points(:, 3) ~= 0); %all the index of the nonzero point intensity in "idx_valid_points_1st_filtered"
candidate_points_1st_filtered =  candidate_points(idx_valid_points_1st_filtered, :); %all the onozero point intensity in candidate_points buffer, e.g. candidate_points (1,:)= the first row of candidate_points

% >>>> check initial seeding
idx_valid_points_2nd_filtered =  find( places_to_check_for_initial_seeding(  ...
                                        sub2ind(    size(places_to_check_for_initial_seeding), ...
                                                    candidate_points_1st_filtered(:, 1), ...
                                                    candidate_points_1st_filtered(:, 2) ) ) );
                                                
candidate_points_2nd_filtered =  candidate_points_1st_filtered(idx_valid_points_2nd_filtered, :); %candidate_points_1st_filtered (1,:)= the first row of candidate_points_1st_filtered
valid_points = ( size(candidate_points_2nd_filtered, 1) >= min_allowed_points_for_initial_seeding); 

% if this center has no hole, then return 1.
return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

%% >>>> Find plane of each patch
function [equations_of_plan_, distances_from_plan_vct] = get_best_fitting_plan_for_points(candidate_points)

%    d = -ax_0-by_0-cz_0
ri_idx = candidate_points(:, 1);    % >>>> row
ci_idx = candidate_points(:, 2);    % >>>> cols
z_data = candidate_points(:, 3);    % >>>> depth
warning('off')
sf = fit( [ri_idx, ci_idx], z_data, 'poly11', 'Robust', 'Bisquare'); %Matlab built-in function
% distances_from_plan_vct is the distance vector of each point to the fitting surface
distances_from_plan_vct = (z_data - feval(sf, ri_idx, ci_idx)).^2; %feval is Matlab built-in function;

%      Linear model Poly11:
%      sf(x,y) = p00 + p10*x + p01*y
%            z = p00 + p10*x + p01*y
%     (jz)   0 = p10*x + p01*y - z + p00
%            0 = A * x + B * y + C * z + D;
%      p00 = - p10*x - p01*y + z
%   [ A, B, C, D] = [ sf.p10, sf.p01, -1, sf.p00];

A = sf.p10;
B = sf.p01;
C = -1;   %(jz) the coefficient of z
D = sf.p00;
% >>>> plane -> 0 = A * x + B * y + C * z + D;
equations_of_plan_{ 1 } = [A, B, C, D, 0]; %(jz) the "0" here is left for future recording the update times, since it could have have more fitting surface functions.

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


