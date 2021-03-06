clear, clc, close all; load('.\colormap.mat');
tic
%% ================== Read Figure  =======================
i = 7;
data_dir = 'TEST';
im_base_name = sprintf( 'TEST%d', i );
data_type = 'bmp';
im_file = sprintf( '%s/%s.%s', data_dir, im_base_name , data_type);
im = double(imread(im_file));     % matlab imread -> unit8
im_d = im(:,:,1); clear im;

%% ================== Parameter ================== 
par             = true;         % using parallel cpu
d_samp          = false;        % using down sampling
d_samp_size     = 150;          % down sampling size
u_samp          = true;         % using up sampling

%% ================== Figure Processing ================== 
[r_, c_] = size(im_d);
if (d_samp)
    numrows = r_ / max(r_,c_) * d_samp_size;
    numcols = c_ / max(r_,c_) * d_samp_size;
    im_d = imresize(im_d, [numrows numcols]);
end

%% ================== generate seed  =======================
% >>>> seeds_data_in
seeds_data_in.initial_seeds_shape                       = ones(4);    % 4 * 4;
seeds_data_in.seed_shape                                = 'square';
seeds_data_in.order_of_generating_seeds                 = 'RasterScan';
seeds_data_in.im_d                                      = im_d;
seeds_data_in.places_to_check_for_initial_seeding       = im_d;
seeds_data_in.min_allowed_points_for_initial_seeding    = sum(sum(seeds_data_in.initial_seeds_shape));     % number of minimum initial seeding

% >>>> Get Initial Seeds
if (par)
    seeds_data = get_initial_seeds_redefine(seeds_data_in, im_base_name);
else
    seeds_data = get_initial_seeds(seeds_data_in,im_base_name);
end

%% ================== plane growing  =======================
fields_data_in.seeds_data{1}.im_d                       = [];
fields_data_in.im_file                                  = im_file;
fields_data_in.im_d                                     = im_d;
fields_data_in.places_need_to_be_checked_for_seeding    = ones(size(fields_data_in.im_d));
fields_data_in.order_of_using_seeds                     = 'ascendorder';
% random, rasterscan_col, raster_col, ascendorder, descendorder, debug
fields_data_in.basic_data.lag_of_update_THR             = 1;
fields_data_in.basic_data.maximum_allowed_dissimilarity = 3;
fields_data_in.seeds_data{1}                            = seeds_data;
weight                                                  = 0.009;
plot_                                                   = true;

[fields_data, ~] = region_growing(fields_data_in, weight, plot_);

%%
if (u_samp)
    im = imresize(fields_data.fields, [r_ c_], 'box');
    index = find(im < 1);
    im(index) = 1;
    im = round(im);
    
    figure('name', 'Original');
    image(uint8(im));colormap(settelments_colormap); axis off;
else
    figure('name', 'Original');
    image(uint8(fields_data.fields));colormap(settelments_colormap); axis off;
    
end

%% ================== over-grwoing correction  =======================
epoch_idx = 1;
basic_penetrating_element_shape = 'line';
show.show_im = 1;
show.fig_idx = 2;
show.settelments_colormap = settelments_colormap;
[fields_data_out] = over_growing_correction(fields_data, basic_penetrating_element_shape, epoch_idx, show);

if (u_samp)
    im = imresize(fields_data_out.fields, [r_ c_]);
    index = find(im <= 1);
    im(index) = 1;
    im = round(im);
    
    figure('name', 'Over Growing Correction');
    image(uint8(im));colormap(settelments_colormap); axis off;
else
    figure('name', 'Over Growing Correction');
    image(uint8(fields_data_out.fields));colormap(settelments_colormap); axis off;
end

%% ================== under-grwowing correction  =======================
fields_data_out1 =  under_growing_correction(fields_data_out);
if (u_samp)
    im = imresize(fields_data_out1.parallel_surface_detection.field_index, [r_ c_]);
    index = find(im <= 1);
    im(index) = 1;
    im = round(im);
    
    figure('name', 'Under Growing Correction');
    image(uint8(im));colormap(settelments_colormap); axis off;
else
    figure('name', 'Under Growing Correction');
    image(uint8(fields_data_out1.parallel_surface_detection.field_index));colormap(settelments_colormap); axis off;
end

%% ================== Save Figure  =======================
saveas(gcf,[data_dir, '\', im_base_name, '_DPD.jpg']);
toc




