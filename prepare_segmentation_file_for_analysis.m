function prepare_segmentation_file_for_analysis(data_dir_name)

%{
INSTRUCTIONS:
This function accepts as input the .tif file with the segmented nuclei / cells, isolates only nuclei / cells that are valid in terms of their size, not
overlapping with the borders of the image, etc. Then, it prepares the data for the next script: calc_morphological_characteristics.

%}

sd_to_size_coeff = 2 * 2.354; % constant - do not change!
last_spacing_filename = 'last_spacing.mat'; % constant - do not change!
n_parallel_workers = 2; % to get info about the availability of workers type 'get_num_of_workers', in the command window.
opengl_mode = 'Hardware';
min_volume = 100; % in cubic microns 
max_volume = 5000;

% checking how many workers are currently active:
warning('off', 'MATLAB:datetime:NonstandardSystemTimeZone');
warning('off', 'MATLAB:DELETE:Permission');
p = gcp('nocreate');
if isempty(p)
    poolsize = 0;
else
    poolsize = p.NumWorkers;
end

% changing the number of active workers - if necessary:
if poolsize ~= n_parallel_workers
    if poolsize ~= 0
        delete(gcp('nocreate'));
    end
    parpool(n_parallel_workers);
end
warning('on', 'MATLAB:datetime:NonstandardSystemTimeZone');
warning('on', 'MATLAB:DELETE:Permission');
%}

% configuring opengl mode:
OS = system_dependent('getos');
if regexpi(OS, 'Windows')
    curr_opengl_mode = opengl('data');
    if curr_opengl_mode.Software && strcmp(opengl_mode, 'Hardware')
        opengl software;
    elseif ~curr_opengl_mode.Software && strcmp(opengl_mode, 'Software')
        opengl hardware;
    end
end



clc;
close all;

last_dir_filename = get_last_dir_filename;

% loading the path of the dir used in the last operation of the script
% (unless it's the first use):
if exist(last_dir_filename, 'file')
    load(last_dir_filename);
else
    last_dir = '.';
    last_spacing = [0,0,0];
end

if ~exist('last_spacing', 'var')
    last_spacing = [0,0,0];
end

if length(last_spacing) ~= 3
    last_spacing = [0,0,0];
end


% asking the user to select files:

datafiles = dir(data_dir_name);

 file_count=0;
 for i=1:length(datafiles)
        t=strfind(datafiles(i).name,'.tif');
        if length(t)~=0
            file_count=file_count+1;
            all_files_name{file_count}=datafiles(i).name;
        end
 end

for fi=1:length(all_files_name)
    files=strcat(data_dir_name,'\',all_files_name{fi});
   
    % validating the response of the user:
    if isempty(files) || isnumeric(files)
        error('No files were selected, aborting script.');
    elseif length(files) == 1
        error('At least two files need to be selected, aborting script.');
    end
         
    % asking the user to set the spacing of the image:
    I.spacing=[0.194,0.194,0.387];
   
    last_spacing = I.spacing;
    save(last_spacing_filename, 'last_spacing');
    
    filename = fullfile (datafiles.folder,datafiles.name);
    %filename = fullfile(directory, files);
    I.info = imfinfo(files);
    I.size = [I.info(1).Height, I.info(1).Width, length(I.info)];
    
    % inferring the data type:
    if I.info(1).MaxSampleValue == intmax('uint16') && I.info(1).MinSampleValue == intmin('uint16')
        I.data_type = 'uint16';
        I.bits = 16;
    elseif I.info(1).MaxSampleValue == intmax('uint8') && I.info(1).MinSampleValue == intmin('uint8')
        I.data_type = 'uint8';
        I.bits = 8;
    elseif I.info(1).MaxSampleValue == intmax('int16') && I.info(1).MinSampleValue == intmin('int16')
        I.data_type = 'int16';
        I.bits = 16;
    elseif I.info(1).MaxSampleValue == intmax('int8') && I.info(1).MinSampleValue == intmin('int8')
        I.data_type = 'int8';
        I.bits = 8;
    else
        disp('Data type could not be identified - contact Tomer.');
        return;
    end
    
    % allocating memory for the image field:
    I.img = zeros(I.size, I.data_type);
    
    % loading the image:
    fprintf('Loading image... ');
    for i = 1 : I.size(3)
        str = [num2str(i), '/', num2str(I.size(3))];
        fprintf(str);
        I.img(:,:,i) = imread(files, 'Index', i);
        fprintf(repmat('\b', 1, length(str)));
    end
    disp(['Voxels = ', mat2str(I.size), ', Microns = ', mat2str(I.size .* I.spacing)]);
    
    % here we binarize the image and label each nucleus / cell, assuming 6 neighbors connectivity:
    [I.img, I.n_conncomp] = bwlabeln(I.img > 0, 6);
    disp(['Labeled binary image and identified ', num2str(I.n_conncomp), ' connected components.']);
    
%     % saving the cleaned image into a mat file:
%     disp('Saving image as a .mat file into the input directory..');
%     save(regexprep(I.info(1).Filename, '\.[^\.]+$', '.mat'), 'I');


% asking the user to set the SD of the kernel and preparing the kernel:
%smooth_sd = inputdlg('Please enter the S.D. of the smoothing kernel to be used (in microns):');

%smooth_sd = str2double(smooth_sd{1});
%smooth_sd = (smooth_sd{1});
%smooth_size = round(smooth_sd ./ I.spacing * sd_to_size_coeff);
smooth_size = round(0.5 ./ I.spacing * sd_to_size_coeff);
pad_size = smooth_size + 1;
ker = fspecial3('gaussian', smooth_size);

% setting the reduction coefficient for the surfaces:
%reduction_coeff = inputdlg('Please enter the coefficient of patch reduction (0, 1]:');

%reduction_coeff = str2double(reduction_coeff{1});

reduction_coeff = 0.1;

% % here we'll identify and isolate the cells:
fprintf('Isolating cells... ');
idcs = find(I.img);
cell_idx = I.img(idcs);
[cell_idx, order] = sort(cell_idx);
idcs = idcs(order);
[x,y,z] = ind2sub(I.size, idcs);
dif = [0; find(diff(cell_idx)); length(cell_idx)];
unq_cell_idx = cell_idx(dif(1:end-1)+1);
clear idcs order;

coords = zeros(I.n_conncomp, 6, 'uint16');
bw = cell(I.n_conncomp,1);
for i = 1 : I.n_conncomp
    str = [num2str(i), '/', num2str(I.n_conncomp)];
    fprintf(str);
    rng = dif(i)+1:dif(i+1);
    coords(i,:) = [min(x(rng)), max(x(rng)), min(y(rng)), max(y(rng)), min(z(rng)), max(z(rng))];
    bw{i,1} = I.img(coords(i,1):coords(i,2), coords(i,3):coords(i,4), coords(i,5):coords(i,6)) == cell_idx(dif(i)+1);
    cc = bwconncomp(bw{i,1});
    bw{i} = imopen(bw{i}, true(3,3,3));
    fprintf(repmat('\b', 1, length(str)));
end
clear dif cell_idx x y z rng;
% I = rmfield(I, 'img');
disp(['Isolated ', num2str(length(bw)), ' cells']);

fprintf('Checking for elements with volume that is outside the specified range... ');
volume = cellfun(@nnz, bw) * prod(I.spacing);
elements_to_remove = volume < min_volume | volume > max_volume;
if nnz(elements_to_remove) > 0
    bw = bw(~elements_to_remove);
    coords(elements_to_remove,:) = [];
    unq_cell_idx = unq_cell_idx(~elements_to_remove,1);
    disp(['identified and removed ', num2str(nnz(elements_to_remove)), ' elements. Current number: ', num2str(length(bw))]);
else
    disp('none were identified.');
end
clear elements_to_remove volume;

fprintf('Removing cells that overlap with the borders of the image... ');
initial_num_cells = size(coords,1);
% non_edge_nuclei = find(~(coords(:,1) == 1 | coords(:,2) == I.size(1) | coords(:,3) == 1 | coords(:,4) == I.size(2) | coords(:,5) == 1 | coords(:,6) == I.size(3)));
% prev_unq_cell_idx = unq_cell_idx;
% unq_cell_idx = unq_cell_idx(non_edge_nuclei);

frame = false(size(I.img));
frame([1,end],:,:) = 1;
frame(:,[1,end],:) = 1;
frame(:,:,[1,end]) = 1;
is_bw_outside = false(length(bw),1);
for i = 1 : size(coords,1)
    curr_region = frame(coords(i,1):coords(i,2), coords(i,3):coords(i,4), coords(i,5):coords(i,6));
    is_bw_outside(i) = any(curr_region(bw{i}));
    if is_bw_outside(i)
        curr_region = I.img(coords(i,1):coords(i,2), coords(i,3):coords(i,4), coords(i,5):coords(i,6));
        curr_region(bw{i}) = 0;
        I.img(coords(i,1):coords(i,2), coords(i,3):coords(i,4), coords(i,5):coords(i,6)) = curr_region;
    end
end
coords(is_bw_outside,:) = [];
bw(is_bw_outside) = [];
disp(['Removed ', num2str(initial_num_cells - length(bw)), ' (', num2str(round((initial_num_cells - length(bw)) / initial_num_cells * 1000) / 10), '%) cells, ' ...
    'current number of cells: ', num2str(length(bw))]);

fprintf('Extracting the surface of each cell/nucleus... ');
surfaces = repmat(struct('faces', [], 'vertices', []), length(bw), 1);
str = '';
is_problematic_nucleus = false(length(bw),1);
for i = 1 : length(bw)
    
    fprintf(repmat('\b', 1, length(str)));
    str = [num2str(i), '/', num2str(length(bw))];
    fprintf(str);
    
    padded_bw = padarray(bw{i,1}, pad_size, 0, 'both');
    smoothed = imfilter(double(padded_bw), ker);
    smoothed(smoothed < 0.001) = 0;
    
    [h, ~] = hist(nonzeros(smoothed(:)), 1000);
    cumulative_sum = cumsum(h(end:-1:1))';
    cumulative_sum = cumulative_sum(end:-1:1);
    
    [~, min_ind] = min(abs(nnz(bw{i,1}) - cumulative_sum));
    iso_value = min_ind / 1000;
    
    bw{i,1} = smoothed(pad_size(1)+1:end-pad_size(1), pad_size(2)+1:end-pad_size(2), pad_size(3)+1:end-pad_size(3)) >= iso_value;
    surface = isosurface(smoothed, iso_value);
    try
        surface.vertices = surface.vertices(:,[2,1,3]);
    catch
        is_problematic_nucleus(i) = true;
        continue;
    end
    
    surface.vertices = (surface.vertices - repmat(pad_size, size(surface.vertices,1), 1) + repmat(double(coords(i,[1,3,5]))-1, size(surface.vertices,1), 1)) .* repmat(I.spacing, size(surface.vertices,1), 1); 
    
    surface = reducepatch(surface, reduction_coeff);
    
    surfaces(i,1) = surface;
    
    if 0
        figure('Color', [0 0 0]); %#ok<*UNRCH>
        patch(surface, 'FaceColor', [1 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        rotate3d on;
        daspect([1,1,1])
        axis vis3d;
        camlight;
        axis off;
    end
end

n_problematic_nuclei = nnz(is_problematic_nucleus);
if n_problematic_nuclei > 0
    coords(is_problematic_nucleus,:) = [];
    bw(is_problematic_nucleus) = [];
    surfaces(is_problematic_nucleus) = [];
    disp(['Removed ', num2str(n_problematic_nuclei), ' nuclei, total remaining nuclei: ', num2str(length(bw))]);
end

% saving the surfaces into a .mat file as a structure array:
spacing = I.spacing; %#ok<*NASGU>
%save(fullfile(datafiles.folder, regexprep(datafiles.name, '\..*$', ' (Surfaces).mat')), 'surfaces', 'bw', 'coords', 'spacing');
 save(regexprep(I.info(1).Filename, '\.[^\.]+$', ' (Surfaces).mat'), 'surfaces', 'bw', 'coords', 'spacing');
%save(fullfile(datafiles.folder, [datafiles.name, '\..*$', ' (Surfaces).mat']), 'surfaces', 'bw', 'coords', 'spacing');
end



end