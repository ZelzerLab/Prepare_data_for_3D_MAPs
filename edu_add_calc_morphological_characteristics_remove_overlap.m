function G = edu_add_calc_morphological_characteristics_03_n_c_ratio(ignore_cells,ignore_edu,...
    data_dir_name, cells_filename,nuclei_filename,edu_filename,output_filename,nuclei_shift)

%{
INSTRUCTIONS:
This function accepts the .mat files of both cells and nuclei (one pair per growth plate zone) that contain the surfaces, from which it calculates most of the
morphological parameters. The names of the input files end with ' (Surfaces).mat', and the output file of each zone is called 'zone (Characteristics).mat'.
After running this function four times (once for each zone), you should rename each of the output files by replacing the word 'zone' with the initials of the
zone itself (i.e. rz, pz, phz, hz), and copy them to the same folder, so you could load all of them together during the run of the next function, which is
'plot_cell_morphology_analysis.m'.
%}

[C, N, E, paths] = load_data_for_calc_morphological_characteristics(ignore_cells, ignore_edu, ...
    data_dir_name, cells_filename,nuclei_filename,edu_filename,nuclei_shift);

%output_filename = inputdlg('Please enter the name of the output file:', '', 1);

if ~ignore_cells
    to_remove = cellfun(@isempty, C.masks);
    C.coords(to_remove,:) = [];
    C.surfaces(to_remove,:) = [];
    C.spacing(to_remove,:) = [];
    C.masks(to_remove,:) = [];
    C.origins(to_remove,:) = [];
    C.edu_signal(to_remove,:)=[];
end

to_remove = cellfun(@isempty, N.masks);
N.coords(to_remove,:) = [];
N.surfaces(to_remove,:) = [];
N.spacing(to_remove,:) = [];
N.masks(to_remove,:) = [];
N.origins(to_remove,:) = [];
N.edu_signal(to_remove,:)=[];

if ~ignore_edu
    to_remove = cellfun(@isempty, E.masks);
    E.coords(to_remove,:) = [];
    E.surfaces(to_remove,:) = [];
    E.spacing(to_remove,:) = [];
    E.masks(to_remove,:) = [];
    E.origins(to_remove,:) = [];
end


if ~ignore_cells
    [C, N] = sort_cells_and_nuclei(C, N, E,ignore_edu);
    C.original_idx = 1:size(C.coords,1);
    C.non_unique_idcs = [];
end

N.original_idx = 1:size(N.coords,1);
N.non_unique_idcs = [];


G.nuc = calc_individual_shape_characteristics(N.surfaces, N.masks, N.origins, N.spacing, N.original_idx, N.non_unique_idcs);
if ~ignore_cells
    G.cel = calc_individual_shape_characteristics(C.surfaces, C.masks, C.origins, C.spacing, C.original_idx, C.non_unique_idcs);
    G.inter = calc_inter_characteristics(G, C.n_overlapping);
end

save_results(paths, G, C, N, ignore_cells, output_filename);


function [C, N,E, paths] = load_data_for_calc_morphological_characteristics(ignore_cells, ignore_edu, ...
    data_dir_name, cells_filename,nuclei_filename,edu_filename,nuclei_shift)

cells_dirname=data_dir_name;
nuclei_dirname=data_dir_name;
edu_dirname=data_dir_name;

if ~ignore_cells
    %[cells_filename, cells_dirname] = get_file_name('* (Surfaces).mat', 'Please choose the .mat file of the MEMBRANES:', 'on');
    %C = struct('bw', [], 'coords', zeros(0,6), 'spacing', zeros(0,3), 'surfaces', []);
    C = load(fullfile(cells_dirname, cells_filename));
    coords=C.coords + uint16([0 0 0 0 nuclei_shift nuclei_shift]);
    spacings = repmat(C.spacing, length(C.bw), 1);
    volume = cellfun(@nnz, C.bw) .* prod(spacings, 2);  
    origins = coords(:,[1,3,5]);
    centroid=calc_centroids(C.bw, origins, spacings);

    
    surfaces={}; 
    for i=1:length(C.surfaces)
        surfaces{i}=C.surfaces(i).vertices;
    end
    
    disp('Time for removing overlap or duplicate cells')
    ia=RemoveBadCells(centroid, volume, surfaces);

    C.masks = C.bw(ia);
    C = rmfield(C, 'bw');
    C=  rmfield(C, 'coords');
    C.spacing=spacings(ia,:);
    C.origins=origins(ia,:);
    C.coords=coords(ia,:);
    C.surfaces=C.surfaces(ia,:);	
    C.edu_signal=zeros(size(C.masks,1),1,'logical');
    paths.cells_filename = cells_filename;
    paths.cells_dirname = cells_dirname;

else
    C = [];
end

%[nuclei_filename, nuclei_dirname] = get_file_name('* (Surfaces).mat', 'Please choose the .mat file of the NUCLEI:', 'on');
%N = struct('bw', [], 'coords', zeros(0,6), 'spacing', zeros(0,3), 'surfaces', []);
N = load(fullfile(nuclei_dirname, nuclei_filename));
%Remove Duplicate or overlap nuclei 

%coords=N.coords + uint16([0 0 0 0 nuclei_shift nuclei_shift]);
%spacings = repmat(N.spacing, length(N.bw), 1);
%origins = N.coords(:,[1,3,5]);

ia=RemoveBadNuclei(N.surfaces);
N.masks = N.bw(ia);
N = rmfield(N, 'bw');
N.coords=N.coords(ia,:);
N.surfaces=N.surfaces(ia,:);
N.edu_signal=zeros(size(N.masks,1),1,'logical');
N.origins = N.coords(:,[1,3,5]);
N.spacing = repmat(N.spacing, length(N.masks), 1);
paths.nuclei_filename = nuclei_filename;
paths.nuclei_dirname = nuclei_dirname;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[edu_filename, edu_dirname] = get_file_name('* (Surfaces).mat', 'Please choose the .mat file of the edu:', 'on');
%E = struct('bw', [], 'coords', zeros(0,6), 'spacing', zeros(0,3), 'surfaces', []);

if ~ignore_edu
    E = load(fullfile(edu_dirname, edu_filename));
    E.masks = E.bw;
    E = rmfield(E, 'bw');
    E.origins = E.coords(:,[1,3,5]);
    E.spacing = repmat(E.spacing, length(E.masks), 1);
    paths.edu_filename = edu_filename;
    paths.edu_dirname = edu_dirname;
else
    E=[];
end

function [C, N] = sort_cells_and_nuclei(C, N, E, ignore_edu)

if ~ignore_edu 
    img_size = max([ E.coords(:,[2,4,6]);  N.coords(:,[2,4,6])    ], [], 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% if overlap between Nuclei and EDU is +ve then label Nuclei with EDU 
    nimg = zeros(img_size, 'uint16');
    for i = 1 : length(N.masks)
        nimg(N.coords(i,1):N.coords(i,2), N.coords(i,3):N.coords(i,4), N.coords(i,5):N.coords(i,6)) =  uint16(N.masks{i}) * uint16(i);
    end
    overlap = zeros(length(E.masks), length(N.masks)); 
    for i = 1 : length(E.masks)
        cellular_region = nimg(E.coords(i,1):E.coords(i,2), E.coords(i,3):E.coords(i,4), E.coords(i,5):E.coords(i,6));
        cellular_vals = unique(nonzeros(cellular_region));
        nucleus_nnz = nnz(E.masks{i});
        for j = 1 : length(cellular_vals)
            overlap(i, cellular_vals(j)) = sum(cellular_region(:)==cellular_vals(j)) / nucleus_nnz;
        end
    end

    %dlmwrite('overlap_E_N.dat',overlap,'\t');
    % Threshold value between nucleus and edu. Change this threshold if you need.
    
    threshold_value=0.1;

    overlap_lg=zeros(size(overlap),'logical');
    for j=1:size(overlap,1)
         [max_val,max_ind]=max(overlap(j,:));
         if max_val>=threshold_value
             overlap_lg(j,max_ind)=1;
         end
    end

    %remove if any column object overlap multiple time with any row object  
    for i=1:size(overlap,2)
         voxel_whose_value_is_1=find(overlap_lg(:,i)==1);
         [max_val,max_ind]=max(overlap(voxel_whose_value_is_1,i));
         for j=1:length(voxel_whose_value_is_1)
             if voxel_whose_value_is_1(j)~=voxel_whose_value_is_1(max_ind)
                 overlap_lg(voxel_whose_value_is_1(j),i)=0;
             end     
         end
    end

    % for i=1:size(b,1)
    %     temp=find(b(i,:));
    %     if length(temp)==1
    %         E.channel(i)=1;
    %     end
    % end

    for i=1:size(overlap_lg,2)
        temp=find(overlap_lg(:,i));
        if length(temp)==1
           N.edu_signal(i)=1;
        end
    end
%dlmwrite('overlap_E_N_lg.dat',overlap_lg,' ');
clear overlap_lg
clear overlap 

end


%%%%%%%%%%%%%%%%% Written  by Tomer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% overlap between Cell and Nuclei
img_size = max([C.coords(:,[2,4,6]); N.coords(:,[2,4,6])    ], [], 1); 
cimg = zeros(img_size, 'uint16');
for i = 1 : length(C.masks)
    cimg(C.coords(i,1):C.coords(i,2), C.coords(i,3):C.coords(i,4), C.coords(i,5):C.coords(i,6)) = uint16(C.masks{i}) * uint16(i);
end
overlap = zeros(length(N.masks), length(C.masks));
for i = 1 : length(N.masks)
    cellular_region = cimg(N.coords(i,1):N.coords(i,2), N.coords(i,3):N.coords(i,4), N.coords(i,5):N.coords(i,6));
    cellular_vals = unique(nonzeros(cellular_region));
    nucleus_nnz = nnz(N.masks{i});
    for j = 1 : length(cellular_vals)
        overlap(i, cellular_vals(j)) =  sum(cellular_region(:)==cellular_vals(j)) / nucleus_nnz;
    end
end


% we say that a nucleus is within a cell if 0.5 of its voxels are within the cells volume:
 
disp('Threshold in percentage for nuclei and cell volume overlap for cross features')
threshold_value=0.5;

overlap_lg=zeros(size(overlap),'logical');
for j=1:size(overlap,1)
     [max_val,max_ind]=max(overlap(j,:));
     if max_val>=threshold_value
         overlap_lg(j,max_ind)=1;
     end
end

%remove if any column object overlap multiple time with any row object  
for i=1:size(overlap,2)
     voxel_whose_value_is_1=find(overlap_lg(:,i)==1);
     [max_val,max_ind]=max(overlap(voxel_whose_value_is_1,i));
     for j=1:length(voxel_whose_value_is_1)
         if voxel_whose_value_is_1(j)~=voxel_whose_value_is_1(max_ind)
             overlap_lg(voxel_whose_value_is_1(j),i)=0;
         end     
     end
end

%dlmwrite('overlap_N_C_lg.dat',overlap_lg,' ');
% checking that each cell has not more than one nucleus inside (we're ignoring dividing cells for now):
if any(sum(overlap_lg, 1)) > 1
    keyboard;
end

% preparing the output, wherein the left column is the index of the nucleus and the right is this of the cell:
[matching(:,1), matching(:,2)] = find(overlap_lg);


% for i = 1 : length(matching(:,1))
%     cimg(N.coords(matching(i,1),1):N.coords(matching(i,1),2), N.coords(matching(i,1),3):N.coords(matching(i,1),4), N.coords(matching(i,1),5):N.coords(matching(i,1),6)) = uint16(N.masks{matching(i,1)}) * uint16(12);
% end

% cleaning the results - some cells overlap with two nuclei
matching_temp = unique(matching(:,2));
for i = 1:length(matching_temp),
    ind_temp = find(matching(:,2) == matching_temp(i));
    
    if length(ind_temp)>1,
        val_temp = [];
        for j = 1:length(ind_temp),
            val_temp = [val_temp;overlap(matching(ind_temp(j),1),matching(ind_temp(j),2))];
        end
        [~,id_max] = max(val_temp);
        matching(setdiff(ind_temp,ind_temp(id_max)),:) = zeros(length(setdiff(ind_temp,ind_temp(id_max))),2);
    end
end

matching(find(matching(:,1) == 0),:) = [];


% cleaning the results - some nuclei overlap with two cells
matching_temp = unique(matching(:,1));
for i = 1:length(matching_temp),
    ind_temp = find(matching(:,1) == matching_temp(i));
    
    if length(ind_temp)>1,
        val_temp = [];
        for j = 1:length(ind_temp),
            val_temp = [val_temp;overlap(matching(ind_temp(j),1),matching(ind_temp(j),2))];
        end
        [~,id_max] = max(val_temp);
        matching(setdiff(ind_temp,ind_temp(id_max)),:) = zeros(length(setdiff(ind_temp,ind_temp(id_max))),2);
    end
end

matching(find(matching(:,1) == 0),:) = [];

non_included_nuclei = setdiff(1:length(N.masks), matching(:,1));
non_included_cells = setdiff(1:length(C.masks), matching(:,2));

prev_N = N;
prev_C = C;


fn = setdiff(fieldnames(N), 'non_unique_idcs');
for i = 1 : length(fn)
    N.(fn{i}) = N.(fn{i})(matching(:,1),:);
    N.(fn{i})(end+1:end+length(non_included_nuclei),:) = prev_N.(fn{i})(non_included_nuclei,:);
    C.(fn{i}) = C.(fn{i})(matching(:,2),:);
    C.(fn{i})(end+1:end+length(non_included_cells),:) = prev_C.(fn{i})(non_included_cells,:);
end

N.n_overlapping = size(matching,1);
C.n_overlapping = size(matching,1);









function save_results(paths, G, C, N, ignore_cells, output_filename) 

% saving a mat file with the results into the directory of the membranes file:
save(fullfile(paths.nuclei_dirname, [output_filename{1}, ' (Characteristics).mat']), 'G', 'C', 'N','-v7.3');

% saving the results into an excel sheet inside the directory of the membranes file:

% if ~ignore_cells
%     writetable(struct2table(G.cel), ...
%         fullfile(paths.cells_dirname, regexprep(output_filename{1}, '\(Surfaces\).mat', '(Characteristics).xlsx')), ...
%         'FileType', 'spreadsheet', 'Sheet', 'Cells');
% 
%     G.inter = rmfield(G.inter, 'n_overlapping');
%     writetable(struct2table(G.inter), ...
%         fullfile(paths.cells_dirname, regexprep(output_filename{1}, '\(Surfaces\).mat', '(Characteristics).xlsx')), ...
%         'FileType', 'spreadsheet', 'Sheet', 'Inter');
% end
% 
% writetable(struct2table(G.nuc), ...
%     fullfile(paths.nuclei_dirname, regexprep(output_filename{1}, '\(Surfaces\).mat', '(Characteristics).xlsx')), ...
%     'FileType', 'spreadsheet', 'Sheet', 'Nuclei');





function inter = calc_inter_characteristics(G, n_overlapping)

if n_overlapping>0
	inter.n_overlapping = n_overlapping;
	inter.volume_ratio(1:n_overlapping,1) = G.cel.volume(1:n_overlapping) ./ G.nuc.volume(1:n_overlapping);
	inter.centroid_shift_physical_coords(1:n_overlapping,:) = G.nuc.centroids(1:n_overlapping,:) - G.cel.centroids(1:n_overlapping,:);
	for i = 1 : n_overlapping
	    inter.centroid_shift_pca_coords(i,:) = inter.centroid_shift_physical_coords(i,:) / G.cel.PCA_coeff{i};
	end
else
	inter.n_overlapping=[];inter.volume_ratio=[];inter.centroid_shift_physical_coords=[];inter.centroid_shift_pca_coords=[];
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





function G = calc_individual_shape_characteristics(surfaces, masks, origins, spacings, idx, non_unique_idcs)

if size(masks,1)>0
	G.index = idx';
	[~, G.surface_area] = calc_surface_volume_and_area(surfaces);
	G.volume = cellfun(@nnz, masks) .* prod(spacings, 2);
	G.sa2vol_ratio = G.surface_area ./ G.volume;
	G.sphericity = (pi^(1/3))*((6*G.volume).^(2/3)) ./ G.surface_area;
	G.centroids = calc_centroids(masks, origins, spacings);
	[G.PCA_coeff, G.PCA_latent] = calc_principal_components(masks, spacings, surfaces);
	[G.ellipsoid_center, G.ellipsoid_radii, G.ellipsoid_evecs, G.ellipsoid_v, G.ellipsoid_chi2, G.PC_range] = approximate_ellipsoid(surfaces);
else
	G.ellipsoid_center=[];G.ellipsoid_radii=[];G.ellipsoid_evecs=[];G.ellipsoid_v=[]; G.ellipsoid_chi2=[];G.PC_range=[];
	G.PCA_coeff=[];G.PCA_latent=[];G.sa2vol_ratio=[];G.sphericity=[]; G.index=[];G.surface_area=[];G.volume=[];G.centroids=[];
end


function [center, radii, evecs, v, chi2, PC_range] = approximate_ellipsoid(surfaces)

center = nan(size(surfaces,1),3);
radii = nan(size(surfaces,1),3);
evecs = cell(size(surfaces,1),1);
v = cell(size(surfaces,1),1);
chi2 = nan(size(surfaces,1),1);

PC_range = nan(size(surfaces,1),3);

for i = 1 : length(surfaces)
    warning off;
    [center(i,:), radii(i,:), evecs{i}, v{i}, chi2(i,:)] = ellipsoid_fit_new(surfaces(i).vertices);

    [~, score, ~] = pca(surfaces(i).vertices, 'NumComponents', 3);
    PC_range(i,:) = range(score);
    
    continue;
    
    mind = min(surfaces(i).vertices);
    maxd = max(surfaces(i).vertices);
    nsteps = 100;
    step = (maxd - mind) / nsteps;
    [ x, y, z ] = meshgrid( linspace( mind(1) - step(1), maxd(1) + step(1), nsteps ), linspace( mind(2) - step(2), maxd(2) + step(2), nsteps ), linspace( mind(3) - step(3), maxd(3) + step(3), nsteps ) );
    
    Ellipsoid = v{i}(1) *x.*x +   v{i}(2) * y.*y + v{i}(3) * z.*z + ...
        2*v{i}(4) *x.*y + 2*v{i}(5)*x.*z + 2*v{i}(6) * y.*z + ...
        2*v{i}(7) *x    + 2*v{i}(8)*y    + 2*v{i}(9) * z;
    
    figure;
    patch( isosurface( x, y, z, Ellipsoid, -v{i}(10) ), 'FaceColor', 'g', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    hold on;
    scatter3(surfaces(i).vertices(:,1), surfaces(i).vertices(:,2), surfaces(i).vertices(:,3), '.b');
    hold off;
    axis image;
    axis vis3d;
    set(gca, 'Position', [0.3 0.3 0.4 0.4]);
    cameratoolbar;
    close;
end
radii = sort(abs(radii), 2, 'descend');
PC_range = sort(PC_range, 2, 'descend');


function [coeff, latent] = calc_principal_components(masks, spacings, surfaces)

% Note: in 'coeff', the first column (i.e. coeff(:,1)) is the orientation vector of PC1.

coeff = cell(length(masks),1);
latent = nan(length(masks),3);

for i = 1 : length(masks)
    [x,y,z] = ind2sub(size(masks{i}), find(masks{i}));
    xyz = [x,y,z] .* repmat(spacings(i,:), length(x), 1);
    xyz = xyz - repmat(mean(xyz,1), size(xyz,1), 1);
    try
        [coeff{i,1}, ~, latent(i,:)] = pca(xyz);
    catch
        keyboard;
    end
    
    if 0
        figure;
%         scatter3(xyz(1:100:end,1), xyz(1:100:end,2), xyz(1:100:end,3), '.');
        s = surfaces(i);
        s.vertices = s.vertices - repmat(mean(s.vertices,1), size(s.vertices,1), 1);
        show_surface(s);
        hold on;
        quiver3(0, 0, 0, coeff{i}(1,1), coeff{i}(2,1), coeff{i}(3,1), 10);
        quiver3(0, 0, 0, coeff{i}(1,2), coeff{i}(2,2), coeff{i}(3,2), 10);
        quiver3(0, 0, 0, coeff{i}(1,3), coeff{i}(2,3), coeff{i}(3,3), 10);
        hold off;
        cameratoolbar;
        axis image;
        axis vis3d;
        axis off;
        close;
    end
end


function centroids = calc_centroids(masks, origins, spacings)

centroids = nan(length(masks),3);
for i = 1 : length(masks)
    [x,y,z] = ind2sub(size(masks{i}), find(masks{i}));
    centroids(i,:) = (mean([x,y,z],1) + double(origins(i,:) - 1)) .* spacings(i,:);
end

function [enclosed_volume, surface_area] = calc_surface_volume_and_area(surfaces)

if isempty(surfaces), return; end
enclosed_volume = zeros(length(surfaces), 1, 'single');
surface_area = zeros(length(surfaces), 1, 'single');
for i = 1 : length(surfaces)
    e1 = surfaces(i).vertices(surfaces(i).faces(:,1),:) - surfaces(i).vertices(surfaces(i).faces(:,2),:);
    e3 = surfaces(i).vertices(surfaces(i).faces(:,3),:) - surfaces(i).vertices(surfaces(i).faces(:,1),:);
    [enclosed_volume(i,1), surface_area(i,1)] = stlVolumeNormals(surfaces(i).vertices', surfaces(i).faces', cross(e1,e3));
end


function goodIndex=RemoveBadNuclei(surfaces)
         goodIndex=[];
         for i=1:length(surfaces)
             if size(surfaces(i).vertices,1)>0
                 goodIndex=[goodIndex;i];
             end
         end
             


function goodcellindex=RemoveBadCells(centroid, volume,surfaces)
    tic
    [~,firstgood]=unique(centroid,'rows');
    %[length(centroid),length(firstgood)]
% spacing=[1, 1, 1];
% delta_spacing=2;
% if the user didn't define 'exclude_boundary_points' we set it to false:
if ~exist('exclude_boundary_points', 'var')
    exclude_boundary_points = false;
end
% calculating the triangulation and the volume of each triangle:
N=centroid(firstgood,:);

TRI = delaunay(N(:,1), N(:,2), N(:,3));
clear neighbor
edges=[];
count=1;
for i = 1 : size(N,1)
    temp=[];
    for j=1:size(TRI,1)
        for k=1:size(TRI,2)
            if TRI(j,k)==i
                temp=[temp,TRI(j,:)];
            end
        end
    end

   
    neighborList=setdiff(unique(temp),i);
    %neighbor(i,1)=length(neighborList{i});
    ids= neighborList;
    
    for j=1:length(ids)
        a=min(ids(j),i);
        b=max(ids(j),i);
        edges(count,:)=[firstgood(a),firstgood(b)];
        count=count+1;
    end
 
end

%size(edges)
[~,ia]=unique(edges,'rows');
edges=edges(ia,:);
%size(edges)

 badcell=[];
 for j=1:size(edges,1)
     
     dist=pdist(centroid(edges(j,:),:));
     if dist<8     
             cell1=surfaces{edges(j,1)}; actualVol(1)=volume( edges(j,1)   );      
             cell2=surfaces{edges(j,2)}; actualVol(2)=volume( edges(j,2)   );  
             combined=[cell1;cell2];
             [~,combV]=convhull(combined);

             volumeFraction= combV/sum(actualVol);
             if volumeFraction<1
                        [sa,sb]=min(actualVol);
                        if sb==1
                            badid=edges(j,1);
                        end
                        if sb==2
                            badid=edges(j,2);
                        end
                        badcell=[badcell;badid];
             end
     end
 end    

  goodcellindex=setdiff(firstgood,badcell);
  toc
  %length(goodcellindex)  
  


        

         


