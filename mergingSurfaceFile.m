


clear all; 

dirname = 'G:\Labs\EliZelzer\sarahru\Analyzed_Data\Nuclei_segmentation_gdf5Mutantandcontrol_data\Newest_XPIWIT\Data\S114\wt_DV\cells\matlab_input\redo\New folder';
datafiles = dir(dirname);

 filecount=0;
 for i=1:length(datafiles)
        t=strfind(datafiles(i).name,'(Surfaces).mat');
        if length(t)~=0
            filecount=filecount+1;
            SurfaceFilename{filecount}=datafiles(i).name;
        end
 end

bw=[];
coords=[];
surfaces=[];

% [filename, dirname] = get_file_name('* (Surfaces).mat', 'Please choose the .mat file of the MEMBRANES:', 'on');
% surfaces = struct('bw', [], 'coords', zeros(0,6), 'spacing', zeros(0,3), 'surfaces', []);
for i = 1 : filecount
    disp(SurfaceFilename{i});
    curr_surfaces= load(fullfile(dirname, SurfaceFilename{i}));
    bw = [bw; curr_surfaces.bw];
    coords = [coords; curr_surfaces.coords];
    surfaces = [surfaces; curr_surfaces.surfaces]; 
end

spacing = [curr_surfaces.spacing];
% or you can write manually
%merge.spacing = [1,1,1]; 

save(fullfile(dirname, SurfaceFilename{i}), 'bw','coords','surfaces','spacing');
   


%     
%    C = struct('bw', [], 'coords', zeros(0,6), 'spacing', zeros(0,3), 'surfaces', []);
%    for i = 1 : length(filename)
%     disp(filename{i});
%     curr_C = load(fullfile(dirname, filename{i}));
%     C.bw = [C.bw; curr_C.bw];
%      C.coords = [C.coords; curr_C.coords];
%        C.spacing = [C.spacing; repmat(curr_C.spacing, length(curr_C.bw), 1)];
%         C.surfaces = [C.surfaces; curr_C.surfaces]; 
%    end  
%     
%  surfaces = C.surfaces;
%    spacing = C.spacing;
%    bw = C.bw;
%    coords = C.coords;
  
    
      