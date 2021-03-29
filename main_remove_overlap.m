
clear all 
close all 


ignore_cells = false;
ignore_edu = true;
%input_directory_name_of_data, all surfaces files of cells, nuclei and edu
%should be one folder 
data_dir_name='F:\Cell Project\Figures for paper\data';
disp('Cell (Surfaces) files start with "cells_" keyword')
disp('Nuclei (Surfaces) files start with "nuclei_" keyword')
disp('Edu files start with "Edu_" keyword') 


datafiles = dir(data_dir_name);
cell_count=0;
nuclei_count=0;
edu_count=0;

if ~ignore_cells
    for i=1:length(datafiles)
        t=strfind(datafiles(i).name,'cells_');
        if length(t)~=0
            cell_count=cell_count+1;
            cells_filename{cell_count}=datafiles(i).name;
        end
    end
end

for i=1:length(datafiles)
	t=strfind(datafiles(i).name,'nuclei_');
	if length(t)~=0
		nuclei_count=nuclei_count+1;
		nuclei_filename{nuclei_count}=datafiles(i).name;
    end
end

if ~ignore_edu
    for i=1:length(datafiles)
        t=strfind(datafiles(i).name,'Edu_');
        if length(t)~=0
            edu_count=edu_count+1;
            edu_filename{edu_count}=datafiles(i).name;
        end
    end
end


disp('...')

if ~ignore_cells
	if cell_count==nuclei_count
		disp('number of cell and nuclei datafiles are correct')
		if ~ignore_edu
		    if edu_count==nuclei_count
		        disp('number of cell, nuclei and edu datafiles are correct')
		    else
		        disp('number of cell, nuclei and edu datafiles are incorrect')
		        return;
		    end
		end 
         else
         	disp('number of cell, nuclei and edu datafiles are incorrect')
            	return;
    	end
end

[cell_count,nuclei_count,edu_count]
disp('If cell, nuclei, and edu files are correct')
disp('press any key to continue')
pause

disp('...')	
disp('reading files')

if ~ignore_cells
    disp('...')
    disp('check the order of cell');
    for i = 1 : length(cells_filename)
            disp(cells_filename{i})
            pos=strfind(cells_filename{i},'_');
            fileNameId{i}=cells_filename{i}(pos(2)+1:pos(3)-1);
    end
end

disp('...')
disp('check the order of nuclei')
for i = 1 : length(nuclei_filename)
    disp(nuclei_filename{i})
    pos=strfind(nuclei_filename{i},'_');
    fileNameId{i}=nuclei_filename{i}(pos(2)+1:pos(3)-1);
end

if ~ignore_edu
    disp('...')
    disp('check the order of edu')
    for i = 1 : length(edu_filename)
        disp(edu_filename{i})
    end
end

flag=1;
while(1)
    disp('...')
    disp('If the order of files are okay then please press 1 to continue and 0 to stop');
    m=input('','s');
    if m=='0'
        flag=0;
        break;
    end
    if m=='1'
        break;
    end
end


% shift 8 in integer z position or 3 micron shift in real coordinate 
% Shift in nuclei are in -Z direction. So inside the script instead of nuclei cell got shifted to +Z direction 
% If you do not want any shift then use 0 
nuclei_shift=0;
    
if flag==1
   for i=1:nuclei_count
            disp(strcat('Start ',num2str(i)))
            output=strcat('c_n_pos',fileNameId{i});
            %output=strcat('nuclei_pos',nuclei_filename{i});
            if ignore_cells
                cells_filename{i}='';
            end
            if ignore_edu    
                edu_filename{i}='';
            end
	     G = edu_add_calc_morphological_characteristics_remove_overlap(ignore_cells,ignore_edu,data_dir_name, ...
                cells_filename{i},nuclei_filename{i},edu_filename{i},{output},nuclei_shift);
   end
end






