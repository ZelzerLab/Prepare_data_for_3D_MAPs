function S = join_surfaces(surfaces)

% joining all cells/nuclei into the same image:
vertices_len = arrayfun(@(x) length(x.vertices), surfaces);
faces_len = arrayfun(@(x) length(x.faces), surfaces);
vertices_cumlen = [0; cumsum(vertices_len(:))];
faces_cumlen = [0; cumsum(faces_len(:))];
 
% updating the indices of the faces in each induvidual surface:
for i = 2 : length(surfaces)
    surfaces(i).faces = surfaces(i).faces + vertices_cumlen(i);
end
 
% allocating memory for the joint surface:
S.vertices = nan(sum(vertices_len),3);
S.faces = nan(sum(faces_len),3);
 
% joining the surfaces:
for i = 1 : length(vertices_cumlen)-1
    S.vertices(vertices_cumlen(i)+1:vertices_cumlen(i+1),:) = surfaces(i).vertices;
    S.faces(faces_cumlen(i)+1:faces_cumlen(i+1),:) = surfaces(i).faces;
end
