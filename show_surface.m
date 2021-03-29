function show_surface(S, colors, axis_excess)

if ~exist('axis_excess', 'var')
    axis_excess = 0.25;
end

if length(S) > 1
    S = join_surfaces(S);
end

if exist('colors', 'var') && isempty(colors)
    clear colors;
end

% prev_opengl = opengl('data');
% if ~strcmp(prev_opengl.HardwareSupportLevel, 'hardware')
%     opengl('Hardware');
% end

f = figure('Color', [1 1 1]);

if isfield(S, 'filename')
    [~, fname] = fileparts(S.filename);
    set(f, 'name', fname);
    S = rmfield(S, 'filename');
end

if exist('colors', 'var')
    all_colors = jet(1000);
    for i = 1 : 3
        scaled_colors(:,i) = interp1(linspace(min(colors), max(colors), 1000), all_colors(:,i), colors); %#ok<AGROW>
    end
    
    hold on;
    for i = 1 : length(S)
        patch(S(i), 'FaceVertexCdata', colors(i) * ones(size(S(i).faces,1),1), 'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    end
    hold off;
    
    colormap jet;
    colorbar;
else
    patch(S, 'FaceColor', [1 1 1]*0.8, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
axis off;
axis image;
axis vis3d;

cameratoolbar;
daspect([1,1,1]);

xl = xlim;
yl = ylim;
zl = zlim;

xlim([xl(1) - range(xl) * axis_excess, xl(2) + range(xl) * axis_excess]);
ylim([yl(1) - range(yl) * axis_excess, yl(2) + range(yl) * axis_excess]);
zlim([zl(1) - range(zl) * axis_excess, zl(2) + range(zl) * axis_excess]);

% if strcmp(prev_opengl.HardwareSupportLevel, 'basic')
%     opengl hardwarebasic;
% elseif strcmp(prev_opengl.HardwareSupportLevel, 'none')
%     opengl software;
% end

