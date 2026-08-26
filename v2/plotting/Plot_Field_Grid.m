function [topo,figure_handle] = Plot_Field_Grid(shape_file,cell_size_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the rotated full simulation grid and exact polygon overlap.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

script_directory = fileparts(mfilename('fullpath'));
v2_directory = fileparts(script_directory);
code_directory = fullfile(v2_directory,'T&C_Code');
addpath(code_directory);

if nargin < 1 || isempty(shape_file)
    shape_file = fullfile(v2_directory,'Spatial inputs', ...
        'Bruce Farm Virginia','applied_yield_shp.shp');
end
if nargin < 2; cell_size_m = []; end

topo = Build_Shapefile_Grid_Topology(shape_file,cell_size_m);
figure_handle = figure('Color','w','Position',[100 100 700 900]);
axes_handle = axes(figure_handle,'Position',[0.16 0.09 0.60 0.82], ...
    'FontSize',10,'XColor','k','YColor','k');
hold(axes_handle,'on');
for c = 1:topo.N_cells
    x = [topo.cell_x_lo(c) topo.cell_x_hi(c) ...
        topo.cell_x_hi(c) topo.cell_x_lo(c)];
    y = [topo.cell_y_lo(c) topo.cell_y_lo(c) ...
        topo.cell_y_hi(c) topo.cell_y_hi(c)];
    patch(axes_handle,x,y,topo.cell_overlap_fraction(c), ...
        'EdgeColor',[0.35 0.35 0.35],'LineWidth',0.7);
end
plot(axes_handle,topo.field_boundary_x_m,topo.field_boundary_y_m, ...
    'k-','LineWidth',1.8);
axis(axes_handle,'equal','tight');
box(axes_handle,'on');
colormap(axes_handle,parula);
clim(axes_handle,[0 1]);
colorbar_handle = colorbar(axes_handle,'Ticks',0:0.2:1);
colorbar_handle.Label.String = 'Cell area inside field [-]';
colorbar_handle.Color = 'k';
colorbar_handle.Label.Color = 'k';
xlabel(axes_handle,'Rotated grid x [m]','Color','k');
ylabel(axes_handle,'Rotated grid y [m]','Color','k');
title(axes_handle,sprintf('%d cells at %.2f m, rotation %.2f degrees', ...
    topo.N_cells,topo.cell_size_m,topo.grid_rotation_degrees), ...
    'Color','k','FontSize',12,'FontWeight','normal');
subtitle(axes_handle,sprintf( ...
    'Field %.0f m^2, simulated %.0f m^2, outside %.0f m^2', ...
    topo.field_area_m2,topo.simulated_area_m2, ...
    topo.grid_calibration.outside_grid_area_m2), ...
    'Color','k','FontSize',10);

[~,shape_name] = fileparts(shape_file);
shape_name = regexprep(shape_name,'[^A-Za-z0-9]+','_');
if isempty(cell_size_m)
    resolution = 'geometry';
else
    resolution = sprintf('cell%gm',cell_size_m);
    resolution = strrep(resolution,'.','p');
end
output_name = sprintf('Field_Grid_%s_%s',shape_name,resolution);
figure_directory = fullfile(v2_directory,'figures');
if ~isfolder(figure_directory); mkdir(figure_directory); end
exportgraphics(figure_handle,fullfile(figure_directory,[output_name '.png']), ...
    'Resolution',300);
savefig(figure_handle,fullfile(figure_directory,[output_name '.fig']));

end
