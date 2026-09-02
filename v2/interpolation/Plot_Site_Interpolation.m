function [result,figure_handle] = Plot_Site_Interpolation(site,points, ...
        attribute,sample_stage,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% One kriged surface over the whole site, drawn the way the notebook drew it.
%%
%% This is the figure plot_kriged_attribute produced: a single ordinary
%% kriging model fitted to every sample inside the selected fields, one grid
%% over their common bounds masked to the union, the field outlines over the
%% top and the measurements on the same colour scale, in Easting and
%% Northing. It is the whole-site view, not the per-polygon simulation grid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

script_directory = fileparts(mfilename('fullpath'));
v2_directory = fileparts(script_directory);

parser = inputParser;
parser.addParameter('StageColumn','Sample Type');
parser.addParameter('XColumn','lng (X)');
parser.addParameter('YColumn','lat (Y)');
parser.addParameter('VariogramModel','exponential');
parser.addParameter('VariogramParameters',struct());
parser.addParameter('GridResolution',200);
parser.addParameter('OutputName','');
parser.addParameter('Visible','off');
parser.parse(varargin{:});
settings = parser.Results;

result = Interpolate_Field_Attribute(points,site,attribute,sample_stage, ...
    'StageColumn',settings.StageColumn,'XColumn',settings.XColumn, ...
    'YColumn',settings.YColumn,'VariogramModel',settings.VariogramModel, ...
    'VariogramParameters',settings.VariogramParameters, ...
    'GridResolution',settings.GridResolution);

figure_handle = figure('Color','w','Position',[100 100 1000 850], ...
    'Visible',settings.Visible);
axes_handle = axes(figure_handle);
axes_handle.Tag = 'SiteInterpolationAxes';
% Set the light palette explicitly; an unstyled axes picks up the session
% theme and exports the masked area as black.
axes_handle.Color = 'w';
axes_handle.XColor = 'k';
axes_handle.YColor = 'k';
hold(axes_handle,'on');
[surface_x,surface_y] = meshgrid(result.grid_x,result.grid_y);
surface(axes_handle,surface_x,surface_y,zeros(size(surface_x)), ...
    result.prediction,'EdgeColor','none','FaceColor','flat', ...
    'Tag','KrigedSurface');
for k = 1:numel(site.feature_polygons)
    outline = site.feature_polygons{k}.Vertices;
    plot(axes_handle,[outline(:,1);outline(1,1)], ...
        [outline(:,2);outline(1,2)],'-','Color',[0.5 0.5 0.5], ...
        'LineWidth',2,'Tag',sprintf('FieldOutline%02d',k));
end
scatter(axes_handle,result.samples.x_grid_m,result.samples.y_grid_m,35, ...
    result.sample_values,'filled','MarkerEdgeColor',[0.2 0.2 0.2], ...
    'LineWidth',0.6,'Tag','SamplePoint');
axis(axes_handle,'equal','tight');
box(axes_handle,'on');
colormap(figure_handle,parula);
colorbar_handle = colorbar(axes_handle);
colorbar_handle.Label.String = attribute;
colorbar_handle.Color = 'k';
colorbar_handle.Label.Color = 'k';
xlabel(axes_handle,'Easting','Color','k');
ylabel(axes_handle,'Northing','Color','k');
title(axes_handle,sprintf('%s at %s',attribute,sample_stage), ...
    'FontWeight','normal','Color','k');
subtitle(axes_handle,sprintf('Ordinary kriging (%s), %d samples', ...
    char(settings.VariogramModel),result.inside_sample_count),'Color','k');

output_name = settings.OutputName;
if isempty(output_name)
    [~,shape_name] = fileparts(site.source_shapefile);
    output_name = sprintf('Site_Interpolation_%s_%s_%s', ...
        regexprep(shape_name,'[^A-Za-z0-9]+','_'), ...
        regexprep(char(attribute),'[^A-Za-z0-9]+','_'), ...
        regexprep(char(sample_stage),'[^A-Za-z0-9]+','_'));
end
figure_directory = fullfile(v2_directory,'figures');
if ~isfolder(figure_directory); mkdir(figure_directory); end
exportgraphics(figure_handle,fullfile(figure_directory, ...
    [output_name '.png']),'Resolution',300);
savefig(figure_handle,fullfile(figure_directory,[output_name '.fig']));
result.output_name = output_name;

end
