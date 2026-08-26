function [cell_size_m,calibration,oriented_polygon] = Calibrate_Field_Grid_Cell_Size(field_polygon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Derive cell width and grid orientation from the polygon convex envelope.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[oriented_polygon,calibration] = ...
    Orient_Field_Polygon_Grid(field_polygon,[]);
cell_size_m = calibration.cell_size_m;
calibration.mode = 'polygon_geometry';

end
