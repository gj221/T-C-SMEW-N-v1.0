function T = Extract_TDep_Dry_Deposition(tdep_dir, lat, lon, years, species)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract_TDep_Dry_Deposition
%%
%% Annual dry deposition at a point from the NADP Total Deposition (TDep)
%% grids, as EPA serves them (gaftp.epa.gov/castnet/tdep/CURRENT_grids,
%% <species>_dw-<year>.zip; TDep MMF 2025.01: 4 km cells, NAD83 / Conus
%% Albers).  The value is that of the 4 km cell holding (lat, lon).  Units
%% are TDep's, kg of the element per ha per year: n_dw kg N (HNO3, NO3 and
%% NH4 particles, NH3, and the gases CMAQ models: NO2, PAN, N2O5, HONO,
%% organic nitrate), s_dw kg S (SO2 + SO4), ca, mg, k, na and cl_dw kg of
%% the ion.  Dry deposition is measured air concentration (CASTNET) times
%% modelled deposition velocity (CMAQ).
%%
%% species default {'n','s','ca','mg','k','na','cl'}.  T has one row per
%% year and a column <species>_dw_kg_ha per species; T.Properties.UserData
%% records the cell (row, column, centre) and the grid's projection.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5 || isempty(species)
    species = {'n','s','ca','mg','k','na','cl'};
end
years = double(years(:));
values = nan(numel(years),numel(species));
cell_rc = [];
scratch = tempname; mkdir(scratch);
cleanup = onCleanup(@() rmdir(scratch,'s'));
for k = 1:numel(species)
    for y = 1:numel(years)
        name = sprintf('%s_dw-%d',species{k},years(y));
        zip_file = fullfile(tdep_dir,[name '.zip']);
        if ~isfile(zip_file)
            error('Extract_TDep_Dry_Deposition:MissingFile','%s is missing.',zip_file);
        end
        files = unzip(zip_file,fullfile(scratch,name));
        tif = files(endsWith(lower(files),'.tif'));
        if numel(tif) ~= 1
            error('Extract_TDep_Dry_Deposition:BadArchive', ...
                '%s holds %d GeoTIFFs, not one.',zip_file,numel(tif));
        end
        [A,R] = readgeoraster(tif{1});
        [x,yw] = projfwd(R.ProjectedCRS,lat,lon);
        [row,col] = worldToDiscrete(R,x,yw);
        if isnan(row) || isnan(col)
            error('Extract_TDep_Dry_Deposition:OutsideGrid', ...
                '(%g, %g) lies outside %s.',lat,lon,zip_file);
        end
        if isempty(cell_rc)
            cell_rc = [row col];
            [cx,cy] = intrinsicToWorld(R,col,row);
            crs_name = char(R.ProjectedCRS.Name);
        elseif ~isequal(cell_rc,[row col])
            error('Extract_TDep_Dry_Deposition:GridChanged', ...
                '%s puts the point in another cell than the first grid.',zip_file);
        end
        v = double(A(row,col));
        if ~(isfinite(v) && v >= 0 && v < 1e6)
            error('Extract_TDep_Dry_Deposition:NoValue', ...
                '%s has no value at (%g, %g).',zip_file,lat,lon);
        end
        values(y,k) = v;
    end
end
T = array2table(values,'VariableNames',strcat(species,'_dw_kg_ha'));
T = addvars(T,years,'Before',1,'NewVariableNames','year');
T.Properties.UserData = struct('cell_row',cell_rc(1),'cell_col',cell_rc(2), ...
    'cell_centre_m',[cx cy],'crs',crs_name,'cell_size_m',R.CellExtentInWorldX, ...
    'point_lat_lon',[lat lon]);
end
