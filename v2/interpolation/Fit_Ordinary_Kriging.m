function model = Fit_Ordinary_Kriging(x_m,y_m,values,variogram_model, ...
        variogram_parameters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assemble the ordinary kriging system for scattered measurements.
%%
%% This is the MATLAB form of pykrige.ok.OrdinaryKriging with an explicit
%% variogram. The unbiasedness row and column are appended to the
%% semivariance matrix and its diagonal is zeroed, so a sample carries no
%% nugget against itself. Left unset, the variogram parameters take the
%% notebook recipe: sill from the population variance, an 80 m range and a
%% nugget of 5% of the variance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4 || isempty(variogram_model); variogram_model = 'exponential'; end
if nargin < 5; variogram_parameters = struct(); end

x_m = reshape(double(x_m),[],1);
y_m = reshape(double(y_m),[],1);
values = reshape(double(values),[],1);
if ~isequal(numel(x_m),numel(y_m),numel(values))
    error('Fit_Ordinary_Kriging:SizeMismatch', ...
        'x_m, y_m and values must hold the same number of samples.');
end
if ~all(isfinite(x_m) & isfinite(y_m) & isfinite(values))
    error('Fit_Ordinary_Kriging:NonFiniteSample', ...
        'Every sample must have a finite coordinate pair and value.');
end
n = numel(values);
if n < 3
    error('Fit_Ordinary_Kriging:TooFewSamples', ...
        'Ordinary kriging needs at least 3 samples, not %d.',n);
end

variogram_model = lower(char(variogram_model));
model_parameters = resolve_parameters( ...
    variogram_model,variogram_parameters,values);

separation_m = hypot(x_m-x_m.',y_m-y_m.');
system_matrix = zeros(n+1);
system_matrix(1:n,1:n) = Variogram_Semivariance( ...
    variogram_model,model_parameters,separation_m);
system_matrix(1:n+2:n*(n+1)) = 0;
system_matrix(n+1,1:n) = 1;
system_matrix(1:n,n+1) = 1;
system_matrix(n+1,n+1) = 0;

model.x_m = x_m;
model.y_m = y_m;
model.values = values;
model.sample_count = n;
model.variogram_model = variogram_model;
model.model_parameters = model_parameters;
model.variogram_parameters = describe_parameters( ...
    variogram_model,model_parameters);
model.system_matrix = system_matrix;
model.exact_tolerance_m = 1e-10;

end

function model_parameters = resolve_parameters( ...
        variogram_model,variogram_parameters,values)
if isnumeric(variogram_parameters)
    model_parameters = double(variogram_parameters(:)).';
    return
end
if ~isstruct(variogram_parameters) || ~isscalar(variogram_parameters)
    error('Fit_Ordinary_Kriging:BadParameters', ...
        'variogram_parameters must be a scalar struct or a parameter vector.');
end

% numpy takes the population variance, so var(...,1) is the matching form.
value_variance = var(values,1);
switch variogram_model
    case 'linear'
        slope = field_or_default(variogram_parameters,'slope', ...
            value_variance/80);
        nugget = field_or_default(variogram_parameters,'nugget', ...
            0.05*value_variance);
        model_parameters = [slope nugget];
    case 'power'
        scale = field_or_default(variogram_parameters,'scale', ...
            value_variance/80);
        exponent = field_or_default(variogram_parameters,'exponent',1);
        nugget = field_or_default(variogram_parameters,'nugget', ...
            0.05*value_variance);
        model_parameters = [scale exponent nugget];
    otherwise
        sill = field_or_default(variogram_parameters,'sill',value_variance);
        range_m = field_or_default(variogram_parameters,'range',80);
        nugget = field_or_default(variogram_parameters,'nugget', ...
            0.05*value_variance);
        if nugget > sill
            error('Fit_Ordinary_Kriging:NuggetAboveSill', ...
                'The nugget must not exceed the sill.');
        end
        % pykrige stores the partial sill, not the sill.
        model_parameters = [sill-nugget range_m nugget];
end
if ~all(isfinite(model_parameters))
    error('Fit_Ordinary_Kriging:NonFiniteParameters', ...
        'The resolved variogram parameters must all be finite.');
end
end

function value = field_or_default(parameters,name,default_value)
if isfield(parameters,name) && ~isempty(parameters.(name))
    value = double(parameters.(name));
else
    value = default_value;
end
end

function described = describe_parameters(variogram_model,model_parameters)
switch variogram_model
    case 'linear'
        described = struct('slope',model_parameters(1), ...
            'nugget',model_parameters(2));
    case 'power'
        described = struct('scale',model_parameters(1), ...
            'exponent',model_parameters(2),'nugget',model_parameters(3));
    otherwise
        described = struct('partial_sill',model_parameters(1), ...
            'range',model_parameters(2),'nugget',model_parameters(3), ...
            'sill',model_parameters(1)+model_parameters(3));
end
end
