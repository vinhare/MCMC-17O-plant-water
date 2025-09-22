function ss = ssfun(par, data)
    % Get full model output
    ymod = model_equisetum3(par, data);

    % Observed values
    yobs = data.ydata;
    sig  = data.sigma;

    % Number of observations
    N = size(yobs, 1);
    nmodel = size(ymod, 1);

    % Subsample N evenly spaced points from model output
    if nmodel < N
        ss = Inf;
        return
    end

%     idx = round(linspace(1, nmodel, N));
%     ymod_sub = ymod(idx, :);
    fpos = data.fractional_positions(:);  % ensure column vector of length N
    x_model = linspace(0, 1, size(ymod, 1))';  % fractional positions in model

    % Interpolate model output to match data fractional positions
    delta18_interp = interp1(x_model, ymod(:,1), fpos, 'linear');
    D17_interp     = interp1(x_model, ymod(:,2), fpos, 'linear');
    ymod_sub = [delta18_interp, D17_interp];
    
    % Residuals
    res = (yobs - ymod_sub) ./ sig;

    % Sum of squares
    ss = sum(res(:).^2);
end

% 
% function ss = ssfun(par, data)
%     % Get model predictions: [?18O, D'17O]
%     ymod = model_equisetum3(par, data);
% 
%     % Extract observed values
%     yobs = data.ydata;     % Matrix: [?18O_obs, D'17O_obs]
%     sig  = data.sigma;     % Matrix: [?18, ?17]
% 
%     % Residuals
%     res = (yobs - ymod) ./ sig;
% 
%     % Sum of squares
%     ss = sum(res(:).^2);
% end
