function prediction_plot_from_n(n, chain, data, predict_equisetum, nsamples)
    % prediction_plot_from_n - ?18O vs ??17O with 95% CI, using n as input
    %
    % Inputs:
    %   chain             - [nsimu x 6] MCMC output matrix: [qk, ea, es, T, d18Omw, m]
    %   data.ydata        - [n x 2] matrix: [?18O, ??17O]
    %   data.sigma        - [n x 2] matrix: observational uncertainties
    %   predict_equisetum - function handle: @(n, par) ? [?18O, ??17O]
    %   nsamples          - number of random samples from the chain (e.g., 1000)

    % Number of segments
    n=round(n); % comment 
    % Preallocate
    delta18_all = NaN(nsamples, n);
    D17_all     = NaN(nsamples, n);
    % Calculate ?17O from each posterior sample
    lambda = 0.528;
     % Sample randomly from MCMC chain
    nsimu = size(chain, 1);
    sample_inds = randi(nsimu, nsamples, 1);

    for i = 1:nsamples
        par = chain(sample_inds(i), :);  % [qk, ea, es, T, d18Omw, n]
        par(:,6)=round(par(:,6));
        out = predict_equisetum(n, par);
        delta18_all(i, :) = out(:,1)';
        D17_all(i, :)     = out(:,2)';
    end
    
     delta17_all = 1000 * ( ...
    exp((D17_all ./ 1e6) + lambda .* log1p(delta18_all ./ 1000)) - 1 );
  

    % Compute 95% bounds across samples, pointwise per segment
    delta18_median = mean(delta18_all, 1);
    delta18_lower  = prctile(delta18_all, 2.5, 1);
    delta18_upper  = prctile(delta18_all, 97.5, 1);

    D17_median = mean(D17_all, 1);
    D17_lower  = prctile(D17_all, 2.5, 1);
    D17_upper  = prctile(D17_all, 97.5, 1);
    
    delta17_median = mean(delta17_all, 1);
    delta17_lower  = prctile(delta17_all, 2.5, 1);
    delta17_upper  = prctile(delta17_all, 97.5, 1);
    
    frac = linspace(0, 1, n); 
    
    % Plotting
    figure; hold on;

    % Confidence region
    fill_area(delta18_lower, D17_lower, delta18_upper, D17_upper, [0.8 0.8 1]);

    % Median prediction line
    plot(delta18_median, D17_median, 'b-', 'LineWidth', 2);

    % Raw data with error bars
    errorbar(data.ydata(:,1), data.ydata(:,2), ...
             data.sigma(:,2), 'ko', ...
             'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 0.75,'MarkerSize', 6, ...
             'CapSize', 3);

    xlabel('\delta^{18}O (?)', 'FontSize', 12);
    ylabel('\Delta''^{17}O (?)', 'FontSize', 12);
    title('\delta^{18}O vs \Delta''^{17}O with 95% Prediction Interval');
    legend({'95% Prediction Interval', 'Median Prediction', 'Observed Data'}, ...
           'Location', 'best');
    grid on;
     figure; hold on;
    fill_area(frac, delta18_lower, frac, delta18_upper, [1 0.8 0.8]);
    plot(frac, delta18_median, 'r-', 'LineWidth', 2);
    if isfield(data, 'fractional_positions') && ~isempty(data.fractional_positions)
    errorbar(data.fractional_positions, data.ydata(:,1), ...
             data.sigma(:,1), 'ko', ...
             'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', ...
             'LineWidth', 0.75, 'MarkerSize', 6, 'CapSize', 3);
    end
    xlabel('Fractional stem length', 'FontSize', 12);
    ylabel('\delta^{18}O (?)', 'FontSize', 12);
    title('\delta^{18}O vs fractional stem length');
    grid on;

    %% --- Plot 3: ?'17O vs fractional length
    figure; hold on;
    fill_area(frac, D17_lower, frac, D17_upper, [0.8 1 0.8]);
    plot(frac, D17_median, 'g-', 'LineWidth', 2);
    if isfield(data, 'fractional_positions') && ~isempty(data.fractional_positions)
    errorbar(data.fractional_positions, data.ydata(:,2), ...
             data.sigma(:,2), 'ko', ...
             'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', ...
             'LineWidth', 0.75, 'MarkerSize', 6, 'CapSize', 3);
    end
    xlabel('Fractional stem length', 'FontSize', 12);
    ylabel('\Delta''^{17}O (?)', 'FontSize', 12);
    title('\Delta''^{17}O vs fractional stem length');
    grid on;
        %% --- Plot 4: d17O vs fractional length    
    figure; hold on;
    fill_area(frac, delta17_lower, frac, delta17_upper, [0.8 0.9 1]);
    plot(frac, delta17_median, 'b-', 'LineWidth', 2);

    % Overlay reconstructed ?17O from data
    if isfield(data, 'fractional_positions') && ~isempty(data.fractional_positions)
        d18_obs = data.ydata(:,1);
        D17_obs = data.ydata(:,2);
        delta17_obs = 1000 * ( ...
        exp((D17_obs ./ 1e6) + lambda .* log1p(d18_obs ./ 1000)) - 1 );

        plot(data.fractional_positions, delta17_obs, 'ko', ...
                 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', ...
                 'LineWidth', 0.75, 'MarkerSize', 6);
    end

    xlabel('Fractional stem length', 'FontSize', 12);
    ylabel('\delta^{17}O (?)', 'FontSize', 12);
    title('\delta^{17}O vs fractional stem length');
    legend({'95% Prediction Interval', 'Median Prediction', 'Observed Data'}, ...
           'Location', 'best');
    grid on;
    
end

function fill_area(xmin, ymin, xmax, ymax, color)
    % fill_area - shade area between (xmin,ymin) and (xmax,ymax)
    fill([xmin, fliplr(xmax)], [ymin, fliplr(ymax)], color, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
