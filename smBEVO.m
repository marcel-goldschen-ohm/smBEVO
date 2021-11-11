function result = smBEVO(data, sigmaX, sigmaY, varargin)
% SMBC Baseline correction for single-molecule data series.
    % An image-based approach to rapidly remove baseline drift/wobble from
    % single-molecule data series (i.e. piecewise continuous series with a 
    % discrete set of amplitude levels).

    %% input parameters
    ip = inputParser;
    % series column data as [y] or [x y]
    addRequired(ip, 'data');
    % image filter sigmas
    addRequired(ip, 'sigmaX', @(x) x > 0);
    addRequired(ip, 'sigmaY', @(x) x > 0);
    % image resolution (pixels per sigma)
    addParameter(ip, 'pixelsPerSigmaX', 4, @(x) x > 0);
    addParameter(ip, 'pixelsPerSigmaY', 4, @(x) x > 0);
    % spline smoothing of level segments (0: no smoothing, >0: more smoothing)
    addParameter(ip, 'smoothing', 1, @(x) x >= 0);
    % minimum level separation
    addParameter(ip, 'minLevelSep', 0, @(x) x >= 0);
    % snake level refinement
    addParameter(ip, 'snakeLevelRefinement', false);
    addParameter(ip, 'maxSnakeIter', 50, @(x) x >= 0);
    addParameter(ip, 'alpha', 1, @(x) x >= 0); % 1st derivative (elasticity)
    addParameter(ip, 'beta', 1, @(x) x >= 0); % 2nd derivative (curvature)
    addParameter(ip, 'gamma', 1, @(x) x >= 0); % scale image gradient forcefield
    addParameter(ip, 'tol', 0.001, @(x) x >= 0);
    % initial guess for baseline
    addParameter(ip, 'ybaseline', zeros(length(data),1), @(x) length(x) == length(data));
    % application handle for aborting this function
    addParameter(ip, 'app', gobjects(0));
    % just return filtered image
    addParameter(ip, 'getFilteredImageOnly', false);
    parse(ip, data, sigmaX, sigmaY, varargin{:});
    
    % x,y data series
    if size(data,2) == 1
        xdata = reshape(1:length(data), [], 1);
        ydata = data;
    elseif size(data,2) == 2
        xdata = data(:,1);
        ydata = data(:,2);
    else
        error("Input 'data' must be [y] or [x y] column format.");
    end
    result.xdata = xdata;
    result.ydata = ydata;
    
    % parameters
    p.sigmaX = sigmaX;
    p.sigmaY = sigmaY;
    p.pixelsPerSigmaX = ip.Results.pixelsPerSigmaX;
    p.pixelsPerSigmaY = ip.Results.pixelsPerSigmaY;
    p.smoothing = ip.Results.smoothing;
    p.minLevelSep = ip.Results.minLevelSep;
    if p.minLevelSep == 0
        p.minLevelSep = 1.75 * p.sigmaY;
    end
    p.snakeLevelRefinement = ip.Results.snakeLevelRefinement;
    p.maxSnakeIter = ip.Results.maxSnakeIter;
    p.alpha = ip.Results.alpha;
    p.beta = ip.Results.beta;
    p.gamma = ip.Results.gamma;
    p.tol = ip.Results.tol;
    result.params = p;
    
    % other input
    ybaseline = reshape(ip.Results.ybaseline, [], 1);
    app = ip.Results.app;
    getFilteredImageOnly = ip.Results.getFilteredImageOnly;
    
    % for convenience
    sx = p.pixelsPerSigmaX;
    sy = p.pixelsPerSigmaY;

    %% data -> image
    ybaselined = ydata - ybaseline;
    xlims = [xdata(1), xdata(end)];
    ylims = [min(ybaselined), max(ybaselined)];
    pixelWidth = p.sigmaX / sx;
    pixelHeight = p.sigmaY / sy;
    nrows = ceil(diff(ylims) / pixelHeight);
    ncols = ceil(diff(xlims) / pixelWidth);
    % pad ylims for desired pixel height (pixel width will be shrunk to fit)
    padY = nrows * pixelHeight - diff(ylims);
    ylims = ylims + [-0.5, 0.5] * padY;
    colXEdges = linspace(xlims(1), xlims(2), ncols + 1);
    rowYEdges = linspace(ylims(1), ylims(2), nrows + 1);
    % pad top and bottom with extra bins
    pixelHeight = diff(rowYEdges(1:2));
    padNumRows = ceil(5 * sy);
    padYEdges = [1:padNumRows] * pixelHeight;
    rowYEdges = [rowYEdges(1) - fliplr(padYEdges), rowYEdges, rowYEdges(end) + padYEdges];
    nrows = numel(rowYEdges) - 1;
    im = histcounts2(xdata, ybaselined, colXEdges, rowYEdges)';
    colXCenters = (colXEdges(1:end-1) + colXEdges(2:end)) / 2;
    rowYCenters = (rowYEdges(1:end-1) + rowYEdges(2:end)) / 2;
    
%     % account for changing ylims with successive images
%     if iter > 1
%         imrow0(iter) = interp1(rowYCenters, 1:nrows, 0);
%     end

%     hold off; imagesc(colXCenters, rowYCenters, im); colormap(gray(256)); axis xy;

    result.colXEdges = colXEdges;
    result.rowYEdges = rowYEdges;
    result.colXCenters = colXCenters;
    result.rowYCenters = rowYCenters;

    %% check for abort
    if ~isempty(app)
        try
            drawnow; % flush event queue
            if app.abort
                return
            end
        catch
        end
    end

    %% filter image
    kx = -ceil(3 * sx):ceil(3 * sx);
    ky = -ceil(3 * sy):ceil(3 * sy);
    fx = exp(-kx.^2 ./ (2 * sx^2));
    fy = exp(-ky.^2 ./ (2 * sy^2)) .* (1 - ky.^2 ./ sy^2) ./ sy^2;
    filterKernel = fx .* fy';
    im = imfilter(im, filterKernel, 'replicate', 'conv');
    im = im / max(abs(im(:))); % scale within [-1, 1]

%     hold off; imagesc(colXCenters, rowYCenters, im); colormap(gray(256)); axis xy;

    result.filterKernel = filterKernel;
    result.im = im;

    if getFilteredImageOnly
        return
    end

    %% check for abort
    if ~isempty(app)
        try
            drawnow; % flush event queue
            if app.abort
                return
            end
        catch
        end
    end

    %% split image idealized traces into segments at jumps
    [~,immaxseq] = max(im, [], 1);
    [~,imminseq] = min(im, [], 1);

    % split idealized trace into segments at jumps
    jumpThreshold = sy / 2;
    colStarts = [1, 1 + find(abs(diff(immaxseq)) > jumpThreshold)];
    colStops = [colStarts(2:end)-1, ncols];
    colStarts2 = [1, 1 + find(abs(diff(imminseq)) > jumpThreshold)];
    colStops2 = [colStarts2(2:end)-1, ncols];
    colStarts = union(colStarts, colStarts2);
    colStops = [colStarts(2:end)-1, ncols];

    % constrain segment ends to segment values slighlty more interior
    % to avoid skew at the edges of the filtered regions
    ncap = ceil(sx / 2);
    for i0 = 1:numel(colStarts)
        cols = colStarts(i0):colStops(i0);
        n = numel(cols);
        if n > 2 * ncap
            immaxseq(cols(1:ncap)) = immaxseq(cols(ncap+1));
            immaxseq(cols(end-ncap+1:end)) = immaxseq(cols(end-ncap));
            imminseq(cols(1:ncap)) = imminseq(cols(ncap+1));
            imminseq(cols(end-ncap+1:end)) = imminseq(cols(end-ncap));
        else
            immaxseq(cols) = mean(immaxseq(cols));
            imminseq(cols) = mean(imminseq(cols));
        end
    end

%     % resplit idealized trace into segments at jumps
%     colStarts = [1, 1 + find(abs(diff(immaxseq)) > jumpThreshold)];
%     colStops = [colStarts(2:end)-1, ncols];
%     colStarts2 = [1, 1 + find(abs(diff(imminseq)) > jumpThreshold)];
%     colStops2 = [colStarts2(2:end)-1, ncols];
%     colStarts = union(colStarts, colStarts2);
%     colStops = [colStarts(2:end)-1, ncols];

    % smooth segments
    if p.smoothing
        try
            for i0 = 1:numel(colStarts)
                cols = colStarts(i0):colStops(i0);
                n = numel(cols);
                pp = splinefit(1:n, immaxseq(cols), ceil(n / (p.smoothing * sx)));
                immaxseq(cols) = ppval(pp, 1:n);
%                 pp = splinefit(1:n, imideal2(cols), ceil(n / (smoothing * sx)));
%                 imideal2(cols) = ppval(pp, 1:n);
            end
        catch err
            disp(err);
            warning("Segments NOT smoothed: Requires 'splinefit' Add-On by Jonas Lundgren.");
        end
    end

%     hold off; imagesc(im); colormap(gray(256)); axis xy; hold on;
%     plot(imideal2, 'm', 'linewidth', 0.5);
%     plot(imideal, 'c', 'linewidth', 2);
%     plot(colStarts, imideal(colStarts), 'o', 'linewidth', 2);
%     plot(colStops, imideal(colStops), 's', 'linewidth', 2);
%     for i = 1:numel(colStarts)
%         cols = colStarts(i):colStops(i);
%         rows = imideal(cols);
%         drows = 2 * pixelsPerSy;
%         x =[cols, fliplr(cols)];
%         y =[rows - drows, fliplr(rows + drows)];
%         patch(x, y, 'c', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%     end

    %% build image levels from ideal segments
    % start with the longest segment
    colWidths = colStops + 1 - colStarts;
    [~,i0] = max(colWidths);
    imlevels = nan(1,ncols);
    cols = colStarts(i0):colStops(i0);
    imlevels(cols) = immaxseq(cols);
    nstarts = numel(colStarts);
    niters = max(i0-1, nstarts-i0);
    minLevelSep_px = p.minLevelSep * (sy / p.sigmaY);
    for j = 1:niters
        a = i0-j;
        b = i0+j;
        if a >= 1
            acols = colStarts(a):colStops(a);
            d = immaxseq(acols(end)) - imlevels(:,acols(end)+1);
            [dmin,ilevel] = min(abs(d));
            if dmin < minLevelSep_px
                % connect to ilevel
                d = imlevels(:,acols(end)+1) - imlevels(ilevel,acols(end)+1);
                imlevels(:,acols) = immaxseq(acols) + d;
            else
                % add new level
                newlevel = imlevels(ilevel,:) + d(ilevel);
                if d(ilevel) >= 0
                    % add above ilevel
                    imlevels = [imlevels(1:ilevel,:); newlevel; imlevels(ilevel+1:end,:)];
                    d = [d(1:ilevel); 0; d(ilevel+1:end)];
                else
                    % add below ilevel
                    imlevels = [imlevels(1:ilevel-1,:); newlevel; imlevels(ilevel:end,:)];
                    d = [d(1:ilevel-1); 0; d(ilevel:end)];
                end
                imlevels(:,acols) = immaxseq(acols) - d;
            end
        end
        if b <= nstarts
            bcols = colStarts(b):colStops(b);
            d = immaxseq(bcols(1)) - imlevels(:,bcols(1)-1);
            [dmin,ilevel] = min(abs(d));
            if dmin < minLevelSep_px
                % connect to ilevel
                d = imlevels(:,bcols(1)-1) - imlevels(ilevel,bcols(1)-1);
                imlevels(:,bcols) = immaxseq(bcols) + d;
            else
                % add new level
                newlevel = imlevels(ilevel,:) + d(ilevel);
                if d(ilevel) >= 0
                    % add above ilevel
                    imlevels = [imlevels(1:ilevel,:); newlevel; imlevels(ilevel+1:end,:)];
                    d = [d(1:ilevel); 0; d(ilevel+1:end)];
                else
                    % add below ilevel
                    imlevels = [imlevels(1:ilevel-1,:); newlevel; imlevels(ilevel:end,:)];
                    d = [d(1:ilevel-1); 0; d(ilevel:end)];
                end
                imlevels(:,bcols) = immaxseq(bcols) - d;
            end
        end
        % smooth levels as we go
        if p.smoothing
            try
                cols = colStarts(max(1,a)):colStops(min(b,nstarts));
                n = numel(cols);
                pp = splinefit(1:n, imlevels(1,cols), ceil(n / (p.smoothing * sx)));
                imlevels(:,cols) = ppval(pp, 1:n) + (imlevels(:,cols(1)) - imlevels(1,cols(1)));
            catch err
                disp(err);
                warning("Segments NOT smoothed: Requires 'splinefit' Add-On by Jonas Lundgren.");
            end
        end
        % check for abort
        if ~isempty(app)
            try
                drawnow; % flush event queue
                if app.abort
                    return
                end
            catch
            end
        end

%         hold off; imagesc(im); colormap(gray(256)); axis xy;
%         hold on; plot(imideal, 'c', 'linewidth', 2);
%         plot((1:ncols)', imlevels', 'm-', 'linewidth', 2);
%         waitforbuttonpress
    end
    nlevels = size(imlevels, 1);

%     hold off; imagesc(im); colormap(gray(256)); axis xy;
%     hold on; plot(imideal, 'c', 'linewidth', 2);
%     plot((1:ncols)', imlevels', 'y-', 'linewidth', 1);
    
    %% refine levels with snakes
    if p.snakeLevelRefinement && p.maxSnakeIter > 0
        N = ncols;
        a = p.gamma * (2 * p.alpha + 6 * p.beta) + 1;
        b = p.gamma * (-p.alpha - 4 * p.beta);
        c = p.gamma * p.beta;
        P = diag(repmat(a,1,N));
        P = P + diag(repmat(b,1,N-1), 1) + diag(   b, -N+1);
        P = P + diag(repmat(b,1,N-1),-1) + diag(   b,  N-1);
        P = P + diag(repmat(c,1,N-2), 2) + diag([c,c],-N+2);
        P = P + diag(repmat(c,1,N-2),-2) + diag([c,c], N-2);
        P(1,:) = 0; 
        P(1,1) = 1;
        P(2,:) = 0; 
        P(2,1:3) = p.gamma*[-p.alpha,1+2*p.alpha,-p.alpha];
        P(end,:) = 0; 
        P(end,end) = 1;
        P(end-1,:) = 0;
        P(end-1,end-2:end) = p.gamma*[-p.alpha,1+2*p.alpha,-p.alpha];
        Pinv = inv(P);

%         hold off; imagesc(im); colormap(gray(256)); axis xy;
%         hold on; plot(imideal, 'c', 'linewidth', 2);
%         plot((1:ncols)', imlevels', 'y-', 'linewidth', 1);

        % external forcefield is image gradient along Y-axis
        im(im < 0) = 0; % ignore valleys in case they are interleaved with peaks
        [~,gy] = imgradientxy(im);
        Fext = gy ./ max(abs(gy(:))); % scale within [-1, 1]

        imlevels = imlevels';
        nlevels = size(imlevels, 2);
        Fy = zeros(N, nlevels);
        changes = zeros(1, p.maxSnakeIter);
        for iter = 1:p.maxSnakeIter
            for i = 1:N
                for j = 1:nlevels
                    try
                        Fy(i,j) = Fext(round(imlevels(i,j)), i);
                    catch
                        % in case level is off of image
                        Fy(i,j) = 0;
                    end
                end
            end
            % AB - 9/7/21 - added convergence check and max move per
            % iteration
            oldlevels = imlevels;
            move = tanh(p.gamma*mean(Fy,2));
            imlevels = Pinv * (imlevels + move);
%             imlevels = P\(imlevels + p.gamma * mean(Fy, 2));
            if p.tol > 0
                max_change = max(abs(oldlevels(:, 1)-imlevels(:,1)));
                changes(iter) = max_change;
                if iter >= 10
                    params = polyfit(1:10, changes(iter-9:iter), 1); 
                    slope = params(1);
                    if abs(slope) < p.tol * sy
                        break
                    end
                end
            end
        end
        imlevels = imlevels';
    
%         plot((1:ncols)', imlevels', 'm-', 'linewidth', 2);
    end
    
    %% baseline and levels in data and image coords
    [~,bi] = min(imlevels(:,1));
    imbaseline = imlevels(bi,:);
    baselineYNodes = interp1(1:nrows, rowYCenters, imbaseline, 'linear', 'extrap');
    ybaseline = ybaseline + interp1(colXCenters', baselineYNodes', xdata, 'linear', 'extrap');
    
    imlevelOffsets = imlevels(:,1) - imlevels(bi,1);
    ylevelOffsets = interp1(1:nrows, rowYCenters, imlevelOffsets' - min(imlevelOffsets) + 1, 'linear', 'extrap');
    ylevelOffsets = ylevelOffsets - ylevelOffsets(bi);
    
%     % account for changing ylims with successive images
%     imbaseline = sum(imbaseline - reshape(imrow0, [], 1), 1);
    
%     hold off; imagesc(result.im); colormap(gray(256)); axis xy;
%     hold on; plot((1:ncols)', (imbaseline + imlevelOffsets)', 'm-', 'linewidth', 2);
    
%     hold off; plot(xdata, ydata - ybaseline);
%     hold on; plot(xdata([1,end]), repmat(ylevelOffsets, 2, 1), 'k--');

    %% idealized sequence in data and image coords
    [~,imstates] = min(abs(imlevels - repmat(immaxseq, nlevels, 1)), [], 1);
%     imideal = zeros(size(imstates));
%     for col = 1:ncols
%         imideal(col) = imlevels(imstates(col),col);
%     end
    ystates = round(interp1(colXCenters', imstates', xdata, 'linear', 'extrap'));
    
%     hold off; imagesc(result.im); colormap(gray(256)); axis xy;
%     hold on; plot((1:ncols)', (imbaseline + imlevelOffsets)', 'm-', 'linewidth', 2);
%     hold on; plot(imideal, 'c-', 'linewidth', 2);
    
%     hold off; plot(xdata, ydata - ybaseline);
%     hold on; plot(xdata([1,end]), repmat(ylevelOffsets, 2, 1), 'k--');
%     hold on; plot(xdata, yideal);
    
    %% output
    result.ybaseline = ybaseline;
    result.ylevelOffsets = ylevelOffsets;
    result.ystates = ystates;
    result.imbaseline = imbaseline;
    result.imlevelOffsets = imlevelOffsets;
    
end