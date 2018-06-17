%%  How to infer connectivity from time series? A basic example. 
%  Patrick Desrosiers, DCClab, 2016-09.

%% Generate surrogate time series
% Time series are matrices of dim (p,n), where p is the population size
% (number of cells or nodes), n is the sample size (number of time steps)
    clear;
    numVariables = 12;      % p 
    numTimeSteps = 1000;    % n
    % Define random time series
    timeSeriesNodes1to4 = sin(12*rand(4,numTimeSteps));
    timeSeriesNodes5to8 = rand(size(timeSeriesNodes1to4));
    timeSeriesNodes9to12 = cos(pi*rand(size(timeSeriesNodes1to4)));
    timeSeries = [timeSeriesNodes1to4 ; timeSeriesNodes5to8 ; ...
        timeSeriesNodes9to12];
    % Introduce causality
    % In the following model 1 -> 2 & 3 & 10, 2 -> 3, 9 -> 11, 10-> 11
    timeSeries(2,2:end) = ...
        2*timeSeries(1,1:end-1)+0.1*randn(1,numTimeSteps-1); 
    timeSeries(3,3:end) = ...
        -timeSeries(2,2:end-1) + timeSeries(1,1:end-2);
    timeSeries(10,2:end) = ...
        2*timeSeries(1,1:end-1)+0.1*randn(1, numTimeSteps-1); 
    timeSeries(11,3:end) = ...
        -timeSeries(10,2:end-1) + timeSeries(9,1:end-2);


%% Visualize time series
    figure(1)
        colormap jet;
        imagesc(timeSeries);
        colorbar;
        xlabel('time step','FontSize', 12);
        ylabel('node id','FontSize', 12);
        set(gca,'TickLength',[0.01, 0.01],'TickDir','in');
        title({'Time series'},'FontSize', 12);

    figure(2)
        h = waterfall(timeSeries);
        title({'Time series'},'FontSize', 12);
        xlabel('time step','FontSize', 12);
        ylabel('node id','FontSize', 12);
        zlabel('amplitude','FontSize', 12);
        set(h, 'FaceColor', 'flat');
        color = rand(numVariables,3);
        set(h, 'FaceVertexCData', color);
        set(h, 'EdgeColor', 'k');
        set(h, 'EdgeAlpha', 0.2)
        set(h, 'FaceAlpha', 0.3)


%% Set parameters for network inference
    % all parameters are saved into a structure
    parameters = struct;
    % alpha = statistical significance level = 1 - (confidence level)
    parameters.alpha = 0.01;
    % maxLag = number of time steps in the future that are supposed to be 
    % influenced by the present 
    parameters.maxLag = 5;

%% Infer connectivity with Granger causality
    G = inferNetworkWithGranger(timeSeries, parameters);

%% Display the inferred network
    figure(3)
        tab1 = uitab('Title','Causal network');
        ax1 = axes(tab1);
        plot(ax1,G);
        set(ax1,'xtick',[],'ytick',[]);

        tab2 = uitab('Title','Weighted adjacency matric');
        ax2 = axes(tab2);
        W = getWeightedAdjacency(G);
        colormap hot;
        imagesc(ax2,W);colorbar;
        %set(ax2,'xtick',[],'ytick',[]);

    disp([newline, 'The inferred causal links are:']);
    disp(G.Edges);