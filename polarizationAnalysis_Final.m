%% Cyst Polarization Analysis (FINAL)
%
% Isaac Breinyn 2022
%
% Code for analyzing the polarization dyanmics of MDCK cysts during electric stimulation.
%
%% Clear Workspace and Load Data

% Clean the workspace and add paths
clear ; close all ; clc ;
addpath('\\latrobe\cohenlabarchive\Gawoon & Isaac\pipeline\functions') ;

% where the data is located
data_dir = '\\latrobe\CohenLabArchive\Gawoon & Isaac\SoRA' ;
cd(data_dir) ;

% create a cell array that stores project directory, cyst names, start and
% end TPs, and junction data
dataInfo = cell(3, 4) ; % rows are projects, columns are [project name, cyst names, [startTP end TP], junction data

dataInfo{1, 1} = '071222_wtMemGlowDAPI' ;
dataInfo{1, 2} = {'series9'} ;
dataInfo{1, 3} = [1 21] ;

dataInfo{2, 1} = '072222_wtMemglowDAPI' ;
dataInfo{2, 2} = {'series2' ; 'series3' ;  'series4' ; 'series4-1' ;  'series5' ; 'series6' ; 'series7' ; 'series8'} ;
dataInfo{2, 3} = [1 19;
    1 23;
    1 18;
    1 14;
    1 17;
    1 21;
    1 18;
    1 20] ;

dataInfo{3, 1} = '072922_MemglowDapiStim' ;
dataInfo{3, 2} = {'series1' ; 'series2' ;  'series3' ; 'series4' ;  'series5' ; 'series6' ; 'series7' ; 'series8' ;'series9'; 'series10' ; 'series11' ; 'series12'} ;
dataInfo{3, 3} = [1 19;
    1 14;
    1 10;
    1 10;
    1 11;
    1 10;
    1 15;
    1 13;
    1 12;
    1 13;
    1 16;
    1 15 ] ;

numCysts = cellfun('size', dataInfo(:, 2), 1) ; % the number of cysts in each project

%% Collect junctional data (manual input)

% check if junction data file exists
if isfile('data.mat')
    load('data.mat') ;
else
    
    for proj = 1:size(dataInfo, 1) % iterate over projects
        
        % storing junction info [before/after, cyst #, junction #, angular position/orientation/length
        junctions = nan(2, 15, 100, 3) ;
        
        for nn = 1:2 % start and end of stim
            for cyst = 1:numCysts(proj)
                
                % load in time series
                fn = fullfile(data_dir, dataInfo{proj, 1}, 'cropped', ['cropped_' dataInfo{proj, 2}{cyst} '_green.tif']) ;
                
                % pull out the specified TP
                image = tiffreadVolume(fn) ;
                image = squeeze(image(:, :, dataInfo{proj, 3}(cyst, nn))) ;
                
                % ----- PERFORM GRAPHICAL INPUT ----- %
                close all
                figure('units','inches','position',[1 3 6 5]); hold on; box on; set(gcf,'color','w') ;
                title({fullfile(dataInfo{proj, 1}, 'cropped', ['cropped_' dataInfo{proj, 2}{cyst} '_green.tif']),'Click on the Center of the Cyst n times, then hit "Enter"'})
                imagesc(image) ; drawnow() ;
                axis equal ;
                [x_ctrs, y_ctrs] = ginput ;
                x_ctr = mean(x_ctrs, 'all', 'omitnan') ;
                y_ctr = mean(y_ctrs, 'all', 'omitnan') ;
                
                close
                figure('units','inches','position',[1 3 6 5]); hold on; box on; set(gcf,'color','w') ;
                title({fullfile(dataInfo{proj, 1}, 'cropped', ['cropped_' dataInfo{proj, 2}{cyst} '_green.tif']),'Click in the junctions of the Cyst (basal --> apical i.e. outside --> in), then hit "Enter"'})
                imagesc(image) ;
                axis equal ;
                [x_jct, y_jct] = ginput ;
                
                % -------------------------- %
                
                for jct = 1:2:length(x_jct)
                    
                    dx = x_jct(jct) - x_jct(jct + 1) ; % x displacement
                    dy = y_jct(jct) - y_jct(jct + 1) ; % y displacement
                    
                    jct_length = sqrt(dx^2 + dy^2) ; % the length of this junction
                    jct_orientation = atan2(dy, dx) ; % the orientation of this junction
                    
                    x_mid = x_jct(jct) + dx/2 ; % the x coord of the junction midpoint
                    y_mid = y_jct(jct) + dy/2 ; % the y coord of the junction midpoint
                    
                    jct_position = [x_mid - x_ctr y_mid - y_ctr] ; % the position of the jct midpoint relative to the cyst center
                    
                    jct_angle = atan2(jct_position(2), jct_position(1)) ; % the angular position of the jct midpoint on the cyst
                    
                    % store the junction angular position, orientation, and length
                    junctions(nn, cyst, jct, 1) = jct_angle ;
                    junctions(nn, cyst, jct, 2) = jct_orientation ;
                    junctions(nn, cyst, jct, 3) = jct_length ;
                    
                end
            end
        end
        % store data
        dataInfo{proj, 4} = junctions ;
    end
    % save the data
    save('data.mat', 'dataInfo') ;
end

%% Some Complicated Plotting
% How does juncitonal orientation and length depend on junction location.

% junctions2 = cat(2, squeeze(junctions(:,1,:,:)),squeeze(junctions(:,2,:,:)),squeeze(junctions(:,3,:,:))) ;
junctions2 = [] ;
for proj = 1:size(dataInfo, 1)
    for cyst = 1:numCysts(proj)
        junctions2 = cat(2, junctions2, squeeze(dataInfo{proj, 4}(:, cyst, :, :))) ;
    end
end


angle_edges = -pi : 0.3 : (pi+0.3) ;
position_edges = -pi : 0.1 : (pi +0.1);
angle_colors = phasemap(length(position_edges)) ;

% make phase color map symmetric
angle_colors(round(length(angle_colors)/2):end, :) = flipud(angle_colors(1:round(length(angle_colors)/2)+1, :)) ;

ctrs = angle_edges(1:length(angle_edges)-1) + mean(diff(angle_edges))/2;              % Calculate Centres

close all
figure('units','inches','position',[1 5 10 4]); set(gcf,'color','w') ;
colormap(angle_colors) ;

sz = 10 ;
px_sz = 325/2000 ; % the um/px
words2 = {'before stim', 'after stim'} ;
% orientations
for nn = 1:2
    
    positions = squeeze(junctions2(nn,:,1)) ;
    orientations = squeeze(junctions2(nn,:,2)) ;
    lengths = squeeze(junctions2(nn,:,3)) ;
    
    [positions, I] = sort(positions(:)) ;
    
    orientations = orientations(I) ;
    lengths = lengths(I) ;
    
    n_good = sum(~isnan(orientations)) ;
    subplot(1, 3, nn)
    hold on; box on;
    title(['Junction Orientations ' words2{nn} ' (n=' num2str(n_good) ')']) ;
    xlim([-pi pi]) ;
    ylim([0 50]) ;
    xlabel('\textbf{Junction Orientation ($\theta$) [rad]}','interpreter','latex')
    ylabel('\textbf{Count}','interpreter','latex') ;
    axis square
    xticks([ -pi -pi/2 0 pi/2 pi ]) % specify x ticks
    xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
    
    counts = nan(length(position_edges), length(angle_edges)-1) ;
    
    for pp = 1:(length(position_edges)-1)
        aoi = orientations(positions >= position_edges(pp) & positions <= position_edges(pp+1)) ;
        counts(pp, :) = histcounts(aoi, angle_edges) ;
    end
    
    total_counts = sum(counts, 1, 'omitnan') ;
    
    ba = bar(ctrs(1:size(counts, 2)),counts,'stacked', 'FaceColor','flat','EdgeColor','none','LineWidth',0.1) ;
    for pp = 1:(length(position_edges)-1)
        for rr = 1:length(ba(pp).CData)
            ba(pp).CData(rr, :) = angle_colors(pp, :) ;
        end
    end
    
    %     phasebar
    
    plot_details = {'b', 'r'} ;
    subplot(1, 3, 3)
    axis equal
    hold on
    box on
    %     plot(positions, orientations, plot_details{nn}, 'LineWidth', 2)
    scatter(positions, orientations, sz, plot_details{nn}, 'filled') ;
    
    xlim([-3.5 3.5]) ;
    ylim([-3.5 3.5]) ;
    title('Junction Correlation') ;
    ylabel('\textbf{Junction Orientation ($\theta$) [rad]}','interpreter','latex') ;
    xlabel('\textbf{Junction Position ($\theta$) [rad]}','interpreter','latex') ;
    xticks([ -pi -pi/2 0 pi/2 pi ]) % specify x ticks
    xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'}) ;
    yticks([ -pi -pi/2 0 pi/2 pi ]) % specify x ticks
    yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
    
end

% legend('Before', 'After', 'Location', 'NorthEastOutside')
figure('units','inches','position',[1 5 10 4]); set(gcf,'color','w') ;
colormap(angle_colors) ;

length_edges = 0:2.5:40 ;
ctrs = length_edges(1:length(length_edges)-1) + mean(diff(length_edges))/2;              % Calculate Centres

% lengths
for nn = 1:2
    
    positions = squeeze(junctions2(nn,:,1)) ;
    orientations = squeeze(junctions2(nn,:,2)) ;
    lengths = squeeze(junctions2(nn,:,3)) ;
    
    [positions, I] = sort(positions(:)) ;
    
    orientations = orientations(I) ;
    lengths = lengths(I) ;
    
    n_good = sum(~isnan(lengths)) ;
    subplot(1, 3, nn)
    hold on; box on;
    axis square
    title(['Junction Lengths ' words2{nn} ' (n=' num2str(n_good) ')']) ;
    xlim([0 30]) ;
    ylim([0 130]) ;
    xlabel('\textbf{Junction Length [$\mu$m]}','interpreter','latex')
    ylabel('\textbf{Count}','interpreter','latex') ;
    
    counts = nan(length(position_edges), length(length_edges)-1) ;
    
    for pp = 1:(length(position_edges)-1)
        aoi = lengths(positions >= position_edges(pp) & positions <= position_edges(pp+1)) ;
        aoi = aoi.*px_sz ;
        counts(pp, :) = histcounts(aoi, length_edges) ;
    end
    
    total_counts = sum(counts, 1, 'omitnan') ;
    
    ba = bar(ctrs(1:size(counts, 2)),counts,'stacked', 'FaceColor','flat','EdgeColor','none','LineWidth',0.1) ;
    for pp = 1:(length(position_edges)-1)
        for rr = 1:length(ba(pp).CData)
            ba(pp).CData(rr, :) = angle_colors(pp, :) ;
        end
    end
    
    phasebar
    
    plot_details = {'b', 'r'} ;
    subplot(1, 3, 3)
    axis square
    hold on
    box on
    %     plot(positions, lengths, plot_details{nn}, 'LineWidth', 2)
    scatter(positions, px_sz*lengths, sz, plot_details{nn}, 'filled') ;
    xlim([-3.5 3.5]) ;
    ylim([0 30]) ;
    title('Junction Correlation') ;
    ylabel('\textbf{Junction Length [$\mu$m]}','interpreter','latex')
    xlabel('\textbf{Junction Position ($\theta$) [rad]}','interpreter','latex') ;
    xticks([ -pi -pi/2 0 pi/2 pi ]) % specify x ticks
    xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
end

% legend('Before', 'After', 'Location', 'NorthEastOutside')

%% Measure Plumpness

% units of junction2 are junction angular position, orientation, and length

meanJunctionLengths = zeros(2,2,3) ; % rows are before and after, column is left and right, 3rd dim is mean/std/N

for nn = 1:2
    
    junctionLengths = squeeze(junctions2(nn,:,3)) ;
    
    % Grab all junction lengths on the left of the cyst
    leftIdx1 = find(junctions2(nn,:,1) < -pi/2 & junctions2(nn,:,1)  > -pi) ;
    leftIdx2 = find(junctions2(nn,:,1) < pi & junctions2(nn,:,1)  > pi/2) ;
    
    leftJunctionsLengths = junctionLengths(vertcat(leftIdx1', leftIdx2')) ;
    
    meanJunctionLengths(nn, 1, 1) = mean(leftJunctionsLengths) ;
    meanJunctionLengths(nn, 1, 2) = std(leftJunctionsLengths) ;
    meanJunctionLengths(nn, 1, 3) = length(leftJunctionsLengths) ;
    % Grab all junction lengths on the right of the cyst
    rightIdx = find(junctions2(nn,:,1) < pi/2 & junctions2(nn,:,1)  > -pi/2) ;
    
    rightJunctionsLengths = junctionLengths(rightIdx) ;
    
    meanJunctionLengths(nn, 2, 1) = mean(rightJunctionsLengths) ;
    meanJunctionLengths(nn, 2, 2) = std(rightJunctionsLengths) ;
    meanJunctionLengths(nn, 2, 3) = length(rightJunctionsLengths) ;
end
