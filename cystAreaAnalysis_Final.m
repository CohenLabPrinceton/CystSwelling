%% Cyst Area Analysis (Final)
%
% Isaac Breinyn & Gawoon Shim 2022
%
% Code for analyzing the swelling dynamics of MDCK cysts during and after
% electric stimulation

%% Clear Workspace and Load Data

% --- Clean the workspace and add paths --- %
clear ; close all ; clc ;
addpath('\\latrobe\cohenlabarchive\Gawoon & Isaac\pipeline\functions') ;

% --- specify where the data is located --- %
%data_dir = 'Z:\Gawoon & Isaac\cyst images' ;
data_dir = '\\latrobe\CohenLabArchive\Gawoon & Isaac\cyst images\';
cd(data_dir) ;

%% Calculate Cyst Geomtric Properties Across Experiments

expts = dir(data_dir) ; % get a list of all the experiments in this folder

% --- Query words for loading in data --- %
query1 = 'wtCysts' ; query2 = 'Y27632'; query3 = '0mA' ;
% ------------------------------------------ %

expts = expts(contains({expts.name},query1) & contains({expts.name},query2) &contains({expts.name},query3) &~contains({expts.name},{'.','..', 'pipeline'}) & ~contains({expts.name},'.mat'));

sizes = nan(length(expts), 75, 200) ; % sizes of the cysts i.e. sum of pixels [experiment, cyst, TP]
for expt = 1:length(expts) % iterate over each experiment
    
    curr_expt = expts(expt).name ; % current experiment name
    cysts = dir(fullfile(data_dir, curr_expt, 'binarized_cysts')) ; % get a list of all the cysts in this folder
    cysts = cysts(~ismember({cysts.name},{'.','..', 'Thumbs.db'}));
    
    for cyst = 1:length(cysts) % iterate over each cyst in this experiment
        
        cyst_name = cysts(cyst).name;  % name of this particular cyst
        info  = fullfile(data_dir, curr_expt, 'binarized_cysts', cyst_name); % full path
        stack = tiffreadVolume(info) ; % load in the binarized stack
        
        for tt = 1:size(stack, 3) % iterate over slices to filter out non-focal cysts and get ellipsicity
            
            slice = stack(:,:,tt) ; % grab this specific slice of the stack
            slice = bwareafilt(logical(slice), 1, 'largest') ; % only keep the largest object in the slice
            
            % Calculate and store geometric properties of the cyst
            rp = regionprops('table',slice, 'MajorAxisLength','MinorAxisLength', 'Area') ; 
            if ~isempty(rp.Area)
                sizes(expt, cyst, tt) = rp.Area ;
                
                a = rp.MajorAxisLength/2 ;
                b = rp.MinorAxisLength/2 ;
                
            end
        end
        disp(['Calculated cyst ' num2str(cyst) '/' num2str(length(cysts)) ' in experiment ' num2str(expt) '/' num2str(length(expts))]) ;
    end
end
sizes(sizes == 0) = NaN ;
save(fullfile(data_dir, [query1 '_' query2 '_' query3 '.mat']), 'sizes') ; 

%% Filter Cyst Geometric Dynamics

frame_rate = 10 ; % min/frame
stim_start = 60 ; % the tp at which stimulation starts
stim_end = 300 ; % the tp at which stimulation ends

norm_sizes = sizes ./ (sizes(:, :, 1)); % normalize all cysts by their initial size
filt_sizes = nan(size(norm_sizes,1)*size(norm_sizes, 2), size(norm_sizes, 3)) ; % filtered cysts (no longer organized by expt)

sizeThreshold = 1.05; % 0.8 for no stim, 1.05 for with stim
gradThreshold = -0.5; % -0.2 for no stim, -0.5 for with stim

goodCount = 0; % how many cysts weren't filtered

geometric_stats = zeros(size(norm_sizes, 1), size(norm_sizes, 2), 4) ; % for each cyst: is it good (0/1)?, initial area, maximal area, time of maximal area

qfilt_norm_sizes = norm_sizes ; % the array where we filter out all cysts that aren't in the top quartile for swelling

for expt = 1:size(norm_sizes, 1) % iterate over experiments
    for cyst = 1:size(norm_sizes, 2) % iterate over cysts
        
        curr_size = squeeze(norm_sizes(expt, cyst, :));
        grad = gradient(curr_size);
        
        % establish a set of thresholding rules
        rule1 = ~isnan(curr_size(1)) ; % check if this cyst actually exists/has an area
        rule2 = isempty(find(grad < gradThreshold, 1)) ; % ensure that the gradient never drops below the threshold
        rule3 = curr_size(stim_end/frame_rate) > sizeThreshold ; % ensure that the cyst actually swells -- that its size at the end of stim is higher than the control case
        
        if rule1 && rule2 && rule3 % check that this cyst satisfies the thresholding rules
            filt_sizes(goodCount+1, :) = squeeze(norm_sizes(expt, cyst, :)) ;
            goodCount = goodCount + 1;
            geometric_stats(expt, cyst, 1) = 1 ;
        else
            geometric_stats(expt, cyst, 1) = 0 ;
        end
        
        geometric_stats(expt, cyst, 2) = sizes(expt, cyst, 1) ;
        [val, idx] = max(sizes(expt, cyst, :), [],'omitnan') ;
        geometric_stats(expt, cyst, 3) = val ;
        
        if idx == 1
            geometric_stats(expt, cyst, 4) = NaN ;
        else
            geometric_stats(expt, cyst, 4) = idx ;
        end
        
    end
end

goodCount = 0 ;
qfilt_sizes = nan(size(filt_sizes)) ;

for expt = 1:size(norm_sizes, 1) % iterate over experiments
    % filter within each experiment to only include the top quartile
    qts = quantile(sort(geometric_stats(expt,:,3)./geometric_stats(expt,:,2)), 3) ;
    for cyst = 1:size(norm_sizes, 2) % iterate over cysts
        if ((geometric_stats(expt,cyst,3)/geometric_stats(expt,cyst,2) > qts(3)) && (geometric_stats(expt,cyst,1)))
            qfilt_sizes(goodCount+1, :) = squeeze(qfilt_norm_sizes(expt,cyst,:)) ;
            goodCount = goodCount + 1;
        end
    end
end

%% Plot all unfiltered cysts

temp = reshape(norm_sizes, [size(norm_sizes, 1)*size(norm_sizes, 2), size(norm_sizes,3)]) ; 

time = (10 : frame_rate : 10*size(temp, 2)) / 60;

close all % close all figures
figure('units', 'inches', 'position', [1 5 5 4]) ; hold on ; box on ; set(gcf, 'color', 'w') ; % create a figure
ylabel('Normalized Area [a.u.]', 'interpreter', 'latex') ; % y axis label
xlabel('Time [h]', 'interpreter', 'latex') ; % x axis label
ylim([0 4]) ; % set y axis limits
xlim([0 8]) ; % set x axis limits
axis square % or axis equal % set axis proportions (depending on your data)

count = 0 ;
for ii = 1:size(temp, 1)
    plot(time, temp(ii, :), '-', 'LineWidth', 1) ; % plot your data
    if length(unique(isnan(temp(ii,:)))) > 1
        count = count + 1 ;
    end
end

xline(stim_start/60,'--k', 'Linewidth', 2, 'label', {'Field = ON'})
xline(stim_end/60,'--k', 'Linewidth', 2,  'label', {'Field = OFF'})
yline(1, '-.k') ;

title(['Unfiltered Cyst Swelling Dynamics (n = ' num2str(count) ')'], 'interpreter', 'latex') ; % title the plot

%% Plot all filtered cysts

temp =filt_sizes; 

time = (10 : frame_rate : 10*size(temp, 2)) / 60;

% close all % close all figures
figure('units', 'inches', 'position', [1 5 5 4]) ; hold on ; box on ; set(gcf, 'color', 'w') ; % create a figure
ylabel('Normalized Area [a.u.]', 'interpreter', 'latex') ; % y axis label
xlabel('Time [h]', 'interpreter', 'latex') ; % x axis label
ylim([0 5]) ; % set y axis limits
xlim([0 8]) ; % set x axis limits
axis square % or axis equal % set axis proportions (depending on your data)

count = 0 ;
for ii = 1:size(temp, 1)
    plot(time, temp(ii, :), '-', 'LineWidth', 1) ; % plot your data
    if length(unique(isnan(temp(ii,:)))) > 1
        count = count + 1 ;
    end
end

xline(stim_start/60,'--k', 'Linewidth', 2, 'label', {'Field = ON'})
xline(stim_end/60,'--k', 'Linewidth', 2,  'label', {'Field = OFF'})
yline(1, '-.k') ;

title(['Filtered Cyst Swelling Dynamics (n = ' num2str(count) ')'], 'interpreter', 'latex') ; % title the plot

%% Plot top quartile filtered cysts

temp = qfilt_sizes; 

time = (10 : frame_rate : 10*size(temp, 2)) / 60;

% close all % close all figures
figure('units', 'inches', 'position', [1 5 5 4]) ; hold on ; box on ; set(gcf, 'color', 'w') ; % create a figure
ylabel('Normalized Area [a.u.]', 'interpreter', 'latex') ; % y axis label
xlabel('Time [h]', 'interpreter', 'latex') ; % x axis label
ylim([0 4]) ; % set y axis limits
xlim([0 8]) ; % set x axis limits
axis square % or axis equal % set axis proportions (depending on your data)

count = 0 ;
for ii = 1:size(temp, 1)
    plot(time, temp(ii, :), '-', 'LineWidth', 1) ; % plot your data
    if length(unique(isnan(temp(ii,:)))) > 1
        count = count + 1 ;
    end
end

xline(stim_start/60,'--k', 'Linewidth', 2, 'label', {'Field = ON'})
xline(stim_end/60,'--k', 'Linewidth', 2,  'label', {'Field = OFF'})
yline(1, '-.k') ;

title(['Top Quartile Filtered Cyst Swelling Dynamics (n = ' num2str(count) ')'], 'interpreter', 'latex') ; % title the plot

%% Plot mean + std

temp = qfilt_sizes; 

time = (10 : frame_rate : 10*size(temp, 2)) / 60;

% close all % close all figures
figure('units', 'inches', 'position', [1 5 5 4]) ; hold on ; box on ; set(gcf, 'color', 'w') ; % create a figure
ylabel('Normalized Area [a.u.]', 'interpreter', 'latex') ; % y axis label
xlabel('Time [h]', 'interpreter', 'latex') ; % x axis label
ylim([0 4]) ; % set y axis limits
xlim([0 8]) ; % set x axis limits
axis square % or axis equal % set axis proportions (depending on your data)
mean_swelling = squeeze(mean(qfilt_sizes, 1, 'omitnan')) ;
std_swelling = squeeze(std(qfilt_sizes, [], 1, 'omitnan')) ;
shadedErrorBar(time, mean_swelling, std_swelling) ;

% count = 0 ;
% for ii = 1:size(temp, 1)
%     plot(time, temp(ii, :), '-', 'LineWidth', 1) ; % plot your data
%     if length(unique(isnan(temp(ii,:)))) > 1
%         count = count + 1 ;
%     end
% end

xline(stim_start/60,'--k', 'Linewidth', 2, 'label', {'Field = ON'})
xline(stim_end/60,'--k', 'Linewidth', 2,  'label', {'Field = OFF'})
yline(1, '-.k') ;

title(['Mean + STD Cyst Swelling Dynamics (n = ' num2str(count) ')'], 'interpreter', 'latex') ; % title the plot

%% Plot Cyst Swelling Dyanmics

max_swelling = max(qfilt_sizes, [], 1) ;
max_swelling(max_swelling == 0) = NaN ;
max_swelling(max_swelling == 1) = NaN ;
max_swelling(isnan(max_swelling)) = [] ;

mean_swelling = squeeze(mean(qfilt_sizes, 1, 'omitnan')) ;
std_swelling = squeeze(std(qfilt_sizes, [], 1, 'omitnan')) ;
sem_swelling = std_swelling/sqrt(goodCount) ;

close all
figure('units','inches','position',[1 5 5 4]); hold on; box on; set(gcf,'color','w') ;
title('Cyst Swelling') ;
ylim([0.8 2.5]) ;
xlim([0 10*length(mean_swelling) / 60]) ;
% xticks(0:2:24) ;
ylabel('Normalized Area [a.u.]') ;
xlabel('Time [hr]') ;
axis square

time = (10 : frame_rate : 10*length(mean_swelling)) / 60;

color = 'm';
plot(time, mean_swelling, color, 'LineWidth', 2);
% shadedErrorBar(time, mean_swelling, std_swelling./sqrt(goodCount), 'lineprops', color) ;
% std_swelling./sqrt(goodCount) if you want SEM

errorbar(time, mean_swelling, std_swelling) ;
xline(stim_start/60,'--k', 'Linewidth', 2, 'label', {'Field = ON'})
xline(stim_end/60,'--k', 'Linewidth', 2,  'label', {'Field = OFF'})
yline(1, '-.k') ;

text(3 ,mean_swelling(34), ['n = ' num2str(goodCount)], 'Color', color) ;
