function run_spatiotemporal_pca_wrapper(folder1, folder2, outputDir, ...
    events, epochWindow, BaselineCorrect, erpWindow, nPermPA, rotationMethods)
% RUN_SPATIOTEMPORAL_PCA_WRAPPER
%   Wrapper for pooled spatiotemporal PCA analysis on EEG .set data.
%
% Usage:
%   run_spatiotemporal_pca_wrapper(folder1, folder2, outputDir, ...
%       events, epochWindow, BaselineCorrect, erpWindow, nPermPA)
%
% Example:
%   run_spatiotemporal_pca_wrapper('example_data/PBO','example_data/THC','results', ...
%       {'S  6','S  7'}, [-0.2 0.7], [-200 0], [300 500], 500);
%
% Inputs:
%   folder1        - path to .SET files for the first group
%   folder2        - path to the second group; '' for a single group
%   outputDir      - directory to save reports/outputs
%   events         - cell array of event codes (default: {'S  6','S  7'})
%   epochWindow    - [min max] epoch window in seconds (default: [-0.2 0.7])
%   BaselineCorrect- [min max] baseline correction window in ms (default: [-200 0])
%   erpWindow      - [min max] ERP component search window in ms (default: [300 500])
%   nPermPA        - number of permutations for parallel analysis (default: 500)
%
% Requirements:
%   - EEGLAB installed and on MATLAB path
%   - pop_loadset, pop_epoch, pop_rmbase available
%
% Author: Matthew Gunn, 2025

%% Handle defaults if not provided
if nargin < 4 || isempty(events)
    error('You must specify at least one event code in "events".');
end
if nargin < 5 || isempty(epochWindow),     epochWindow = [-0.2 0.7]; end
if nargin < 6 || isempty(BaselineCorrect), BaselineCorrect = [-200 0]; end
if nargin < 7 || isempty(erpWindow),       erpWindow = [300 500]; end
if nargin < 8 || isempty(nPermPA),         nPermPA = 500; end
if nargin < 9 || isempty(rotationMethods), rotationMethods = {{'promax','Power',4}}; end

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

validateattributes(nPermPA,{'numeric'},{'scalar','integer','positive'});
folders={folder1}; groupNames={'Group1'};
if ~isempty(folder2), folders{2}=folder2; groupNames{2}='Group2'; end
for k=1:numel(folders)
    if ~isfolder(folders{k}), error('PCA:InputFolder','Input folder does not exist: %s',folders{k}); end
end
if isempty(which('pop_loadset')), error('PCA:EEGLAB','Initialize EEGLAB before running this workflow.'); end
% Use MATLAB Statistics PCA even if an EEGLAB plugin shadows it, and restore
% the caller's path after the analysis.
oldPath=path; pathCleanup=onCleanup(@() path(oldPath)); %#ok<NASGU>
statsPath=fullfile(matlabroot,'toolbox','stats','stats');
if ~isfile(fullfile(statsPath,'pca.m')), error('PCA:Toolbox','Statistics and Machine Learning Toolbox is required.'); end
addpath(statsPath,'-begin');
initialFigures=findall(groot,'Type','figure');

%% =============================
% Run pooled analysis (all conditions together)
% =============================
[missing_files, mismatch_files, labels] = run_spatiotemporal_pca_core( ...
    folders, groupNames, events, epochWindow, BaselineCorrect, erpWindow, nPermPA, rotationMethods, outputDir);



save_list(fullfile(outputDir,'missing_epochs_all.txt'), missing_files);
save_list(fullfile(outputDir,'mismatch_length_all.txt'), mismatch_files);
save(fullfile(outputDir,'analysis_metadata.mat'),'folders','labels','events','epochWindow','BaselineCorrect','erpWindow','nPermPA','rotationMethods');
newFigures=setdiff(findall(groot,'Type','figure'),initialFigures);
for k=1:numel(newFigures)
    savefig(newFigures(k),fullfile(outputDir,sprintf('figure_%02d.fig',k)));
    exportgraphics(newFigures(k),fullfile(outputDir,sprintf('figure_%02d.png',k)));
end

end
% =============================
% Helper: Get files in folder
% =============================
    function flist = files_in_folder(folderpath)
        flist = dir(fullfile(folderpath, '*.set'));
    end

% =============================
% Helper: Save list to .txt
% =============================
    function save_list(filename, list)
        fid = fopen(filename, 'w');
        if isempty(list)
            fprintf(fid, 'None\n');
        else
            fprintf(fid, '%s\n', list{:});
        end
        fclose(fid);
    end

%% =============================
% Main PCA function (pooled)
% =============================
    function [missing_files, mismatch_files, labels] = run_spatiotemporal_pca_core( ...
            folders, groupNames, events, epochWindow, BaselineCorrect, erpWindow, nPermPA, rotationMethods, outputDir)


        all_data = [];
        labels = {};  % condition labels for each subject × condition ERP
        missing_files = {};
        mismatch_files = {};
        times_ref = [];

        % ----- Load & preprocess all folders -----
        for g = 1:length(folders)
            flist = files_in_folder(folders{g});
            for i = 1:length(flist)% This is so only half the sample is grabbed.
                EEG = pop_loadset('filename', flist(i).name, 'filepath', flist(i).folder);

                % Skip if no matching events
                eventTypes = {EEG.event.type};
                if ~any(ismember(eventTypes, events))
                    warning('No matching events in %s', flist(i).name);
                    missing_files{end+1} = flist(i).name;
                    continue;
                end

                % Epoching
                EEG = pop_epoch(EEG, events, [epochWindow(1) epochWindow(2)], 'epochinfo', 'yes');
                if EEG.trials == 0
                    warning('No epochs for %s', flist(i).name);
                    missing_files{end+1} = flist(i).name;
                    continue;
                end

                % Baseline correction
                EEG = pop_rmbase(EEG, BaselineCorrect);

                % Save time reference
                if isempty(times_ref), times_ref = EEG.times; chanlocs_ref = EEG.chanlocs; end

                % ERP for subject
                erp = mean(EEG.data, 3); % channels × timepoints

                % Length match check
                if isempty(all_data)
                    all_data(:,:,1) = erp;
                else
                    if isequal(size(erp),size(all_data(:,:,1))) && ...
                            isequal(EEG.times,times_ref) && ...
                            isequal({EEG.chanlocs.labels},{chanlocs_ref.labels})
                        all_data(:,:,end+1) = erp;
                    else
                        warning('Skipping %s (time mismatch)', flist(i).name);
                        mismatch_files{end+1} = flist(i).name;
                        continue;
                    end
                end

                % Save condition label
                labels{end+1} = groupNames{g};
            end
        end

        if isempty(all_data)
            error('PCA:NoData','No usable datasets for pooled PCA; check files and event codes.');
        end
        % PCA path conflicts are handled by the wrapper using matlabroot;
        % its cleanup restores the caller's path even after errors.

        % =============================
        % Step 1: Temporal PCA Setup
        % =============================
        %rotationMethods = {...
        % {'varimax'}, ...
        % {'quartimax'}, ...
        % {'equamax'}, ...
        % {'parsimax'}, ...
        %{'promax','Power',4} ...
        % };
        if nargin < 8 || isempty(rotationMethods)
            rotationMethods = {{'promax','Power',4}}; % default
        end

        [nChan, nTime, nSubj] = size(all_data);
        if nSubj<3, error('PCA:Subjects','At least three usable datasets are required.'); end
        reshaped_data = reshape(permute(all_data, [3,1,2]), nSubj*nChan, nTime);




        % =============================
        % Step 1a: Decide factor numbers via Parallel Analysis
        % =============================


        % Temporal Parallel Analysis
        fprintf('Running Parallel Analysis for Temporal PCA...\n');
        nTemporalFactors = parallel_analysis(reshaped_data, nPermPA);

        fprintf('Suggested number of Temporal Factors = %d\n', nTemporalFactors);

        % Spatial Parallel Analysis
        fprintf('Running Parallel Analysis for Spatial PCA...\n');
        % Here we just test on subject × channel matrix for each temporal factor
        % (will be re-run inside run_spatiotemporal_pca, but this gives global guidance)
        % Spatial Parallel Analysis (use subject × channel ERP matrix)
        fprintf('Running Parallel Analysis for Spatial PCA...\n');
        erpMatrix = reshape(all_data, nChan, nTime * nSubj)';  % subjects*trials × channels
        nSpatialFactors = parallel_analysis(erpMatrix, nPermPA);
        nSpatialFactors=min([nSpatialFactors,nSubj-1,nChan]);
        if nTemporalFactors<1 || nSpatialFactors<1
            error('PCA:NoFactors','Parallel analysis retained no factors. No decomposition was forced.');
        end

        fprintf('Suggested number of Spatial Factors = %d\n', nSpatialFactors);
        % === Reporting block ===
        fprintf('\n--- Parallel Analysis Results ---\n');
        fprintf('Temporal factors chosen = %d\n', nTemporalFactors);
        fprintf('Spatial factors chosen  = %d\n', nSpatialFactors);
        fprintf('These numbers are retained because real eigenvalues > random eigenvalues in parallel analysis.\n');
        fprintf('--------------------------------\n\n');

        % Save decision for reproducibility
        fid = fopen(fullfile(outputDir,'parallel_analysis_report.txt'),'w');
        fprintf(fid,'Temporal factors chosen = %d\n', nTemporalFactors);
        fprintf(fid,'Spatial factors chosen  = %d\n', nSpatialFactors);
        fprintf(fid,'(Criteria: real eigenvalue > mean random eigenvalue, %d permutations)\n', nPermPA);
        fclose(fid);


        % =============================
        % Step 2: Pre-allocate scree structures & run unrotated PCA once
        % =============================
        scree.RotationNames = cell(1, length(rotationMethods));
        scree.VarExplained = zeros(length(rotationMethods), nTemporalFactors);
        scree.VarCumulative = zeros(length(rotationMethods), nTemporalFactors);
        screeSpatial.VarExplained = cell(length(rotationMethods), nTemporalFactors);
        screeSpatial.VarCumulative = cell(length(rotationMethods), nTemporalFactors);



        [coeff_time_unrot, score_time_unrot, latent_time] = pca(reshaped_data, ...
            'Centered', true, 'NumComponents', nTemporalFactors);

        unrot_var = 100 * latent_time(:) / sum(latent_time);
        unrot_cum = cumsum(unrot_var);

        % =============================
        % Step 3: Loop through each rotation method
        % =============================
        %for r = 1:length(rotationMethods)
        %    rotArgs = rotationMethods{r};
        for r = 1:length(rotationMethods)
            rotArgs = rotationMethods{r};
            rotName = upper(rotArgs{1});
            scree.RotationNames{r} = rotName;
            screeSpatial.RotationNames{r} = rotName;

            % Apply rotation
            if size(coeff_time_unrot,2)==1
                coeff_time_rot=coeff_time_unrot; rotation=1;
            else
                [coeff_time_rot,rotation] = rotatefactors(coeff_time_unrot, 'Method', rotArgs{:});
            end
            coeff_time_rot_scaled = (coeff_time_unrot .* sqrt(latent_time(1:size(coeff_time_rot,2)))') * rotation;

            % Rotated scores & variance explained
            rotated_scores = score_time_unrot / rotation';
            var_explained_rot = 100 * var(rotated_scores) / sum(var(rotated_scores));
            scree.VarExplained(r,:) = var_explained_rot(1:nTemporalFactors);
            scree.VarCumulative(r,:) = cumsum(var_explained_rot(1:nTemporalFactors));


            % -----------------------------------------
            % Temporal PCA plot: Raw (top) vs covariance-scaled (bottom)
            % -----------------------------------------
            figure('Name', sprintf('Temporal PCA (%s) Raw (top) vs covariance-scaled (bottom)', ...
                rotName));
            colors = lines(size(coeff_time_rot,2));

            % Raw loadings
            subplot(2,1,1);
            hold on;
            for f = 1:size(coeff_time_rot,2)
                plot(times_ref, coeff_time_rot(:,f), '-', 'Color', colors(f,:), 'LineWidth', 1.5);
            end

            condition_label = '.';

            xlabel('Time (ms)'); ylabel('Raw Loading');
            title(sprintf('%s Temporal PCA (%s) - Raw\nnTemporalFactors = %d (Parallel Analysis)', ...
                condition_label, rotName, nTemporalFactors));
            legend(arrayfun(@(x) sprintf('F%d',x), 1:size(coeff_time_rot,2), 'UniformOutput', false), ...
                'Location', 'best');
            grid on; hold off;

            % covariance-scaled loadings
            subplot(2,1,2);
            hold on;
            for f = 1:size(coeff_time_rot_scaled,2)
                plot(times_ref, coeff_time_rot_scaled(:,f), '-', 'Color', colors(f,:), 'LineWidth', 1.5);
            end
            xlabel('Time (ms)'); ylabel('Scaled Loading');
            title(sprintf('%s Temporal PCA (%s) - covariance-scaled\nnTemporalFactors = %d (Parallel Analysis)', ...
                condition_label, rotName, nTemporalFactors));
            legend(arrayfun(@(x) sprintf('F%d',x), 1:size(coeff_time_rot_scaled,2), 'UniformOutput', false), ...
                'Location', 'best');
            grid on; hold off;

            % =============================
            % Step 4: Spatial PCA & Detailed Per-TF Plots
            % =============================
            erpShadeColor = [0.85 0.85 1]; % ERP window shade
            erpShadeAlpha = 0.35;

            for f = 1:size(coeff_time_rot,2)
                var_explained_temporal = var_explained_rot(f);

                [~, peak_idx] = max(abs(coeff_time_rot_scaled(:,f)));
                peak_latency = times_ref(peak_idx);

                % Spatial PCA
                factor_scores = rotated_scores(:,f);
                subj_factor_scores = reshape(factor_scores, nSubj, nChan);
                [coeff_space, score_space, latent_space] = pca(subj_factor_scores, ...
                    'Centered', true, 'NumComponents', nSpatialFactors);

                if size(coeff_space,2)==1
                    coeff_space_rot=coeff_space; spatialRotation=1;
                else
                    [coeff_space_rot,spatialRotation] = rotatefactors(coeff_space, 'Method', rotArgs{:});
                end
                coeff_space_rot_scaled = (coeff_space .* sqrt(latent_space(1:size(coeff_space_rot,2)))') * spatialRotation;
                save(fullfile(outputDir,sprintf('pca_%s_TF%02d.mat',rotName,f)), ...
                    'coeff_time_unrot','score_time_unrot','coeff_time_rot','coeff_time_rot_scaled', ...
                    'rotated_scores','rotation','coeff_space_rot','coeff_space_rot_scaled', ...
                    'latent_time','latent_space','times_ref','chanlocs_ref','labels');

                variance_explained_space = 100 * latent_space / sum(latent_space);
                cumulative_variance_space = cumsum(variance_explained_space);

                % Store spatial scree data
                screeSpatial.VarExplained{r,f} = variance_explained_space;
                screeSpatial.VarCumulative{r,f} = cumulative_variance_space;

                %if showChanLabels
                elecOpt = 'labels';
                % else
                %     elecOpt = 'on';
                % end

                raw_min = min(coeff_space_rot(:));
                raw_max = max(coeff_space_rot(:));
                scaled_min = min(coeff_space_rot_scaled(:));
                scaled_max = max(coeff_space_rot_scaled(:));

                figure('Name', sprintf('%s TF%d (%s) Detailed Overview', ...
                    condition_label, f, rotName), 'Position', [100, 100, 1400, 800]);

                % Temporal Raw
                subplot(3,3,1);
                hold on;
                patch([erpWindow(1) erpWindow(2) erpWindow(2) erpWindow(1)], ...
                    [min(coeff_time_rot(:,f)) min(coeff_time_rot(:,f)) ...
                    max(coeff_time_rot(:,f)) max(coeff_time_rot(:,f))], ...
                    erpShadeColor, 'FaceAlpha', erpShadeAlpha, 'EdgeColor', 'none');
                plot(times_ref, coeff_time_rot(:,f), 'b-', 'LineWidth', 1.5);
                plot([peak_latency peak_latency], ylim, 'r--', 'LineWidth', 1);
                hold off;
                xlabel('Time (ms)'); ylabel('Raw Loading');
                title(sprintf('TF%d Raw (%s)\n%.1f%% var (Peak: %d ms)\nFactors chosen: %d temporal, %d spatial', ...
                    f, rotName, var_explained_temporal, round(peak_latency), nTemporalFactors, nSpatialFactors));

                grid on;

                % Temporal Scaled
                subplot(3,3,2);
                hold on;
                patch([erpWindow(1) erpWindow(2) erpWindow(2) erpWindow(1)], ...
                    [min(coeff_time_rot_scaled(:,f)) min(coeff_time_rot_scaled(:,f)) ...
                    max(coeff_time_rot_scaled(:,f)) max(coeff_time_rot_scaled(:,f))], ...
                    erpShadeColor, 'FaceAlpha', erpShadeAlpha, 'EdgeColor', 'none');
                plot(times_ref, coeff_time_rot_scaled(:,f), 'g-', 'LineWidth', 1.5);
                plot([peak_latency peak_latency], ylim, 'r--', 'LineWidth', 1);
                hold off;
                xlabel('Time (ms)'); ylabel('Scaled Loading');
                title(sprintf('TF%d Scaled (%s)\n%.1f%% var (Peak: %d ms)\nFactors chosen: %d temporal, %d spatial', ...
                    f, rotName, var_explained_temporal, round(peak_latency), nTemporalFactors, nSpatialFactors));

                grid on;

                % Spatial Topoplots (Raw)
                nCompsToPlot = size(coeff_space_rot,2);
                for sf = 1:nCompsToPlot
                    subplot(3, nCompsToPlot, nCompsToPlot + sf);
                    topoplot(coeff_space_rot(:,sf), chanlocs_ref, ...
                        'maplimits', [raw_min raw_max], 'electrodes', elecOpt, ...
                        'style', 'both', 'numcontour', 6);
                    title(sprintf('Raw C%d (%s)\n%.1f%%', ...
                        sf, rotName, variance_explained_space(sf)));
                end

                % Spatial Topoplots (Scaled)
                for sf = 1:nCompsToPlot
                    subplot(3, nCompsToPlot, 2*nCompsToPlot + sf);
                    topoplot(coeff_space_rot_scaled(:,sf), chanlocs_ref, ...
                        'maplimits', [scaled_min scaled_max], 'electrodes', elecOpt, ...
                        'style', 'both', 'numcontour', 6);
                    title(sprintf('Scaled C%d (%s)\n%.1f%%', ...
                        sf, rotName, variance_explained_space(sf)));
                end

                % Bar chart of spatial variance
                subplot(3,3,3);
                bar(variance_explained_space, 'FaceColor', [0.2 0.6 0.8]);
                hold on; plot(cumulative_variance_space, 'r-o', 'LineWidth', 1.2);
                hold off;
                xlabel('Spatial Component'); ylabel('% Variance Explained');
                legend({'Individual', 'Cumulative'}, 'Location', 'best');
                title(sprintf('Spatial PCA Variance (%s)', rotName));
                grid on;
            end
        end

        % =============================
        % Step 5: Combined Temporal + Spatial Scree Summary
        % =============================
        maxCompsToShow = min(10, nTemporalFactors);
        colorsRotations = lines(length(rotationMethods));

        figure('Name', sprintf('%s - Temporal & Spatial Scree Summary', condition_label), ...
            'Position', [100 100 1200 800]);

        % Temporal variance per component
        subplot(2,1,1);
        hold on;
        for r = 1:length(rotationMethods)
            plot(1:maxCompsToShow, scree.VarExplained(r,1:maxCompsToShow), '-o', ...
                'Color', colorsRotations(r,:), 'LineWidth', 1.5, ...
                'MarkerFaceColor', colorsRotations(r,:));
        end
        plot(1:maxCompsToShow, unrot_var(1:maxCompsToShow), 'k--', 'LineWidth', 2);
        yline(80, 'k:', '80% Threshold');
        xlabel('Component Number'); ylabel('% Variance Explained');
        title(sprintf('%s Temporal PCA Scree\nnTemporalFactors = %d (Parallel Analysis)', ...
            condition_label, nTemporalFactors));

        legend([scree.RotationNames, {'Unrotated'}], 'Location', 'best');
        grid on;

        % Spatial variance averaged across TFs
        maxCompsToShow=min([maxCompsToShow,nChan,nSubj-1]);
        subplot(2,1,2);
        hold on;
        for r = 1:length(rotationMethods)
            avgSpatialVar = zeros(1, maxCompsToShow);
            for f = 1:nTemporalFactors
                avgSpatialVar = avgSpatialVar + screeSpatial.VarExplained{r,f}(1:maxCompsToShow)';
            end
            avgSpatialVar = avgSpatialVar / nTemporalFactors;
            plot(1:maxCompsToShow, avgSpatialVar, '-o', ...
                'Color', colorsRotations(r,:), 'LineWidth', 1.5, ...
                'MarkerFaceColor', colorsRotations(r,:));
        end
        yline(80, 'k:', '80% Threshold');
        xlabel('Spatial Component Number'); ylabel('Avg % Variance Explained');
        title(sprintf('%s Spatial PCA Scree (Mean Across TFs)\nnSpatialFactors = %d (Parallel Analysis)', ...
            condition_label, nSpatialFactors));

        legend(screeSpatial.RotationNames, 'Location', 'best');
        grid on;

    end


    function nFactors = parallel_analysis(dataMatrix, nPerm)
        % dataMatrix: observations × variables (e.g., subjects*channels × time)
        % nPerm: number of random permutations
        % Output: number of components to retain

        % --- Step 1. Compute real eigenvalues ---
        [~,~,latent] = pca(dataMatrix,'Centered',true);
        realEig = latent(:);

        % --- Step 2. Random eigenvalues ---
        randEigMat = zeros(nPerm, length(realEig));
        for i = 1:nPerm
            randData = zeros(size(dataMatrix));
            for col=1:size(dataMatrix,2)
                randData(:,col)=dataMatrix(randperm(size(dataMatrix,1)),col);
            end
            [~,~,latent_rand] = pca(randData,'Centered',true);
            randEigMat(i,:) = latent_rand(:);
        end

        randMeanEig = mean(randEigMat,1);

        % --- Step 3. Decide retained factors ---
        nFactors = sum(realEig > randMeanEig(:));

        % --- Step 4. Plot real vs. random eigenvalues ---
        figure('Name','Parallel Analysis Scree','Position',[100 100 800 600]);
        hold on;
        plot(1:length(realEig), realEig, 'b-o','LineWidth',2,'DisplayName','Real Data');
        plot(1:length(randMeanEig), randMeanEig, 'r--s','LineWidth',2,'DisplayName','Random Data (mean)');
        xline(nFactors,'k:','LineWidth',2,'DisplayName',sprintf('Cutoff at %d',nFactors));
        hold off;
        xlabel('Component Number'); ylabel('Eigenvalue');
        title(sprintf('Parallel Analysis: Retained %d Factors', nFactors));
        legend('Location','northeast');
        grid on;
    end

