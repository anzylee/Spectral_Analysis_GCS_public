clear all; close all; clc; warning('off') 

path_fig = ".\figures_ZnWn";

exec_mscohere = 1; exec_sig_freq = exec_mscohere;
exec_boxplot = 0;
exec_violinplot = 0;
Fbf_factor = 999; % 999 if no threshold for frequency 
co_thresholds_avg = [];

ds = tabularTextDatastore('.\4_14_21 Full data set.csv');
[full_data_set, info] = read(ds);

ds_inf = spreadsheetDatastore('.\inflection_points.xlsx');
[inflection_data, info_inf] = read(ds_inf);

ind_W = inflection_data.variable == "W";
ind_Z = inflection_data.variable == "Z";


for channel_type = 1:5   % 1:5
    channel_dir = ".\DATASET_by_reach\SC0"...
        +num2str(channel_type);
    site_names_dir = dir(channel_dir);
    site_names = strings(1, length(site_names_dir)-3+1);

    for ii = 3:length(site_names_dir)               % skipping '.' and '..' folders
        site_name = site_names_dir(ii).name;
        site_names(ii-2) = site_name;
    end

    cum_locs_1 = []; cum_ang_pks_1 = [];
    cum_locs_2 = []; cum_ang_pks_2 = [];
    cum_locs_3 = []; cum_ang_pks_3 = [];

    for ii = 1:length(site_names) %1:length(site_names)
        site_name = site_names(ii);
        char_site_name = char(site_name);
        site_num = char_site_name(6:end);
        table_dir = channel_dir + '\' + site_name + "\Tables";
        table_all_dir = dir(table_dir);
        table_all = strings(1, length(table_all_dir)-3+1);

        for jj = 3:length(table_all_dir)            % skipping '.' and '..' folders
            table_name = table_all_dir(jj).name;
            table_all(jj-2) = table_name;
        end

        stage_num = 1;
        water_stages = [];
        %% 

        for table = table_all
            if length(split(table, 'ft')) == 2
                water_stage = split(table, 'ft');
                water_stages = [water_stages, str2num(replace(water_stage(1), 'p', '.'))];
            end
        end
        
        [water_stages_sort, water_stages_sort_ind] = sort(water_stages);
        table_all(1:length(water_stages)) = table_all(water_stages_sort_ind);

        for table = table_all % table_all
            char_table = char(table);

            if length(char_table) > 21

                if char_table(end-20:end) == 'WD_analysis_table.csv'
                    stage_num

                    table_split = table.split("_");
                    water_stage_p = table_split(1);
                    water_stage_pt = water_stage_p.replace("p", ".");

                    site_stage = site_name + "_" +num2str(stage_num);

                    table_file = table_dir + '\' + table;
                    datastore = tabularTextDatastore(table_file);
                    geoseries = read(datastore);    % in US cumstomary units
                    
                    %% normalization
                    index_comid = find(full_data_set.comid == str2num(site_num));
                    bf_width = full_data_set.BF_width_ft(index_comid);
                    geoseries.W_n = geoseries.W/bf_width;
                    geoseries.Z_n = geoseries.Z/bf_width;

                    W = geoseries.W_n;
                    Z = geoseries.Z_n;

                    Fs = 1 / geoseries.dist_down(2);

                    %Cxy_rand = zeros(1000, length(W));

                    n = 1000; Cxy_rand = [];
                    for imc = 1:n
                        W_rand_tmp = rand(size(W));
                        W_rand_tmp_mse = sum((W_rand_tmp-mean(W_rand_tmp)).^2);
                        W_mse = sum((W-mean(W)).^2);
                        W_rand = W_rand_tmp*sqrt((W_mse)/(W_rand_tmp_mse));

                        Z_rand_tmp = rand(size(Z));
                        Z_rand_tmp_mse = sum((Z_rand_tmp-mean(Z_rand_tmp)).^2);
                        Z_mse = sum((Z-mean(Z)).^2);
                        Z_rand = Z_rand_tmp*sqrt((Z_mse)/(Z_rand_tmp_mse));

                        [Cxy,F] = mscohere(Z_rand,W_rand,[],[],[],Fs);
                        Cxy_rand(imc,:) = Cxy;
                    end

                    for ifreq = 1:size(Cxy_rand,2) % 1:size(Cxy_rand,2)
                        x = Cxy_rand(:,ifreq);

                        % empirical CDF
                        figure(512)
                        [Fi,xi] = ecdf(x);
                        stairs(xi,Fi,'r');

                        xj = xi(2:end);
                        Fj = (Fi(1:end-1)+Fi(2:end))/2;
                        %hold on
                        %plot(xj,Fj,'b.', xj,Fj,'b-');
                        %hold off
                        %legend({'ECDF' 'Breakpoints' 'Piecewise Linear Estimate'},'location','NW');
                        xj = [xj(1)-Fj(1)*(xj(2)-xj(1))/((Fj(2)-Fj(1)));
                              xj;
                              xj(n)+(1-Fj(n))*((xj(n)-xj(n-1))/(Fj(n)-Fj(n-1)))];
                        Fj = [0; Fj; 1];
                        %hold on
                        %plot(xj,Fj,'b-');

                        %figure;
                        %stairs(Fi,[xi(2:end); xi(end)],'r');
                        %hold on
                        %plot(Fj,xj,'b-');
                        %hold off
                        %legend({'ECDF' 'Piecewise Linear Estimate'},'location','NW');
                        close(512)
                        Finv = @(u) interp1(Fj,xj,u,'linear','extrap');
                        co_threshold(ifreq) = Finv(0.99);
                    end

                    co_threshold_avg = mean(co_threshold);
                    co_thresholds_avg = [co_thresholds_avg; co_threshold_avg];


                    if exec_mscohere == 1

                        cohere_threshold = co_threshold_avg;

                        %% mscohere
                        Fbf = 1 / bf_width;
                        Fs = 1 / geoseries.dist_down(2);

                        [Cxy,F] = mscohere(Z,W,[],[],[],Fs);
                        [Pxy,F] = cpsd(Z,W,[],[],[],Fs);

                        f41 = figure; 
                        subplot(2,1,1)

                        % Peak analysis
                        
                        findpeaks(Cxy,F,'MinPeakHeight', cohere_threshold)
                        [pks, locs] = findpeaks(Cxy,F, 'MinPeakHeight',cohere_threshold);
                        
                        if length(gca().Children) == 2
                            obj1 = gca().Children(1); % current figure handle
                            obj1.Color = [0.7656    0.2812    0.1602];
                            obj1.MarkerFaceColor = [0.7656    0.2812    0.1602];
                            obj2 = gca().Children(2); % current figure handle
                            obj2.Color = [0.7656    0.2812    0.1602];
                        else
                            obj = gca().Children; % current figure handle
                            obj.Color = [0.7656    0.2812    0.1602];
                        end
                        hold on;
                        plot([Fbf, Fbf], [0, 1], 'k--'); 
                        legend(site_name, 'Peaks', "Bankfull width = "+num2str(round(bf_width, 3))+" ft")

                        xlabel('Frequency'); ylabel('Magnitude-squared coherence')
                        title("SC0" + num2str(channel_type) + ", Magnitude-Squared Coherence at " + water_stage_pt)
                        legend(site_name)
                        ylim([0.5, 1])
                        xlim([-inf inf])

                        % Angle
                        subplot(2,1,2)
                        angle_Pxy = angle(Pxy);
                        plot(F,angle_Pxy); hold on;
                        cur_ax = gca; XL = cur_ax.XLim;
                        
                        ang_pks = zeros(length(locs),1);
                        for kk = 1:length(locs)
                            ind_peak = find(F==locs(kk));
                            ang_pks(kk) = angle_Pxy(ind_peak);
                        end
                        plot(locs, ang_pks+0.4, 'v', ...
                            'MarkerEdgeColor', [0.7656    0.2812    0.1602], ...
                            'MarkerFaceColor', [0.7656    0.2812    0.1602])
                        xlim(XL);
                        xlabel('Frequency'); ylabel('Phase of the cross spectrum')
                        
                        saveas(f41, path_fig+"\mscohere\SC0"+num2str(channel_type)+"_"+site_stage,'png')
                        savefig(f41, path_fig+"\mscohere\fig\SC0"+num2str(channel_type)+"_"+site_stage)

                        %pause;
                        close(f41)
                        
                        %% plot
                        %f51 = figure(100); 
                        %plot(locs, ang_pks, '*'); hold on;

                        %figure(101)
                        %f52 = plot(locs, pks, '*'); hold on;
                        %pause;

                        if stage_num == 1
                            cum_locs_1 = [cum_locs_1; locs];
                            cum_ang_pks_1 = [cum_ang_pks_1; ang_pks];
                        elseif stage_num == 2
                            cum_locs_2 = [cum_locs_2; locs];
                            cum_ang_pks_2 = [cum_ang_pks_2; ang_pks];
                        else
                            cum_locs_3 = [cum_locs_3; locs];
                            cum_ang_pks_3 = [cum_ang_pks_3; ang_pks];
                        end

                    end

                    stage_num = stage_num + 1;
                end

            end

        end

    end    
    
    %% significant frequency analysis
    if exec_sig_freq == 1
        for kk = 1:3
            kk_str = num2str(kk);
            eval("cum_locs = cum_locs_"+kk_str+";")
            eval("cum_ang_pks = cum_ang_pks_"+kk_str+";")
    
            fc1 = figure;
            histogram(cum_locs, BinMethod ='fd')
            title("Significant frequency histogram - Channel type "+num2str(channel_type))
        
            saveas(fc1, path_fig+"\sig_freq\SC0"+num2str(channel_type)+"_freq_"+kk_str,'png')
            close(fc1) 
            save(path_fig+"\sig_freq\SC0"+num2str(channel_type)+"_freq_"+kk_str, "cum_locs_"+kk_str, '-mat')
            
            fc2 = figure;
            histogram(cum_ang_pks, BinMethod='fd')
            title("Angle histogram - Channel type "+num2str(channel_type))
        
            saveas(fc2, path_fig+"\sig_freq\SC0"+num2str(channel_type)+"_angle_"+kk_str,'png')
            close(fc2)
            save(path_fig+"\sig_freq\SC0"+num2str(channel_type)+"_angle_"+kk_str, "cum_ang_pks_"+kk_str, '-mat')
        end
    end
end




% Chi-square distribution: https://www.mathworks.com/help/stats/chi-square-distribution.html#mw_1b8489a4-79c3-484c-94cd-6b186c41b728
% Chi-square cumulative distribution function: p = chi2cdf(x,nu) https://www.mathworks.com/help/stats/chi2cdf.html
% Chi-square inverse cumulative distribution function: x = chi2inv(p,nu) https://www.mathworks.com/help/stats/chi2inv.html
