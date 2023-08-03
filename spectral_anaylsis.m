clear all; close all; clc; warning('off') 

path_fig = ".\figures_ZnWn";

exec_variogram = 0;
exec_pmtm = 1;
exec_cpsd = 0;
exec_mscohere = 0; exec_sig_freq = exec_mscohere;
cohere_threshold = 0.5;
exec_boxplot = 1;
exec_violinplot = 1;
Fbf_factor = 1; % 999 if no threshold for frequency 

ds = tabularTextDatastore('.\4_14_21 Full data set.csv');
[full_data_set, info] = read(ds);

ds_inf = spreadsheetDatastore('.\inflection_points.xlsx');
[inflection_data, info_inf] = read(ds_inf);

ind_W = inflection_data.variable == "W";
ind_Z = inflection_data.variable == "Z";


for channel_type = 1:5 % 1:5
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

    tmp = num2str(cohere_threshold); cohere_threshold_p = replace(tmp, '.', 'p');

    for ii = 1:length(site_names)
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

        for table = table_all
            char_table = char(table);

            if length(char_table) > 21

                if char_table(end-20:end) == 'WD_analysis_table.csv'
                    %stage_num

                    table_split = table.split("_");
                    water_stage_p = table_split(1);
                    water_stage_pt = water_stage_p.replace("p", ".");
                    water_stage_m = str2double(water_stage_pt.replace("ft",""))* 0.3048;
                    %water_stage_pt_m = water_stage_pt.replace("ft", "m");

                    site_stage = site_name + "_" +num2str(stage_num);

                    table_file = table_dir + '\' + table;
                    datastore = tabularTextDatastore(table_file);
                    geoseries = read(datastore);    % in US cumstomary units
                    
                    %% normalization
                    index_comid = find(full_data_set.comid == str2num(site_num));
                    bf_width = full_data_set.BF_width_ft(index_comid) * 0.3048; % ft to meter
                    wet_width = mean(geoseries.W);

                    %% when figure_ZnWn

                    geoseries.W_n = geoseries.W/wet_width;
                    geoseries.Z_n = geoseries.Z/wet_width;

                    geoseries.X_n = geoseries.dist_down * 0.3048; % ft to meter


                    W = geoseries.W_n;
                    Z = geoseries.Z_n;
                    dist = geoseries.X_n;

                    %if stage_num == 1
                    %    channel_type
                    %    site_num
                    %    fprintf('Average Baseflow Width = %6.6f\n',mean(geoseries.W))
                    %    fprintf('Average Baseflow Depth = %6.6f\n',mean(geoseries.Z))
                    %elseif stage_num == 3
                    %     channel_type
                    %     site_num
                    %     fprintf('Average FP Width = %6.6f\n',mean(geoseries.W))
                    %     fprintf('Average FP Depth = %6.6f\n',mean(geoseries.Z))
                    %end

                    %% Variogram & fitting for fractal dimension (Spectral slope)

                    % 1. Variogram f(h) is the average squared difference 
                    %                      between the values at points
                    %                      separated by distance h.

                    if exec_variogram == 1
                        f11 = figure; 
                        d = variogram(dist, Z, ...
                            'plotit',true,'nrbins',50);
                        title("Isotropic variogram for Z at "+ num2str(round(water_stage_m, 2))+" m")
                        legend(site_name, 'Location', 'southeast')
                        saveas(f11, path_fig+"\variogram\SC0"+num2str(channel_type)+"_Z_n_"+site_stage,'png')

                        f12 = figure;
                        a0 = 150; % initial value: range 
                        c0 = 0.1; % initial value: sill 
                        [a,c,n] = variogramfit(d.distance,d.val,a0,c0,[],...
                           'solver','fminsearchbnd',...
                           'nugget',0,...
                           'plotit',true);
                        saveas(f12, path_fig+"\variogram\SC0"+num2str(channel_type)+"_Z_n_fit_"+site_stage,'png')

                        f13 = figure; 
                        d = variogram(dist, W, ...
                            'plotit',true,'nrbins',50);
                        title("Isotropic variogram for W at "+ num2str(round(water_stage_m, 2))+" m")
                        legend(site_name, 'Location', 'southeast')
                        saveas(f13, path_fig+"\variogram\SC0"+num2str(channel_type)+"_W_n_"+site_stage,'png')
    
                        f14 = figure;
                        a0 = 200; % initial value: range 
                        c0 = 100; % initial value: sill 
                        [a,c,n] = variogramfit(d.distance,d.val,a0,c0,[],...
                           'solver','fminsearchbnd',...
                           'nugget',0,...
                           'plotit',true);
                        saveas(f14, path_fig+"\variogram\SC0"+num2str(channel_type)+"_W_n_fit_"+site_stage,'png')

                        %pause;
                        close(f11, f12, f13, f14);
                    end
                      

                    %% Spectral analysis & Significance test
                    if exec_pmtm == 1
                        Fbf_orig = 1 / bf_width ;
                        Fbf = Fbf_orig * Fbf_factor;
                        Fs = 1 / dist(2);
                        ind_comid = inflection_data.COMID==site_name;
                        ind_stage = inflection_data.stage==stage_num;
                        
                        %%
                        f21 = figure;                        
                        [Pxx_Z,F] = pmtm(Z, [], [], Fs);

                        % linear fitting
                        XX = log10(F); YY = log10(Pxx_Z);
                        
                        %ind_var = ind_stage(inflection_data.variable(ind_stage)=="Z");
                        ind_var = find(ind_stage.*ind_comid.*ind_Z);
                        
                        %% inflection points for Z

                        first_ceiled = inflection_data.first_ceiled(ind_var) - log10(0.3048);
                        second_ceiled = inflection_data.second_ceiled(ind_var) -log10(0.3048);


                        ind_inf1 = find(XX-(first_ceiled) <= 0); ind_inf1 = ind_inf1(end);
                        ind_inf2 = find(XX-(second_ceiled) <= 0); ind_inf2 = ind_inf2(end);
                        
                        ind_end = length(XX);
                        %% defining cut-off freq
%                         if log10(Fbf) > inflection_data.second_ceiled(ind_var) && log10(Fbf) < max(XX)
%                             ind_end = find(XX-log10(Fbf) < 0); ind_end = ind_end(end);
%                         else 
%                             ind_end = length(XX);
%                         end
% 
%                         if ind_end == ind_inf2
%                             ind_end = length(XX);
%                         end
% 
%                         if Fbf_factor == 999
%                             ind_end = length(XX);
%                         end

                        XX1 = XX(ind_inf1:ind_inf2); YY1 = YY(ind_inf1:ind_inf2); f1 = fit(XX1, YY1, 'poly1');
                        XX2 = XX(ind_inf2:ind_end); YY2 = YY(ind_inf2:ind_end); f2 = fit(XX2, YY2, 'poly1');
                        
                        % PSD
                        plot(XX, YY); hold on;
                        obj = gca().Children; % current figure handle
                        obj.Color = [1 0.411764705882353 0.16078431372549];
                        
                        % Fitting curves
                        %plot(XX1, f1.p1*XX1 + f1.p2, '--'); 
                        plot(XX2, f2.p1*XX2 + f2.p2, '--'); 

                        plot(log10([Fbf_orig, Fbf_orig]), log10([min(Pxx_Z); max(Pxx_Z)]), 'k--'); 
                        %plot(log10([Fbf, Fbf]), log10([min(Pxx_Z); max(Pxx_Z)]), 'k--');

                        % Slope breaks
                        %plot([XX(ind_inf1), XX(ind_inf2), XX(ind_end)], ...
                        %    [YY(ind_inf1), YY(ind_inf2), YY(ind_end)],'*')
                        plot(XX(ind_inf2), YY(ind_inf2),'*')

                        legend('PSD', 'Fitting curve', ...
                            "Freq. for bankfull width", ...
                            "Slope break")

                        %ax = gca;
                        %ax.XScale = 'log';
                        %ax.YScale = 'log';
                        xlabel('Log_{10}(Frequency)'); ylabel('Log_{10}(Power/frequency)')
                        title("Multitaper PSD Estimate for Z at "+ num2str(round(water_stage_m, 2))+" m")
                        saveas(f21, path_fig+"\pmtm\SC0"+num2str(channel_type)+"_"+site_stage+"_Z_n",'png')
                        savefig(f21, path_fig+"\pmtm\fig\SC0"+num2str(channel_type)+"_"+site_stage+"_Z_n")
                        inflection_data.Fbf(ind_var) = Fbf;
                        inflection_data.slope1(ind_var) = f1.p1; inflection_data.slope2(ind_var) = f2.p1;

                        %%
                        f22 = figure;
                        [Pxx_W, F] = pmtm(W, [], [], Fs);

                        % linear fitting
                        XX = log10(F); YY = log10(Pxx_W);
                        
                        ind_var = find(ind_stage.*ind_comid.*ind_W);

                        %% inflection points for W

                        first_ceiled = inflection_data.first_ceiled(ind_var) - log10(0.3048);
                        second_ceiled = inflection_data.second_ceiled(ind_var) -log10(0.3048);

                        ind_end = length(XX);

                        %% defining cut-off freq
%                         if log10(Fbf) > inflection_data.second_ceiled(ind_var) && log10(Fbf) < max(XX)
%                             ind_end = find(XX-log10(Fbf) < 0); ind_end = ind_end(end);
%                         else 
%                             ind_end = length(XX);
%                         end
%                         
%                         if ind_end < ind_inf2 + 10
%                             ind_end = length(XX);
%                         end
% 
%                         if Fbf_factor == 999
%                             ind_end = length(XX);
%                         end

                        XX1 = XX(ind_inf1:ind_inf2); YY1 = YY(ind_inf1:ind_inf2); f1 = fit(XX1, YY1, 'poly1');
                        XX2 = XX(ind_inf2:ind_end); YY2 = YY(ind_inf2:ind_end); f2 = fit(XX2, YY2, 'poly1');

                        plot(XX, YY); hold on;

                        %plot(XX1, f1.p1*XX1 + f1.p2, '--'); 
                        plot(XX2, f2.p1*XX2 + f2.p2, '--'); 
                        
                        plot(log10([Fbf_orig, Fbf_orig]), log10([min(Pxx_W), max(Pxx_W)]), 'k--'); 
                        %plot(log10([Fbf, Fbf]), log10([min(Pxx_W), max(Pxx_W)]), 'k--'); 
                        
                        plot(XX(ind_inf2), YY(ind_inf2), '*')

                        legend('PSD', 'Fitting curve', ...
                            "Freq. for bankfull width", ...
                            "Slope break")

                        %ax = gca;
                        %ax.XScale = 'log';
                        %ax.YScale = 'log';
                        xlabel('Log_{10}(Frequency)'); ylabel('Log_{10}(Power/frequency)')
                        title("Multitaper PSD Estimate for W at "+ num2str(round(water_stage_m, 2))+" m")
                        saveas(f22, path_fig+"\pmtm\SC0"+num2str(channel_type)+"_"+site_stage+"_W_n",'png')
                        savefig(f22, path_fig+"\pmtm\fig\SC0"+num2str(channel_type)+"_"+site_stage+"_W_n")
                        inflection_data.Fbf(ind_var) = Fbf;
                        inflection_data.slope1(ind_var) = f1.p1; inflection_data.slope2(ind_var) = f2.p1;

                        %pause;
                        close(f21, f22)
                    end

                    %% Magnitude-Squared Coherence and phase of the cross spectrum
                    if exec_cpsd == 1
                        Fs = 1 / geoseries.dist_down(2);
                        f31 = figure;
                        [Pxy,F] = cpsd(Z,W,[],[],[],Fs);
                        plot(F, Pxy)
                        legend(site_name)
                        saveas(f31, path_fig+"\cpsd\SC0"+num2str(channel_type)+"_"+site_stage,'png')

                        %pause;
                        close(f31);
                    end

                    if exec_mscohere == 1

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
                        legend(site_name, 'Peaks', "Bankfull width = "+num2str(round(bf_width, 2))+" m")

                        xlabel('Frequency'); ylabel('Magnitude-squared coherence')
                        title("SC0" + num2str(channel_type) + ", Magnitude-Squared Coherence at "+ num2str(round(water_stage_m, 2))+" m")
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
                        
                        saveas(f41, path_fig+"\mscohere\SC0"+num2str(channel_type)+"_"+cohere_threshold_p+"_"+site_stage,'png')
                        savefig(f41, path_fig+"\mscohere\fig\SC0"+num2str(channel_type)+"_"+cohere_threshold_p+"_"+site_stage)

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
        
            saveas(fc1, path_fig+"\sig_freq\SC0"+num2str(channel_type)+"_freq_"+cohere_threshold_p+"_"+kk_str,'png')
            close(fc1) 
            save(path_fig+"\sig_freq\SC0"+num2str(channel_type)+"_freq_"+cohere_threshold_p+"_"+kk_str, "cum_locs_"+kk_str, '-mat')
       
            % weighted significant frequency - weighted by its coherence ??
            
            fc2 = figure;
            histogram(cum_ang_pks, BinMethod='fd')
            title("Angle histogram - Channel type "+num2str(channel_type))
        
            saveas(fc2, path_fig+"\sig_freq\SC0"+num2str(channel_type)+"_angle_"+cohere_threshold_p+"_"+kk_str,'png')
            close(fc2)
            save(path_fig+"\sig_freq\SC0"+num2str(channel_type)+"_angle_"+cohere_threshold_p+"_"+kk_str, "cum_ang_pks_"+kk_str, '-mat')
        end
    end

end



stage_num = 2;

ind_stage = inflection_data.stage==stage_num;
ind_W2 =  find(ind_W.*ind_stage);
ind_Z2 = find(ind_Z.*ind_stage);




writetable(inflection_data,'.\inflection_points_slope.xlsx')


% Chi-square distribution: https://www.mathworks.com/help/stats/chi-square-distribution.html#mw_1b8489a4-79c3-484c-94cd-6b186c41b728
% Chi-square cumulative distribution function: p = chi2cdf(x,nu) https://www.mathworks.com/help/stats/chi2cdf.html
% Chi-square inverse cumulative distribution function: x = chi2inv(p,nu) https://www.mathworks.com/help/stats/chi2inv.html
