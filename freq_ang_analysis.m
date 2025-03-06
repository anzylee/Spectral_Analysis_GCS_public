clear all; close all; clc

% coherence thresholds for each river type: 
% [0.5000, 0.5254, 0.5408, 0.5369, 0.5406]

exec_dot = 1;
exec_cip = 0;

num_sites = [6, 5, 7, 8, 8];
colors = ["Winter", "Autumn", "Summer"];
colors = [[0, 0.45, 0.74];
    [0.85, 0.33, 0.1];
    [0.93, 0.69, 0.13]];
colors = [[255, 70, 70];
    [20, 210, 255];
    [100, 255, 50]]/256;
colors = [[0 0.4470 0.7410];
    [0.8500 0.3250 0.0980];
    [0.9290 0.6940 0.1250]];
markers = ["square", "diamond", "o"];
scolor_stage = ["#0072BD", "#D95319", "#EDB120"];
Fbf = 1./[59.6560, 384.2256, 46.5278, 58.1230, 68.9276];

unit_freq = 'feet';
med_ang_all = []; std_ang_all = [];
med_freq_all = []; std_freq_all = []; 
perc = 95; perc_freq_all = []; perc_bf_all = [];

if matches(unit_freq, 'feet')
    Fbf = Fbf/0.3048;
end

path_fig = ".\figures_ZnWn";

ang_all_all = []; freq_all_all = []; water_stage_all_all = [];

for channel_type = 1:5
    fdot = figure;
    fhist = figure;
    %fcon = figure;
    %fbox = figure;
    fcip = figure;
    fcip_cust = figure;
    len = 0;
    ang_all = []; freq_all = []; water_stage_all = [];

    for water_stage = 1:3
        %len_stage = 0;
        ang = importdata(path_fig+"\sig_freq\SC0"+num2str(channel_type)+"_angle_"+num2str(water_stage));
        ind_angpiover = find(ang>pi);
        ang(ind_angpiover) = ang(ind_angpiover) - 2*pi;
        ind_angpiunder = find(ang<-pi);
        ang(ind_angpiunder) = ang(ind_angpiunder) + 2*pi;

        freq = importdata(path_fig+"\sig_freq\SC0"+num2str(channel_type)+"_freq_"+num2str(water_stage));
        if matches(unit_freq, 'feet')
            freq = freq/0.3048;
            unit_freq_plt = "m";
        else
            unit_freq_plt = "ft";
        end
        figure(fdot);
        plot(freq, ang, '*'); hold on;
        xlabel("Frequency"); ylabel("Phase"); ylim([-3.14 3.14]);

        len_stage = length(freq);
        len = len + length(freq);

        fprintf("Channel type = "+num2str(channel_type)+ ...
        ", Water stage = "+ num2str(water_stage)+...
        ", Significant frequency per site = "+ num2str(len_stage/num_sites(channel_type))+"\n");
        
        color_stage = colors(water_stage, 1:3);

%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Confidence level plot
        clear freq_bin
        if exec_dot == 1
            if exec_cip == 1
                for ii = 1:length(freq)
                    figure(fhist)
                    %h = histogram(freq, 'BinMethod', 'sturges');
                    h = histogram(freq, 20);
                    ind_bin = find(freq(ii) < h.BinEdges);
                    freq_bin(ii) = h.BinEdges(ind_bin(1));
                end
                
                figure(fcip)
                cip = confidenceIntervalPlot(freq_bin, ang);
                cip.NumSteps = 10;
        
                title("95% Confidence Interval Plot, Sine Curve with Random Noise");
                subtitle("Alpha = 0.05");
        
                fprintf("Press Enter: \n")
                pause(5)
                
                
                XX = [0, cip.CenterXData];
                YY = [0, cip.CenterYData];
                YUB = [0, cip.UpperBoundData];
                YLB = [0, cip.LowerBoundData];
        
        
                XX_sp = linspace(min(XX), max(XX), 100);
                YY_sp = interp1(XX, YY, XX_sp, 'makima'); % interp1 -> spline, makima
                YUB_sp = interp1(XX, YUB, XX_sp, 'makima');
                YLB_sp = interp1(XX, YLB, XX_sp, 'makima');
        
                xconf = [XX_sp XX_sp(end:-1:1)];
                yconf = [YUB_sp YLB_sp(end:-1:1)];
                
                figure(fcip_cust)
                p = fill(xconf,yconf,'red');
                p.FaceColor = color_stage;      
                p.EdgeColor = 'none';   
                p.FaceAlpha = 0.1;
        
                hold on;
                cline = plot(XX_sp, YY_sp, '--');
                cline.Color = color_stage;
                %plot(cip.XData, cip.YData, markers(water_stage), 'MarkerSize', 8, 'MarkerFaceColor', color_stage,'MarkerEdgeColor', color_stage);     
            end
            
            figure(fcip_cust)
            plot(freq, ang, markers(water_stage), ...
            'MarkerSize', 8, 'MarkerFaceColor', color_stage,'MarkerEdgeColor', color_stage); hold on;
            
            freq_all = cat(1, freq_all, freq);
            ang_all = cat(1, ang_all, ang);
            water_stage_all = cat(1, water_stage_all, water_stage*ones(length(freq),1));
            med_ang_all = cat(1, med_ang_all, median(ang));
            std_ang_all = cat(1, std_ang_all, std(ang));
            med_freq_all = cat(1, med_freq_all, median(freq));
            std_freq_all = cat(1, std_freq_all, std(freq));
            perc_freq_all = cat(1, perc_freq_all, prctile(freq, perc));
            perc_bf_all = cat(1, perc_bf_all, Fbf(channel_type)/prctile(freq, perc)*100);

            if water_stage == 3
                plot([Fbf(channel_type), Fbf(channel_type)], [-pi, pi], 'k-'); 
                plot([Fbf(channel_type)*5, Fbf(channel_type)*5], [-pi, pi], 'k--', 'LineWidth',1); 
                xlim_org = get(gca, 'xlim');
                plot(xlim_org, zeros(size(xlim_org)), 'k-')
            end
            
            xlabel("Significant frequency ("+unit_freq_plt+"^{-1})");
            ylabel('Phase of cross spectrum (rad)');
            acip = gca; acip.FontSize = 18;
%             xlim([0, 0.2])
            ylim([-3.14 3.14]);

        end
    end
    
    figure(fdot)
    legend("Base", "Bankfull", "Floodplain")
    saveas(fdot, path_fig+"\sig_freq\SC0"+num2str(channel_type)+"_hist.png")

    figure
    sh = scatterhist(freq_all, ang_all,'Group',water_stage_all, ...
        'Marker', markers, 'MarkerSize',[8,8,8], 'Kernel','on'); hold on;
    for sh_ind = 1:3
        sh(1).Children(4-sh_ind).MarkerFaceColor = scolor_stage(sh_ind);
    end
    xlim([0, 0.6]); 
    ylim([-3.14 3.14]);
%     plot([Fbf(channel_type), Fbf(channel_type)], [-pi, pi], 'k-'); 
%     plot([Fbf(channel_type)*5, Fbf(channel_type)*5], [-pi, pi], 'k--', 'LineWidth',1); 
    xlim_org = get(gca, 'xlim');
    plot(xlim_org, zeros(size(xlim_org)), 'k-')

    xlabel("Significant frequency ("+unit_freq_plt+"^{-1})", 'VerticalAlignment','middle');
    ylabel('Phase (rad)', 'VerticalAlignment','middle');
    legend('Baseflow', 'Bankfull', 'Flood')
    set(gca, 'FontSize', 16)

    saveas(gcf, path_fig+"\sig_freq\SC0"+num2str(channel_type)+"_scatterhist.png")
    saveas(gcf, path_fig+"\sig_freq\SC0"+num2str(channel_type)+"_scatterhist.pdf")
    
    figure(fcip_cust)
    if channel_type == 1
        legend("Base", "Bankfull", "Flood")
    end
    saveas(fcip_cust, path_fig+"\sig_freq\SC0"+num2str(channel_type)+"_cip.emf")
    savefig(fcip_cust, path_fig+"\sig_freq\SC0"+num2str(channel_type)+"_cip")
    close all;

    fprintf("Channel type = "+num2str(channel_type)+ ...
        ", Significant frequency per site = "+ num2str(len/num_sites(channel_type))+"\n");


    freq_all_all = cat(1, freq_all_all, freq_all);
    ang_all_all = cat(1, ang_all_all, ang_all);
    water_stage_all_all = cat(1, water_stage_all_all, water_stage_all);

end
% median(ang_all_all)
% median(ang_all_all(water_stage_all_all==1))
% median(ang_all_all(water_stage_all_all==2))
% median(ang_all_all(water_stage_all_all==3))
% prctile(ang_all_all(water_stage_all_all==1), 100-perc)
% prctile(ang_all_all(water_stage_all_all==2), 100-perc)
% prctile(ang_all_all(water_stage_all_all==3), 100-perc)
% prctile(ang_all_all(water_stage_all_all==1), perc)
% prctile(ang_all_all(water_stage_all_all==2), perc)
% prctile(ang_all_all(water_stage_all_all==3), perc)

fprintf('std by water stage')
std(ang_all_all(water_stage_all_all==1))
std(ang_all_all(water_stage_all_all==2))
std(ang_all_all(water_stage_all_all==3))

fprintf('std, median, percentile all')
std(ang_all_all)
median(freq_all_all)
prctile(freq_all_all, perc)


figure
sh = scatterhist(freq_all_all, ang_all_all, 'Group', water_stage_all_all, ...
    'Marker', markers, 'MarkerSize',[4, 3, 3], 'Kernel','on'); hold on;
for sh_ind = 1:3
    sh(1).Children(4-sh_ind).MarkerFaceColor = scolor_stage(sh_ind);
end

% xlim([0, 0.6]); 
ylim([-3.14 3.14]);
%     plot([Fbf(channel_type), Fbf(channel_type)], [-pi, pi], 'k-'); 
%     plot([Fbf(channel_type)*5, Fbf(channel_type)*5], [-pi, pi], 'k--', 'LineWidth',1); 
xlim_org = get(gca, 'xlim');
plot(xlim_org, zeros(size(xlim_org)), 'k-')


xlabel("Significant frequency ("+unit_freq_plt+"^{-1})", 'VerticalAlignment','middle');
ylabel('Phase (rad)', 'VerticalAlignment','middle');
legend('Baseflow', 'Bankfull', 'Flood')
set(gca, 'FontSize', 16)


saveas(gcf, path_fig+"\sig_freq\SC0"+num2str(channel_type)+"_scatterhist_all.png")
saveas(gcf, path_fig+"\sig_freq\SC0"+num2str(channel_type)+"_scatterhist_all.pdf")