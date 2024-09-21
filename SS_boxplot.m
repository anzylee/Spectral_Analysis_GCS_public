

stage_name = ["Baseflow", "Bankfull", "Flood"];
ds_inf = spreadsheetDatastore('.\regression_analysis.xlsx');
[var, info] = read(ds_inf);
color_scatter = ["#0072BD","#D95319","#EDB120","#7E2F8E","#4DBEEE","#A2142F"];
        
indep_vars = ["Channel_Type"]';

indep_var_names = ["Channel type"]';

dep_vars = ["Baseflow_W_SS"; "Bankfull_W_SS";
    "Floodplain_W_SS"; "Baseflow_Z_SS";
    "Bankfull_Z_SS"; "Floodplain_Z_SS"]';

dep_var_names = ["Spectral slope of Wn"; "Spectral slope of Wn";
    "Spectral slope of Wn"; "Spectral slope of Zn";
    "Spectral slope Zn"; "Spectral slope of Zn"]';

legend_all = ["Spectral slope of baseflow width"; "Power-law fitting";
    "Spectral slope of bankfull width"; "Power-law fitting";
    "Spectral slope of floodplain width"; "Power-law fitting";
    "Spectral slope of baseflow elevation"; "Power-law fitting";
    "Spectral slope of bankfull elevation"; "Power-law fitting";
    "Spectral slope of floodplain elevation"; "Power-law fitting"]';

stage_names = ["Baseflow", "Bankfull", "Flood"];

fitting_method = ["power1"];
ind_indep = 1;
ind_st = 0;

for indep_var = indep_vars 

    fig_num = 1;
    ind_dep = 1;
    coeff_dep = 0;

    
    for dep_var = dep_vars 

        
        eval("X = var."+indep_var+";"); 
        eval("Y = var."+dep_var+";");
        
        [XY, ind_remove] = rmmissing([X, Y]);
        X = XY(:,1); Y = XY(:,2);

        figure
        hc = boxchart(X, Y, 'MarkerStyle','none'); % group by index
        hold on

        % overlay the scatter plots
        for n=1:max(unique(X))        
            hs = scatter(ones(sum(X==n),1) + n-1, Y(X==n),"filled",'jitter','on','JitterAmount',0.1);
            hs.MarkerFaceAlpha = 0.5;        
        end
        
        xlabel(indep_var_names(ind_indep))
        ylabel(dep_var_names(ind_dep))
        ind_dep = ind_dep + 1;

        xticks([1, 2, 3, 4, 5])
        xticklabels({'1', '2', '3', '4', '5'})

        
        ax = gca;
        set(ax,'Fontsize',19)
        title(stage_names(mod(ind_st, 3)+1), "FontSize", 23)
        set(ax.XLabel, 'FontWeight', 'bold')
        set(ax.YLabel, 'FontWeight', 'bold')

        ind_st = ind_st + 1;
    end
    
    ind_indep = ind_indep + 1;
end
