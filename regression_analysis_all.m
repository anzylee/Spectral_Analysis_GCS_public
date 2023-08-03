clear all; close all; clc


stage_name = ["Baseflow", "Bankfull", "Flood"];
ds_inf = spreadsheetDatastore('.\regression_analysis.xlsx');
[var, info] = read(ds_inf);
color_scatter = ["#0072BD","#D95319","#EDB120","#7E2F8E","#4DBEEE","#A2142F"];
        
indep_vars = ["Stream_order"; "Catchment_area_km2";
     "Valley_confinement"; "Channel_slope";
    "CV_bf_d"; "CV_bf_W";
    "Baseflow_W"; "Baseflow_d";
    "Bankfull_W"; "Bankfull_d";
    "Floodplain_W"; "Floodplain_d";
    "Baseflow_Wd"; "Bankfull_Wd"; "Floodplain_Wd"]';

indep_var_names = ["Stream order"; "Catchment area (km^2)";
    "Valley confinement (m)"; "Channel slope"
    "Coefficient of variation of bankfull depth"; "Coefficient of variation of bankfull width"
    "Baseflow width (m)"; "Baseflow depth (m)";
    "Bankfull width (m)"; "Bankfull depth (m)";
    "Floodplain width (m)"; "Floodplain depth (m)";    
    "Baseflow W/d"; "Bankfull W/d"; "Floodplain W/d";]';

dep_vars = ["Baseflow_W_SS"; "Bankfull_W_SS";
    "Floodplain_W_SS"; "Baseflow_Z_SS";
    "Bankfull_Z_SS"; "Floodplain_Z_SS"]';

dep_var_names = ["Spectral slope of baseflow width"; "Spectral slope of bankfull width";
    "Spectral slope of floodplain width"; "Spectral slope of baseflow elevation";
    "Spectral slope of bankfull elevation"; "Spectral slope of floodplain elevation"]';

legend_all = ["Spectral slope of baseflow width"; "Power-law fitting";
    "Spectral slope of bankfull width"; "Power-law fitting";
    "Spectral slope of floodplain width"; "Power-law fitting";
    "Spectral slope of baseflow elevation"; "Power-law fitting";
    "Spectral slope of bankfull elevation"; "Power-law fitting";
    "Spectral slope of floodplain elevation"; "Power-law fitting"]';

fitting_method = ["power1"];
ind_indep = 1;

for indep_var = indep_vars 

    fig_num = 1;
    ind_dep = 1;
    coeff_dep = 0;
            
    for dep_var = dep_vars 
        

        eval("X = var."+indep_var+";"); 
        eval("Y = var."+dep_var+";");
        
        [XY, ind_remove] = rmmissing([X, Y]);
        X = XY(:,1); Y = XY(:,2);
        
        % Original data
        figure(fig_num)
        s = scatter(X, Y, 'o', 'filled', 'MarkerFaceColor', color_scatter(ind_dep), ...
            'MarkerEdgeColor',color_scatter(ind_dep)); hold on; 
        scolor = s.MarkerFaceColor;

        % Fitting curve
        if fitting_method == "power1"
            [curvefit,gof,output] = fit(X, Y, 'power1');
    
            XX = linspace(min(X), max(X), 30);
            YY = curvefit(XX);
            plot(XX, YY, '--', "Color", scolor)
            
            %title(indep_var_names(ind_indep), "FontSize", 19)
            ax = gca;
            set(ax, "FontSize", 14)
    
            % Annotation
            a = num2str(round(curvefit.a,2));
            b = num2str(round(curvefit.b,2));
            R2 = num2str(round(gof.rsquare,2));
    
            dim = [0.52, 0.85-coeff_dep*0.08, 0.1, 0.1];
    
            str = "y = "+a+"x^{"+b+"}, R^2 = "+R2;
            t = annotation("textbox", dim, ...
                'LineStyle', "none",'String', str, ...
                'FitBoxToText','on', ...
                "FontSize", 15, "Color", scolor);

        end

        ind_dep = ind_dep + 1;
        coeff_dep = coeff_dep + 1;

        if coeff_dep == 3
            coeff_dep = coeff_dep + 4;
        end

    end
    
    xlabel(indep_var_names(ind_indep))
    ylabel('Spectral slopes')

    %legend(legend_all)

    ax = gca;
    xlim_orig = get(ax, "XLim");
    set(ax, "XLim", [xlim_orig(1), round(xlim_orig(2)*1.2,1)])
    set(ax, "YLim", [0 3])

    %pause(100);

    savefig('.\regression_analysis\'+indep_var)
    saveas(gcf, '.\regression_analysis\'+indep_var+'.emf')
    close all;


    ind_indep = ind_indep + 1;
end