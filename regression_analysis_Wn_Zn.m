clear all; close all; clc


stage_name = ["Baseflow", "Bankfull", "Flood"];
stages =  ["Base", "Bankfull", "Floodplain"] ;
scolor_stage = ["#0072BD", "#D95319", "#EDB120"];
ind_ab = 1;
method = "all"; %"by_stage"; "all";
color = "by_stage";

for variable = ["Zn"]

    ds_inf = spreadsheetDatastore('.\ZnWn\slope_'+variable+'.xlsx');
    [slope, info] = read(ds_inf);

    if method == "by_stage"
        
        ind_st = 1;
        for stage =  stages
         
            ind = slope.water_stage == stage;
            X = slope.X(ind);
            Y = slope.Y(ind);


            [curvefit,gof,output] = fit(X, Y, 'power1');

            a = num2str(round(curvefit.a,2));
            b = num2str(round(curvefit.b,2));
            R2 = num2str(round(gof.rsquare,2));

            XX = linspace(min(X), max(X), 30);
            YY = curvefit(XX);

            figure(); 
            s = scatter(X, Y, 'o', 'filled', 'MarkerFaceColor',scolor_stage(ind_st)); hold on; 
%             plot(X,Y, 'o', 'MarkerFaceColor',"#0072BD"); hold on;
            plot(XX, YY, 'k--', "LineWidth", 3)
            
            ax = gca;
            set(ax, "FontSize", 14)
            xlabel(stage_name(ind_st)+" width (m)", "FontSize", 17)
            ylabel("Spectral slope of "+variable, "FontSize", 17)
            str = "y = "+num2str(a)+"x^{"+num2str(b)+"}, R^2 = "+num2str(R2);
            dim = [0.44, 0.15, 0.1, 0.1];
            ylim([0,3])

            title(stage_name(ind_st), "FontSize", 23)
            
            t = annotation("textbox", dim, 'LineStyle', "none",'String', str, ...
                'FitBoxToText','on', "FontSize", 17);
    
            saveas(gcf,".\slope_"+variable+"_"+stage+".emf")
    
            ind_ab = ind_ab + 1;
            ind_st = ind_st + 1;
        end
    else
        X = slope.X; Y= slope.Y;

        % Original data
        figure
        if color == "by_stage"
            for stage = stages
                ind = slope.water_stage == stage;
                %[X, ind_sort_X] = sort(slope.X(ind));
                X = slope.X(ind);
                Y = slope.Y(ind);
                s = scatter(X, Y, 'o', 'filled'); hold on; 
                scolor = [0, 0, 0];
            end

        else
            s = scatter(X, Y, 'o', 'filled'); hold on; 
            scolor = s.CData;
        end

        X = slope.X; Y= slope.Y;
        if variable == "Zn"
            [curvefit,gof,output] = fit(X, Y, 'power1');

            a = num2str(round(curvefit.a,2));
            b = num2str(round(curvefit.b,2));
            R2 = num2str(round(gof.rsquare,2));
        else
            [curvefit,gof,output] = fit(X, Y, 'poly1');

            a = num2str(round(curvefit.p1,2));
            b = num2str(round(curvefit.p2,2));
            R2 = num2str(round(gof.rsquare,2));
            
        end

        XX = linspace(min(X), max(X), 30);
        YY = curvefit(XX);
        plot(XX, YY,  'k--', "LineWidth", 3)
        
        %title(indep_var_names(ind_indep), "FontSize", 19)
        ax = gca;
        set(ax, "FontSize", 14)


        if variable == "Zn"
            xlabel("Wetted width (m)", "FontSize", 17)
            ylabel("Spectral slope of Zn", "FontSize", 17)
            str = "y = "+a+"x^{"+b+"}, R^2 = "+R2;
            dim = [0.44, 0.15, 0.1, 0.1];
            ylim([0,3])
        else
            xlabel("Spectral slope of Zn", "FontSize", 17)
            ylabel("Spectral slope of Wn", "FontSize", 17)
            str = "y = "+a+"x + "+b+ ", R^2 = "+R2;
            dim = [0.41, 0.15, 0.1, 0.1];
            ylim([0, 3])
        end
        
        legend("Baseflow", "Bankfull", "Flood", "Fitting curve")
        t = annotation("textbox", dim, 'LineStyle', "none",'String', str, ...
            'FitBoxToText','on', "FontSize", 17);
        
        title("All stages", "FontSize", 23)

        saveas(gcf,".\slope_"+variable+".emf")

        ind_ab = ind_ab + 1;
        
    end

end