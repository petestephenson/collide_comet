
%%
[E_threshold,cross_sections,lines]=import_cross_sections('xs_set','old');
%Importing the new ajello xsection for 4PG emission from CO for comparison 
onedrive = 'C:\Users\pps18\OneDrive - Imperial College London\';
new_CO_4PG=csvread([one_drive, "ajello_2019_4PG_CO_xs.csv"],1,0);
new_CO_4PG=[new_CO_4PG(:,1),sum(new_CO_4PG(:,2:3),2)*0.2786e-18]; %Factor of 0.2786 is the percentage of the v'=0 emission that
%enters the 1653 band (see Beegle 1999)
%%
i_plot_CO2_CO_1657=1;
i_plot_CO2_CO_1335=0;


if i_plot_CO2_CO_1657
    figure;
    plot(cross_sections{2}{2}(1,:),cross_sections{2}{2}(5,:),'DisplayName', 'CO2 CI1657 and 4PG','LineWidth',2);
    hold on
    plot(cross_sections{2}{3}(1,:),cross_sections{2}{3}(5,:),'DisplayName', 'CO 4PG Beegle','LineWidth',2);
    hold on
    plot(new_CO_4PG(:,1),new_CO_4PG(:,2),'DisplayName', 'CO 4PG Ajello 2019','LineWidth',2)
    legend('location','southeast');
    set(gca, 'XLim',[6,8e2]);
    set(gca,'Ylim',[1e-19,inf]);
    set(gca, {'XScale', 'YScale'},{'log','log'});
    grid on
    xlabel('Energy$/eV$','Interpreter','latex');
    ylabel('Cross Section $/cm^{-2}$','Interpreter','latex');
    set(gca,'FontSize',15);
end


if i_plot_CO2_CO_1335
    figure;
    plot(cross_sections{2}{2}(1,:),cross_sections{2}{2}(6,:),'DisplayName', 'CO2 CII1335','LineWidth',2);
    hold on
    plot(cross_sections{2}{3}(1,:),cross_sections{2}{3}(4,:),'DisplayName', 'CO CII1335','LineWidth',2);
    legend('location','southeast');
    set(gca, 'XLim',[10,8e2]);
    set(gca, {'XScale', 'YScale'},{'log','log'});
    grid on
    xlabel('Energy$/eV$','Interpreter','latex');
    ylabel('Cross Section $/cm^{-2}$','Interpreter','latex');
    set(gca,'FontSize',15);
end
%% Grouped by Species

lines=[1026,1304,1356,1561,1657,1335];
lineName={'LyB','OI1304','OI1356','CI1561','CI1657 and 4PG','CII1335'};
speciesName={'H20','CO2','CO','O2'};
figure('Renderer', 'painters', 'Position', [10 10 900 700]);
colorset=get(gca,'ColorOrder');
for i_species=1:4
    subplot(2,2,i_species);
    nLines=size(cross_sections{2}{i_species},1)-1;
    
    ix_lines=arrayfun(@(x) find(x==lines,1,'first'),cross_sections{3}{i_species});
    for i_line=1:nLines
        plot(cross_sections{2}{i_species}(1,:),cross_sections{2}{i_species}(i_line+1,:),...
            'DisplayName',lineName{ix_lines(i_line)},'LineWidth', 2,'color',colorset(ix_lines(i_line),:));
        hold on
    end
    title(speciesName(i_species));
    set(gca, {'XScale', 'YScale'},{'log','log'});
    set(gca, {'XLim','Fontsize'},{[10,5e2],15});
    if i_species<4
        set(gca,'YLim',[1e-20,6e-18]);
    end
    xlabel('Energy$/eV$','Interpreter','latex');
    ylabel('Cross Section $/cm^{-2}$','Interpreter','latex');
    grid on
    legend('location','southeast');
end
%% Ajello 2019 XSections
E_ajello=[30,100];
XS_Ajello=cell(4,1);
ix_ajello_line=[2,3,5];
XS_Ajello{1}=nan(3,2);
XS_Ajello{2}=[1.61,6.03;0.51,1.87;nan,nan]*1e-19;
XS_Ajello{3}=[0.7,8.8;0.4,1.8;15.3,8.19]*1e-19;
XS_Ajello{4}=nan(3,2);

%%  Grouped by line
lines=[1026,1304,1356,1561,1657,1335];
lineName={'LyB','OI1304','OI1356','CI1561','CI1657 and 4PG','CII1335'};
speciesName={'H2O','CO2','CO','O2'};
figure('Renderer', 'painters', 'Position', [10 10 1600 800]);
colorset=get(gca,'ColorOrder');
nSpecies=4;
plot_lines=[1,2,3,5,6]; %Excludes CI1561
for i_plot_line=1:5
    i_line=plot_lines(i_plot_line);
    subplot(2,3,i_plot_line+1);

    ix_lines=cellfun(@(x) [boolean(0),ismember(x,lines(i_line))],cross_sections{3},'UniformOutput',0); 
    %Adding zero to front of matrix so that the energy is not plotted
    for i_species=1:nSpecies
        if any(ix_lines{i_species})
            plot(cross_sections{2}{i_species}(1,:),cross_sections{2}{i_species}(ix_lines{i_species},:),...
                'DisplayName',speciesName{i_species},'LineWidth', 2,'color',colorset(i_species,:));
            hold on
        end
%         [ajello_log,ajello_pos]=ismember(i_line,ix_ajello_line);
%         if ismember(i_species,[2,3])
%         if ajello_log
%             ax0=scatter(E_ajello,XS_Ajello{i_species}(ajello_pos,:),'DisplayName',[speciesName{i_species},' Ajello+19'],...
%                 'MarkerEdgeColor',colorset(i_species,:),'MarkerFaceColor',colorset(i_species,:),'Marker','*',...
%                 'SizeData',150,'LineWidth',1.2,'HandleVisibility','off'); 
%             hold on
%         end
%         end
    end
    title(lineName(i_line));
    set(gca, {'XScale', 'YScale'},{'log','log'});
    set(gca, {'XLim','Fontsize'},{[8,5e2],15});
%     if i_species<4
%         set(gca,'YLim',[1e-20,6e-18]);
%     end
    xlabel('Energy$/eV$','Interpreter','latex');
    ylabel('Cross Section $/cm^{-2}$','Interpreter','latex');
    grid on
    legend('location','southeast');
end
