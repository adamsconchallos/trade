%---------------------------------------------------------------
%This program constructs all figures and some auxilliary results
%---------------------------------------------------------------
%Preliminary calculations
clear all
close all
clc

versions={'Main','Low_Sig','Medlow_Sig','Medhigh_Sig','High_Sig','Novar_Sig'};
projectRoot = fileparts(mfilename('fullpath'));
origDir = pwd;
figDir = fullfile(projectRoot,'Figures');
vers = versions{1};
% Copy key datasets and helper code into Figures/ for self-contained
% plotting. These temporary files are removed at the end of the script.
copyfile(fullfile(projectRoot,'Data','Main','*.mat'),figDir)
copyfile(fullfile(projectRoot,'Programs','*.m'),figDir)
copyfile(fullfile(projectRoot,'Results','Optimal tariffs',vers,'*.mat'),figDir)
copyfile(fullfile(projectRoot,'Results','Trade wars',vers,'*.mat'),figDir)
cd(figDir)
mycalculations

% Figure 1 uses OPTIMALTARIFFBAS.mat to display optimal tariffs without
% lobbying. Tariffs are reshaped to an exporter-by-sector matrix, merged
% with elasticity data (SIGMA) for ranking, and plotted for each importer.
%Constructing Figure 1: Optimal tariffs without lobbying
load OPTIMALTARIFFBAS
for j=1:N
OPTIMALTARIFFj=OPTIMALTARIFFBAS(:,j);
TEMP=reshape(OPTIMALTARIFFj,N-1,S);   % (N-1) exporters × S sectors
TEMP=[TEMP(1:j-1,:);zeros(1,S);TEMP(j:end,:)];
TEMP=[100*TEMP;SIGMA'];
TEMP=sortrows(TEMP',8)';
TEMP=[TEMP;1:1:S];
if j<=3
    subplot(4,3,6+j)
    position=get(gca,'position');
    position=position+[0.01 0 -0.02 0];
    subplot('position',position)
    if j~=1
    scatter(TEMP(9,:),TEMP(1,:),10,'k','o')
    hold on
    end
    if j~=2
    scatter(TEMP(9,:),TEMP(2,:),10,'k','+')
    hold on
    end
    if j~=3
    scatter(TEMP(9,:),TEMP(3,:),10,'k','*')
    hold on
    end
    if j~=4
    scatter(TEMP(9,:),TEMP(4,:),10,'k','x')
    hold on
    end
    if j~=5
    scatter(TEMP(9,:),TEMP(5,:),10,'k','s')
    hold on
    end
    if j~=6
    scatter(TEMP(9,:),TEMP(6,:),10,'k','d')
    hold on
    end
    if j~=7
    scatter(TEMP(9,:),TEMP(7,:),10,'k','^')
    end
    Names={'Brazil','China','EU','India','Japan','RoW','US'};
    axis([1 S 0 100])
    set(gca,'fontsize',7,'xtick',1:4:S);
    title(['Optimal tariffs ' Names{j}])
    xlabel('Industry rank (lowest sigma to highest sigma)');
    ylabel('Optimal tariff in %');
elseif j>3 && j<7
    subplot(4,3,6+j)
    position=get(gca,'position');
    position=position+[0.01 -0.04 -0.02 0];
    subplot('position',position)
    if j~=1
    scatter(TEMP(9,:),TEMP(1,:),10,'k','o')
    hold on
    end
    if j~=2
    scatter(TEMP(9,:),TEMP(2,:),10,'k','+')
    hold on
    end
    if j~=3
    scatter(TEMP(9,:),TEMP(3,:),10,'k','*')
    hold on
    end
    if j~=4
    scatter(TEMP(9,:),TEMP(4,:),10,'k','x')
    hold on
    end
    if j~=5
    scatter(TEMP(9,:),TEMP(5,:),10,'k','s')
    hold on
    end
    if j~=6
    scatter(TEMP(9,:),TEMP(6,:),10,'k','d')
    hold on
    end
    if j~=7
    scatter(TEMP(9,:),TEMP(7,:),10,'k','^')
    end
    Names={'Brazil','China','EU','India','Japan','ROW','US'};
    axis([1 S 0 100])
    set(gca,'fontsize',7,'xtick',1:4:S);
    title(['Optimal tariffs ' Names{j}])
    xlabel('Industry rank (lowest sigma to highest sigma)');
    ylabel('Optimal tariff in %');
elseif j==7
    subplot(4,3,[1 6])
    position=get(gca,'position');
    position=position+[0.2 0.04 -0.4 0];
    subplot('position',position)
    if j~=1
    scatter(TEMP(9,:),TEMP(1,:),10,'k','o')
    hold on
    end
    if j~=2
    scatter(TEMP(9,:),TEMP(2,:),10,'k','+')
    hold on
    end
    if j~=3
    scatter(TEMP(9,:),TEMP(3,:),10,'k','*')
    hold on
    end
    if j~=4
    scatter(TEMP(9,:),TEMP(4,:),10,'k','x')
    hold on
    end
    if j~=5
    scatter(TEMP(9,:),TEMP(5,:),10,'k','s')
    hold on
    end
    if j~=6
    scatter(TEMP(9,:),TEMP(6,:),10,'k','d')
    hold on
    end
    if j~=7
    scatter(TEMP(9,:),TEMP(7,:),10,'k','^')
    end
    Names={'Brazil','China','EU','India','Japan','RoW','US'};
    axis([1 S 0 100])
    set(gca,'fontsize',7,'xtick',1:2:S);
    title(['Optimal tariffs ' Names{j}])
    xlabel('Industry rank (lowest sigma to highest sigma)');
    ylabel('Optimal tariff in %');
    Names(j)=[];
    legend1=legend(Names);
    set(legend1,'FontSize',7,'location','southwest');
    M=findobj(legend1,'type','patch');
    set(M,'MarkerSize',sqrt(10))
end
end
saveas (gcf,'figure1','fig');
close

% Figure 2 compares observed target tariffs with model-implied MFN optimal
% tariffs under lobbying. Sorting by tariff level highlights how the model
% tracks the data across sectors.
%Constructing Figure 2: Optimal tariffs with lobbying
for j=1:N
TEMP=sortrows(100*[TARGETTARIFF(:,j),MFNOPTIMALTARIFFPOL(:,j)],1);
TEMP=[TEMP,[1:1:S]'];
if j<=3
    subplot(4,3,6+j)
    position=get(gca,'position');
    position=position+[0.01 0 -0.02 0];
    subplot('position',position)
    scatter(TEMP(:,3),TEMP(:,1),10,'k','*')
    hold on
    scatter(TEMP(:,3),TEMP(:,2),10,'k','o')
    hold on
    Names={'Brazil','China','EU','India','Japan','RoW','US'};
    axis([1 S 0 250])
    set(gca,'fontsize',7,'xtick',1:4:S);
    title(['Optimal tariffs ' Names{j}])
    xlabel('Industry rank (lowest tariff to highest tariff)');
    ylabel('Optimal tariff in %');
elseif j>3 && j<7
    subplot(4,3,6+j)
    position=get(gca,'position');
    position=position+[0.01 -0.04 -0.02 0];
    subplot('position',position)
    scatter(TEMP(:,3),TEMP(:,1),10,'k','*')
    hold on
    scatter(TEMP(:,3),TEMP(:,2),10,'k','o')
    hold on
    Names={'Brazil','China','EU','India','Japan','RoW','US'};
    if j==5
        axis([1 S 0 800])
    elseif j~=5
        axis([1 S 0 250])
    end
    set(gca,'fontsize',7,'xtick',1:4:S);
    title(['Optimal tariffs ' Names{j}])
    xlabel('Industry rank (lowest tariff to highest tariff)');
    ylabel('Optimal tariff in %');
elseif j==7
    subplot(4,3,[1 6])
    position=get(gca,'position');
    position=position+[0.2 0.04 -0.4 0];
    subplot('position',position)
    scatter(TEMP(:,3),TEMP(:,1),10,'k','*')
    hold on
    scatter(TEMP(:,3),TEMP(:,2),10,'k','o')
    hold on
    Names={'Brazil','China','EU','India','Japan','RoW','US'};
    axis([1 S 0 250])
    set(gca,'fontsize',7,'xtick',1:4:S);
    title(['Optimal tariffs ' Names{j}])
    xlabel('Industry rank (lowest tariff to highest tariff)');
    ylabel('Optimal tariff in %');
    legend1=legend({'Data' 'Model'});
    set(legend1,'FontSize',7,'location','northwest');
    M=findobj(legend1,'type','patch');
    set(M,'MarkerSize',sqrt(10))   
end
end
saveas (gcf,'figure2','fig');
close

% Figure 3 visualizes Nash equilibrium tariffs without lobbying. Bilateral
% Nash tariffs are reshaped to N×S matrices and sorted by elasticity.
%Constructing Figure 3: Nash tariffs without lobbying
load NASHTARIFFBASs
for j=1:N
TEMP=reshape(NASHTARIFFBASs(:,j,:),N,S);   % N partners × S sectors for importer j
TEMP=[100*TEMP;SIGMA'];
TEMP=sortrows(TEMP',8)';
TEMP=[TEMP;1:1:S];
if j<=3
    subplot(4,3,6+j)
    position=get(gca,'position');
    position=position+[0.01 0 -0.02 0];
    subplot('position',position)
    if j~=1
    scatter(TEMP(9,:),TEMP(1,:),10,'k','o')
    hold on
    end
    if j~=2
    scatter(TEMP(9,:),TEMP(2,:),10,'k','+')
    hold on
    end
    if j~=3
    scatter(TEMP(9,:),TEMP(3,:),10,'k','*')
    hold on
    end
    if j~=4
    scatter(TEMP(9,:),TEMP(4,:),10,'k','x')
    hold on
    end
    if j~=5
    scatter(TEMP(9,:),TEMP(5,:),10,'k','s')
    hold on
    end
    if j~=6
    scatter(TEMP(9,:),TEMP(6,:),10,'k','d')
    hold on
    end
    if j~=7
    scatter(TEMP(9,:),TEMP(7,:),10,'k','^')
    end
    Names={'Brazil','China','EU','India','Japan','RoW','US'};
    axis([1 S 0 100])
    set(gca,'fontsize',7,'xtick',1:4:S);
    title(['Nash tariffs ' Names{j}])
    xlabel('Industry rank (lowest sigma to highest sigma)');
    ylabel('Nash tariff in %');
elseif j>3 && j<7
    subplot(4,3,6+j)
    position=get(gca,'position');
    position=position+[0.01 -0.04 -0.02 0];
    subplot('position',position)
    if j~=1
    scatter(TEMP(9,:),TEMP(1,:),10,'k','o')
    hold on
    end
    if j~=2
    scatter(TEMP(9,:),TEMP(2,:),10,'k','+')
    hold on
    end
    if j~=3
    scatter(TEMP(9,:),TEMP(3,:),10,'k','*')
    hold on
    end
    if j~=4
    scatter(TEMP(9,:),TEMP(4,:),10,'k','x')
    hold on
    end
    if j~=5
    scatter(TEMP(9,:),TEMP(5,:),10,'k','s')
    hold on
    end
    if j~=6
    scatter(TEMP(9,:),TEMP(6,:),10,'k','d')
    hold on
    end
    if j~=7
    scatter(TEMP(9,:),TEMP(7,:),10,'k','^')
    end
    Names={'Brazil','China','EU','India','Japan','ROW','US'};
    axis([1 S 0 100])
    set(gca,'fontsize',7,'xtick',1:4:S);
    title(['Nash tariffs ' Names{j}])
    xlabel('Industry rank (lowest sigma to highest sigma)');
    ylabel('Nash tariff in %');
elseif j==7
    subplot(4,3,[1 6])
    position=get(gca,'position');
    position=position+[0.2 0.04 -0.4 0];
    subplot('position',position)
    if j~=1
    scatter(TEMP(9,:),TEMP(1,:),10,'k','o')
    hold on
    end
    if j~=2
    scatter(TEMP(9,:),TEMP(2,:),10,'k','+')
    hold on
    end
    if j~=3
    scatter(TEMP(9,:),TEMP(3,:),10,'k','*')
    hold on
    end
    if j~=4
    scatter(TEMP(9,:),TEMP(4,:),10,'k','x')
    hold on
    end
    if j~=5
    scatter(TEMP(9,:),TEMP(5,:),10,'k','s')
    hold on
    end
    if j~=6
    scatter(TEMP(9,:),TEMP(6,:),10,'k','d')
    hold on
    end
    if j~=7
    scatter(TEMP(9,:),TEMP(7,:),10,'k','^')
    end
    Names={'Brazil','China','EU','India','Japan','RoW','US'};
    axis([1 S 0 100])
    set(gca,'fontsize',7,'xtick',1:2:S);
    title(['Nash tariffs ' Names{j}])
    xlabel('Industry rank (lowest sigma to highest sigma)');
    ylabel('Nash tariff in %');
    Names(j)=[];
    legend1=legend(Names);
    set(legend1,'FontSize',7,'location','southwest');
    M=findobj(legend1,'type','patch');
    set(M,'MarkerSize',sqrt(10))
end
end
saveas (gcf,'figure3','fig');
close

% Figure 4 depicts Nash equilibrium tariffs when lobbying is present. We
% adjust Japan's column to match bounds used in Figure 2, derive MFN
% averages, and compare them with observed target tariffs.
%Constructing Figure 4: Nash tariffs with lobbying
load MFNNASHTARIFFPOLs
MFNNASHTARIFFPOLCLEANs=MFNNASHTARIFFPOLs;
j=5;
UB=max(TARGETTARIFF(:,j)+0.03,2.25);
MFNOPTIMALTARIFFj=mymfnoptimaltariffj(j,MFNNASHTARIFFPOLs,LAMBDAPOL,0,UB,TARGETTARIFF(:,j)); %Recomputing for Japan with same UB as in Figure 2
TEMP=repmat(reshape(MFNOPTIMALTARIFFj,[1 1 S]),[N-1,1,1]); % replicate S×1 vector across exporters
MFNNASHTARIFFPOLCLEANs(:,j,:)=[TEMP(1:j-1,:,:);zeros(1,1,S);TEMP(j:end,:,:)];
save('MFNNASHTARIFFPOLCLEANs','MFNNASHTARIFFPOLCLEANs')

MFNNASHTARIFFPOLCLEAN=zeros(S,N);
for j=1:N
    TEMP=MFNNASHTARIFFPOLCLEANs(:,j,:);
    TEMP(j,:,:)=[];
    MFNNASHTARIFFPOLCLEAN(:,j)=reshape(mean(TEMP,1),S,1); % average across exporters -> S×1
end

for j=1:N
TEMP=sortrows(100*[TARGETTARIFF(:,j),MFNNASHTARIFFPOLCLEAN(:,j)],1);
TEMP=[TEMP,[1:1:S]'];
if j<=3
    subplot(4,3,6+j)
    position=get(gca,'position');
    position=position+[0.01 0 -0.02 0];
    subplot('position',position)
    scatter(TEMP(:,3),TEMP(:,1),10,'k','*')
    hold on
    scatter(TEMP(:,3),TEMP(:,2),10,'k','o')
    hold on
    Names={'Brazil','China','EU','India','Japan','RoW','US'};
    axis([1 S 0 250])
    set(gca,'fontsize',7,'xtick',1:4:S);
    title(['Nash tariffs ' Names{j}])
    xlabel('Industry rank (lowest tariff to highest tariff)');
    ylabel('Nash tariff in %');
elseif j>3 && j<7
    subplot(4,3,6+j)
    position=get(gca,'position');
    position=position+[0.01 -0.04 -0.02 0];
    subplot('position',position)
    scatter(TEMP(:,3),TEMP(:,1),10,'k','*')
    hold on
    scatter(TEMP(:,3),TEMP(:,2),10,'k','o')
    hold on
    Names={'Brazil','China','EU','India','Japan','RoW','US'};
    if j==5
        axis([1 S 0 800])
    elseif j~=5
        axis([1 S 0 250])
    end
    set(gca,'fontsize',7,'xtick',1:4:S);
    title(['Nash tariffs ' Names{j}])
    xlabel('Industry rank (lowest tariff to highest tariff)');
    ylabel('Nash tariff in %');
elseif j==7
    subplot(4,3,[1 6])
    position=get(gca,'position');
    position=position+[0.2 0.04 -0.4 0];
    subplot('position',position)
    scatter(TEMP(:,3),TEMP(:,1),10,'k','*')
    hold on
    scatter(TEMP(:,3),TEMP(:,2),10,'k','o')
    hold on
    Names={'Brazil','China','EU','India','Japan','RoW','US'};
    axis([1 S 0 250])
    set(gca,'fontsize',7,'xtick',1:4:S);
    title(['Nash tariffs ' Names{j}])
    xlabel('Industry rank (lowest tariff to highest tariff)');
    ylabel('Nash tariff in %');
    legend1=legend({'Data' 'Model'});
    set(legend1,'FontSize',7,'location','northwest');
    M=findobj(legend1,'type','patch');
    set(M,'MarkerSize',sqrt(10))   
end
end
saveas (gcf,'figure4','fig');
close

% Figure 5 examines cooperative tariffs when talks begin from the Nash
% equilibrium. Bilateral cooperative tariffs are reshaped and ranked by
% elasticity to illustrate negotiated outcomes.
%Constructing Figure 5: Cooperative tariffs starting at Nash tariffs
cd(origDir)
copyfile(fullfile(projectRoot,'Results','Trade talks',vers,'Nash_Bas','*.mat'),figDir)
cd(figDir)
load UNRESTRICTEDCOOPERATIVETARIFFBASs
COOPERATIVETARIFFs=UNRESTRICTEDCOOPERATIVETARIFFBASs;

for j=1:N
TEMP=reshape(COOPERATIVETARIFFs(:,j,:),N,S); % N partners × S sectors
TEMP=[100*TEMP;SIGMA'];
TEMP=sortrows(TEMP',8)';
TEMP=[TEMP;1:1:S];
if j<=3
    subplot(4,3,6+j)
    position=get(gca,'position');
    position=position+[0.01 0 -0.02 0];
    subplot('position',position)
    if j~=1
    scatter(TEMP(9,:),TEMP(1,:),10,'k','o')
    hold on
    end
    if j~=2
    scatter(TEMP(9,:),TEMP(2,:),10,'k','+')
    hold on
    end
    if j~=3
    scatter(TEMP(9,:),TEMP(3,:),10,'k','*')
    hold on
    end
    if j~=4
    scatter(TEMP(9,:),TEMP(4,:),10,'k','x')
    hold on
    end
    if j~=5
    scatter(TEMP(9,:),TEMP(5,:),10,'k','s')
    hold on
    end
    if j~=6
    scatter(TEMP(9,:),TEMP(6,:),10,'k','d')
    hold on
    end
    if j~=7
    scatter(TEMP(9,:),TEMP(7,:),10,'k','^')
    end
    Names={'Brazil','China','EU','India','Japan','RoW','US'};
    axis([1 S -50 50])
    set(gca,'fontsize',7,'xtick',1:4:S);
    title(['Cooperative tariffs ' Names{j}])
    xlabel('Industry rank (lowest sigma to highest sigma)');
    ylabel('Cooperative tariff in %');
elseif j>3 && j<7
    subplot(4,3,6+j)
    position=get(gca,'position');
    position=position+[0.01 -0.04 -0.02 0];
    subplot('position',position)
    if j~=1
    scatter(TEMP(9,:),TEMP(1,:),10,'k','o')
    hold on
    end
    if j~=2
    scatter(TEMP(9,:),TEMP(2,:),10,'k','+')
    hold on
    end
    if j~=3
    scatter(TEMP(9,:),TEMP(3,:),10,'k','*')
    hold on
    end
    if j~=4
    scatter(TEMP(9,:),TEMP(4,:),10,'k','x')
    hold on
    end
    if j~=5
    scatter(TEMP(9,:),TEMP(5,:),10,'k','s')
    hold on
    end
    if j~=6
    scatter(TEMP(9,:),TEMP(6,:),10,'k','d')
    hold on
    end
    if j~=7
    scatter(TEMP(9,:),TEMP(7,:),10,'k','^')
    end
    Names={'Brazil','China','EU','India','Japan','ROW','US'};
    axis([1 S -50 50])
    set(gca,'fontsize',7,'xtick',1:4:S);
    title(['Cooperative tariffs ' Names{j}])
    xlabel('Industry rank (lowest sigma to highest sigma)');
    ylabel('Cooperative tariff in %');
elseif j==7
    subplot(4,3,[1 6])
    position=get(gca,'position');
    position=position+[0.2 0.04 -0.4 0];
    subplot('position',position)
    if j~=1
    scatter(TEMP(9,:),TEMP(1,:),10,'k','o')
    hold on
    end
    if j~=2
    scatter(TEMP(9,:),TEMP(2,:),10,'k','+')
    hold on
    end
    if j~=3
    scatter(TEMP(9,:),TEMP(3,:),10,'k','*')
    hold on
    end
    if j~=4
    scatter(TEMP(9,:),TEMP(4,:),10,'k','x')
    hold on
    end
    if j~=5
    scatter(TEMP(9,:),TEMP(5,:),10,'k','s')
    hold on
    end
    if j~=6
    scatter(TEMP(9,:),TEMP(6,:),10,'k','d')
    hold on
    end
    if j~=7
    scatter(TEMP(9,:),TEMP(7,:),10,'k','^')
    end
    Names={'Brazil','China','EU','India','Japan','RoW','US'};
    axis([1 S -50 50])
    set(gca,'fontsize',7,'xtick',1:2:S);
    title(['Cooperative tariffs ' Names{j}])
    xlabel('Industry rank (lowest sigma to highest sigma)');
    ylabel('Cooperative tariff in %');
    Names(j)=[];
    legend1=legend(Names);
    set(legend1,'FontSize',7,'location','southeast');
    M=findobj(legend1,'type','patch');
    set(M,'MarkerSize',sqrt(10))
end
end
saveas (gcf,'figure5','fig');
close

% Figure 6 repeats the cooperative-tariff exercise but starts from factual
% tariffs. Results from the 'Fact' scenario are reshaped and ranked by
% elasticity for comparison.
%Constructing Figure 6: Cooperative tariffs starting at factual tariffs

cd(origDir)
copyfile(fullfile(projectRoot,'Results','Trade talks',vers,'Fact','*.mat'),figDir)
cd(figDir)
load UNRESTRICTEDCOOPERATIVETARIFFBASs
COOPERATIVETARIFFs=UNRESTRICTEDCOOPERATIVETARIFFBASs;

for j=1:N
TEMP=reshape(COOPERATIVETARIFFs(:,j,:),N,S); % N partners × S sectors
TEMP=[100*TEMP;SIGMA'];
TEMP=sortrows(TEMP',8)';
TEMP=[TEMP;1:1:S];
if j<=3
    subplot(4,3,6+j)
    position=get(gca,'position');
    position=position+[0.01 0 -0.02 0];
    subplot('position',position)
    if j~=1
    scatter(TEMP(9,:),TEMP(1,:),10,'k','o')
    hold on
    end
    if j~=2
    scatter(TEMP(9,:),TEMP(2,:),10,'k','+')
    hold on
    end
    if j~=3
    scatter(TEMP(9,:),TEMP(3,:),10,'k','*')
    hold on
    end
    if j~=4
    scatter(TEMP(9,:),TEMP(4,:),10,'k','x')
    hold on
    end
    if j~=5
    scatter(TEMP(9,:),TEMP(5,:),10,'k','s')
    hold on
    end
    if j~=6
    scatter(TEMP(9,:),TEMP(6,:),10,'k','d')
    hold on
    end
    if j~=7
    scatter(TEMP(9,:),TEMP(7,:),10,'k','^')
    end
    Names={'Brazil','China','EU','India','Japan','RoW','US'};
    axis([1 S -50 50])
    set(gca,'fontsize',7,'xtick',1:4:S);
    title(['Cooperative tariffs ' Names{j}])
    xlabel('Industry rank (lowest sigma to highest sigma)');
    ylabel('Cooperative tariff in %');
elseif j>3 && j<7
    subplot(4,3,6+j)
    position=get(gca,'position');
    position=position+[0.01 -0.04 -0.02 0];
    subplot('position',position)
    if j~=1
    scatter(TEMP(9,:),TEMP(1,:),10,'k','o')
    hold on
    end
    if j~=2
    scatter(TEMP(9,:),TEMP(2,:),10,'k','+')
    hold on
    end
    if j~=3
    scatter(TEMP(9,:),TEMP(3,:),10,'k','*')
    hold on
    end
    if j~=4
    scatter(TEMP(9,:),TEMP(4,:),10,'k','x')
    hold on
    end
    if j~=5
    scatter(TEMP(9,:),TEMP(5,:),10,'k','s')
    hold on
    end
    if j~=6
    scatter(TEMP(9,:),TEMP(6,:),10,'k','d')
    hold on
    end
    if j~=7
    scatter(TEMP(9,:),TEMP(7,:),10,'k','^')
    end
    Names={'Brazil','China','EU','India','Japan','ROW','US'};
    axis([1 S -50 50])
    set(gca,'fontsize',7,'xtick',1:4:S);
    title(['Cooperative tariffs ' Names{j}])
    xlabel('Industry rank (lowest sigma to highest sigma)');
    ylabel('Cooperative tariff in %');
elseif j==7
    subplot(4,3,[1 6])
    position=get(gca,'position');
    position=position+[0.2 0.04 -0.4 0];
    subplot('position',position)
    if j~=1
    scatter(TEMP(9,:),TEMP(1,:),10,'k','o')
    hold on
    end
    if j~=2
    scatter(TEMP(9,:),TEMP(2,:),10,'k','+')
    hold on
    end
    if j~=3
    scatter(TEMP(9,:),TEMP(3,:),10,'k','*')
    hold on
    end
    if j~=4
    scatter(TEMP(9,:),TEMP(4,:),10,'k','x')
    hold on
    end
    if j~=5
    scatter(TEMP(9,:),TEMP(5,:),10,'k','s')
    hold on
    end
    if j~=6
    scatter(TEMP(9,:),TEMP(6,:),10,'k','d')
    hold on
    end
    if j~=7
    scatter(TEMP(9,:),TEMP(7,:),10,'k','^')
    end
    Names={'Brazil','China','EU','India','Japan','RoW','US'};
    axis([1 S -50 50])
    set(gca,'fontsize',7,'xtick',1:2:S);
    title(['Cooperative tariffs ' Names{j}])
    xlabel('Industry rank (lowest sigma to highest sigma)');
    ylabel('Cooperative tariff in %');
    Names(j)=[];
    legend1=legend(Names);
    set(legend1,'FontSize',7,'location','southeast');
    M=findobj(legend1,'type','patch');
    set(M,'MarkerSize',sqrt(10))
end
end
saveas (gcf,'figure6','fig');
close

% Figure 7 shows cooperative tariffs when negotiations begin from free
% trade. We load the "Free" scenario results, reshape the tariff cube and
% rank sectors by elasticity to show outcomes from a zero‐tariff baseline.
%Constructing Figure 7: Cooperative tariffs starting at free trade

cd(origDir)
copyfile(fullfile(projectRoot,'Results','Trade talks',vers,'Free','*.mat'),figDir)
cd(figDir)
load UNRESTRICTEDCOOPERATIVETARIFFBASs
COOPERATIVETARIFFs=UNRESTRICTEDCOOPERATIVETARIFFBASs;

for j=1:N
TEMP=reshape(COOPERATIVETARIFFs(:,j,:),N,S); % N partners × S sectors
TEMP=[100*TEMP;SIGMA'];
TEMP=sortrows(TEMP',8)';
TEMP=[TEMP;1:1:S];
if j<=3
    subplot(4,3,6+j)
    position=get(gca,'position');
    position=position+[0.01 0 -0.02 0];
    subplot('position',position)
    if j~=1
    scatter(TEMP(9,:),TEMP(1,:),10,'k','o')
    hold on
    end
    if j~=2
    scatter(TEMP(9,:),TEMP(2,:),10,'k','+')
    hold on
    end
    if j~=3
    scatter(TEMP(9,:),TEMP(3,:),10,'k','*')
    hold on
    end
    if j~=4
    scatter(TEMP(9,:),TEMP(4,:),10,'k','x')
    hold on
    end
    if j~=5
    scatter(TEMP(9,:),TEMP(5,:),10,'k','s')
    hold on
    end
    if j~=6
    scatter(TEMP(9,:),TEMP(6,:),10,'k','d')
    hold on
    end
    if j~=7
    scatter(TEMP(9,:),TEMP(7,:),10,'k','^')
    end
    Names={'Brazil','China','EU','India','Japan','RoW','US'};
    axis([1 S -50 50])
    set(gca,'fontsize',7,'xtick',1:4:S);
    title(['Cooperative tariffs ' Names{j}])
    xlabel('Industry rank (lowest sigma to highest sigma)');
    ylabel('Cooperative tariff in %');
elseif j>3 && j<7
    subplot(4,3,6+j)
    position=get(gca,'position');
    position=position+[0.01 -0.04 -0.02 0];
    subplot('position',position)
    if j~=1
    scatter(TEMP(9,:),TEMP(1,:),10,'k','o')
    hold on
    end
    if j~=2
    scatter(TEMP(9,:),TEMP(2,:),10,'k','+')
    hold on
    end
    if j~=3
    scatter(TEMP(9,:),TEMP(3,:),10,'k','*')
    hold on
    end
    if j~=4
    scatter(TEMP(9,:),TEMP(4,:),10,'k','x')
    hold on
    end
    if j~=5
    scatter(TEMP(9,:),TEMP(5,:),10,'k','s')
    hold on
    end
    if j~=6
    scatter(TEMP(9,:),TEMP(6,:),10,'k','d')
    hold on
    end
    if j~=7
    scatter(TEMP(9,:),TEMP(7,:),10,'k','^')
    end
    Names={'Brazil','China','EU','India','Japan','ROW','US'};
    axis([1 S -50 50])
    set(gca,'fontsize',7,'xtick',1:4:S);
    title(['Cooperative tariffs ' Names{j}])
    xlabel('Industry rank (lowest sigma to highest sigma)');
    ylabel('Cooperative tariff in %');
elseif j==7
    subplot(4,3,[1 6])
    position=get(gca,'position');
    position=position+[0.2 0.04 -0.4 0];
    subplot('position',position)
    if j~=1
    scatter(TEMP(9,:),TEMP(1,:),10,'k','o')
    hold on
    end
    if j~=2
    scatter(TEMP(9,:),TEMP(2,:),10,'k','+')
    hold on
    end
    if j~=3
    scatter(TEMP(9,:),TEMP(3,:),10,'k','*')
    hold on
    end
    if j~=4
    scatter(TEMP(9,:),TEMP(4,:),10,'k','x')
    hold on
    end
    if j~=5
    scatter(TEMP(9,:),TEMP(5,:),10,'k','s')
    hold on
    end
    if j~=6
    scatter(TEMP(9,:),TEMP(6,:),10,'k','d')
    hold on
    end
    if j~=7
    scatter(TEMP(9,:),TEMP(7,:),10,'k','^')
    end
    Names={'Brazil','China','EU','India','Japan','RoW','US'};
    axis([1 S -50 50])
    set(gca,'fontsize',7,'xtick',1:2:S);
    title(['Cooperative tariffs ' Names{j}])
    xlabel('Industry rank (lowest sigma to highest sigma)');
    ylabel('Cooperative tariff in %');
    Names(j)=[];
    legend1=legend(Names);
    set(legend1,'FontSize',7,'location','southeast');
    M=findobj(legend1,'type','patch');
    set(M,'MarkerSize',sqrt(10))
end
end
saveas (gcf,'figure7','fig');
close

% Figure 8 explores liberalization starting from MFN Nash tariffs. We load the Nash baseline, convert trade and tariff cubes to (N*S)×N matrices for the GE solver, then simulate tariff cuts among major economies.
%Constructing Figure 8: Liberalization scenarios starting at MFN Nash tariffs
load MFNNASHTARIFFBASs
LAMBDA=LAMBDABAS;
NASHTARIFFs=MFNNASHTARIFFBASs;

[GOVERNMENTWELFAREHAT WELFAREHAT WAGEHAT TRADECs LOBBYWELFAREHAT EXPENDITUREHAT]=mycounterfactuals(NASHTARIFFs,zeros(N,1),LAMBDA);
TRADE=reshape(permute(TRADECs,[1 3 2]),[N*S,N,1]); % trade flows cube -> (N*S)×N
TARIFF=reshape(permute(NASHTARIFFs,[1 3 2]),[N*S,N,1]);  % Nash tariffs cube -> (N*S)×N
save('DATA','SIGMA','TARIFF','TRADE','LAMBDAPOL','LAMBDABAS')
mycalculations

GSUM=[];
WSUM=[];
WAGESUM=[];
PROFITSSUM=[];
for r=1:-0.05:0 %Simulating non-MFN tariff reduction between US, EU, and Japan starting at MFN Nash tariffs
    TARIFFCs=TARIFFs;
    for j=[3 5 7]
        for i=[3 5 7]
            TARIFFCs(i,j,:)=r*TARIFFCs(i,j,:);
        end
    end
    [GOVERNMENTWELFAREHAT WELFAREHAT WAGEHAT TRADECs LOBBYWELFAREHAT EXPENDITUREHAT]=mycounterfactuals(TARIFFCs,zeros(N,1),LAMBDA);
    GSUM=[GSUM;[100*(r-1) 100*(GOVERNMENTWELFAREHAT'-1)]];
    WSUM=[WSUM;[100*(r-1) 100*(WELFAREHAT'-1)]];
    PROFITS=sum((1./SIGMAs_N_1_S).*sum(TRADEs,2),3);
    PROFITSC=sum((1./SIGMAs_N_1_S).*sum(TRADECs,2),3);
    PROFITSHAT=PROFITSC./PROFITS;
    PROFITSHATADJ=PROFITSHAT./WAGEHAT;
    PROFITSSUM=[PROFITSSUM;[100*(r-1) 100*(PROFITSHATADJ'-1)]];
    WAGEHAT=WAGEHAT-mean(WAGEHAT)+1;
    WAGESUM=[WAGESUM;[100*(r-1) 100*(WAGEHAT'-1)]];   
end

GSUM=[GSUM(:,1) mean([GSUM(:,4) GSUM(:,6) GSUM(:,8)],2) mean([GSUM(:,2) GSUM(:,3) GSUM(:,5)  GSUM(:,7)],2)];
PROFITSSUM=[PROFITSSUM(:,1) mean([PROFITSSUM(:,4) PROFITSSUM(:,6) PROFITSSUM(:,8)],2) mean([PROFITSSUM(:,2) PROFITSSUM(:,3) PROFITSSUM(:,5)  PROFITSSUM(:,7)],2)];
WAGESUM=[WAGESUM(:,1) mean([WAGESUM(:,4) WAGESUM(:,6) WAGESUM(:,8)],2) mean([WAGESUM(:,2) WAGESUM(:,3) WAGESUM(:,5)  WAGESUM(:,7)],2)];

subplot(2,2,1)
plot(-GSUM(:,1),GSUM(:,2),'kx','MarkerSize',sqrt(10))
hold on
plot(-GSUM(:,1),GSUM(:,3),'ko','MarkerSize',sqrt(10))
axis([0 100 -0.5 1.5])
set(gca,'fontsize',7,'xtick',0:20:100);
title('Non-MFN liberalization among EU, Japan, and US')
xlabel('Tariff cut relative to Nash in %');
ylabel('Average welfare change in %');
legend1=legend({'Liberalizing countries' 'Other countries'});
set(legend1,'FontSize',7,'location','northwest');

subplot(2,2,2)
plot(-WAGESUM(:,1),WAGESUM(:,2),'kx','MarkerSize',sqrt(10))
hold on
plot(-WAGESUM(:,1),WAGESUM(:,3),'ko','MarkerSize',sqrt(10))
axis([0 100 -1 1.5])
set(gca,'fontsize',7,'xtick',0:20:100);
title('Non-MFN liberalization among EU, Japan, and US')
xlabel('Tariff cut relative to Nash in %');
ylabel('Average wage change in %');
legend1=legend({'Liberalizing countries' 'Other countries'});
set(legend1,'FontSize',7,'location','northwest');

GSUM=[];
WSUM=[];
WAGESUM=[];
PROFITSSUM=[];
for r=1:-0.05:0 %Simulating MFN tariff reduction between US, EU, and Japan starting at MFN Nash tariffs
    TARIFFCs=TARIFFs;
    for j=[3 5 7]
        TARIFFCs(:,j,:)=r*TARIFFCs(:,j,:);
    end
    [GOVERNMENTWELFAREHAT WELFAREHAT WAGEHAT TRADECs LOBBYWELFAREHAT EXPENDITUREHAT]=mycounterfactuals(TARIFFCs,zeros(N,1),LAMBDA);
    GSUM=[GSUM;[100*(r-1) 100*(GOVERNMENTWELFAREHAT'-1)]];
    WSUM=[WSUM;[100*(r-1) 100*(WELFAREHAT'-1)]];
    PROFITS=sum((1./SIGMAs_N_1_S).*sum(TRADEs,2),3);
    PROFITSC=sum((1./SIGMAs_N_1_S).*sum(TRADECs,2),3);
    PROFITSHAT=PROFITSC./PROFITS;
    PROFITSHATADJ=PROFITSHAT./WAGEHAT;
    PROFITSSUM=[PROFITSSUM;[100*(r-1) 100*(PROFITSHATADJ'-1)]];
    WAGEHAT=WAGEHAT-mean(WAGEHAT)+1;
    WAGESUM=[WAGESUM;[100*(r-1) 100*(WAGEHAT'-1)]];  
end

GSUM=[GSUM(:,1) mean([GSUM(:,4) GSUM(:,6) GSUM(:,8)],2) mean([GSUM(:,2) GSUM(:,3) GSUM(:,5)  GSUM(:,7)],2)];
PROFITSSUM=[PROFITSSUM(:,1) mean([PROFITSSUM(:,4) PROFITSSUM(:,6) PROFITSSUM(:,8)],2) mean([PROFITSSUM(:,2) PROFITSSUM(:,3) PROFITSSUM(:,5)  PROFITSSUM(:,7)],2)];
WAGESUM=[WAGESUM(:,1) mean([WAGESUM(:,4) WAGESUM(:,6) WAGESUM(:,8)],2) mean([WAGESUM(:,2) WAGESUM(:,3) WAGESUM(:,5)  WAGESUM(:,7)],2)];

subplot(2,2,3)
plot(-GSUM(:,1),GSUM(:,2),'kx','MarkerSize',sqrt(10))
hold on
plot(-GSUM(:,1),GSUM(:,3),'ko','MarkerSize',sqrt(10))
axis([0 100 -1 4])
set(gca,'fontsize',7,'xtick',0:20:100);
title('MFN liberalization among EU, Japan, and US')
xlabel('Tariff cut relative to Nash in %');
ylabel('Average welfare change in %');
legend1=legend({'Liberalizing countries' 'Other countries'});
set(legend1,'FontSize',7,'location','northwest');

subplot(2,2,4)
plot(-WAGESUM(:,1),WAGESUM(:,2),'kx','MarkerSize',sqrt(10))
hold on
plot(-WAGESUM(:,1),WAGESUM(:,3),'ko','MarkerSize',sqrt(10))
axis([0 100 -15 15])
set(gca,'fontsize',7,'xtick',0:20:100);
title('MFN liberalization among EU, Japan, and US')
xlabel('Tariff cut relative to Nash in %');
ylabel('Average wage change in %');
legend1=legend({'Liberalizing countries' 'Other countries'});
set(legend1,'FontSize',7,'location','northwest');

saveas (gcf,'figure8','fig');
close

%Cleaning up
% Remove the temporary data (.mat) and code (.m) files copied to the
% Figures/ directory above. They were needed only for plot generation and
% are deleted to avoid leaving redundant artifacts.
delete(fullfile(figDir,'*.mat'))  % temporary data files
delete(fullfile(figDir,'*.m'))    % temporary helper scripts
cd(origDir)

clear all
close all
clc

%This is checked and correct
