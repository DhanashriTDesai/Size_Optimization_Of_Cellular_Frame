close all; clear all; clc; % Housekeeping commands
global nCellH nCellV Nrows Ncolumns BoundaryNodeID CentralNodeID nnodes nelem nx ny elemcon L LCell LCellH LCellV E extloaddof1 extloaddof2 dispbcdof ExtLoadDir ExtLoad tStrut bStrut elemconUpdated tStrutUpdated Linclined maxstress maxstressId

nCellH = 5; nCellV = 5; % No of cells in horizontal and vertical  direction
Nrows = nCellV+1; Ncolumns = nCellH+1; % no of rows and columns in the grid

tmin=0.003; tmax=0.006;
tStrut = (tmin+tmax)/2; bStrut = 0.002; % width (layer height) of strut in meters
RhoRel = 0.5; LCell=6.83*(tStrut)/RhoRel;
LCellH = LCell/2; LCellV = LCellH; Linclined = sqrt(LCellH^2+LCellV^2);

[BoundaryNodeID,CentralNodeID,nnodes,nelem,nx,ny,elemcon,L] = ...
GroundStructure(nCellH,nCellV,Nrows,Ncolumns,LCell,LCellH,LCellV); % Load data from ground structure

tStrut = tStrut*ones(nelem,1); % thickness (layer width) of strut
A = tStrut*bStrut; I = bStrut*tStrut.^3/12; % c/s and area MI of of struts
E = 150E6*ones(nelem,1); % YM of material in N/m2
t_lb=tmin*ones(nelem,1); t_ub=tmax*ones(nelem,1);
t0 = (tmin+tmax)/2 * ones(nelem,1);
Vmax=0.9*sum(bStrut*(tmin+tmax)/2*ones(nelem,1).*L'); % max volume in volume constraint

%% Traction + Disp BCs

% Case1 - Compression
ExtLoad = -100; ExtLoadDir = 2; % Extrenal load and direction
FixedNodes=[1:6]'; k=1; % Displacement BC DOF
for i=1:length(FixedNodes)
    for j=1:3
        dispbcdof(k,1) = 3*(FixedNodes(i)-1) + j;  k = k+1;    
    end
end
LoadNodes1=[31:36]'; k=1;
for i=1:length(LoadNodes1)
    extloaddof1(k,1)=3*(LoadNodes1(i)-1)+ExtLoadDir; k=k+1;
end


% Case2 - Cantilever
% ExtLoad = -100; ExtLoadDir = 2; % Extrenal load and direction
% FixedNodes=[1:11:111]'; k=1; % Displacement BC DOF
% for i=1:length(FixedNodes)
%     for j=1:3
%         dispbcdof(k,1) = 3*(FixedNodes(i)-1) + j;  k = k+1;    
%     end
% end
% LoadNodes1=[11]'; k=1; ExtLoadDir = 2; % External load DOF
% for i=1:length(LoadNodes1)
%     extloaddof1(k,1)=3*(LoadNodes1(i)-1)+ExtLoadDir; k=k+1;
% end


% Case3 - SSB
% ExtLoad = -100; ExtLoadDir = 2; % Extrenal load and direction
% dispbcdof = [1 2 31 32]';
% LoadNodes1=[116]'; k=1; ExtLoadDir = 2; % External load DOF
% for i=1:length(LoadNodes1)
%     extloaddof1(k,1)=3*(LoadNodes1(i)-1)+ExtLoadDir; k=k+1;
% end


% Case4 - Combined loading (left edge fixed)
% ExtLoad = -100; ExtLoadDir = 2; % Extrenal load and direction
% FixedNodes=[1:11:111]'; k=1; % Displacement BC DOF
% for i=1:length(FixedNodes)
%     for j=1:3
%         dispbcdof(k,1) = 3*(FixedNodes(i)-1) + j;  k = k+1;    
%     end
% end
% extloaddof1=3*(66-1)+1;
% extloaddof2=3*(66-1)+2;


% Pure traction BCs
% ExtLoad = -100; ExtLoadDir = 2; % Extrenal load and direction
% LoadNodes1=[1:6:31]'; k=1; % External load DOF
% for i=1:length(LoadNodes1)
%     extloaddof1(k,1)=3*(LoadNodes1(i)-1)+ExtLoadDir; k=k+1;
% end
% LoadNodes2=[6:6:36]'; k=1; % External load DOF
% for i=1:length(LoadNodes2)
%     extloaddof2(k,1)=3*(LoadNodes2(i)-1)+ExtLoadDir; k=k+1;
% end

%% Fmincon algorithm
options = optimoptions('fmincon','Algorithm','interior-point','EnableFeasibilityMode',true,'Display',...
    'iter','MaxIter',1000,'MaxFunEvals',10000,'GradObj','on','SpecifyObjectiveGradient',...
    true,'SpecifyConstraintGradient',true,'DerivativeCheck','off','PlotFcns',@optimplotfval,...
    'TolX',1e-10,'TolFun',1e-10);

[Thickness,fval,exitflag,output] = fmincon(@Frame_obj_grad,t0,L*bStrut,Vmax,[],[],t_lb,t_ub,[],options);

j=1;
for i=1:nelem
    if Thickness(i)>=tmin
        elemconUpdated(j,:)=elemcon(i,:);
        tStrutUpdated(j)=Thickness(i); j=j+1;
    end
end
tStrutUpdated=tStrutUpdated';

%% Plotting
figure(1); hold on;
for i = 1:size(elemconUpdated,1)
    node1 = elemconUpdated(i,1);  node2 = elemconUpdated(i,2); % nodes of element
    if tStrutUpdated(i)>=tmin && tStrutUpdated(i)<0.004
        plot([nx(node1) nx(node2)], [ny(node1) ny(node2)],'-b','linewidth',tStrutUpdated(i)*300); hold on;
    elseif tStrutUpdated(i)>=0.004 && tStrutUpdated(i)<0.005
        plot([nx(node1) nx(node2)], [ny(node1) ny(node2)],'-g','linewidth',tStrutUpdated(i)*300); hold on;
    else
        plot([nx(node1) nx(node2)], [ny(node1) ny(node2)],'-r','linewidth',tStrutUpdated(i)*300); hold on;
    end
end
disp(['Volume of the size optimized frame in cm3: ' num2str(sum(bStrut*tStrutUpdated.*L')*1E6)]);
axis equal; hold off;

% close all;
