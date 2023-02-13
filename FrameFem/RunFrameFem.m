close all; clear all; clc; % Housekeeping commands
global nCellH nCellV Nrows Ncolumns BoundaryNodeID CentralNodeID nnodes nelem nx ny elemcon L LCell LCellH LCellV E extloaddof1  extloaddof2 dispbcdof ExtLoadDir ExtLoad tStrut bStrut elemconUpdated tStrutUpdated
nCellH = 10; nCellV = 10; % No of cells in horizontal and vertical  direction
Nrows = nCellV+1; Ncolumns = nCellH+1; % no of rows and columns in the grid

tmin=0.003; tmax=0.006;
tStrut = (tmin+tmax)/2; bStrut = 0.002; % width (layer height) of strut in meters
RhoRel = 0.5; LCell=6.83*tStrut/RhoRel;
LCellH = LCell/2; LCellV = LCellH;

[BoundaryNodeID,CentralNodeID,nnodes,nelem,nx,ny,elemcon,L] = ...
GroundStructure(nCellH,nCellV,Nrows,Ncolumns,LCell,LCellH,LCellV); % Load data from ground structure

tStrut = tStrut*ones(nelem,1); % thickness (layer width) of strut
A = tStrut*bStrut; I = bStrut*tStrut.^3/12; % c/s and area MI of of struts
E = 150E6*ones(nelem,1); % YM of material in N/m2

%% Traction + Disp BCs

% Case1 - Compression
ExtLoad = -100; ExtLoadDir = 2; % Extrenal load and direction
FixedNodes=[1:11]'; k=1; % Displacement BC DOF
for i=1:length(FixedNodes)
    for j=1:3
        dispbcdof(k,1) = 3*(FixedNodes(i)-1) + j;  k = k+1;    
    end
end
LoadNodes1=[111:121]'; k=1;
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

[SE,gdisp,gstiffness,reactions,nodalforces,elementforces,stresses] = ... % Run frame fem code
    BernoulliFrameFem(nnodes,nelem,A,E,I,L,nx,ny,elemcon,dispbcdof,extloaddof1,extloaddof2,ExtLoad);
% disp('Displacements at Nodes in mm'); [(1:3*nnodes)'  gdisp*1000] % Display solutions
% disp('Reactions at Supports in kN'); [dispbcdof reactions/1000]
% disp('Stresses in elements in MPa'); [(1:nelem)' abs(stresses(:,1))]
maxval = max(stresses(:,1)); idx = find(stresses == maxval);
disp(['Maximum stress ' num2str(max(stresses(:,1))) ' MPa in elements:']); disp(['element' '  node1' '  node2']); disp([idx elemcon(idx,:)]);
disp(['SE of frame in Joules: ' num2str(SE)]);
disp(['Volume of the uniform frame in cm3: ' num2str(sum(tStrut*bStrut.*L')*1E6)]);
PlotConfigurations(nnodes,nelem,elemcon,nx,ny,gdisp); % Plot the results

%% Benchmarking
% Parameter = 'Load'; ParameterVector = (1:1:100);
Parameter = 'YM'; ParameterVector = 1E6:1E6:200E6;
% Parameter = 'CellWallLength'; ParameterVector = 0.1:0.001:0.5; % ParameterVector is relative density
% Parameter = 'CellWallThickness'; ParameterVector = 0.1:0.001:0.5;
% Parameter = 'CellWallWidth'; ParameterVector = 0.001:0.0001:0.005;

switch Parameter            
    case {'Load','YM'}
        [BoundaryNodeID,CentralNodeID,nnodes,nelem,nx,ny,elemcon,L] = ...
            GroundStructure(nCellH,nCellV,Nrows,Ncolumns,LCell,LCellH,LCellV); % Load data from ground structure

        if Parameter==convertCharsToStrings('Load')
            ExtLoad = ParameterVector;
            for i = 1:length(ParameterVector)
                [SE,gdisp,gstiffness,reactions,nodalforces,elementforces] = ... % Run frame fem code
                    BernoulliFrameFem(nnodes,nelem,A,E,I,L,nx,ny,elemcon,dispbcdof,extloaddof1,extloaddof2,ExtLoad(i));
                SE_FEM(i) = SE; 
                SE_Analytical(i) = nCellH*nCellV*(ExtLoad(i)^2 * ( (3*LCell/(4*A(1)*E(1))) + (LCell^3/(48*sqrt(2)*E(1)*I(1))) ));
            end
        else
            for i = 1:length(ParameterVector)
                E=ParameterVector(i)*ones(nelem,1);
                [SE,gdisp,gstiffness,reactions,nodalforces,elementforces] = ... % Run frame fem code
                    BernoulliFrameFem(nnodes,nelem,A,E,I,L,nx,ny,elemcon,dispbcdof,extloaddof1,extloaddof2,ExtLoad);
                SE_FEM(i) = SE; 
                SE_Analytical(i) = nCellH*nCellV*(ExtLoad^2 * ( (3*LCell/(4*A(1)*E(1))) + (LCell^3/(48*sqrt(2)*E(1)*I(1))) ));
            end
        end

    case {'CellWallLength','CellWallThickness','CellWallWidth'}
        if Parameter==convertCharsToStrings('CellWallLength')
            LCell=6.83*tStrut(1)./ParameterVector;
            LCellH = LCell/2; LCellV = LCellH;
            ParameterVector=LCell;
            for i = 1:length(ParameterVector)
                [BoundaryNodeID,CentralNodeID,nnodes,nelem,nx,ny,elemcon,L] = GroundStructure(nCellH,nCellV,Nrows,Ncolumns,LCell(i),LCellH(i),LCellV(i));
                [SE,gdisp,gstiffness,reactions,nodalforces,elementforces] = ... % Run frame fem code
                    BernoulliFrameFem(nnodes,nelem,A,E,I,L,nx,ny,elemcon,dispbcdof,extloaddof1,extloaddof2,ExtLoad);
                SE_FEM(i) = SE; 
                SE_Analytical(i) = nCellH*nCellV*(ExtLoad^2 * ( (3*LCell(i)/(4*A(1)*E(1))) + (LCell(i)^3/(48*sqrt(2)*E(1)*I(1))) ));
            end
        else
            if Parameter==convertCharsToStrings('CellWallThickness')
                tStrut=LCell*ParameterVector/6.83;
                ParameterVector=tStrut;
                for i = 1:length(ParameterVector)
                    tStrut = ParameterVector(i)*ones(nelem,1); % thickness (layer width) of strut
                    A = tStrut*bStrut; I = bStrut*tStrut.^3/12; % c/s and area MI of of struts
                    [BoundaryNodeID,CentralNodeID,nnodes,nelem,nx,ny,elemcon,L] = GroundStructure(nCellH,nCellV,Nrows,Ncolumns,LCell,LCellH,LCellV);
                    [SE,gdisp,gstiffness,reactions,nodalforces,elementforces] = ... % Run frame fem code
                        BernoulliFrameFem(nnodes,nelem,A,E,I,L,nx,ny,elemcon,dispbcdof,extloaddof1,extloaddof2,ExtLoad);
                    SE_FEM(i) = SE; 
                    SE_Analytical(i) = nCellH*nCellV*(ExtLoad^2 * ( (3*LCell/(4*A(1)*E(1))) + (LCell^3/(48*sqrt(2)*E(1)*I(1))) ));
                end
            else
                bStrut = ParameterVector;
                for i = 1:length(ParameterVector)
                    A = tStrut*bStrut(i); I = bStrut(i)*tStrut.^3/12; % c/s and area MI of of struts
                    [BoundaryNodeID,CentralNodeID,nnodes,nelem,nx,ny,elemcon,L] = GroundStructure(nCellH,nCellV,Nrows,Ncolumns,LCell,LCellH,LCellV);
                    [SE,gdisp,gstiffness,reactions,nodalforces,elementforces] = ... % Run frame fem code
                        BernoulliFrameFem(nnodes,nelem,A,E,I,L,nx,ny,elemcon,dispbcdof,extloaddof1,extloaddof2,ExtLoad);
                    SE_FEM(i) = SE; 
                    SE_Analytical(i) = nCellH*nCellV*(ExtLoad^2 * ( (3*LCell/(4*A(1)*E(1))) + (LCell^3/(48*sqrt(2)*E(1)*I(1))) ));
                end
            end
        end
end

figure();
plot(abs(ParameterVector),abs(SE_FEM),LineWidth=2); hold on;
plot(abs(ParameterVector),abs(SE_Analytical),LineWidth=2); hold off;
NameTheGraphModified(convertCharsToStrings(Parameter),'$$\mathrm{SE}\;\textbf{[J]}$$',[],...
            2,'$$\;\mathrm{\textbf{FEM}}$$','$$\;\mathrm{\textbf{Analytical}}$$',[],[],[],[],[],[],'northeast');


