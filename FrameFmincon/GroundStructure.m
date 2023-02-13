% Matlab function to create ground structure with diamond cell
% Written by Dhanashri Desai 

function [BoundaryNodeID,CentralNodeID,nnodes,nelem,nx,ny,elemcon,L] = ...
GroundStructure(nCellH,nCellV,Nrows,Ncolumns,LCell,LCellH,LCellV)

nelem = 0; nnodes = 0;

for i = 1:Nrows % Co-ordinates of boundary nodes of square unit cell
    for j = 1:Ncolumns
        nnodes = nnodes + 1;
        BoundaryNodeID(i,j) = nnodes;
        nx(BoundaryNodeID(i,j)) = LCell * (j-1);
        ny(BoundaryNodeID(i,j)) = LCell * (i-1);
    end
end

for i = 1:nCellV % Co-ordinates of central nodes of square unit cell
    for j = 1:nCellH
        nnodes = nnodes+1;
        CentralNodeID(i,j) = nnodes;
        nx(CentralNodeID(i,j)) = (2*j-1)*LCellH;
        ny(CentralNodeID(i,j)) = (2*i-1)*LCellV;
    end
end
    
for i = 1:Nrows % Connect horizontal elements
    for j = 1:(Ncolumns-1)
        nelem = nelem + 1;
        elemcon(nelem,1) = BoundaryNodeID(i,j);
        elemcon(nelem,2) = BoundaryNodeID(i,j+1);
    end
end

for j = 1:Ncolumns % Connect vertical elements
    for i = 1:(Nrows-1)
        nelem = nelem + 1;
        elemcon(nelem,1) = BoundaryNodeID(i,j);
        elemcon(nelem,2) = BoundaryNodeID(i+1,j);
    end
end

for i=1:nCellV % Connect inclined elements
    for j=1:nCellH
        nelem = nelem + 1;
        elemcon(nelem,1) = CentralNodeID(i,j);
        elemcon(nelem,2) = BoundaryNodeID(i,j);

        nelem = nelem + 1;
        elemcon(nelem,1) = CentralNodeID(i,j);
        elemcon(nelem,2) = BoundaryNodeID(i,j+1);

        nelem = nelem + 1;
        elemcon(nelem,1) = CentralNodeID(i,j);
        elemcon(nelem,2) = BoundaryNodeID(i+1,j);

        nelem = nelem + 1;
        elemcon(nelem,1) = CentralNodeID(i,j);
        elemcon(nelem,2) = BoundaryNodeID(i+1,j+1);
    end
end

for i = 1:nelem % Calculate lengths of elements
    node1 = elemcon(i,1); node2 = elemcon(i,2);
    L(i) = sqrt((nx(node1)-nx(node2))^2 + (ny(node1)-ny(node2))^2);
end

figure(1); % Plotting
for i = 1:nelem
    node1 = elemcon(i,1);  node2 = elemcon(i,2); 
    plot([nx(node1) nx(node2)],[ny(node1) ny(node2)],'-k'); hold on; grid on;
%     plot([nx(node1) nx(node2)],[ny(node1) ny(node2)],'ko'); hold on;
end

for i = 1:nnodes
    text(nx(i),ny(i),num2str(i),'Fontsize',14); hold on;
    plot(nx(i),ny(i),'bo',MarkerSize=4); hold on;
end
axis equal;
 
% if nelem==2*nnodes-3
%     disp('Stiff frame');
% else
%     disp('Non-stiff frame');
% end

end
