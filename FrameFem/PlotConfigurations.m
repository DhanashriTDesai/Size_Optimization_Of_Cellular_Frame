% Matlab function to plot the output
% Written by Dhanashri Desai 

function PlotConfigurations(nnodes,nelem,elemcon,nx,ny,gdisp)
figure(1); hold on; grid on;
% plotting nodal displacements
for i = 1:nelem
    node1 = elemcon(i,1);  node2 = elemcon(i,2); % nodes of element
    plot([nx(node1) nx(node2)],[ny(node1) ny(node2)],'-k'); hold on; % plot line for undeformed element
    scalefactor = 1;
    dx1 = gdisp(3*node1-2)*scalefactor; dx2 = gdisp(3*node2-2)*scalefactor;
    dy1 = gdisp(3*node1-1)*scalefactor; dy2 = gdisp(3*node2-1)*scalefactor; 
    plot([nx(node1)+dx1 nx(node2)+dx2],[ny(node1)+dy1 ny(node2)+dy2],'-r'); hold on; % plot line for deformed element
    plot([nx(node1) nx(node2)],[ny(node1) ny(node2)],'ko'); % plot cicle for node at undeformed position
    plot([nx(node1)+dx1 nx(node2)+dx2],[ny(node1)+dy1 ny(node2)+dy2],'ro'); % plot cicle for node at deformed position
end
hold on
for i = 1:nnodes
    text(nx(i),ny(i),num2str(i),'Fontsize',14);
end
legend('$$\;\mathrm{\textbf{Undeformed config.}}$$', ...
    '$$\;\mathrm{\textbf{Deformed config.}}$$', ...
    'FontSize',14,'Location','northeast','Interpreter','latex');
axis equal
end