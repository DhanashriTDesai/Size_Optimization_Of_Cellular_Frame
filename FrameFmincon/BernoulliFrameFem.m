% Matlab function to run 2D Bernoulli FEM code
% Written by Dhanashri Desai 

function [SE,gdisp,gstiffness,reactions,nodalforces,elementforces,stresses] = ...
    BernoulliFrameFem(nnodes,nelem,A,E,I,L,nx,ny,elemcon,dispbcdof,extloaddof1,extloaddof2,ExtLoad)

% Initialize global arrays to zero
gdisp = zeros(3*nnodes,1); % global displacement vector
gforce = zeros(3*nnodes,1); % global force vector
gstiffness = zeros(3*nnodes,3*nnodes); % global stiffness matrix

gforce(extloaddof1) = ExtLoad;
gforce(extloaddof2) = ExtLoad;
gdisp(dispbcdof) = 0;

for i = 1:nelem % gstiffness and gforce assembly
    node1 = elemcon(i,1); node2 = elemcon(i,2); % nodes of element
    
    id1 = 3*(node1-1)+1; id2 = id1+1; id3 = id2+1;  % dofs corresponding to element
    id4 = 3*(node2-1)+1; id5 = id4+1; id6 = id5+1;   
    elemdof = [id1 id2 id3 id4 id5 id6]; 

    a=E(i)*A(i)/L(i); b=12*E(i)*I(i)/L(i)^3; c=6*E(i)*I(i)/L(i)^2;
    d=4*E(i)*I(i)/L(i); e=2*E(i)*I(i)/L(i);

    localstiffness = [a 0 0 -a 0 0 ;
                      0 b c 0 b c ;
                      0 c d 0 -c e ;
                     -a 0 0 a 0 0 ;
                      0 b -c 0 b -c ;
                      0 c e 0 -c d]; 

    l = (nx(node2)-nx(node1))/L(i); % l=cos(theta) m=sin(theta)
    m = (ny(node2)-ny(node1))/L(i); 
    
    Lmat = [l m 0 0 0 0 ; % transformation matrix
           -m l 0 0 0 0 ;
            0 0 1 0 0 0 ;
            0 0 0 l m 0 ;
            0 0 0 -m l 0 ;
            0 0 0 0 0 1];

    elemstiffness = Lmat' * localstiffness * Lmat; % element stiffness matrix
    gstiffness(elemdof,elemdof) = gstiffness(elemdof,elemdof) + elemstiffness; % assembly
end

alldof = (1:3*nnodes)'; % total dofs
activedof = setdiff(alldof,dispbcdof); % active dofs i.e. one not subjected to displacement BC

gdisp(activedof) = gstiffness(activedof,activedof) \ gforce(activedof); % solving linear system
nodalforces = gstiffness * gdisp; % forces at nodes = KU
reactions = nodalforces(dispbcdof); % reactions at dispbc nodes

SE=0.5*gdisp'*gstiffness*gdisp; % SE of entire structure

for i = 1:nelem % Calculation of internal forces and stresses in elements
    node1 = elemcon(i,1); node2 = elemcon(i,2); % nodes of element
    
    id1 = 3*(node1-1)+1; id2 = id1+1; id3 = id2+1;  % dofs corresponding to element
    id4 = 3*(node2-1)+1; id5 = id4+1; id6 = id5+1;   
    elemdof = [id1 id2 id3 id4 id5 id6]; 

    a=E(i)*A(i)/L(i); b=12*E(i)*I(i)/L(i)^3; c=6*E(i)*I(i)/L(i)^2;
    d=4*E(i)*I(i)/L(i); e=2*E(i)*I(i)/L(i);

    localstiffness = [a 0 0 -a 0 0 ;
                      0 b c 0 b c ;
                      0 c d 0 -c e ;
                     -a 0 0 a 0 0 ;
                      0 b -c 0 b -c ;
                      0 c e 0 -c d]; 

    l = (nx(node2)-nx(node1))/L(i); % l=cos(theta) m=sin(theta)
    m = (ny(node2)-ny(node1))/L(i); 
    
    Lmat = [l m 0 0 0 0 ; % transformation matrix
           -m l 0 0 0 0 ;
            0 0 1 0 0 0 ;
            0 0 0 l m 0 ;
            0 0 0 -m l 0 ;
            0 0 0 0 0 1];

    elemstiffness = Lmat' * localstiffness * Lmat; % element stiffness matrix

    Uelem = gdisp(elemdof);
    Pint = elemstiffness*Uelem;
    elementforces(i,1) = Pint(1); elementforces(i,2) = Pint(2);
    elementforces(i,3) = Pint(3); elementforces(i,4) = Pint(4);
    elementforces(i,5) = Pint(5); elementforces(i,6) = Pint(6);
    stresses(i,1) = abs(elementforces(i,1))/A(i)*1E-6; % stress in element
    
    nxd1 = nx(node1) + Uelem(1); nyd1 = ny(node1) + Uelem(2);
    nxd2 = nx(node2) + Uelem(4); nyd2 = ny(node2) + Uelem(5);
    Ld = sqrt( (nxd2-nxd1)^2 + (nyd2-nyd1)^2 );
    if Ld>L(i) % member in tension
        elementforces(i,7)=1; stresses(i,2)=1;
    elseif Ld<L(i) % member in compression
        elementforces(i,7)=-1; stresses(i,2)=-1;
    else % undeformed member
        elementforces(i,7)=0; stresses(i,2)=0;
    end
end
stresses(:,1)=round(stresses(:,1),4); % round upto 4 decimal places

end