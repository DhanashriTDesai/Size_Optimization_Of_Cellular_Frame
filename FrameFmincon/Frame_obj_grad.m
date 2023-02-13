function [obj,grad] = Frame_obj_grad(t0)

global nCellH nCellV Nrows Ncolumns nodeID nnodes nelem nx ny elemcon L E extloaddof1 extloaddof2 dispbcdof ExtLoadDir ExtLoad bStrut maxstress maxstressId

A = t0*bStrut; I = bStrut*t0.^3/12; % c/s and area MI of of struts
[SE,gdisp,gstiffness,reactions,nodalforces,elementforces,stresses] = ... % Run frame fem code
    BernoulliFrameFem(nnodes,nelem,A,E,I,L,nx,ny,elemcon,dispbcdof,extloaddof1,extloaddof2,ExtLoad);

obj = 0.5 * gdisp' * gstiffness * gdisp;

for i = 1:nelem
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
    elemdisp = [gdisp(id1) ; gdisp(id2) ; gdisp(id3) ; ... % element displacement vector
                gdisp(id4) ; gdisp(id5) ; gdisp(id6)];
    SEelem = elemdisp' * elemstiffness * elemdisp; % SE of element 
    SEelemV(i) = SEelem / (A(i)*L(i)); % Volumetric SE of element

    grad(i,1) = - elemdisp' * elemstiffness * elemdisp; % gradient of obj fun
end
maxstress = max(stresses(:,1)); maxstressId = find(stresses == maxstress);
% disp(['Maximum stress ' num2str(max(stresses(:,1))) ' MPa in elements:']); disp(['element' '  node1' '  node2']); disp([idx elemcon(idx,:)]);
% disp(['SE of frame in Joules: ' num2str(SE)]);

end