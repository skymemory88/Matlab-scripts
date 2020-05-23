function idx=latmod(params,r_ij)
% computes an index for the matrices of dipolar interactions 

L=params.L;

% the coordinates of rij give 3 indices, and we transform them into 3
% integer nonnegative indices
x=2*(r_ij(1)+L+0.5)+1;
y=2*(r_ij(2)+L+0.5)+1;
z=4*(r_ij(3)+L+0.75)+1;

% we now transform these 3 indices into a single one
% maxX=2*(2*(L+0.5))+1; % 2*(r_ij_max+L+0.5)+1
maxY=2*(2*(L+0.5))+1; % 2*(r_ij_max+L+0.5)+1
maxZ=4*(2*(L+0.75))+1; % 4*(r_ij_max+L+0.5)+1


idx=((x-1)*(2*maxY+1)+(y-1))*(4*maxZ+1)+z;

end