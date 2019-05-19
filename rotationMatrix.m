function [Q]=rotationMatrix(v,d)
% Function for evaluation of the rotation matrix
% Given vector of direction v and length 1 is rotated to be vector (0,0,1)
% Given vector d of length 1, orthogonal to v is rotated to (1,0,0)
v=v(:); % premena na stlpcovy vektor
cx=v(3)/sqrt(v(2)^2+v(3)^2); % kosinus rotacie okolo osi x
sx=v(2)/sqrt(v(2)^2+v(3)^2); % sinus rotacie okolo osi x 
Qx=[1 0 0;0 cx -sx; 0 sx cx]; % matica rotacie okolo osi x
vx=Qx*v; 
cy=vx(3)/sqrt(vx(1)^2+vx(3)^2);  % kosinus rotacie okolo osi y
sy=vx(1)/sqrt(vx(1)^2+vx(3)^2);  % sinus rotacie okolo osi y
Qy=[cy 0 -sy;0 1 0; sy 0 cy]; % matica rotacie okolo osi y
Q=Qy*Qx; % vysledna matica rotacie suradnicoveho systemu  
 dr=Q*d;
Qz=[dr(1) dr(2) 0; -dr(2) dr(1) 0;0 0 1]; % rotacia okolo osi z
Q=Qz*Q; % vysledna rotacia
