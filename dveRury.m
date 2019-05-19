% Prienik dvoch lubovolnych rur
clear all; close all;
epsd=1.e-13;
%A1=input('Bod A1 na osi 1. rury: '); 
A1=[-5 2 15]; A1=A1(:) ;               % stlpcove vektory 1 2 3    1 0 0
%B1=input('Bod B1 na osi 1. rury: '); 
B1=[8 2 -11]; B1=B1(:) ;               % stlpcove vektory 2 2 1    0 2 0
% s1=input('Smerovy vektor s1 osi 1. rury: '); 
% s1=[1 0 -2];  
s1=B1-A1; % smerovy vektor 1. rury
s1=s1(:); s1=s1/norm(s1); % jednotkovy vektor
%A2=input('Bod A2 na osi 2. rury: '); 
A2=[-1 -3 -13 ]; A2=A2(:);   % 2 0 -1     0 0 5
%B2=input('Bod B2 na osi 2. rury: '); 
B2=[6 4 15]; B2=B2(:);   % 3 1 3    3 1 1      4 2 -3    4 2 5
s2=B2-A2; % smerovy vektor 2. rury
% s2=input('Smerovy vektor s2 osi 2. rury: ');
% s2=[1 1 2]; 
 s2=s2(:); s2=s2/norm(s2); % jednotkovy vektor
r1=input('Polomer 1. rury: ');
r2=input('Polomer 2. rury: ');

vymena=0;
if r2>r1, % vymena rur, na konci vymenime naspat
  disp('Pocas rotacii a vypoctov prieniku vymenime poradie rur.');
  vymena=1;
  A=A1; A1=A2; A2=A;
  B=B1; B1=B2; B2=B;
  s=s1; s1=s2; s2=s;
  r=r1; r1=r2; r2=r;	
end

n=90;
lambda=2; 

s12=s1'*s2;  % s12 je skalarny sucin vektorov s1 a s2
uhol=acos(abs(s12))*180/pi;
disp(['Uhol osi je priblizne ',num2str(uhol,4),' stupnov.'])

if (rank([s1 s2])==1), 
	disp('Osi su rovnobezne alebo splyvaju!');
	disp('Tento pripad je nezaujimavy, rezy su rovnobezne osiam')
	
	
else
	disp('Osi su mimobezne alebo roznobezne!');	
   a=[1 -s12; -s12 1]; % matica systemu na urcenie parametrov priecky predpokladanych mimo/roznobeziek
   b=[(A2-A1)'*s1; (A1-A2)'*s2];  % prava strana 
   tuv=inv(a)*b; % parametre krajnych bodov priecky
   P1=A1+tuv(1)*s1; % Krajny bod z priamky p1 - osi rury 1
   P2=A2+tuv(2)*s2; % Krajny bod z priamky p2 - osi rury 2
   d=norm(P2-P1);
   disp(['Vzdialenost osi je ',num2str(d),'.'])
%plot3([A1(1),P1(1)],[A1(2),P1(2)],[A1(3),P1(3)],'b'); % usecka na prvej priamke
   plot3([A1(1),P1(1)],[A1(2),P1(2)],[A1(3),P1(3)],'r'); % usecka na prvej priamke
%   plot3([A1(1),2*P1(1)-A1(1)],[A1(2),2*P1(2)-A1(2)],[A1(3),2*P1(3)-A1(3)],'r'); % usecka na prvej priamke
   hold on;
   plot3(A1(1),A1(2),A1(3),'*'); plot3(B1(1),B1(2),B1(3),'*'); plot3(A2(1),A2(2),A2(3),'*'); plot3(B2(1),B2(2),B2(3),'*');
   
%plot3([A2(1),P2(1)],[A2(2),P2(2)],[A2(3),P2(3)],'g'); % usecka na druhej priamke
   plot3([A2(1),P2(1)],[A2(2),P2(2)],[A2(3),P2(3)],'r'); % usecka na druhej priamke
%   plot3([A2(1),2*P2(1)-A2(1)],[A2(2),2*P2(2)-A2(2)],[A2(3),2*P2(3)-A2(3)],'r'); % usecka na druhej priamke
   plot3([P1(1),P2(1)],[P1(2),P2(2)],[P1(3),P2(3)],'r'); % priecka
   if d>epsd,
     disp('Osi su mimobezne!')	
      d1=P2-P1;  d1=d1/norm(d1); % jednotkovy vektor priecky od osi rury 1
      d2=-d1; % jednotkovy vektor priecky od osi rury 2
  else   
     disp('Osi su prakticky roznobezne!')	
     d1=cross(s1,s2);
     d2=-d1;
  end;   
      w1=cross(d1,s1);  % druhy jednotkovy vektor kolmy na os rury 1      
      w2=cross(d2,s2);  % druhy jednotkovy vektor kolmy na os rury 2      
      phi=0:(2*pi)/n:(2*pi); % chystame sa kreslit kruznicu
      bk11=A1*ones(size(phi))+r1*(d1*cos(phi)+w1*sin(phi)); % kruznica okolo bodu A1 s polomerom r1
%      bk12=bk11+lambda*(P1-A1)*ones(size(phi)); % kruznica oproti
      bk12=bk11+(B1-A1)*ones(size(phi)); % kruznica oproti
      bk21=A2*ones(size(phi))+r2*(d2*cos(phi)+w2*sin(phi)); % kruznica okolo bodu A2 s polomerom r2
 %    bk22=bk21+lambda*(P2-A2)*ones(size(phi)); % kruznica oproti
      bk22=bk21+(B2-A2)*ones(size(phi)); % kruznica oproti
      plot3(bk11(1,:),bk11(2,:),bk11(3,:),'m'); % hranicna kruznica okolo bodu A1
      plot3(bk12(1,:),bk12(2,:),bk12(3,:),'m');% hranicna kruznica okolo bodu oproti bodu A1
      plot3(bk21(1,:),bk21(2,:),bk21(3,:),'g'); % hranicna kruznica okolo bodu A2
      plot3(bk22(1,:),bk22(2,:),bk22(3,:),'g');% hranicna kruznica okolo bodu oproti bodu A2
      for k=1:n+1,
      	plot3([bk11(1,k),bk12(1,k)],[bk11(2,k),bk12(2,k)],[bk11(3,k),bk12(3,k)],'m'); % plast prvej rury
      	plot3([bk21(1,k),bk22(1,k)],[bk21(2,k),bk22(2,k)],[bk21(3,k),bk22(3,k)],'g'); % plast druhej rury
      end
      
%Q=rotationMatrixZ(s1);   % 
%s1r=Q*s1, s2r=Q*s2, d1r=Q*d1,
%Qz=[d1r(1) d1r(2) 0; -d1r(2) d1r(1) 0;0 0 1];   % rotacia okolo osi z    [x y;-y x]*[d1 d2]'=[1 0]'   d1*d1-(-d2)*d2=1   -d2*d1+d1*d2=0
%s1rr=Qz*s1r, s2rr=Qz*s2r, d1rr=Qz*d1r,
%Qr=Qz*Q; % vysledna matica rotacii

disp('Matica rotacii:')
Qr=rotationMatrix(s1,d1)
disp(' ')
A1r=Qr*(A1-P1);
B1r=Qr*(B1-P1);
disp(A1r)
disp(B1r)
disp(' ')
A2r=Qr*(A2-P1);
B2r=Qr*(B2-P1);
disp(A2r)
disp(B2r)

%[s1 s2 d1]
%S=Qr*[s1 s2 d1]
%Qr'*S

%Q*v/norm(v)

disp(['Po rotacii pomocou Qr sa krajny bod priecky z osi 1 posuva do bodu (0,0,0) =rotacia okolo P1  a  krajny bod z osi 2 do bodu (',num2str(d,4),',0,0).'])

   if (d>(r1+r2)),
      disp('Rury sa nepretinaju!');
   elseif (d==(r1+r2)),
       xd=Qr'*[r1;0;0]+P1;
      disp(['Rury sa dotykaju v bode ',num2str(xd')]);
      plot3(xd(1),xd(2),xd(3),'*k');
   else
      disp('Rury sa pretinaju!');   	
      s2r=Qr*s2;
      otoc1=0;
      t1rmin=A1r(3);
      t1rmax=B1r(3);
      figure(2); 
      plot([t1rmin t1rmin],[-pi*r1 pi*r1],'b');
      hold on;
      plot([t1rmax t1rmax],[-pi*r1 pi*r1],'r');
      if (t1rmin>t1rmax),
      	otoc1=1;
      	tt=t1min;
      	t1rmin=t1rmax;
      	t1rmax=tt;
      end;	
      otoc2=0;
      t2rmin=A2r(2)/s2r(2);
      t2rmax=B2r(2)/s2r(2);
      figure(3); 
      plot([t2rmin t2rmin],[0 2*pi*r2],'b');
      hold on;
      plot([t2rmax t2rmax],[0 2*pi*r2],'r');
      if (t2rmin>t2rmax),
      	otoc2=1;
      	tt=t2rmin;
      	t2rmin=t2rmax;
      	t2rmax=tt;
      end;	

      s=[0;s2r(2);s2r(3)]; % vektor kolmy na os 2. rury
%      vo=[0;-s(3);s(2)];    % 2. vektor kolmy na os 2. rury   - nepotrebujeme ho :)
      if (d+r2<=r1)   % dve krivky - jedna "vpredu" jedna "vzadu"
         disp('Prienik pozostava z 2 kriviek.')
          np=length(phi); % uhol phi sa meni od 0 do 2*pi
          rura2a(:,2)=r2*phi;
          rura2b(:,2)=rura2a(:,2);
          for k=1:np, % pocitame budy prieniku na sietke uhlov
             xx=d+r2*cos(phi(k));
             y1=sqrt(r1^2-xx^2);   % predna krivka
             y2=-y1;                         % zadna krivka
             pom=r2*sin(phi(k));
             t1=(y1-s(3)*pom)/s(2);  % parameter vpredu
             t2=(y2-s(3)*pom)/s(2);	  % parameter vzadu
             z1=t1*s(3)-s(2)*pom;   % predne z
             z2=t2*s(3)-s(2)*pom;   % zadne z
             xp1(:,k)=[xx;y1;z1];    % matica bodov prednej krivky
             xp2(:,k)=[xx;y2;z2];    % matica bodov zadnej krivky
             rura1a(k,1)=z1;
             rura1a(k,2)=atan(y1/xx); 
             rura1b(k,1)=z2;
             rura1b(k,2)=atan(y2/xx);
             if (xx<0),
                 if (y1>0), 
             	rura1a(k,2)=rura1a(k,2)+pi; 
                 else
                 	rura1a(k,2)=rura1a(k,2)-pi;
                 end;	
                 if (y2>0),	 	
             	rura1b(k,2)=rura1b(k,2)+pi; 
                 else	
             	rura1b(k,2)=rura1b(k,2)-pi;
                 end;
             end;    
             rura2a(k,1)=t1;
             rura2b(k,1)=t2;             
          end
          xp1=Qr'*xp1; xp1=xp1+[P1(1)*ones(1,np);P1(2)*ones(1,np);P1(3)*ones(1,np)];
          xp2=Qr'*xp2; xp2=xp2+[P1(1)*ones(1,np);P1(2)*ones(1,np);P1(3)*ones(1,np)];
          figure(1);
          plot3(xp1(1,:),xp1(2,:),xp1(3,:),'k'); 	
          plot3(xp2(1,:),xp2(2,:),xp2(3,:),'r');
          figure(2);
          plot(rura1a(:,1),rura1a(:,2),'b'); 	
          plot(rura1b(:,1),rura1b(:,2),'r');
          tmin=min([min(rura1a(:,1)),min(rura1b(:,1))]);
          tmax=max([max(rura1a(:,1)),max(rura1b(:,1))]);
          if t1rmin<tmin, tmin=t1rmin; end;
          if t1rmax>tmax, tmax=t1rmax; end;
          plot([tmin tmax],[-pi*r1 -pi*r1],'k');
          plot([tmin tmax],[pi*r1 pi*r1],'k');
          axis('equal');
          axis([1.05*tmin-0.05*tmax 1.05*tmax-0.05*tmin -3.5*r1  3.5*r1]);
          title('Rozvinutie krivky na 1. ruru'); 	
          figure(3);
          plot(rura2a(:,1),rura2a(:,2),'b'); 	
          plot(rura2b(:,1),rura2b(:,2),'r');
          tmin=min([min(rura2a(:,1)),min(rura2b(:,1))]);
          tmax=max([max(rura2a(:,1)),max(rura2b(:,1))]);
          if t2rmin<tmin, tmin=t2rmin; end;
          if t2rmax>tmax, tmax=t2rmax; end;
          plot([tmin tmax],[0 0],'k');
          plot([tmin tmax],[2*pi*r2 2*pi*r2],'k');
          axis('equal');
          axis([1.05*tmin-0.05*tmax 1.05*tmax-0.05*tmin -0.314*r2  6.594*r2]);
          title('Rozvinutie krivky na 2. ruru'); 	
      else % jedina krivka
         disp('Prienik tvori 1 krivka.')
         phimax=acos((r1-d)/r2); %    phimax <= phi <= 2*pi-phimax 
         disp(['phimax=',num2str(180*phimax/pi)])
         phid=phimax:(2*(pi-phimax)/n):(2*pi-phimax);
         np=length(phid); % uhol phi sa meni od 0 do 2*pi
          rura2a(:,2)=r2*phid;
          rura2b(:,2)=rura2a(:,2);
         for k=1:np, % pocitame budy prieniku na sietke uhlov
             xx=d+r2*cos(phid(k));
             y1=sqrt(r1^2-xx^2);   % predna krivka
             y2=-y1;                         % zadna krivka
             pom=r2*sin(phid(k));
             t1=(y1-s(3)*pom)/s(2);  % parameter vpredu
             t2=(y2-s(3)*pom)/s(2);	  % parameter vzadu
             z1=t1*s(3)-s(2)*pom;   % predne z
             z2=t2*s(3)-s(2)*pom;   % zadne z
             xp1(:,k)=[xx;y1;z1];    % matica bodov prednej krivky
             xp2(:,k)=[xx;y2;z2];    % matica bodov zadnej krivky
             rura1a(k,1)=z1;
             rura1a(k,2)=atan(y1/xx); 
             rura1b(k,1)=z2;
             rura1b(k,2)=atan(y2/xx);
             if (xx<0),
                 if (y1>0), 
             	rura1a(k,2)=rura1a(k,2)+pi; 
                 else
                 	rura1a(k,2)=rura1a(k,2)-pi;
                 end;	
                 if (y2>0),	 	
             	rura1b(k,2)=rura1b(k,2)+pi; 
                 else	
             	rura1b(k,2)=rura1b(k,2)-pi;
                 end;
             end;    
             rura2a(k,1)=t1;
             rura2b(k,1)=t2;             
          end
          xp1=Qr'*xp1; xp1=xp1+[P1(1)*ones(1,np);P1(2)*ones(1,np);P1(3)*ones(1,np)];
          xp2=Qr'*xp2; xp2=xp2+[P1(1)*ones(1,np);P1(2)*ones(1,np);P1(3)*ones(1,np)];
          figure(1);
          plot3(xp1(1,:),xp1(2,:),xp1(3,:),'k'); 	
          plot3(xp2(1,:),xp2(2,:),xp2(3,:),'r'); 
          figure(2);
          plot(rura1a(:,1),rura1a(:,2),'b'); 	
          plot(rura1b(:,1),rura1b(:,2),'r');
          tmin=min([min(rura1a(:,1)),min(rura1b(:,1))]);
          tmax=max([max(rura1a(:,1)),max(rura1b(:,1))]);
          if t1rmin<tmin, tmin=t1rmin; end;
          if t1rmax>tmax, tmax=t1rmax; end;
          plot([tmin tmax],[-pi*r1 -pi*r1],'k');
          plot([tmin tmax],[pi*r1 pi*r1],'k');
          axis('equal');
          axis([1.05*tmin-0.05*tmax 1.05*tmax-0.05*tmin -3.5*r1  3.5*r1]);
          title('Rozvinutie krivky na 1. ruru'); 	
          figure(3);
          plot(rura2a(:,1),rura2a(:,2),'b'); 	
          plot(rura2b(:,1),rura2b(:,2),'r');
          tmin=min([min(rura2a(:,1)),min(rura2b(:,1))]);
          tmax=max([max(rura2a(:,1)),max(rura2b(:,1))]);
          if t2rmin<tmin, tmin=t2rmin; end;
          if t2rmax>tmax, tmax=t2rmax; end;
          plot([tmin tmax],[0 0],'k');
          plot([tmin tmax],[2*pi*r2 2*pi*r2],'k');
          axis('equal');
          axis([1.05*tmin-0.05*tmax 1.05*tmax-0.05*tmin -0.314*r2  6.594*r2]);  
          title('Rozvinutie krivky na 2. ruru'); 	
      end
      
      
   end;
end;

% rotacia 
%v=input('Zadaj vektor v: '); v=v(:);
