function res = variants()
L=1000; E=100e3; nu=1/3; G=E/(2*(1+nu));
b=L/40; a=2*b; t=L/400;
A=a*b-(a-2*t)*(b-2*t);
Iz=(b*a^3-(b-2*t)*(a-2*t)^3)/12;
Iy=(a*b^3-(a-2*t)*(b-2*t)^3)/12;
am=a-t; bm=b-t; Om=am*bm; J=4*Om^2/(2*(am+bm)/t);
q0x=5; u_load=L/20;

c_nod = [0 0 0; L/6 0 L/3; L/3 0 2*L/3; L/2 0 L; 5*L/6 0 L; 5*L/6 0 L;
         7*L/6 0 L; 3*L/2 0 L; 3*L/2 L/4 L; 3*L/2 L/2 L; 3*L/2 L/2 2*L/3;
         3*L/2 L/2 L/3; 3*L/2 L/2 0; 3*L/2 -L/4 L; 3*L/2 -L/2 L;
         3*L/2 -L/2 2*L/3; 3*L/2 -L/2 L/3; 3*L/2 -L/2 0];
con = [1 2;2 3;3 4;4 5;6 7;7 8;8 9;9 10;8 14;14 15;10 11;11 12;12 13;15 16;16 17;17 18];
n = 18*6;
vref = cell(16,1); vref(1:3)={[-1;0;0]}; vref(4:10)={[0;0;1]}; vref(11:16)={[0;-1;0]};
Le=zeros(16,1); Te=zeros(12,12,16); cg=zeros(16,12);
for i=1:16
  nd=con(i,:);
  [Le(i), Te(:,:,i)] = gv(c_nod(nd(1),:), c_nod(nd(2),:), vref{i});
  cg(i,:) = [6*(nd(1)-1)+(1:6), 6*(nd(2)-1)+(1:6)];
end
K=zeros(n); Kg=zeros(12,12,16);
for i=1:16
  Ke = ke(E,G,Le(i),A,Iy,Iz,J);
  Kge = Te(:,:,i)'*Ke*Te(:,:,i);
  K(cg(i,:),cg(i,:)) = K(cg(i,:),cg(i,:)) + Kge;
  Kg(:,:,i)=Kge;
end
F=zeros(n,1);
for e=1:16
  nA=con(e,1); nB=con(e,2);
  if ismember(e,[1 2 3 4 5 6]), qA=0; qB=0; else qA=q0x*c_nod(nA,3)/L; qB=q0x*c_nod(nB,3)/L; end
  fl = cde(Te(1:3,1:3,e), Le(e), [qA;0;0], [qB;0;0]);
  F(cg(e,:)) = F(cg(e,:)) + Te(:,:,e)'*fl;
end

casos = {
 'A: rotula 3 trasl, apoyos empotrados 1/13/18, u3=50'
 'B: sin rotula'
 'C: rotula 3 trasl, n1 empotrado, n13/n18 empotrados, u3=-50'
 'D: rotula, columna cargada tambien'
 'E: rotula, apoyos n1/n13/n18 empotrados + cargas nodales P1,P2'
 'F: rotula, n1 empotrado, n13/n18 solo traslacion'
 };
for c=1:size(casos,1)
  fprintf('\n=== %s\n', casos{c,1});
  escl = []; maes = []; Fn = F;
  switch c
    case 1
      escl = 31:33; maes = 25:27;
    case 2
      escl = []; maes = [];
    case 3
      escl = 31:33; maes = 25:27; sgn=-1;
    case 4
      escl = 31:33; maes = 25:27;
      for e=1:4
        nA=con(e,1); nB=con(e,2);
        fl = cde(Te(1:3,1:3,e), Le(e), [q0x*c_nod(nA,3)/L;0;0], [q0x*c_nod(nB,3)/L;0;0]);
        Fn(cg(e,:)) = Fn(cg(e,:)) + Te(:,:,e)'*fl;
      end
    case 5
      escl = 31:33; maes = 25:27;
      Fn(6*(8-1)+3) = Fn(6*(8-1)+3) - 10000;
      Fn(6*(14-1)+3) = Fn(6*(14-1)+3) - 5000;
    case 6
      escl = 31:33; maes = 25:27;
  end
  sgn = 1; if c==3, sgn=-1; end
  if c==6
    sup = [1:3, 73:78, 103:108];
  else
    sup = [1:6, 73:78, 103:108];
  end
  imp = 13;
  libres = setdiff(1:n, [sup, escl, imp]);
  T = zeros(n, numel(libres)); T(libres,:) = eye(numel(libres));
  for k=1:numel(escl), T(escl(k),:) = T(maes(k),:); end
  U0 = zeros(n,1); U0(imp) = sgn*u_load;
  q = (T'*K*T) \ (T'*(Fn - K*U0));
  U = T*q + U0;
  R = K*U - Fn;
  [um, iu] = max(abs(U));
  fprintf('  max|U| = %.3f mm (nodo %d, comp %d)\n', um, ceil(iu/6), mod(iu-1,6)+1);
  fprintf('  U(n5)=%s\n', mat2str(U(25:27).',5));
  fprintf('  FR1=%s\n', mat2str(R(1:6).',5));
  fprintf('  FR2=%s\n', mat2str(R(73:78).',5));
  if ~isempty(escl)
    fprintf('  Frot=%s\n', mat2str(R(escl).',5));
  end
  fprintf('  resid other=%.2e\n', max(abs(R(setdiff(libres,maes)))));
end

end

function [Ke] = ke(E,G,Le,A,Iy,Iz,J)
  ka=E*A/Le; kt=G*J/Le;
  az=12*E*Iz/Le^3; bz=6*E*Iz/Le^2; cz=4*E*Iz/Le; dz=2*E*Iz/Le;
  ay=12*E*Iy/Le^3; by=6*E*Iy/Le^2; cy=4*E*Iy/Le; dy=2*E*Iy/Le;
  Ke=zeros(12);
  Ke([1 7],[1 7])=[ka -ka;-ka ka];
  Ke([4 10],[4 10])=[kt -kt;-kt kt];
  Ke([2 6 8 12],[2 6 8 12])=[az bz -az bz; bz cz -bz dz; -az -bz az -bz; bz dz -bz cz];
  Ke([3 5 9 11],[3 5 9 11])=[ay -by -ay -by; -by cy by dy; -ay by ay by; -by dy by cy];
end
function f = cde(Re,Le,w1,w2)
  p=[Re*w1, Re*w2]; p1=p(:,1); dp=p(:,2)-p(:,1);
  f=[Le*(p1(1)/2+dp(1)/6); Le*(p1(2)/2+dp(2)*3/20); Le*(p1(3)/2+dp(3)*3/20); 0;
     -Le^2*(p1(3)/12+dp(3)/30); Le^2*(p1(2)/12+dp(2)/30);
     Le*(p1(1)/2+dp(1)/3); Le*(p1(2)/2+dp(2)*7/20); Le*(p1(3)/2+dp(3)*7/20); 0;
     Le^2*(p1(3)/12+dp(3)/20); Le^2*(-p1(2)/12-dp(2)/20)];
end
function [Le,Te] = gv(n1,n2,v)
  n1=n1(:); n2=n2(:); d=n2-n1; Le=norm(d); nx=d/Le;
  nz=cross(nx,v(:)); nz=nz/norm(nz); ny=cross(nz,nx);
  Re=[nx';ny';nz']; R0=zeros(3);
  Te=[Re R0 R0 R0; R0 Re R0 R0; R0 R0 Re R0; R0 R0 R0 Re];
end
