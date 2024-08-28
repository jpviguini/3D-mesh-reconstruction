close all;
clear;
clc;
tic

 %Leitura do arquivo
  [V,F] = leitura_obj('pasta_destino/face_rosto2.obj');  #'pasta_destino/face_rosto.obj'
  
  %Numero de vértices e de faces

  nverts = size(V,1);
  nfaces = size(F,1);

##%------------------------------ Contas Discretas -------------------------------
##
##
##
## 
##
  S=[];
  Face=[];
  L=[];
  
  %Definindo a estrela, a quantidade de ligação de cada ponto e os triângulos que 
  %o ponto pertence

  for i=1:nverts
      count3=1; %conta a quantidade de vertices na estrela
      count3a=1;%conta a quantidade de faces na estrela
      for j=1:nfaces
        if (i == F(j,1))
          S(i,count3)=F(j,2);
          Face(i,count3a)=j;
          count3=count3+1;
          S(i,count3)=F(j,3);
          count3=count3+1;
          count3a=count3a+1;
        elseif (i == F(j,2))
          S(i,count3)=F(j,3);
          Face(i,count3a)=j;
          count3=count3+1;
          S(i,count3)=F(j,1);
          count3=count3+1;
          count3a=count3a+1;
        elseif (i == F(j,3))
          S(i,count3)=F(j,1);
          Face(i,count3a)=j;
          count3=count3+1;  
          S(i,count3)=F(j,2);
          count3=count3+1;
          count3a=count3a+1;
        end
        L(i,1)=count3-1; %numero de pontos na estrela
        LL(i,1)=count3a-1; %numero de triangulos na estrela
      end  
  end
  
  
  Sc=[];
  for i=1:size(S,1)
    X=S(i,1:L(i))';
    Xc=unique(X,"rows");
    for j=1:size(Xc,1)
      Sc(i,j)=Xc(j,1);
    endfor
  endfor


for i=1:nverts
  count=0;  
  for j=1:size(Sc,2)
    if(Sc(i,j)!=0)
    count=count+1;
    endif
  endfor
  L(i,1)=count;
endfor

SS=zeros(nverts,max(L));
idx=zeros(nverts,1);

  for i=1:nverts
      temp=1;
      
      while(idx(i)==0)
      count33=1;
      m=Sc(i,temp); 
      for k=1:L(i)
        for j=1:nfaces
          if (i == F(j,1) && m == F(j,2))
            SS(i,count33)=F(j,2);
            if (F(j,3) == SS(i,1))
              break
            else
              SS(i,count33+1)=F(j,3);
              m=SS(i,count33+1);
              count33=count33+1;
              break
            endif
            end
          if (i == F(j,2)  && m == F(j,3))
            SS(i,count33)=F(j,3);
            if (F(j,1) == SS(i,1))
              break
            else
              SS(i,count33+1)=F(j,1);
              m=SS(i,count33+1);
              count33=count33+1;
              break
            endif
            end
          if (i == F(j,3)  && m == F(j,1))
            SS(i,count33)=F(j,1);
            if (F(j,2) == SS(i,1))
              break
            else
              SS(i,count33+1)=F(j,2);
              m=SS(i,count33+1);
              count33=count33+1;
              break  
            endif
          end
        end 
       if (SS(i,1)==0)
          m=Sc(i,2);
       endif 
      Lx(i)=count33;
    end
    
    if(L(i)==Lx(i))
      idx(i)=1;
    else
      temp=temp+1;
    end
    endwhile       
  end
  
  printf('Parte 1');

  %Descobrindo vetor normal no vértice

        %Normal da Face
  function [vn] = vetnorm(A,B,C)
    vn=cross((B-A),(C-A));
    vn=vn/norm(vn);  
  endfunction 

         %Área da Face

  function [af] = areaface(A,B,C)
    af=1/2*norm(cross(B-A,C-A));  
  endfunction 
  
 for i=1:nverts
    for j=1:LL(i)
      areafaces1(i,j)=areaface(V(F(Face(i,j),1),:),V(F(Face(i,j),2),:),V(F(Face(i,j),3),:));
    endfor
  endfor
  
  b=zeros(nverts,3);
  c=zeros(nverts,3);
  Nvt=zeros(nverts,3);
  
  for i=1:nverts
    sum33=[0 0 0];
    for j=1:LL(i)
      sum33 = sum33 + areafaces1(i,j) * vetnorm(V(F(Face(i,j),1),:),V(F(Face(i,j),2),:),V(F(Face(i,j),3),:));
    end  

    Nvt(i,:)=sum33/norm(sum33);
    
    xn=Nvt(i,:)';
    [Q R]=qr(xn);
    x1=Q(:,1);
    if (x1'*xn<0)
      Q=-Q;
    end
    x1=Q(:,1);
    b(i,:)=Q(:,2)';
    c(i,:)=Q(:,3)';
  end  

function [ang] = angulo(A,B,C)
  vet = dot(B-A,C-A);
  co = vet/(norm(B-A)*norm(C-A));
  ang = acos(co);
endfunction

Amixed=zeros(nverts,1);

for i=1:nverts
  for j=1:LL(i)
    Amixed(i)= Amixed(i) + 1/3 * areafaces1(i,j);
  endfor  
endfor

%Curvatura Gaussiana


Kd1=[];
Kd2=[];

for i=1:nverts
  sumK = 0.0;
    for j = 1:LL(i)
    if (F(Face(i,j),1) == i)
      sumK = sumK + angulo(V(F(Face(i,j),1),:),V(F(Face(i,j),2),:),V(F(Face(i,j),3),:));  
    endif
    if (F(Face(i,j),2) == i)
      sumK = sumK + angulo(V(F(Face(i,j),2),:),V(F(Face(i,j),1),:),V(F(Face(i,j),3),:));  
    endif
    if (F(Face(i,j),3) == i)
      sumK = sumK + angulo(V(F(Face(i,j),3),:),V(F(Face(i,j),2),:),V(F(Face(i,j),1),:));  
    endif    
  endfor
  
  % Curvatura Gaussiana 1: Sulivan
    
  Kd1(i,:)=2*pi-sumK;
  
  %Curvatura Gaussiana 2: K= (2pi-soma angulos)/A_mixed : Meyer
  
  Kd2(i,:)=Kd1(i,:)/Amixed(i,:);
  
##  endif  
endfor


  printf('Parte 2');


% Curvatura Média

sumHn3=zeros(nverts,1);

for i=1:nverts
  vp=zeros(2,L(i));
  countp=1;
  vp(1,countp)= 0;
  vp(2,countp)= 0;
  for j=1:LL(i)
    
    if (F(Face(i,j),1) == i)
      vj=F(Face(i,j),2);
      vp(2,countp)=F(Face(i,j),3);
      
      for k=1:nverts
        for l=1:LL(k)
          if ((F(Face(k,l),2) == i) && (F(Face(k,l),1) == vj) && (vp(1,countp) == 0) )
            vp(1,countp) = F(Face(k,l),3);
            eij=V(vj,:)-V(i,:);
            eij2=eij/norm(eij);
            nijk=vetnorm(V(i,:),V(vj,:),V(vp(1,countp),:));
            njil=vetnorm(V(vj,:),V(i,:),V(vp(2,countp),:));
            phi=atan2(dot(eij2,cross(nijk,njil)),dot(nijk,njil));
            sumHn3(i,:)=sumHn3(i,:) + phi*norm(eij);
          elseif ((F(Face(k,l),3) == i) && (F(Face(k,l),2) == vj) && (vp(1,countp) == 0) )
            vp(1,countp) = F(Face(k,l),1);
            eij=V(vj,:)-V(i,:);
            eij2=eij/norm(eij);
            nijk=vetnorm(V(i,:),V(vj,:),V(vp(1,countp),:));
            njil=vetnorm(V(vj,:),V(i,:),V(vp(2,countp),:));
            phi=atan2(dot(eij2,cross(nijk,njil)),dot(nijk,njil));
            sumHn3(i,:)=sumHn3(i,:) + phi*norm(eij);
          elseif ((F(Face(k,l),1) == i) && (F(Face(k,l),3) == vj) && (vp(1,countp) == 0) )
            vp(1,countp) = F(Face(k,l),2);
            eij=V(vj,:)-V(i,:);
            eij2=eij/norm(eij);
            nijk=vetnorm(V(i,:),V(vj,:),V(vp(1,countp),:));
            njil=vetnorm(V(vj,:),V(i,:),V(vp(2,countp),:));
            phi=atan2(dot(eij2,cross(nijk,njil)),dot(nijk,njil));
            sumHn3(i,:)=sumHn3(i,:) + phi*norm(eij);
          else
           continue;
          endif
        endfor
      endfor
      countp=countp+1;
    endif
      
     if (F(Face(i,j),2) == i)
      vj=F(Face(i,j),3);
      vp(2,countp)=F(Face(i,j),1);
      
      for k=1:nverts
        for l=1:LL(k)
          if ((F(Face(k,l),2) == i) && (F(Face(k,l),1) == vj) && (vp(1,countp) == 0) )
            vp(1,countp) = F(Face(k,l),3);
            eij=V(vj,:)-V(i,:);
            eij2=eij/norm(eij);
            nijk=vetnorm(V(i,:),V(vj,:),V(vp(1,countp),:));
            njil=vetnorm(V(vj,:),V(i,:),V(vp(2,countp),:));
            phi=atan2(dot(eij2,cross(nijk,njil)),dot(nijk,njil));
            sumHn3(i,:)=sumHn3(i,:) + phi*norm(eij);
          elseif ((F(Face(k,l),3) == i) && (F(Face(k,l),2) == vj) && (vp(1,countp) == 0) )
            vp(1,countp) = F(Face(k,l),1);
            eij=V(vj,:)-V(i,:);
            eij2=eij/norm(eij);
            nijk=vetnorm(V(i,:),V(vj,:),V(vp(1,countp),:));
            njil=vetnorm(V(vj,:),V(i,:),V(vp(2,countp),:));
            phi=atan2(dot(eij2,cross(nijk,njil)),dot(nijk,njil));
            sumHn3(i,:)=sumHn3(i,:) + phi*norm(eij);
          elseif ((F(Face(k,l),1) == i) && (F(Face(k,l),3) == vj) && (vp(1,countp) == 0) )
            vp(1,countp) = F(Face(k,l),2);
            eij=V(vj,:)-V(i,:);
            eij2=eij/norm(eij);
            nijk=vetnorm(V(i,:),V(vj,:),V(vp(1,countp),:));
            njil=vetnorm(V(vj,:),V(i,:),V(vp(2,countp),:));
            phi=atan2(dot(eij2,cross(nijk,njil)),dot(nijk,njil));
            sumHn3(i,:)=sumHn3(i,:) + phi*norm(eij);
          else
            continue;
          endif
        endfor
      endfor
      countp=countp+1;
    endif
    if (F(Face(i,j),3) == i)
      vj=F(Face(i,j),1);
      vp(2,countp)=F(Face(i,j),2);
      
      for k=1:nverts
        for l=1:LL(k)
          if ((F(Face(k,l),2) == i) && (F(Face(k,l),1) == vj) && (vp(1,countp) == 0) )
            vp(1,countp) = F(Face(k,l),3);
            eij=V(vj,:)-V(i,:);
            eij2=eij/norm(eij);
            nijk=vetnorm(V(i,:),V(vj,:),V(vp(1,countp),:));
            njil=vetnorm(V(vj,:),V(i,:),V(vp(2,countp),:));
            phi=atan2(dot(eij2,cross(nijk,njil)),dot(nijk,njil));
            sumHn3(i,:)=sumHn3(i,:) + phi*norm(eij);
          elseif ((F(Face(k,l),3) == i) && (F(Face(k,l),2) == vj) && (vp(1,countp) == 0) )
            vp(1,countp) = F(Face(k,l),1);
            eij=V(vj,:)-V(i,:);
            eij2=eij/norm(eij);
            nijk=vetnorm(V(i,:),V(vj,:),V(vp(1,countp),:));
            njil=vetnorm(V(vj,:),V(i,:),V(vp(2,countp),:));
            phi=atan2(dot(eij2,cross(nijk,njil)),dot(nijk,njil));
            sumHn3(i,:)=sumHn3(i,:) + phi*norm(eij);
          elseif ((F(Face(k,l),1) == i) && (F(Face(k,l),3) == vj) && (vp(1,countp) == 0) )
            vp(1,countp) = F(Face(k,l),2);
            eij=V(vj,:)-V(i,:);
            eij2=eij/norm(eij);
            nijk=vetnorm(V(i,:),V(vj,:),V(vp(1,countp),:));
            njil=vetnorm(V(vj,:),V(i,:),V(vp(2,countp),:));
            phi=atan2(dot(eij2,cross(nijk,njil)),dot(nijk,njil));
            sumHn3(i,:)=sumHn3(i,:) + phi*norm(eij);
          else
            continue;
          endif
        endfor
      endfor
      countp=countp+1;
    endif
  endfor
  
  
  % Curvatura Media 2: Crane
  
  Hd3(i,:) = (1/4*sumHn3(i,:))/Amixed(i,:);

endfor

printf('Parte 4');

% Curvaturas Principais 3

for i=1:nverts
  k14(i,:)=Hd3(i,:)+sqrt((Hd3(i,:))^2-(Kd2(i,:)));
  k24(i,:)=Hd3(i,:)-sqrt((Hd3(i,:))^2-(Kd2(i,:)));
  if ((Hd3(i,:))^2-(Kd2(i,:))>=1e-10)
    k13(i,:)=k14(i,:);
    k23(i,:)=k24(i,:);
  else
    k13(i,:)=Hd3(i,:);
    k23(i,:)=Hd3(i,:);
  endif
  
endfor


  Positivo3=[];
  Negativo3=[];
  count5551=1;
  count5651=1;
  for i=1:nverts  
    if((Hd3(i,:))^2-(Kd2(i,:))>1e-10)
      Positivo3(count5551,:)=V(i,:);
      count5551=count5551+1;
    end
    if ((Hd3(i,:))^2-(Kd2(i,:))<-1e-10)
      Negativo3(count5651,:)=V(i,:);
      count5651=count5651+1;
    end
  endfor
 
    count5721=1;
  for i=1:nfaces
    if ((((Hd3(F(i,1),1))^2-(Kd2(F(i,1),1)))>1e-10 && ((Hd3(F(i,2),1))^2-(Kd2(F(i,2),1)))<-1e-10))
      t2=((Hd3(F(i,1),1))^2-(Kd2(F(i,1),1)))/(((Hd3(F(i,1),1))^2-(Kd2(F(i,1),1)))-((Hd3(F(i,2),1))^2-(Kd2(F(i,2),1))));
      MMw1(count5721,:)=V(F(i,1),:)+t2*(V(F(i,2),:)-V(F(i,1),:));
      PMM1(count5721,1)=F(i,1);
      PMM1(count5721,2)=F(i,2);
      count5721=count5721+1;
    end
    if ((((Hd3(F(i,2),1))^2-(Kd2(F(i,2),1)))>1e-10 && ((Hd3(F(i,3),1))^2-(Kd2(F(i,3),1)))<-1e-10))
      t2=((Hd3(F(i,2),1))^2-(Kd2(F(i,2),1)))/(((Hd3(F(i,2),1))^2-(Kd2(F(i,2),1)))-((Hd3(F(i,3),1))^2-(Kd2(F(i,3),1))));
      MMw1(count5721,:)=V(F(i,2),:)+t2*(V(F(i,3),:)-V(F(i,2),:));
      PMM1(count5721,1)=F(i,2);
      PMM1(count5721,2)=F(i,3);
      count5721=count5721+1;
    end
    if ((((Hd3(F(i,3),1))^2-(Kd2(F(i,3),1)))>1e-10 && ((Hd3(F(i,1),1))^2-(Kd2(F(i,1),1)))<-1e-10) )
      t2=((Hd3(F(i,3),1))^2-(Kd2(F(i,3),1)))/(((Hd3(F(i,3),1))^2-(Kd2(F(i,3),1)))-((Hd3(F(i,1),1))^2-(Kd2(F(i,1),1))));
      MMw1(count5721,:)=V(F(i,3),:)+t2*(V(F(i,1),:)-V(F(i,3),:));
      PMM1(count5721,1)=F(i,3);
      PMM1(count5721,2)=F(i,1);
      count5721=count5721+1;
    endif
  end
  
  
  
##MM3w1 = unique(MMw1,'rows');
  
  printf('Parte 5');
  
  
##  MM3 = unique(MM,'rows');
##  
  


kn=[];
d=[];
for i=1:nverts
  for j=L(i):-1:1
    kn(i,j)=2*dot(V(i,:)-V(SS(i,j),:),Nvt(i,:))/(norm(V(i,:)-V(SS(i,j),:)))^2;
    d(i,3*(j-1)+1:3*(j-1)+3)=((V(SS(i,j),:)-V(i,:))-((dot(V(SS(i,j),:)-V(i,:),Nvt(i,:)))*Nvt(i,:)))/norm((V(SS(i,j),:)-V(i,:))-((dot(V(SS(i,j),:)-V(i,:),Nvt(i,:)))*Nvt(i,:)));
  endfor
endfor


ujvj=zeros(nverts,2*max(L)+2);
  for i=1:nverts
    count5=1;
    for j=1:3:3*L(i)
      ujvj(i,count5)=d(i,j)*b(i,1)+d(i,j+1)*b(i,2)+d(i,j+2)*b(i,3);
      count5=count5+1;
      ujvj(i,count5)=d(i,j)*c(i,1)+d(i,j+1)*c(i,2)+d(i,j+2)*c(i,3);
      count5=count5+1;
    end   
  end

  P=zeros(nverts,3);
  for i=1:nverts
    count6=1;
    A=zeros(L(i)+1,3);
    B=zeros(L(i)+1,1);
    for j=1:2:2*L(i)
      A(count6,:)=[(ujvj(i,j))^2 ujvj(i,j)*ujvj(i,j+1) (ujvj(i,j+1))^2];
      B(count6,:)=[kn(i,j-(count6-1))];
      count6=count6+1;
    end
    A(count6,:)=[1 0 1];
    B(count6,:)=[-2*Hd3(i)];
    P(i,:)=(A\B)';
  end
  
  for i=1:nverts
    BB(1,1)=P(i,1);
    BB(1,2)=P(i,2);
    BB(2,1)=P(i,2);
    BB(2,2)=P(i,3);
    [Vb,D] = eig(BB);
    AV(i,1)=Vb(1,1);
    AV(i,2)=Vb(2,1);
    AV(i,3)=Vb(1,2);
    AV(i,4)=Vb(2,2);
    k1D(i,1)=D(1,1);
    k2D(i,1)=D(2,2);
  end
  
   for i=1:nverts
    ttt=AV(i,1)*b(i,:)+AV(i,2)*c(i,:);
    qqq=AV(i,3)*b(i,:)+AV(i,4)*c(i,:);
    AV2(i,1)=ttt(1,1);
    AV2(i,2)=ttt(1,2);
    AV2(i,3)=ttt(1,3);
    AV2(i,4)=qqq(1,1);
    AV2(i,5)=qqq(1,2);
    AV2(i,6)=qqq(1,3);
  end
  
   %Função que encontra o cruzamento entre dois segmentos
  
  function [flag,x] = EncontraCruzamento(p1,p2,p,d)
    
    tol = 1e-10;  
    flag = 0; % 0 se não se cruzam; 1 caso contrario.
    
    M = [p2(1)-p1(1) -d(1); p2(2)-p1(2) -d(2); p2(3)-p1(3) -d(3)];
    b = [p(1)-p1(1); p(2)-p1(2); p(3)-p1(3)];
    x = (M\b)';
    if(x(1) > 0+tol && x(1) < 1-tol && x(2) > 0+tol && x(2) < 1-tol )
    flag = 1;
    end
  
endfunction

  
cruzat11=zeros(1,nverts); %u1 do cruzamento de t1 com a estrela
cruzat12=zeros(1,nverts); %u2 do cruzamento de t1 com a estrela
cruzamt11=zeros(1,nverts); %u1 do cruzamento de t2 com a estrela
cruzamt12=zeros(1,nverts); %u2 do cruzamento de t2 com a estrela

cruzat21=zeros(1,nverts); %u1 do cruzamento de -t1 com a estrela
cruzat22=zeros(1,nverts); %u2 do cruzamento de -t1 com a estrela
cruzamt21=zeros(1,nverts); %u1 do cruzamento de -t2 com a estrela
cruzamt22=zeros(1,nverts); %u2 do cruzamento de -t2 com a estrela

xcruz1=zeros(nverts,3);
xcruz2=zeros(nverts,3);

for i=1:nverts
 for jj=L(i)-1:-1:1
      [cruz1,xc1]=EncontraCruzamento(V(SS(i,jj),:),V(SS(i,jj+1),:),V(i,:),AV2(i,1:3));
      if(cruz1==1)
        cruzat11(i)=SS(i,jj);
        cruzat12(i)=SS(i,jj+1);
        tt1(i)=xc1(1,1); %t>0  do cruzamento de t1 com a estrela
        %ponto do cruzamento
        xcruz1(i,:)=V(cruzat11(i),:)+tt1(i)*(V(cruzat12(i),:)-V(cruzat11(i),:));
      end
      [cruz2,xc2]=EncontraCruzamento(V(SS(i,jj),:),V(SS(i,jj+1),:),V(i,:),AV2(i,4:6));
      if(cruz2==1)
        cruzamt11(i)=SS(i,jj);
        cruzamt12(i)=SS(i,jj+1);
        tt2(i)=xc2(1,1);  %t>0  do cruzamento de t2 com a estrela
        %ponto do cruzamento
        xcruz2(i,:)=V(cruzamt11(i),:)+tt2(i)*(V(cruzamt12(i),:)-V(cruzamt11(i),:));
      end
    end
    
    %verificação se o cruzamento for entre u6 e u1
    if(cruzat11(i)==0)
      [cruz1,xc1]=EncontraCruzamento(V(SS(i,L(i)),:),V(SS(i,1),:),V(i,:),AV2(i,1:3));
      if(cruz1==1)
        cruzat11(i)=SS(i,L(i));
        cruzat12(i)=SS(i,1);
        tt1(i)=xc1(1,1); %t>0  do cruzamento de t1 com a estrela
        %ponto do cruzamento
        xcruz1(i,:)=V(cruzat11(i),:)+tt1(i)*(V(cruzat12(i),:)-V(cruzat11(i),:));
      end
    end
    if(cruzamt11(i)==0)
      [cruz2,xc2]=EncontraCruzamento(V(SS(i,L(i)),:),V(SS(i,1),:),V(i,:),AV2(i,4:6));
      if(cruz2==1)
        cruzamt11(i)=SS(i,L(i));
        cruzamt12(i)=SS(i,1);
        tt2(i)=xc2(1,1);
        xcruz2(i,:)=V(cruzamt11(i),:)+tt2(i)*(V(cruzamt12(i),:)-V(cruzamt11(i),:));
      end
    end
  end
  
  xcruz21=zeros(nverts,3);
  xcruz22=zeros(nverts,3);
  for i=1:nverts
 for jj=L(i)-1:-1:1
      [cruz1,xc1]=EncontraCruzamento(V(SS(i,jj),:),V(SS(i,jj+1),:),V(i,:),-AV2(i,1:3));
      if(cruz1==1)
        cruzat21(i)=SS(i,jj);
        cruzat22(i)=SS(i,jj+1);
        tt21(i)=xc1(1,1); %t>0  do cruzamento de -t1 com a estrela
        %ponto do cruzamento
        xcruz21(i,:)=V(cruzat21(i),:)+tt21(i)*(V(cruzat22(i),:)-V(cruzat21(i),:));
      end
      [cruz2,xc2]=EncontraCruzamento(V(SS(i,jj),:),V(SS(i,jj+1),:),V(i,:),-AV2(i,4:6));
      if(cruz2==1)
        cruzamt21(i)=SS(i,jj);
        cruzamt22(i)=SS(i,jj+1);
        tt22(i)=xc2(1,1);  %t>0  do cruzamento de -t2 com a estrela
        %ponto do cruzamento
        xcruz223(i,:)=[ujvj(i,jj) ujvj(i,jj+1)]+tt22(i)*([ujvj(i,jj+2) ujvj(i,jj+3)]-[ujvj(i,jj) ujvj(i,jj+1)]);
        xcruz22(i,:)=V(cruzamt21(i),:)+tt22(i)*(V(cruzamt22(i),:)-V(cruzamt21(i),:));
      end
    end
    
    %verificação se o cruzamento for entre u6 e u1
    if(cruzat21(i)==0)
      [cruz1,xc1]=EncontraCruzamento(V(SS(i,L(i)),:),V(SS(i,1),:),V(i,:),-AV2(i,1:3));
      if(cruz1==1)
        cruzat21(i)=SS(i,L(i));
        cruzat22(i)=SS(i,1);
        tt21(i)=xc1(1,1);
        xcruz21(i,:)=V(cruzat21(i),:)+tt21(i)*(V(cruzat22(i),:)-V(cruzat21(i),:));
      end
    end
    if(cruzamt21(i)==0)
      [cruz2,xc2]=EncontraCruzamento(V(SS(i,L(i)),:),V(SS(i,1),:),V(i,:),-AV2(i,4:6));
      if(cruz2==1)
        cruzamt21(i)=SS(i,L(i));
        cruzamt22(i)=SS(i,1);
        tt22(i)=xc2(1,1);
        xcruz22(i,:)=V(cruzamt21(i),:)+tt22(i)*(V(cruzamt22(i),:)-V(cruzamt21(i),:));
      end
    end
  end
##  
##  for i=1:nverts
##    if (cruzat11(i) !=0)
##      if(L(cruzat11(i)) < 6)
##        cruzat11(i)=0;
##      end
##       for k=1:L(i)
##          if(L(Sc(i,k))< 6)
##            cruzat11(i)=0;
##            break
##          end
##       endfor
##     end
##    if (cruzat12(i) !=0)
##      if(L(cruzat12(i)) < 6)
##        cruzat12(i)=0;
##      end
##      for k=1:L(i)
##          if(L(Sc(i,k))< 6)
##            cruzat12(i)=0;
##            break
##          end
##       endfor
##    end
##    if (cruzamt11(i) !=0)
##      if(L(cruzamt11(i)) < 6)
##        cruzamt11(i)=0;
##      end
##      for k=1:L(i)
##          if(L(Sc(i,k))< 6)
##            cruzamt11(i)=0;
##            break
##          end
##       endfor
##    end
##    if (cruzamt12(i) !=0)
##      if(L(cruzamt12(i)) < 6)
##        cruzamt12(i)=0;
##      end
##      for k=1:L(i)
##          if(L(Sc(i,k))< 6)
##            cruzamt12(i)=0;
##            break
##          end
##       endfor
##    end
##    if (cruzat21(i) !=0)
##      if(L(cruzat21(i)) < 6)
##        cruzat21(i)=0;
##      end
##      for k=1:L(i)
##          if(L(Sc(i,k))< 6)
##            cruzat21(i)=0;
##            break
##          end
##       endfor
##    end
##    if (cruzat22(i) !=0)
##      if(L(cruzat22(i)) < 6)
##        cruzat22(i)=0;
##      end
##      for k=1:L(i)
##          if(L(Sc(i,k))< 6)
##            cruzat22(i)=0;
##            break
##          end
##       endfor
##    end
##    if (cruzamt21(i) !=0)
##      if(L(cruzamt21(i)) < 6)
##        cruzamt21(i)=0;
##      end
##      for k=1:L(i)
##          if(L(Sc(i,k))< 6)
##            cruzamt21(i)=0;
##            break
##          end
##       endfor
##    end
##    if (cruzamt22(i) !=0)
##      if(L(cruzamt22(i)) < 6)
##        cruzamt22(i)=0;
##      end
##      for k=1:L(i)
##          if(L(Sc(i,k))< 6)
##            cruzamt22(i)=0;
##            break
##          end
##       endfor
##    end
##  end

  
  printf('Parte 6');
  
  kxcruz1=zeros(nverts,1);
  kxcruz2=zeros(nverts,1);
  kxcruz21=zeros(nverts,1);
  kxcruz22=zeros(nverts,1);
  
  for i=1:nverts
    if(cruzat11(i) == 0 || cruzat12(i) == 0)
      continue
    else
      kxcruz1(i,1)=k13(cruzat11(i))+tt1(i)*(k13(cruzat12(i))-k13(cruzat11(i)));
    end
    if(cruzamt11(i) == 0 || cruzamt12(i) == 0)
      continue
    else  
      kxcruz2(i,1)=k23(cruzamt11(i))+tt2(i)*(k23(cruzamt12(i))-k23(cruzamt11(i)));
    end
    if(cruzat21(i) == 0 || cruzat22(i) == 0)
      continue
    else      
      kxcruz21(i,1)=k13(cruzat21(i))+tt21(i)*(k13(cruzat22(i))-k13(cruzat21(i)));
    end
    if(cruzamt21(i) == 0 || cruzamt22(i) == 0)
      continue
    else
      kxcruz22(i,1)=k23(cruzamt21(i))+tt22(i)*(k23(cruzamt22(i))-k23(cruzamt21(i)));
    endif
  endfor
  
  Ddk1e1=zeros(nverts,1);
  Ddk2e2=zeros(nverts,1);
  Positivo=[];
  Negativo=[];
  Positivo2=[];
  Negativo2=[];
  count55=1;
  count56=1;
  count555=1;
  count565=1;
  for i=1:nverts
    if(abs(xcruz1(i,:)-xcruz21(i,:)) < 1e-10 || abs(xcruz2(i,:)-xcruz22(i,:)) < 1e-10)
      continue
    else
      Ddk1e1(i,1)=(kxcruz1(i,1)-kxcruz21(i,1))/(dot(xcruz1(i,:)-xcruz21(i,:),[AV2(i,1) AV2(i,2) AV2(i,3)]/norm([AV2(i,1) AV2(i,2) AV2(i,3)])));
      Ddk2e2(i,1)=(kxcruz2(i,1)-kxcruz22(i,1))/(dot(xcruz2(i,:)-xcruz22(i,:),[AV2(i,4) AV2(i,5) AV2(i,6)]/norm([AV2(i,4) AV2(i,5) AV2(i,6)])));
    end
    if(Ddk1e1(i,1)>1e-10)
      Positivo(count55,:)=V(i,:);
      count55=count55+1;
    end
    if (Ddk1e1(i,1)<-1e-10)
      Negativo(count56,:)=V(i,:);
      count56=count56+1;
    end
    if(Ddk2e2(i,1)>1e-10)
      Positivo2(count555,:)=V(i,:);
      count555=count555+1;
    end
    if (Ddk2e2(i,1)<-1e-10)
      Negativo2(count565,:)=V(i,:);
      count565=count565+1;
    end
  endfor
##  Ddk1e1
    count57=1;
    count572=1;
  for i=1:nfaces
    if ((Ddk1e1(F(i,1),1)>1e-10 && Ddk1e1(F(i,2),1)<-1e-10))
      t=Ddk1e1(F(i,1),1)/(Ddk1e1(F(i,1),1)-Ddk1e1(F(i,2),1));
      MM(count57,:)=V(F(i,1),:)+t*(V(F(i,2),:)-V(F(i,1),:));
      PMM(count57,1)=F(i,1);
      PMM(count57,2)=F(i,2);
      count57=count57+1;
    end
    if ((Ddk1e1(F(i,2),1)>1e-10 && Ddk1e1(F(i,3),1)<-1e-10))
      t=Ddk1e1(F(i,2),1)/(Ddk1e1(F(i,2),1)-Ddk1e1(F(i,3),1));
      MM(count57,:)=V(F(i,2),:)+t*(V(F(i,3),:)-V(F(i,2),:));
      PMM(count57,1)=F(i,2);
      PMM(count57,2)=F(i,3);
      count57=count57+1;
    end
    if ((Ddk1e1(F(i,3),1)>1e-10 && Ddk1e1(F(i,1),1)<-1e-10))
      t=Ddk1e1(F(i,3),1)/(Ddk1e1(F(i,3),1)-Ddk1e1(F(i,1),1));
      MM(count57,:)=V(F(i,3),:)+t*(V(F(i,1),:)-V(F(i,3),:));
      PMM(count57,1)=F(i,3);
      PMM(count57,2)=F(i,1);
      count57=count57+1;
    endif
  end
  
  for i=1:nfaces
    if ((Ddk2e2(F(i,1),1)>1e-10 && Ddk2e2(F(i,2),1)<-1e-10))
      t2=Ddk2e2(F(i,1),1)/(Ddk2e2(F(i,1),1)-Ddk2e2(F(i,2),1));
      MMw(count572,:)=V(F(i,1),:)+t2*(V(F(i,2),:)-V(F(i,1),:));
      PMMw(count572,1)=F(i,1);
      PMMw(count572,2)=F(i,2);
      count572=count572+1;
    end
    if ((Ddk2e2(F(i,2),1)>1e-10 && Ddk2e2(F(i,3),1)<-1e-10))
      t2=Ddk2e2(F(i,2),1)/(Ddk2e2(F(i,2),1)-Ddk2e2(F(i,3),1));
      MMw(count572,:)=V(F(i,2),:)+t2*(V(F(i,3),:)-V(F(i,2),:));
      PMMw(count572,1)=F(i,2);
      PMMw(count572,2)=F(i,3);
      count572=count572+1;
    end
    if ((Ddk2e2(F(i,3),1)>1e-10 && Ddk2e2(F(i,1),1)<-1e-10))
      t2=Ddk2e2(F(i,3),1)/(Ddk2e2(F(i,3),1)-Ddk2e2(F(i,1),1));
      MMw(count572,:)=V(F(i,3),:)+t2*(V(F(i,1),:)-V(F(i,3),:));
      PMMw(count572,1)=F(i,3);
      PMMw(count572,2)=F(i,1);
      count572=count572+1;
    endif
  end
  
##  MMMMM=MMw1;
  uuwx=0;
  for i=1:size(PMM1,1)
    for j=i:size(PMM1,1)
      if (PMM1(i,1)==PMM1(j,1) || PMM1(i,1)==PMM1(j,2))
        MMw1(i,:)=V(PMM1(i,1),:);
##        umb(uuwx+1,:)=V(PMM1(i,1),:);
##        uuwx=uuwx+1;
      endif
      if (PMM1(i,2)==PMM1(j,1) ||PMM1(i,2)==PMM1(j,2))
        MMw1(i,:)=V(PMM1(i,2),:);
##        umb(uuwx+1,:)=V(PMM1(i,2),:);
##        uuwx=uuwx+1;
      endif
    endfor
  endfor  
  

  
  for i=1:size(MM)
    for j=1:size(MMw1)
      if((PMM(i,1)==PMM1(j,1) & PMM(i,2)==PMM1(j,2)) || (PMM(i,1)==PMM1(j,2) & PMM(i,2)==PMM1(j,1)))
        MM(i,:)=MMw1(j,:);
      end
    endfor
  endfor
  
  for i=1:size(MMw)
    for j=1:size(MMw1)
      if((PMMw(i,1)==PMM1(j,1) & PMMw(i,2)==PMM1(j,2)) || (PMMw(i,1)==PMM1(j,2) & PMMw(i,2)==PMM1(j,1)))
        MMw(i,:)=MMw1(j,:);
      end
    endfor
  endfor
  

  
uuw=0;
uuwx=0;
for i=1:size(MM,1)
    for j=1:size(MMw,1)
        if((MM(i,1)==MMw(j,1) && MM(i,2)==MMw(j,2)))
          posumbw(uuw+1,:)=MM(i,:);
          uuw=uuw+1;
          for k=1:size(MMw1)
            if(posumbw(uuw,:)==MMw1(k,:))
              umb(uuwx+1,:)=posumbw(uuw,:);
              uuwx=uuwx+1;
            end
          endfor
        end
    end
endfor

  umb=unique(umb,"rows");
  aa=0;
  for i=1:nverts
    for j=1:size(umb,1)
      if (V(i,:)==umb(j,:))
        Vumb(aa+1,:)=i;
        aa=aa+1;
      endif
    endfor
  endfor

MM3=[MM PMM];
MM3w=[MMw PMMw];
  
ct=1;
ct2=1;
zero=[];
zero2=[];
##for i=1:nverts
##  if(abs(Ddk1e1(i))<=1e-10 && L(i)==6)
##    zero(ct,1)=i;
##    ct=ct+1;
##  end
##  if(abs(Ddk2e2(i))<=1e-10 && L(i)==6)
##    zero2(ct2,1)=i;
##    ct2=ct2+1;
##  end
##endfor

for i=1:nverts
  if(abs(Ddk1e1(i))<=1e-10)
    zero(ct,1)=i;
    ct=ct+1;
  end
  if(abs(Ddk2e2(i))<=1e-10)
    zero2(ct2,1)=i;
    ct2=ct2+1;
  end
endfor

guardar=[];
for i=1:size(zero,1)
  ccccc=1;
  guardar(i,ccccc)=zero(i,1);
  for j=1:L(i)
    for k=1:size(zero,1)
      if (SS(zero(i),j)==zero(k))
        guardar(i,ccccc+1)=zero(k);
        ccccc=ccccc+1;
        break
        endif
    endfor
    if (ccccc>1)
      break
    endif
  endfor
  idx=0;
  while (idx!=1)
    ccc=ccccc; 
   if (j+1<=L(i))  
    for k=1:size(zero,1)
      if (SS(zero(i),j+1)==zero(k))
      guardar(i,ccccc+1)=zero(k);
      ccccc=ccccc+1;
      j=j+1;
      break
      endif
    endfor
   endif
   if(ccc==ccccc)
      idx=1;
    else
      ccc=ccccc;
    endif
  endwhile
  nguardar(i,1)=ccccc;
end

guardar2=[];
for i=1:size(zero2,1)
  ccccc=1;
  guardar2(i,ccccc)=zero2(i,1);
  for j=1:L(i)
    for k=1:size(zero2,1)
      if (SS(zero2(i),j)==zero2(k))
        guardar2(i,ccccc+1)=zero2(k);
        ccccc=ccccc+1;
        break
      endif
    endfor
    if (ccccc>1)
      break
    endif
  endfor
  idx=0;
  while (idx!=1)
    ccc=ccccc;
    if (j+1<=L(i))  
     for k=1:size(zero2,1)
      if (SS(zero2(i),j+1)==zero2(k))
      guardar2(i,ccccc+1)=zero2(k);
      ccccc=ccccc+1;
      j=j+1;
      break
      endif
     endfor
    endif
    if(ccc==ccccc)
      idx=1;
    else
      ccc=ccccc;
    endif
  endwhile
  nguardar2(i,1)=ccccc;
end
  
  printf('Parte 7');

function [A] = deletarColuna(Matrix,coluna)
  
  if (coluna == 1)
    A = Matrix(:,2:size(Matrix)(2));
    return;
  end
  
  if (coluna == size(Matrix)(2))
    A = Matrix(:,1:size(Matrix)(2)-1);
    return;
  end
  
  A1 = Matrix(:,1:coluna-1);
  A2 = Matrix(:,coluna+1:size(Matrix)(2));

  A = [A1  A2];
    
end
  
function [A] = deletarLinha(Matrix,linha)
  
  if (linha == 1)
    A = Matrix(2:size(Matrix)(1),:);
    return;
  end
  
  if (linha == size(Matrix)(1))
    A = Matrix(1:size(Matrix)(1)-1,:);
    return;
  end
  
  A1 = Matrix(1:linha-1,:);
  A2 = Matrix(linha+1:size(Matrix)(1),:);

  A = [A1 ; A2];
    
end

##for i=1:size(MM3)(1)
##  nn=size(MM3)(1);
##  if (i>=nn)
##    break
##  elseif ((abs(MM3(i,1)-MM3(i+1,1))<10e-6) && (abs(MM3(i,2)-MM3(i+1,2))<10e-6))
##      MM3=deletarLinha(MM3,i+1);
##    endif
##end

##for i=1:size(MM3)(1)
##  nn=size(MM3)(1);
##  if (i>=nn)
##    break
##  elseif ((L(MM3(i,4))< 6) || (L(MM3(i,5))< 6))
##      MM3=deletarLinha(MM3,i);
##    endif
##end

##for i=1:size(MM3w)(1)
##  nnw=size(MM3w)(1);
##  if (i>=nnw)
##    break
##  elseif ((abs(MM3w(i,1)-MM3w(i+1,1))<10e-6) && (abs(MM3w(i,2)-MM3w(i+1,2))<10e-6))
##      MM3w=deletarLinha(MM3w,i+1);
##    endif
##end

##for i=1:size(MM3w)(1)
##  nnw=size(MM3w)(1);
##  if (i>=nnw)
##    break
##  elseif ((L(MM3w(i,4))< 6) || (L(MM3w(i,5))< 6))
##      MM3w=deletarLinha(MM3w,i);
##    endif
##end

MM2w=[];
jafoiw=zeros(size(MM3w)(1),1);
kkw=1;
for iiw=1:size(MM3w)(1)
  if (jafoiw(iiw)==0)
    MM2w(1,3*kkw-2:3*kkw)=MM3w(iiw,1:3);
    PMM2w(1,2*kkw-1:2*kkw)=MM3w(iiw,4:5);
ccw=1;
llw=iiw;
jjjw=0;

jafoiw(iiw)=1;
while (llw != jjjw)
  jjjw=llw;
  for i=1:nfaces
    if (MM3w(llw,4)==F(i,1) && MM3w(llw,5)==F(i,2))
      for k=1:size(MM3w)(1)
        if(jafoiw(k)==0)
        if ((MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3)) || (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1)) || (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2)) || (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3)))
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          llw=k;
          jafoiw(k)=1;
          ccw=ccw+1;
          MM2w(ccw,3*kkw-2:3*kkw)=MM3w(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MM3w(llw,4)==F(i,1) && MM3w(llw,5)==F(i,2))
      for k=1:size(MM3w)(1)
        if(jafoiw(k)==0)
        if ((MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2)) || (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1)) || (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1)) || (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3)))
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          llw=k;
          jafoiw(k)=1;
          ccw=ccw+1;
          MM2w(ccw,3*kkw-2:3*kkw)=MM3w(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MM3w(llw,4)==F(i,2) && MM3w(llw,5)==F(i,3))
      for k=1:size(MM3w)(1)
        if(jafoiw(k)==0)
        if ((MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1)) || (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2)) || (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3)) || (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1)))
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          llw=k;
          jafoiw(k)=1;
          ccw=ccw+1;
          MM2w(ccw,3*kkw-2:3*kkw)=MM3w(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MM3w(llw,4)==F(i,2) && MM3w(llw,5)==F(i,3))
      for k=1:size(MM3w)(1)
        if(jafoiw(k)==0)
        if ((MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3)) || (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2)) || (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2)) || (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1)))
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          llw=k;
          jafoiw(k)=1;
          ccw=ccw+1;
          MM2w(ccw,3*kkw-2:3*kkw)=MM3w(k,1:3);
          break
      endif
      end
    endfor
  end
  if (MM3w(llw,4)==F(i,3) && MM3w(llw,5)==F(i,2))
      for k=1:size(MM3w)(1)
        if(jafoiw(k)==0)
        if ((MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1)) || (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2)) || (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3)) || (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1)))
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          llw=k;
          jafoiw(k)=1;
          ccw=ccw+1;
          MM2w(ccw,3*kkw-2:3*kkw)=MM3w(k,1:3);
          break
      endif
      end
    endfor
  end
  if (MM3w(llw,4)==F(i,2) && MM3w(llw,5)==F(i,1))
      for k=1:size(MM3w)(1)
        if(jafoiw(k)==0)
        if ((MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3)) || (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1)) || (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2)) || (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3)))
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          llw=k;
          jafoiw(k)=1;
          ccw=ccw+1;
          MM2w(ccw,3*kkw-2:3*kkw)=MM3w(k,1:3);
          break
      endif
      end
    endfor
  end

      if (MM3w(llw,4)==F(i,3) && MM3w(llw,5)==F(i,1))
      for k=1:size(MM3w)(1)
        if(jafoiw(k)==0)
        if ((MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2)) || (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3)) || (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1)) || (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2)))
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          llw=k;
          jafoiw(k)=1;
          ccw=ccw+1;
          MM2w(ccw,3*kkw-2:3*kkw)=MM3w(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MM3w(llw,4)==F(i,3) && MM3w(llw,5)==F(i,1))
      for k=1:size(MM3w)(1)
        if(jafoiw(k)==0)
        if ((MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3)) || (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1)) || (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2)) || (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3)))
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          llw=k;
          jafoiw(k)=1;
          ccw=ccw+1;
          MM2w(ccw,3*kkw-2:3*kkw)=MM3(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MM3w(llw,4)==F(i,1) && MM3w(llw,5)==F(i,3))
      for k=1:size(MM3w)(1)
        if(jafoiw(k)==0)
        if ((MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2)) || (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3)) || (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1)) || (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2)))
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          llw=k;
          jafoiw(k)=1;
          ccw=ccw+1;
          MM2w(ccw,3*kkw-2:3*kkw)=MM3w(k,1:3);
          break
      endif
      end
    endfor
  end
  if (MM3w(llw,4)==F(i,1) && MM3w(llw,5)==F(i,3))
      for k=1:size(MM3w)(1)
        if(jafoiw(k)==0)
        if ((MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2)) || (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1)) || (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1)) || (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3)))
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,2))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,2);
          endif
          if (MM3w(k,4)==F(i,3) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,3);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          if (MM3w(k,4)==F(i,2) && MM3w(k,5)==F(i,1))
            PMM2w(ccw+1,2*kkw-1)=F(i,2);
            PMM2w(ccw+1,2*kkw)=F(i,1);
          endif
          if (MM3w(k,4)==F(i,1) && MM3w(k,5)==F(i,3))
            PMM2w(ccw+1,2*kkw-1)=F(i,1);
            PMM2w(ccw+1,2*kkw)=F(i,3);
          endif
          llw=k;
          jafoiw(k)=1;
          ccw=ccw+1;
          MM2w(ccw,3*kkw-2:3*kkw)=MM3w(k,1:3);
          break
      endif
      end
    endfor
    end

endfor
end
pontocurvasw(kkw,1)=ccw;
kkw=kkw+1;
endif
endfor

MM2=[];
jafoi=zeros(size(MM3)(1),1);
kk=1;
for ii=1:size(MM3)(1)
  if (jafoi(ii)==0)
    MM2(1,3*kk-2:3*kk)=MM3(ii,1:3);
    PMM2(1,2*kk-1:2*kk)=MM3(ii,4:5);
cc=1;
ll=ii;
jjj=0;

jafoi(ii)=1;
while (ll != jjj)
  jjj=ll;
  for i=1:nfaces
    if (MM3(ll,4)==F(i,1) && MM3(ll,5)==F(i,2))
      for k=1:size(MM3)(1)
        if(jafoi(k)==0)
        if ((MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3)) || (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1)) || (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2)) || (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3)))
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          ll=k;
          jafoi(k)=1;
          cc=cc+1;
          MM2(cc,3*kk-2:3*kk)=MM3(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MM3(ll,4)==F(i,1) && MM3(ll,5)==F(i,2))
      for k=1:size(MM3)(1)
        if(jafoi(k)==0)
        if ((MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2)) || (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1)) || (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1)) || (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3)))
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          ll=k;
          jafoi(k)=1;
          cc=cc+1;
          MM2(cc,3*kk-2:3*kk)=MM3(k,1:3);
          break
      endif
      end
    endfor
    end

      if (MM3(ll,4)==F(i,2) && MM3(ll,5)==F(i,3))
      for k=1:size(MM3)(1)
        if(jafoi(k)==0)
        if ((MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1)) || (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2)) || (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3)) || (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1)))
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          ll=k;
          jafoi(k)=1;
          cc=cc+1;
          MM2(cc,3*kk-2:3*kk)=MM3(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MM3(ll,4)==F(i,2) && MM3(ll,5)==F(i,3))
      for k=1:size(MM3)(1)
        if(jafoi(k)==0)
        if ((MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3)) || (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2)) || (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2)) || (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1)))
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          ll=k;
          jafoi(k)=1;
          cc=cc+1;
          MM2(cc,3*kk-2:3*kk)=MM3(k,1:3);
          break
      endif
      end
    endfor
  end
  if (MM3(ll,4)==F(i,3) && MM3(ll,5)==F(i,2))
      for k=1:size(MM3)(1)
        if(jafoi(k)==0)
        if ((MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1)) || (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2)) || (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3)) || (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1)))
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          ll=k;
          jafoi(k)=1;
          cc=cc+1;
          MM2(cc,3*kk-2:3*kk)=MM3(k,1:3);
          break
      endif
      end
    endfor
  end
  if (MM3(ll,4)==F(i,2) && MM3(ll,5)==F(i,1))
      for k=1:size(MM3)(1)
        if(jafoi(k)==0)
        if ((MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3)) || (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1)) || (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2)) || (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3)))
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          ll=k;
          jafoi(k)=1;
          cc=cc+1;
          MM2(cc,3*kk-2:3*kk)=MM3(k,1:3);
          break
      endif
      end
    endfor
  end

      if (MM3(ll,4)==F(i,3) && MM3(ll,5)==F(i,1))
      for k=1:size(MM3)(1)
        if(jafoi(k)==0)
        if ((MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2)) || (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3)) || (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1)) || (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2)))
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          ll=k;
          jafoi(k)=1;
          cc=cc+1;
          MM2(cc,3*kk-2:3*kk)=MM3(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MM3(ll,4)==F(i,3) && MM3(ll,5)==F(i,1))
      for k=1:size(MM3)(1)
        if(jafoi(k)==0)
        if ((MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3)) || (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1)) || (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2)) || (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3)))
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          ll=k;
          jafoi(k)=1;
          cc=cc+1;
          MM2(cc,3*kk-2:3*kk)=MM3(k,1:3);
          break
      endif
      end
    endfor
    end
      if (MM3(ll,4)==F(i,1) && MM3(ll,5)==F(i,3))
      for k=1:size(MM3)(1)
        if(jafoi(k)==0)
        if ((MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2)) || (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3)) || (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1)) || (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2)))
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          ll=k;
          jafoi(k)=1;
          cc=cc+1;
          MM2(cc,3*kk-2:3*kk)=MM3(k,1:3);
          break
      endif
      end
    endfor
  end
  if (MM3(ll,4)==F(i,1) && MM3(ll,5)==F(i,3))
      for k=1:size(MM3)(1)
        if(jafoi(k)==0)
        if ((MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2)) || (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1)) || (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1)) || (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3)))
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,2))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,2);
          endif
          if (MM3(k,4)==F(i,3) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,3);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          if (MM3(k,4)==F(i,2) && MM3(k,5)==F(i,1))
            PMM2(cc+1,2*kk-1)=F(i,2);
            PMM2(cc+1,2*kk)=F(i,1);
          endif
          if (MM3(k,4)==F(i,1) && MM3(k,5)==F(i,3))
            PMM2(cc+1,2*kk-1)=F(i,1);
            PMM2(cc+1,2*kk)=F(i,3);
          endif
          ll=k;
          jafoi(k)=1;
          cc=cc+1;
          MM2(cc,3*kk-2:3*kk)=MM3(k,1:3);
          break
      endif
      end
    endfor
    end

endfor
end
pontocurvas(kk,1)=cc;
kk=kk+1;
endif
endfor


###### Vermelha #####

idx2=0;
while(idx2!=1)
pontocurvasax=pontocurvas;
#1
idx=0;
while(idx!=1)
pontocurvasa=pontocurvas;
for tt=1:2:size(PMM2,2)
  ppt=size(PMM2,2);
  if (tt>ppt)
    break
  endif
for pp=tt+2:2:size(PMM2,2)
  pptz=size(PMM2,2);
  if (pp>pptz)
    break
  endif
for i=1:nfaces
  if((PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2) && PMM2(1,pp)==F(i,2) && PMM2(1,pp+1)==F(i,3)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2) && PMM2(1,pp+1)==F(i,2) && PMM2(1,pp)==F(i,3)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2) && PMM2(1,pp)==F(i,1) && PMM2(1,pp+1)==F(i,3)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2) && PMM2(1,pp+1)==F(i,1) && PMM2(1,pp)==F(i,3)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(pontocurvas((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(pontocurvas((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
  elseif((PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3) && PMM2(1,pp)==F(i,3) && PMM2(1,pp+1)==F(i,1)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3) && PMM2(1,pp+1)==F(i,3) && PMM2(1,pp)==F(i,1)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3) && PMM2(1,pp)==F(i,2) && PMM2(1,pp+1)==F(i,1)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3) && PMM2(1,pp+1)==F(i,2) && PMM2(1,pp)==F(i,1)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(pontocurvas((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(pontocurvas((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
  elseif((PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1) && PMM2(1,pp)==F(i,1) && PMM2(1,pp)==F(i,2)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1) && PMM2(1,pp+1)==F(i,1) && PMM2(1,pp)==F(i,2)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1) && PMM2(1,pp)==F(i,3) && PMM2(1,pp)==F(i,2)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1) && PMM2(1,pp+1)==F(i,3) && PMM2(1,pp)==F(i,2)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(pontocurvas((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(pontocurvas((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
  elseif((PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1) && PMM2(1,pp)==F(i,1) && PMM2(1,pp+1)==F(i,3)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1) && PMM2(1,pp+1)==F(i,1) && PMM2(1,pp)==F(i,3)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1) && PMM2(1,pp)==F(i,2) && PMM2(1,pp+1)==F(i,3)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1) && PMM2(1,pp+1)==F(i,2) && PMM2(1,pp)==F(i,3)))
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(pontocurvas((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(pontocurvas((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
  elseif((PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2) && PMM2(1,pp)==F(i,2) && PMM2(1,pp+1)==F(i,1)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2) && PMM2(1,pp+1)==F(i,2) && PMM2(1,pp)==F(i,1)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2) && PMM2(1,pp)==F(i,3) && PMM2(1,pp+1)==F(i,1)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2) && PMM2(1,pp+1)==F(i,3) && PMM2(1,pp)==F(i,1)))
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(pontocurvas((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(pontocurvas((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
  elseif((PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3) && PMM2(1,pp)==F(i,3) && PMM2(1,pp+1)==F(i,2)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3) && PMM2(1,pp+1)==F(i,3) && PMM2(1,pp)==F(i,2)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3) && PMM2(1,pp)==F(i,1) && PMM2(1,pp+1)==F(i,2)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3) && PMM2(1,pp+1)==F(i,1) && PMM2(1,pp)==F(i,2)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(pontocurvas((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(pontocurvas((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
endif
endfor
endfor
endfor
if(size(pontocurvasa,1)==size(pontocurvas,1))
idx=1;
end
endwhile

printf('Vermelha 1');

#2
idx=0;
while(idx!=1)
pontocurvasa=pontocurvas;
for tt=1:2:size(PMM2,2)
  ppt=size(PMM2,2);
  if (tt>ppt)
    break
  endif
for pp=tt+2:2:size(PMM2,2)
  pptz=size(PMM2,2);
  if (pp>pptz)
    break
  endif
for i=1:nfaces
  if((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt+1)==F(i,2) && PMM2(1,tt)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt+1)==F(i,1) && PMM2(1,tt)==F(i,3)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
  elseif((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt+1)==F(i,3) && PMM2(1,tt)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt+1)==F(i,2) && PMM2(1,tt)==F(i,1)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
  elseif((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt+1)==F(i,1) && PMM2(1,tt)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt+1)==F(i,3) && PMM2(1,tt)==F(i,2)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
  elseif((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt+1)==F(i,1) && PMM2(1,tt)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt+1)==F(i,2) && PMM2(1,tt)==F(i,3)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
  elseif((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt+1)==F(i,2) && PMM2(1,tt)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt+1)==F(i,3) && PMM2(1,tt)==F(i,1)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
  elseif((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt+1)==F(i,3) && PMM2(1,tt)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt+1)==F(i,1) && PMM2(1,tt)==F(i,2)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  MM2=deletarColuna(MM2,3*(pp+1)/2-5);
  PMM2=deletarColuna(PMM2,pp-2);
  PMM2=deletarColuna(PMM2,pp-2);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2-1);
  break
endif
endfor
endfor
endfor
if(size(pontocurvasa,1)==size(pontocurvas,1))
idx=1;
end
endwhile

printf('Vermelha 2');

#3
idx=0;
while(idx!=1)
pontocurvasa=pontocurvas;
for pp=1:2:size(PMM2,2)
  ppt=size(PMM2,2);
  if (pp>ppt)
    break
  endif
for tt=pp+2:2:size(PMM2,2)
  pptz=size(PMM2,2);
  if (tt>pptz)
    break
  endif
for i=1:nfaces
  if((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt+1)==F(i,2) && PMM2(1,tt)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt+1)==F(i,1) && PMM2(1,tt)==F(i,3)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  PMM2=deletarColuna(PMM2,pp);
  PMM2=deletarColuna(PMM2,pp);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  PMM2=deletarColuna(PMM2,tt-2);
  PMM2=deletarColuna(PMM2,tt-2);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2-1);
  break
  elseif((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt+1)==F(i,3) && PMM2(1,tt)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt+1)==F(i,2) && PMM2(1,tt)==F(i,1)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  PMM2=deletarColuna(PMM2,pp);
  PMM2=deletarColuna(PMM2,pp);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  PMM2=deletarColuna(PMM2,tt-2);
  PMM2=deletarColuna(PMM2,tt-2);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2-1);
  break
  elseif((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt+1)==F(i,1) && PMM2(1,tt)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt+1)==F(i,3) && PMM2(1,tt)==F(i,2)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  PMM2=deletarColuna(PMM2,pp);
  PMM2=deletarColuna(PMM2,pp);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  PMM2=deletarColuna(PMM2,tt-2);
  PMM2=deletarColuna(PMM2,tt-2);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2-1);
  break
  elseif((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt+1)==F(i,1) && PMM2(1,tt)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,2) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,1) && PMM2(1,tt+1)==F(i,2) && PMM2(1,tt)==F(i,3)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  PMM2=deletarColuna(PMM2,pp);
  PMM2=deletarColuna(PMM2,pp);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  PMM2=deletarColuna(PMM2,tt-2);
  PMM2=deletarColuna(PMM2,tt-2);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2-1);
  break
  elseif((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt+1)==F(i,2) && PMM2(1,tt)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,3) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,2) && PMM2(1,tt+1)==F(i,3) && PMM2(1,tt)==F(i,1)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  PMM2=deletarColuna(PMM2,pp);
  PMM2=deletarColuna(PMM2,pp);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  PMM2=deletarColuna(PMM2,tt-2);
  PMM2=deletarColuna(PMM2,tt-2);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2-1);
  break
  elseif((PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt+1)==F(i,3) && PMM2(1,tt)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2)) || (PMM2(pontocurvas((pp+1)/2,1),pp)==F(i,1) && PMM2(pontocurvas((pp+1)/2,1),pp+1)==F(i,3) && PMM2(1,tt+1)==F(i,1) && PMM2(1,tt)==F(i,2)) )
  MM2(1:pontocurvas((pp+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:pontocurvas((pp+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((pp+1)/2,1),pp:(pp+1));
  PMM2(pontocurvas((pp+1)/2,1)+1:pontocurvas((pp+1)/2,1)+pontocurvas((tt+1)/2,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+pontocurvas((pp+1)/2,1);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  MM2=deletarColuna(MM2,3*(pp+1)/2-2);
  PMM2=deletarColuna(PMM2,pp);
  PMM2=deletarColuna(PMM2,pp);
  pontocurvas=deletarLinha(pontocurvas,(pp+1)/2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  MM2=deletarColuna(MM2,3*(tt+1)/2-5);
  PMM2=deletarColuna(PMM2,tt-2);
  PMM2=deletarColuna(PMM2,tt-2);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2-1);
  break
endif
endfor
endfor
endfor
if(size(pontocurvasa,1)==size(pontocurvas,1))
idx=1;
end
endwhile

printf('Vermelha 3');

#4
idx=0;
while(idx!=1)
mx=0;
pontocurvasa=pontocurvas;
for tt=1:2:size(PMM2,2)
for pp=1:size(zero,1)
for i=1:nfaces
  if((PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,2) && zero(pp,1)==F(i,3)))
  MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-2:end)=V(guardar(pp,1:nguardar(pp,1)),:);
  PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1)=guardar(pp,1:nguardar(pp,1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end)=guardar(pp,1:nguardar(pp,1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,3) && zero(pp,1)==F(i,1)))
  MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-2:end)=V(guardar(pp,1:nguardar(pp,1)),:);
  PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1)=guardar(pp,1:nguardar(pp,1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end)=guardar(pp,1:nguardar(pp,1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,1) && zero(pp,1)==F(i,2)))
  MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-2:end)=V(guardar(pp,1:nguardar(pp,1)),:);
  PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1)=guardar(pp,1:nguardar(pp,1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end)=guardar(pp,1:nguardar(pp,1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,1) && zero(pp,1)==F(i,3)))
  MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-2:end)=V(guardar(pp,1:nguardar(pp,1)),:);
  PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1)=guardar(pp,1:nguardar(pp,1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end)=guardar(pp,1:nguardar(pp,1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,2) && zero(pp,1)==F(i,1)))
  MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-2:end)=V(guardar(pp,1:nguardar(pp,1)),:);
  PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1)=guardar(pp,1:nguardar(pp,1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end)=guardar(pp,1:nguardar(pp,1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,3) && zero(pp,1)==F(i,2)))
  MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-2:end)=V(guardar(pp,1:nguardar(pp,1)),:);
  PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1)=guardar(pp,1:nguardar(pp,1));
  PMM2(pontocurvas((tt+1)/2,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end)=guardar(pp,1:nguardar(pp,1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
endif
endfor
endfor
endfor
if(mx==0)
idx=1;
end
endwhile

printf('Vermelha 4');

##5
idx=0;
while(idx!=1)
mx=0;
pontocurvasa=pontocurvas;
for tt=1:2:size(PMM2,2)
for pp=1:size(zero,1)
for i=1:nfaces
  if((PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2) && zero(pp,1)==F(i,3)) )
  MM2(1:nguardar(pp,1),end+1:end+3)=V(guardar(pp,nguardar(pp,1):-1:1),:);
  MM2(nguardar(pp,1)+1:nguardar(pp,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(nguardar(pp,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3) && zero(pp,1)==F(i,1)))
  MM2(1:nguardar(pp,1),end+1:end+3)=V(guardar(pp,nguardar(pp,1):-1:1),:);
  MM2(nguardar(pp,1)+1:nguardar(pp,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(nguardar(pp,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1) && zero(pp,1)==F(i,2)))
  MM2(1:nguardar(pp,1),end+1:end+3)=V(guardar(pp,nguardar(pp,1):-1:1),:);
  MM2(nguardar(pp,1)+1:nguardar(pp,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(nguardar(pp,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1) && zero(pp,1)==F(i,3)))
  MM2(1:nguardar(pp,1),end+1:end+3)=V(guardar(pp,nguardar(pp,1):-1:1),:);
  MM2(nguardar(pp,1)+1:nguardar(pp,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(nguardar(pp,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2) && zero(pp,1)==F(i,1)))
  MM2(1:nguardar(pp,1),end+1:end+3)=V(guardar(pp,nguardar(pp,1):-1:1),:);
  MM2(nguardar(pp,1)+1:nguardar(pp,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(nguardar(pp,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3) && zero(pp,1)==F(i,2)) )
  MM2(1:nguardar(pp,1),end+1:end+3)=V(guardar(pp,nguardar(pp,1):-1:1),:);
  MM2(nguardar(pp,1)+1:nguardar(pp,1)+pontocurvas((tt+1)/2,1),end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(1:nguardar(pp,1),end+1)=guardar(pp,nguardar(pp,1):-1:1);
  PMM2(nguardar(pp,1)+1:pontocurvas((tt+1)/2,1)+nguardar(pp,1),end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
  pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+nguardar(pp,1);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  MM2=deletarColuna(MM2,3*(tt+1)/2-2);
  PMM2=deletarColuna(PMM2,tt);
  PMM2=deletarColuna(PMM2,tt);
  pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
  mx=1;
  break
endif
endfor
endfor
endfor
if(mx==0)
idx=1;
end
endwhile

printf('Vermelha 5');

if(size(pontocurvasax,1)==size(pontocurvas,1))
idx2=1;
end
endwhile

for tt=1:2:size(PMM2,2)
  if (PMM2(1,tt)==PMM2(1,tt+1))
    for pp=1:size(Vumb,1)
      for i=1:nfaces
        if((PMM2(1,tt)==F(i,1) && Vumb(pp,1)==F(i,2))||(PMM2(1,tt)==F(i,2) && Vumb(pp,1)==F(i,3))||(PMM2(1,tt)==F(i,3) && Vumb(pp,1)==F(i,1))||(PMM2(1,tt)==F(i,1) && Vumb(pp,1)==F(i,3))||(PMM2(1,tt)==F(i,3) && Vumb(pp,1)==F(i,2))||(PMM2(1,tt)==F(i,2) && Vumb(pp,1)==F(i,1)))
          MM2(1,end+1:end+3)=V(Vumb(pp,1),:);
          MM2(2:pontocurvas((tt+1)/2,1)+1,end-2:end)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
          PMM2(1,end+1)=Vumb(pp,1);
          PMM2(1,end+1)=Vumb(pp,1);
          PMM2(2:pontocurvas((tt+1)/2,1)+1,end-1:end)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
          pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
          MM2=deletarColuna(MM2,3*(tt+1)/2-2);
          MM2=deletarColuna(MM2,3*(tt+1)/2-2);
          MM2=deletarColuna(MM2,3*(tt+1)/2-2);
          PMM2=deletarColuna(PMM2,tt);
          PMM2=deletarColuna(PMM2,tt);
          pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
          break
        endif
        if((PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,1) && Vumb(pp,1)==F(i,2))||(PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,2) && Vumb(pp,1)==F(i,3))||(PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,3) && Vumb(pp,1)==F(i,1))||(PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,1) && Vumb(pp,1)==F(i,3))||(PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,3) && Vumb(pp,1)==F(i,2))||(PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,2) && Vumb(pp,1)==F(i,1)))
          MM2(1:pontocurvas((tt+1)/2,1),end+1:end+3)=MM2(1:pontocurvas((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
          MM2(1+pontocurvas((tt+1)/2,1),end-2:end)=V(Vumb(pp,1),:);
          PMM2(1:pontocurvas((tt+1)/2,1),end+1:end+2)=PMM2(1:pontocurvas((tt+1)/2,1),tt:(tt+1));
          PMM2(pontocurvas((tt+1)/2,1)+1,end-1)=Vumb(pp,1);
          PMM2(pontocurvas((tt+1)/2,1)+1,end)=Vumb(pp,1);
          pontocurvas(end+1,1)=pontocurvas((tt+1)/2,1)+1;
          MM2=deletarColuna(MM2,3*(tt+1)/2-2);
          MM2=deletarColuna(MM2,3*(tt+1)/2-2);
          MM2=deletarColuna(MM2,3*(tt+1)/2-2);
          PMM2=deletarColuna(PMM2,tt);
          PMM2=deletarColuna(PMM2,tt);
          pontocurvas=deletarLinha(pontocurvas,(tt+1)/2);
          break
        endif
      endfor
    endfor
  endif
endfor

for tt=1:2:size(PMM2,2)
for i=1:nfaces
##if (pontocurvasw((tt+1)/2,1)>4)
if((PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,3)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,3)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,3)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,3)) )
  MM2(pontocurvas((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvas((tt+1)/2,1)=pontocurvas((tt+1)/2,1)+1;
  break
  elseif((PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,1)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,1)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,1)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,1)) )
  MM2(pontocurvas((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvas((tt+1)/2,1)=pontocurvas((tt+1)/2,1)+1;
  break
  elseif((PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,2)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,2)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,2)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,2)) )
  MM2(pontocurvas((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvas((tt+1)/2,1)=pontocurvas((tt+1)/2,1)+1;
  break
  elseif((PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,3)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,3)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,3)) || (PMM2(1,tt)==F(i,2) && PMM2(1,tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,3)) )
  MM2(pontocurvas((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvas((tt+1)/2,1)=pontocurvas((tt+1)/2,1)+1;
  break
  elseif((PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,1)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,1)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,1)) || (PMM2(1,tt)==F(i,3) && PMM2(1,tt+1)==F(i,2) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,1)) )
  MM2(pontocurvas((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvas((tt+1)/2,1)=pontocurvas((tt+1)/2,1)+1;
  break
  elseif((PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,2)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,2)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,2)) || (PMM2(1,tt)==F(i,1) && PMM2(1,tt+1)==F(i,3) && PMM2(pontocurvas((tt+1)/2,1),tt+1)==F(i,1) && PMM2(pontocurvas((tt+1)/2,1),tt)==F(i,2)) )
  MM2(pontocurvas((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvas((tt+1)/2,1)=pontocurvas((tt+1)/2,1)+1;
  break
endif
##endif
endfor
endfor

printf('Vermelha fecho');

####### AZUL ######

idx2=0;
while(idx2!=1)
pontocurvasaxw=pontocurvasw;
#1
idx=0;
while(idx!=1)
pontocurvasaw=pontocurvasw;
for tt=1:2:size(PMM2w,2)
  ppt=size(PMM2w,2);
  if (tt>ppt)
    break
  endif
for pp=tt+2:2:size(PMM2w,2)
  pptz=size(PMM2w,2);
  if (pp>pptz)
    break
  endif
for i=1:nfaces
  if((PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,pp)==F(i,2) && PMM2w(1,pp+1)==F(i,3)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,pp+1)==F(i,2) && PMM2w(1,pp)==F(i,3)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,pp)==F(i,1) && PMM2w(1,pp+1)==F(i,3)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,pp+1)==F(i,1) && PMM2w(1,pp)==F(i,3)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(pontocurvasw((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(pontocurvasw((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
  elseif((PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,pp)==F(i,3) && PMM2w(1,pp+1)==F(i,1)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,pp+1)==F(i,3) && PMM2w(1,pp)==F(i,1)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,pp)==F(i,2) && PMM2w(1,pp+1)==F(i,1)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,pp+1)==F(i,2) && PMM2w(1,pp)==F(i,1)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(pontocurvasw((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(pontocurvasw((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
  elseif((PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,pp)==F(i,1) && PMM2w(1,pp)==F(i,2)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,pp+1)==F(i,1) && PMM2w(1,pp)==F(i,2)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,pp)==F(i,3) && PMM2w(1,pp)==F(i,2)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,pp+1)==F(i,3) && PMM2w(1,pp)==F(i,2)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(pontocurvasw((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(pontocurvasw((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
  elseif((PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,pp)==F(i,1) && PMM2w(1,pp+1)==F(i,3)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,pp+1)==F(i,1) && PMM2w(1,pp)==F(i,3)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,pp)==F(i,2) && PMM2w(1,pp+1)==F(i,3)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,pp+1)==F(i,2) && PMM2w(1,pp)==F(i,3)))
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(pontocurvasw((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(pontocurvasw((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
  elseif((PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,pp)==F(i,2) && PMM2w(1,pp+1)==F(i,1)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,pp+1)==F(i,2) && PMM2w(1,pp)==F(i,1)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,pp)==F(i,3) && PMM2w(1,pp+1)==F(i,1)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,pp+1)==F(i,3) && PMM2w(1,pp)==F(i,1)))
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(pontocurvasw((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(pontocurvasw((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
  elseif((PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,pp)==F(i,3) && PMM2w(1,pp+1)==F(i,2)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,pp+1)==F(i,3) && PMM2w(1,pp)==F(i,2)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,pp)==F(i,1) && PMM2w(1,pp+1)==F(i,2)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,pp+1)==F(i,1) && PMM2w(1,pp)==F(i,2)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(pontocurvasw((pp+1)/2,1):-1:1,3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(pontocurvasw((pp+1)/2,1):-1:1,pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
endif
endfor
endfor
endfor
if(size(pontocurvasaw,1)==size(pontocurvasw,1))
idx=1;
end
endwhile

printf('Azul 1');

#2
idx=0;
while(idx!=1)
pontocurvasaw=pontocurvasw;
for tt=1:2:size(PMM2w,2)
  ppt=size(PMM2w,2);
  if (tt>ppt)
    break
  endif
for pp=tt+2:2:size(PMM2w,2)
  pptz=size(PMM2w,2);
  if (pp>pptz)
    break
  endif
for i=1:nfaces
  if((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,tt)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,tt)==F(i,3)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
  elseif((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,tt)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,tt)==F(i,1)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
  elseif((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,tt)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,tt)==F(i,2)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
  elseif((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,tt)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,tt)==F(i,3)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
  elseif((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,tt)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,tt)==F(i,1)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
  elseif((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,tt)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,tt)==F(i,2)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-5);
  PMM2w=deletarColuna(PMM2w,pp-2);
  PMM2w=deletarColuna(PMM2w,pp-2);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2-1);
  break
endif
endfor
endfor
endfor
if(size(pontocurvasaw,1)==size(pontocurvasw,1))
idx=1;
end
endwhile

printf('Azul 2');


idx=0;
while(idx!=1)
pontocurvasaw=pontocurvasw;
for pp=1:2:size(PMM2w,2)
  ppt=size(PMM2w,2);
  if (pp>ppt)
    break
  endif
for tt=pp+2:2:size(PMM2w,2)
  pptz=size(PMM2w,2);
  if (tt>pptz)
    break
  endif
for i=1:nfaces
  if((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,tt)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,tt)==F(i,3)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  PMM2w=deletarColuna(PMM2w,pp);
  PMM2w=deletarColuna(PMM2w,pp);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  PMM2w=deletarColuna(PMM2w,tt-2);
  PMM2w=deletarColuna(PMM2w,tt-2);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2-1);
  break
  elseif((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,tt)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,tt)==F(i,1)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  PMM2w=deletarColuna(PMM2w,pp);
  PMM2w=deletarColuna(PMM2w,pp);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  PMM2w=deletarColuna(PMM2w,tt-2);
  PMM2w=deletarColuna(PMM2w,tt-2);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2-1);
  break
  elseif((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,tt)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,tt)==F(i,2)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  PMM2w=deletarColuna(PMM2w,pp);
  PMM2w=deletarColuna(PMM2w,pp);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  PMM2w=deletarColuna(PMM2w,tt-2);
  PMM2w=deletarColuna(PMM2w,tt-2);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2-1);
  break
  elseif((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,tt)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,2) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,1) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,tt)==F(i,3)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  PMM2w=deletarColuna(PMM2w,pp);
  PMM2w=deletarColuna(PMM2w,pp);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  PMM2w=deletarColuna(PMM2w,tt-2);
  PMM2w=deletarColuna(PMM2w,tt-2);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2-1);
  break
  elseif((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt+1)==F(i,2) && PMM2w(1,tt)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,3) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,2) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,tt)==F(i,1)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  PMM2w=deletarColuna(PMM2w,pp);
  PMM2w=deletarColuna(PMM2w,pp);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  PMM2w=deletarColuna(PMM2w,tt-2);
  PMM2w=deletarColuna(PMM2w,tt-2);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2-1);
  break
  elseif((PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt+1)==F(i,3) && PMM2w(1,tt)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2)) || (PMM2w(pontocurvasw((pp+1)/2,1),pp)==F(i,1) && PMM2w(pontocurvasw((pp+1)/2,1),pp+1)==F(i,3) && PMM2w(1,tt+1)==F(i,1) && PMM2w(1,tt)==F(i,2)) )
  MM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((pp+1)/2,1),3*(pp+1)/2-2:3*(pp+1)/2);
  MM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:pontocurvasw((pp+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((pp+1)/2,1),pp:(pp+1));
  PMM2w(pontocurvasw((pp+1)/2,1)+1:pontocurvasw((pp+1)/2,1)+pontocurvasw((tt+1)/2,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+pontocurvasw((pp+1)/2,1);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(pp+1)/2-2);
  PMM2w=deletarColuna(PMM2w,pp);
  PMM2w=deletarColuna(PMM2w,pp);
  pontocurvasw=deletarLinha(pontocurvasw,(pp+1)/2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-5);
  PMM2w=deletarColuna(PMM2w,tt-2);
  PMM2w=deletarColuna(PMM2w,tt-2);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2-1);
  break
endif
endfor
endfor
endfor
if(size(pontocurvasaw,1)==size(pontocurvasw,1))
idx=1;
end
endwhile

printf('Azul 3');

#4
idx=0;
while(idx!=1)
mx=0;
for tt=1:2:size(PMM2w,2)
for pp=1:size(zero2,1)
for i=1:nfaces
  if((PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,2) && zero2(pp,1)==F(i,3)))
  MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-2:end)=V(guardar2(pp,1:nguardar2(pp,1)),:);
  PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1)=guardar2(pp,1:nguardar2(pp,1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end)=guardar2(pp,1:nguardar2(pp,1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,3) && zero2(pp,1)==F(i,1)))
  MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-2:end)=V(guardar2(pp,1:nguardar2(pp,1)),:);
  PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1)=guardar2(pp,1:nguardar2(pp,1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end)=guardar2(pp,1:nguardar2(pp,1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,1) && zero2(pp,1)==F(i,2)))
  MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-2:end)=V(guardar2(pp,1:nguardar2(pp,1)),:);
  PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1)=guardar2(pp,1:nguardar2(pp,1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end)=guardar2(pp,1:nguardar2(pp,1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,1) && zero2(pp,1)==F(i,3)))
  MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-2:end)=V(guardar2(pp,1:nguardar2(pp,1)),:);
  PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1)=guardar2(pp,1:nguardar2(pp,1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end)=guardar2(pp,1:nguardar2(pp,1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,2) && zero2(pp,1)==F(i,1)))
  MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-2:end)=V(guardar2(pp,1:nguardar2(pp,1)),:);
  PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1)=guardar2(pp,1:nguardar2(pp,1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end)=guardar2(pp,1:nguardar2(pp,1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,3) && zero2(pp,1)==F(i,2)))
  MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  MM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-2:end)=V(guardar2(pp,1:nguardar2(pp,1)),:);
  PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1)=guardar2(pp,1:nguardar2(pp,1));
  PMM2w(pontocurvasw((tt+1)/2,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end)=guardar2(pp,1:nguardar2(pp,1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
endif
endfor
endfor
endfor
if(mx==0)
idx=1;
end
endwhile

printf('Azul 4');

##5
idx=0;
while(idx!=1)
mx=0;
for tt=1:2:size(PMM2w,2)
for pp=1:size(zero2,1)
for i=1:nfaces
  if((PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2) && zero2(pp,1)==F(i,3)) )
  MM2w(1:nguardar2(pp,1),end+1:end+3)=V(guardar2(pp,nguardar2(pp,1):-1:1),:);
  MM2w(nguardar2(pp,1)+1:nguardar2(pp,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(nguardar2(pp,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3) && zero2(pp,1)==F(i,1)))
  MM2w(1:nguardar2(pp,1),end+1:end+3)=V(guardar2(pp,nguardar2(pp,1):-1:1),:);
  MM2w(nguardar2(pp,1)+1:nguardar2(pp,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(nguardar2(pp,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1) && zero2(pp,1)==F(i,2)))
  MM2w(1:nguardar2(pp,1),end+1:end+3)=V(guardar2(pp,nguardar2(pp,1):-1:1),:);
  MM2w(nguardar2(pp,1)+1:nguardar2(pp,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(nguardar2(pp,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1) && zero2(pp,1)==F(i,3)))
  MM2w(1:nguardar2(pp,1),end+1:end+3)=V(guardar2(pp,nguardar2(pp,1):-1:1),:);
  MM2w(nguardar2(pp,1)+1:nguardar2(pp,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(nguardar2(pp,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2) && zero2(pp,1)==F(i,1)))
  MM2w(1:nguardar2(pp,1),end+1:end+3)=V(guardar2(pp,nguardar2(pp,1):-1:1),:);
  MM2w(nguardar2(pp,1)+1:nguardar2(pp,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(nguardar2(pp,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
  elseif((PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3) && zero2(pp,1)==F(i,2)) )
  MM2w(1:nguardar2(pp,1),end+1:end+3)=V(guardar2(pp,nguardar2(pp,1):-1:1),:);
  MM2w(nguardar2(pp,1)+1:nguardar2(pp,1)+pontocurvasw((tt+1)/2,1),end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(1:nguardar2(pp,1),end+1)=guardar2(pp,nguardar2(pp,1):-1:1);
  PMM2w(nguardar2(pp,1)+1:pontocurvasw((tt+1)/2,1)+nguardar2(pp,1),end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
  pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+nguardar2(pp,1);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
  PMM2w=deletarColuna(PMM2w,tt);
  PMM2w=deletarColuna(PMM2w,tt);
  pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
  mx=1;
  break
endif
endfor
endfor
endfor
if(mx==0)
idx=1;
end
endwhile

printf('Azul 5');

if(size(pontocurvasaxw,1)==size(pontocurvasw,1))
idx2=1;
end
endwhile

for tt=1:2:size(PMM2w,2)
  if (PMM2w(1,tt)==PMM2w(1,tt+1))
    for pp=1:size(Vumb,1)
      for i=1:nfaces
        if((PMM2w(1,tt)==F(i,1) && Vumb(pp,1)==F(i,2))||(PMM2w(1,tt)==F(i,2) && Vumb(pp,1)==F(i,3))||(PMM2w(1,tt)==F(i,3) && Vumb(pp,1)==F(i,1))||(PMM2w(1,tt)==F(i,1) && Vumb(pp,1)==F(i,3))||(PMM2w(1,tt)==F(i,3) && Vumb(pp,1)==F(i,2))||(PMM2w(1,tt)==F(i,2) && Vumb(pp,1)==F(i,1)))
          MM2w(1,end+1:end+3)=V(Vumb(pp,1),:);
          MM2w(2:pontocurvasw((tt+1)/2,1)+1,end-2:end)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
          PMM2w(1,end+1)=Vumb(pp,1);
          PMM2w(1,end+1)=Vumb(pp,1);
          PMM2w(2:pontocurvasw((tt+1)/2,1)+1,end-1:end)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
          pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
          MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
          MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
          MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
          PMM2w=deletarColuna(PMM2w,tt);
          PMM2w=deletarColuna(PMM2w,tt);
          pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
          break
        endif
        if((PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,1) && Vumb(pp,1)==F(i,2))||(PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,2) && Vumb(pp,1)==F(i,3))||(PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,3) && Vumb(pp,1)==F(i,1))||(PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,1) && Vumb(pp,1)==F(i,3))||(PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,3) && Vumb(pp,1)==F(i,2))||(PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,2) && Vumb(pp,1)==F(i,1)))
          MM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+3)=MM2w(1:pontocurvasw((tt+1)/2,1),3*(tt+1)/2-2:3*(tt+1)/2);
          MM2w(1+pontocurvasw((tt+1)/2,1),end-2:end)=V(Vumb(pp,1),:);
          PMM2w(1:pontocurvasw((tt+1)/2,1),end+1:end+2)=PMM2w(1:pontocurvasw((tt+1)/2,1),tt:(tt+1));
          PMM2w(pontocurvasw((tt+1)/2,1)+1,end-1)=Vumb(pp,1);
          PMM2w(pontocurvasw((tt+1)/2,1)+1,end)=Vumb(pp,1);
          pontocurvasw(end+1,1)=pontocurvasw((tt+1)/2,1)+1;
          MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
          MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
          MM2w=deletarColuna(MM2w,3*(tt+1)/2-2);
          PMM2w=deletarColuna(PMM2w,tt);
          PMM2w=deletarColuna(PMM2w,tt);
          pontocurvasw=deletarLinha(pontocurvasw,(tt+1)/2);
          break
        endif
      endfor
    endfor
  endif
endfor

for tt=1:2:size(PMM2w,2)
for i=1:nfaces
if((PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,3)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,3)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,3)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,3)) )
  MM2w(pontocurvasw((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2w(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvasw((tt+1)/2,1)=pontocurvasw((tt+1)/2,1)+1;
  break
  elseif((PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,1)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,1)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,1)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,1)) )
  MM2w(pontocurvasw((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2w(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvasw((tt+1)/2,1)=pontocurvasw((tt+1)/2,1)+1;
  break
  elseif((PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,2)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,2)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,2)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,2)) )
  MM2w(pontocurvasw((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2w(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvasw((tt+1)/2,1)=pontocurvasw((tt+1)/2,1)+1;
  break
  elseif((PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,3)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,3)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,3)) || (PMM2w(1,tt)==F(i,2) && PMM2w(1,tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,3)) )
  MM2w(pontocurvasw((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2w(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvasw((tt+1)/2,1)=pontocurvasw((tt+1)/2,1)+1;
  break
  elseif((PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,1)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,1)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,1)) || (PMM2w(1,tt)==F(i,3) && PMM2w(1,tt+1)==F(i,2) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,1)) )
  MM2w(pontocurvasw((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2w(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvasw((tt+1)/2,1)=pontocurvasw((tt+1)/2,1)+1;
  break
  elseif((PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,2)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,2)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,2)) || (PMM2w(1,tt)==F(i,1) && PMM2w(1,tt+1)==F(i,3) && PMM2w(pontocurvasw((tt+1)/2,1),tt+1)==F(i,1) && PMM2w(pontocurvasw((tt+1)/2,1),tt)==F(i,2)) )
  MM2w(pontocurvasw((tt+1)/2,1)+1,3*(tt+1)/2-2:3*(tt+1)/2)=MM2w(1,3*(tt+1)/2-2:3*(tt+1)/2);
  pontocurvasw((tt+1)/2,1)=pontocurvasw((tt+1)/2,1)+1;
  break
endif
##endif
endfor
endfor

printf('Azul fecho');

printf('Parte 8');

con6=0;
ncon6=0;
for tt=1:3:size(MM2,2)
  if (pontocurvas((tt+2)/3,1)>1 && MM2(1,tt)==MM2(pontocurvas((tt+2)/3,1),tt) && MM2(1,tt+1)==MM2(pontocurvas((tt+2)/3,1),tt+1))
    circulo(con6+1,1)=tt;
    con6=con6+1;
  else
    ncirculo(ncon6+1,1)=tt;
    ncon6=ncon6+1;
  endif
endfor
con6w=0;
ncon6w=0;
for tt=1:3:size(MM2w,2)
  if (pontocurvasw((tt+2)/3,1)>1 && MM2w(1,tt)==MM2w(pontocurvasw((tt+2)/3,1),tt) && MM2w(1,tt+1)==MM2w(pontocurvasw((tt+2)/3,1),tt+1))
    circulow(con6w+1,1)=tt;
    con6w=con6w+1;
  else
    ncirculow(ncon6w+1,1)=tt;
    ncon6w=ncon6w+1;
  endif
endfor


printf('Parte 9');


   %Plotagem dos pontos para triangulacao
  x = V(:,1);
  y = V(:,2);
  z = V(:,3);

##  Resultado da implementacao octave
  hold on;
##  tetramesh (Tdo,V);
##  trimesh (F,x,y,z);
##  plot3(Positivo(:,1), Positivo(:,2), Positivo(:,3),'gx','LineWidth',[2]);
##  plot3(Negativo(:,1), Negativo(:,2), Negativo(:,3),'mx','LineWidth',[2]);
##  plot3(posumbw(:,1), posumbw(:,2), posumbw(:,3), 'cx','LineWidth',[4]);
##  plot3(umb(:,1), umb(:,2), umb(:,3),'kx','LineWidth',[5]);
  for j=1:size(pontocurvas)(1)
    plot3(MM2(1:pontocurvas(j,1),3*j-2), MM2(1:pontocurvas(j,1),3*j-1), MM2(1:pontocurvas(j,1),3*j),'b-','LineWidth',[2]);
  endfor
  for j=1:size(pontocurvasw)(1)
    plot3(MM2w(1:pontocurvasw(j,1),3*j-2), MM2w(1:pontocurvasw(j,1),3*j-1), MM2w(1:pontocurvasw(j,1),3*j),'r-','LineWidth',[2]);
endfor

##  plot3(V(zero,1), V(zero,2), V(zero,3), 'k*','LineWidth',[5]);
##    plot3(V(zero2,1),V(zero2,2),V(zero2,3), 'y*','LineWidth',[5]);

##    plot(umb3(:,1), umb3(:,2), 'b-*','LineWidth',[5]);
##    plot(curvpar(:,1), curvpar(:,2),'g*','LineWidth',[5]);

##  plot3(MMMMM(:,1), MMMMM(:,2), MMMMM(:,3),'g*','LineWidth',[5]);
##  plot3(V(385,1), V(385,2), V(385,3),'k*','LineWidth',[5]);
##  plot3(V(383,1), V(383,2), V(383,3),'k*','LineWidth',[5]);
##  plot3(V(259,1), V(259,2), V(259,3),'m*','LineWidth',[5]);
##  plot3(V(443,1), V(443,2), V(443,3),'m*','LineWidth',[5]);
##  plot3(V(287,1), V(287,2), V(287,3),'m*','LineWidth',[5]);
##  plot3(MM5(:,1), MM5(:,2), MM5(:,3),'g-*','LineWidth',[5]);
##  plot3(MM6(:,1), MM6(:,2), MM6(:,3),'g-*','LineWidth',[5]);
##  plot3(MM7(:,1), MM7(:,2), MM7(:,3),'g-*','LineWidth',[5]);
##  plot3(MM8(:,1), MM8(:,2), MM8(:,3),'g-*','LineWidth',[5]);
##  plot3(MM9(:,1), MM9(:,2), MM9(:,3),'g-*','LineWidth',[5]);
##  plot3(MM10(:,1), MM10(:,2), MM10(:,3),'g-*','LineWidth',[5]);
##  plot3(MM11(:,1), MM11(:,2), MM11(:,3),'g-*','LineWidth',[5]);
##  plot3(MM12(:,1), MM12(:,2), MM12(:,3),'g-*','LineWidth',[5]);
##  plot3(MM13(:,1), MM13(:,2), MM13(:,3),'g-*','LineWidth',[5]);
##  plot3(MM14(:,1), MM14(:,2), MM14(:,3),'g-*','LineWidth',[5]);
##  plot3(MM15(:,1), MM15(:,2), MM15(:,3),'g-*','LineWidth',[5]);
##  plot3(V(98,1), V(98,2), V(98,3),'k*','LineWidth',[5]);
##  plot3(V(109,1), V(109,2), V(109,3),'k*','LineWidth',[5]);
## end
##  end
##
##  triplot (Tdo,x,y);
##  title ("Delaunay Octave ");
axis square
  hold off;
 toc
 
 
