close all;
clear;
clc;
tic
 %Leitura do arquivo
  [V,F] = leitura_obj('pasta_destino/face_rosto2.obj'); #faceGerada   #pasta_destino/face_rosto.obj
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
        L(i,1)=count3-1;
        LL(i,1)=count3a-1;
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
  
 Positivo=[];
 Negativo=[];
 count55=1;
 count56=1;
  
  for i=1:nverts
    if(Kd2(i,1)>1e-6)
      Positivo(count55,:)=V(i,:);
      count55=count55+1;
    end
    if (Kd2(i,1)<-1e-6)
      Negativo(count56,:)=V(i,:);
      count56=count56+1;
    end
  endfor
  
    count57=1;
    
    for i=1:nfaces
    if ((Kd2(F(i,1),1)>1e-6 && Kd2(F(i,2),1)<-1e-6))
      t=Kd2(F(i,1),1)/(Kd2(F(i,1),1)-Kd2(F(i,2),1));
      MM(count57,:)=V(F(i,1),:)+t*(V(F(i,2),:)-V(F(i,1),:));
      PMM(count57,1)=F(i,1);
      PMM(count57,2)=F(i,2);
      count57=count57+1;
    end
    if ((Kd2(F(i,2),1)>1e-6 && Kd2(F(i,3),1)<-1e-6))
      t=Kd2(F(i,2),1)/(Kd2(F(i,2),1)-Kd2(F(i,3),1));
      MM(count57,:)=V(F(i,2),:)+t*(V(F(i,3),:)-V(F(i,2),:));
      PMM(count57,1)=F(i,2);
      PMM(count57,2)=F(i,3);
      count57=count57+1;
    end
    if ((Kd2(F(i,3),1)>1e-6 && Kd2(F(i,1),1)<-1e-6))
      t=Kd2(F(i,3),1)/(Kd2(F(i,3),1)-Kd2(F(i,1),1));
      MM(count57,:)=V(F(i,3),:)+t*(V(F(i,1),:)-V(F(i,3),:));
      PMM(count57,1)=F(i,3);
      PMM(count57,2)=F(i,1);
      count57=count57+1;
    endif
  end
  
  
  printf('Parte 3');
  
  MM3=[MM PMM];
##  MM3 = unique(MM3,'rows')
  
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
##  elseif ((MM3(i,4)==MM3(i+1,5)) && (MM3(i,5)==MM3(i+1,4)))
##      MM3=deletarLinha(MM3,i+1);
##    endif
##end


##for i=1:size(MM3)(1)
##  nn=size(MM3)(1);
##  for j=1:nn
##    nnn=size(MM3)(1);
##    if (j>nnn)
##      break
##    elseif ((L(MM3(j,4))<6) || (L(MM3(j,5))<6))
##      MM3=deletarLinha(MM3,j);
##    endif
##  endfor
##  
##end

    
##  for i=1:nfaces
##    if ((Kd2(F(i,1),1)>0.01 & Kd2(F(i,2),1)<-0.01) || (Kd2(F(i,1),1)<-0.01 & Kd2(F(i,2),1)>0.01))
##      t=Kd2(F(i,1),1)/(Kd2(F(i,1),1)-Kd2(F(i,2),1));
##      MM(count57,:)=V(F(i,1),:)+t*(V(F(i,2),:)-V(F(i,1),:));
##      PMM(count57,1)=F(i,1);
##      PMM(count57,2)=F(i,2);
##      count57=count57+1;
##    end
##    if ((Kd2(F(i,2),1)>0.01 & Kd2(F(i,3),1)<-0.01) || (Kd2(F(i,2),1)<-0.01 & Kd2(F(i,3),1)>0.01))
##      t=Kd2(F(i,2),1)/(Kd2(F(i,2),1)-Kd2(F(i,3),1));
##      MM(count57,:)=V(F(i,2),:)+t*(V(F(i,3),:)-V(F(i,2),:));
##      PMM(count57,1)=F(i,2);
##      PMM(count57,2)=F(i,3);
##      count57=count57+1;
##    end
##    if ((Kd2(F(i,3),1)>0.01 & Kd2(F(i,1),1)<-0.01) || (Kd2(F(i,3),1)<-0.01 & Kd2(F(i,1),1)>0.01))
##      t=Kd2(F(i,3),1)/(Kd2(F(i,3),1)-Kd2(F(i,1),1));
##      MM(count57,:)=V(F(i,3),:)+t*(V(F(i,1),:)-V(F(i,3),:));
##      PMM(count57,1)=F(i,3);
##      PMM(count57,2)=F(i,1);
##      count57=count57+1;
##    endif
##  end
##  
##
##  
##  printf('Parte 3');
##  
##  MM3=[MM PMM];
##  MM3 = unique(MM3,'rows');
##
##function [A] = deletarColuna(Matrix,coluna)
##  
##  if (coluna == 1)
##    A = Matrix(:,2:size(Matrix)(2));
##    return;
##  end
##  
##  if (coluna == size(Matrix)(2))
##    A = Matrix(:,1:size(Matrix)(2)-1);
##    return;
##  end
##  
##  A1 = Matrix(:,1:coluna-1);
##  A2 = Matrix(:,coluna+1:size(Matrix)(2));
##
##  A = [A1  A2];
##    
##end
##  
##function [A] = deletarLinha(Matrix,linha)
##  
##  if (linha == 1)
##    A = Matrix(2:size(Matrix)(1),:);
##    return;
##  end
##  
##  if (linha == size(Matrix)(1))
##    A = Matrix(1:size(Matrix)(1)-1,:);
##    return;
##  end
##  
##  A1 = Matrix(1:linha-1,:);
##  A2 = Matrix(linha+1:size(Matrix)(1),:);
##
##  A = [A1 ; A2];
##    
##end
##
##for i=1:size(MM3)(1)
##  nn=size(MM3)(1);
##  if (i>=nn)
##    break
##  elseif ((abs(MM3(i,1)-MM3(i+1,1))<10e-6) && (abs(MM3(i,2)-MM3(i+1,2))<10e-6))
##      MM3=deletarLinha(MM3,i+1);
##    endif
##end
##
##
##for i=1:size(MM3)(1)
##  nn=size(MM3)(1);
##  for j=1:nn
##    nnn=size(MM3)(1);
##    if (j>nnn)
##      break
##    elseif ((L(MM3(j,4))<4) || (L(MM3(j,5))<4))
##      MM3=deletarLinha(MM3,j);
##    endif
##  endfor
##  
##end

printf('Parte 4');

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
for i=1:nfaces
if((PMM2(1,2*kk-1)==F(i,1) && PMM2(1,2*kk)==F(i,2) && PMM2(cc,2*kk-1)==F(i,3)) || (PMM2(1,2*kk-1)==F(i,1) && PMM2(1,2*kk)==F(i,2) && PMM2(cc,2*kk)==F(i,3)) )
  MM2(cc+1,3*kk-2:3*kk)=MM2(1,3*kk-2:3*kk);
  cc=cc+1;
  break
  elseif((PMM2(1,2*kk-1)==F(i,2) && PMM2(1,2*kk)==F(i,3) && PMM2(cc,2*kk-1)==F(i,1)) || (PMM2(1,2*kk-1)==F(i,2) && PMM2(1,2*kk)==F(i,3) && PMM2(cc,2*kk)==F(i,1)) )
  MM2(cc+1,3*kk-2:3*kk)=MM2(1,3*kk-2:3*kk);
  cc=cc+1;
  break
  elseif((PMM2(1,2*kk-1)==F(i,3) && PMM2(1,2*kk)==F(i,1) && PMM2(cc,2*kk-1)==F(i,2)) || (PMM2(1,2*kk-1)==F(i,3) && PMM2(1,2*kk)==F(i,1) && PMM2(cc,2*kk)==F(i,2)) )
  MM2(cc+1,3*kk-2:3*kk)=MM2(1,3*kk-2:3*kk);
  cc=cc+1;
  break
  elseif((PMM2(1,2*kk-1)==F(i,2) && PMM2(1,2*kk)==F(i,1) && PMM2(cc,2*kk-1)==F(i,3)) || (PMM2(1,2*kk-1)==F(i,2) && PMM2(1,2*kk)==F(i,1) && PMM2(cc,2*kk)==F(i,3)) )
  MM2(cc+1,3*kk-2:3*kk)=MM2(1,3*kk-2:3*kk);
  cc=cc+1;
  break
  elseif((PMM2(1,2*kk-1)==F(i,3) && PMM2(1,2*kk)==F(i,2) && PMM2(cc,2*kk-1)==F(i,1)) || (PMM2(1,2*kk-1)==F(i,3) && PMM2(1,2*kk)==F(i,2) && PMM2(cc,2*kk)==F(i,1)) )
  MM2(cc+1,3*kk-2:3*kk)=MM2(1,3*kk-2:3*kk);
  cc=cc+1;
  break
  elseif((PMM2(1,2*kk-1)==F(i,1) && PMM2(1,2*kk)==F(i,3) && PMM2(cc,2*kk-1)==F(i,2)) || (PMM2(1,2*kk-1)==F(i,1) && PMM2(1,2*kk)==F(i,3) && PMM2(cc,2*kk)==F(i,2)) )
  MM2(cc+1,3*kk-2:3*kk)=MM2(1,3*kk-2:3*kk);
  cc=cc+1;
  break
endif
endfor
pontocurvas(kk,1)=cc;
kk=kk+1;
endif
endfor

printf('Parte 5');

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



   %Plotagem dos pontos para triangulacao
  x = V(:,1);
  y = V(:,2);
  z = V(:,3);

##  Resultado da implementacao octave
  hold on;
##  tetramesh (Tdo,V);
##  trimesh (F,x,y,z);
##  plot3(Positivo(:,1), Positivo(:,2), Positivo(:,3),'bp');
##  plot3(Negativo(:,1), Negativo(:,2), Negativo(:,3),'rp');
##  plot3(MM(:,1), MM(:,2), MM(:,3),'kp');
  for j=1:size(pontocurvas)(1)
    plot3(MM2(1:pontocurvas(j,1),3*j-2), MM2(1:pontocurvas(j,1),3*j-1), MM2(1:pontocurvas(j,1),3*j),'b-','LineWidth',[3]);
##    pause()
  endfor
##plot3(V(214,1), V(214,2), V(214,3),'k*','LineWidth',[5]);
##plot3(V(193,1), V(193,2), V(193,3),'k*','LineWidth',[5]);
##plot3(V(148,1), V(148,2), V(148,3),'k*','LineWidth',[5]);
##plot3(V(309,1), V(309,2), V(309,3),'k*','LineWidth',[5]);
##plot3(V(362,1), V(362,2), V(362,3),'k*','LineWidth',[5]);
##plot3(V(378,1), V(378,2), V(378,3),'k*','LineWidth',[5]);
##  for j=1:size(pontocurvasw)(1)
##    plot(MM2w(1:pontocurvasw(j,1),3*j-2), MM2w(1:pontocurvasw(j,1),3*j-1),'r-*');
####    pause()
##endfor

##    plot(umb1(:,1), umb1(:,2), 'k-*','LineWidth',[5]);
##    plot(umb2(:,1), umb2(:,2), 'r-*','LineWidth',[5]);
##    plot(umb3(:,1), umb3(:,2), 'b-*','LineWidth',[5]);
##    plot(curvpar(:,1), curvpar(:,2),'g*','LineWidth',[5]);

##  plot3(MM4(:,1), MM4(:,2), MM4(:,3),'g-*','LineWidth',[5]);
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
  

 tempo=toc
