function FociAnalysis()

    
param_fit=50; %% Minimo de puntos para considerar el fit exponencial valido. 
list_1 = dir('*xy*');% todas las xy carpetas
 
%parfor?
for p=1:1:length(list_1)
    
    name_root=list_1(p).name; 
    name_folder=strcat(name_root,'/morphometrics');
    cd(name_folder)
    list_2 = dir('*pill_MESH.mat'); %pill mesh file
    name_morph=list_2(1).name;
    load(name_morph);

    

    %%%%% Inicia aca

    ind=zeros(Nframe,Ncell);
    cont=zeros(Ncell, Nframe);
    limites=zeros(Ncell,2);

    for N=1:1:Ncell
        for fr=1:1:Nframe
            allCN = [frame(fr).object.cellID]; % comma separated list expansion 
            ind = find(allCN == N);
            aux = isempty(ind);
            if aux==1
                cont(N,fr)=-10;    
            end
        end
    end  



    for j=1:1:Ncell
        
         maximo=Nframe;
         min=1;
         temp=find(cont(j,:)==0);
         if isempty(temp)==1
             continue
         else
         loc=find(diff(temp)~=1);
         
         if isempty(loc)==1
            limites(j,1)=temp(1);
            limites(j,2)=temp(end);
            continue
         else
            delta=[temp(1),temp(loc),temp(loc(end)+1),temp(end)];
            dif_d=diff(delta);
            pos=find(dif_d==max(dif_d));
            limites(j,1)=delta(pos(1));
            limites(j,2)=delta(pos(1)+1);
          end
         end
    end

            %%%% Termina aca. 
            
    %Spot detection and Cell cycle analysis
    salida=fun_anal4N(Ncell,frame,name_morph,limites,param_fit);

    %Saving results
    name=name_morph(1:end-4); 
    arch_nombre=strcat(name,'_OUTPUT_FOCI.mat');
    arch_nombre=strcat('cell_cycle/',arch_nombre);
    save(arch_nombre,'salida');

 
     
 end       
end

