function FociAnalysis(dirname,paramFit)

%%%%% Encuentra el frame inicial y final por celula con la
    %%%%% mayor longitud in-interrumpida de frames. Genera la matrix
    %%%%% limites. Cada fila es una celula. Esto podria ser una
    %%%%% funcion?



%param_fit=50; %% Minimo de puntos para considerar el fit exponencial valido. 
%list_1 = dir('*xy*');% todas las xy carpetas

dirname=pwd;
dirname=fixDir(dirname);
contents = dir([dirname,'xy*']); %List all xy folders
num_dir_tmp = numel(contents);
nxy = [];
num_xy = 0;

%% Creates a struct for fluor1 structures and the tif stack

dirnamelist=cell(1,num_dir_tmp);
for i = 1:num_dir_tmp
    if (contents(i).isdir) && (numel(contents(i).name) > 2)
    num_xy = num_xy+1;
    nxy = [nxy, str2double(contents(i).name(3:end))];
    dirnamelist{i} = [dirname,contents(i).name,filesep,'morphometrics'];
    end
end


%% Lopps through all XY folders

for p= 1:num_xy
    
    %load morphometrics file
    morphoFiles = dir([dirnamelist{p},filesep,'*pill_MESH.mat']);
    morphoFile= [morphoFiles.folder,filesep, morphoFiles.name];
    load(morphoFile);
   
    %set variables
    ind=zeros(Nframe,Ncell);
    cont=zeros(Ncell, Nframe);
    limits=zeros(Ncell,2);

    for N=1:Ncell
        for fr=1:Nframe
            allCN = [frame(fr).object.cellID]; % comma separated list expansion 
            ind = find(allCN == N);
            aux = isempty(ind);
            if aux==1
                cont(N,fr)=-10;    
            end
        end
    end  


    for j=1:Ncell
        
         maximo=Nframe;
         min=1;
         
         temp=find(cont(j,:)==0);
         if isempty(temp)==1
             continue
         else
         loc=find(diff(temp)~=1);
         
         if isempty(loc)==1
            limits(j,1)=temp(1);
            limits(j,2)=temp(end);
            continue
         else
            delta=[temp(1),temp(loc),temp(loc(end)+1),temp(end)];
            dif_d=diff(delta);
            pos=find(dif_d==max(dif_d));
            limits(j,1)=delta(pos(1));
            limits(j,2)=delta(pos(1)+1);
          end
         end
    end

   
    
    %Spot detection and Cell cycle analysis        
    CCResults=fun_anal4N(Ncell,frame,morphoFile,limits,paramFit);

    %Save results 
    dataname=[dirname,filesep,'cell_cycle',filesep,morphoFile(1:end-22),'CC_RESULTS.mat']; 
    save(dataname,'CCResults');

 
     
 end       
end
