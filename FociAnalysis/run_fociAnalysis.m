function run_fociAnalysis(dirname,paramFit)

%%%%% Encuentra el frame inicial y final por celula con la
    %%%%% mayor longitud in-interrumpida de frames. Genera la matrix
    %%%%% limites. Cada fila es una celula. Esto podria ser una
    %%%%% funcion?

%param_fit=50; %% Minimo de puntos para considerar el fit exponencial valido. 
%list_1 = dir('*xy*');% todas las xy carpetas

% CC_RESULTS:
% First Row: cell length.
% 2: Fitting: Doubling time and Rsqr. 
% 3: differences frame by frame of pole up
% 4: differences frame by frame of pole down\
% 5: Diego's foci counting 
% 6: Wavelet foci counting 
% 7: B times:1foci,2foci,4foci,8foci,16foci,32foci
% 8: C Times:1foci,2foci,4foci,8foci,16foci,32foci
% 9: Synchonization of foci counitng. Check
%                    first coumn
%                     -1: starts Ctime:4 
%                     -2: starts Btime:2 
%                     -3: starts Ctime:2   
%                     -4: starts Btime:1 

%clearvars 
%dirname=pwd;
dirname=fixDir(dirname); %once i join this code i wont need it
contents = dir([dirname,'xy*']); %List all xy folders
num_dir_tmp = numel(contents);
num_xy = 0;


%% Creates a struct for fluor1 structures and the tif stack

dirnamelist=cell(1,num_dir_tmp);
dirnamelist2=cell(1,num_dir_tmp);

for i = 1:num_dir_tmp
    if (contents(i).isdir) && (numel(contents(i).name) > 2)
    num_xy = num_xy+1;
    dirnamelist{i} = [dirname,contents(i).name,filesep,'morphometrics'];
    dirnamelist2{i} = [dirname,contents(i).name,filesep,'fluor1'];
    end
end

%% Results folder

fociresfolder=[dirname,'foci_analysis',filesep,'cc_results',filesep];
mkdir(fociresfolder);

%% Lopps through all XY folders

for p= 1:num_xy
    
    %load morphometrics file
    morphoFiles = dir([dirnamelist{p},filesep,'*pill_MESH.mat']);
    morphoFile= [morphoFiles.folder,filesep, morphoFiles.name];
    load(morphoFile);
    
    %get stack name for fluorescence data
    stacknames = dir([dirnamelist2{p},filesep,'*.tif']);
    stackname= [stacknames.folder,filesep, stacknames.name];
    
    %set variables
    ind=zeros(Nframe,Ncell);
    cont=zeros(Ncell, Nframe);
    limits=zeros(Ncell,2);

    
    %We maybe can use the table that generates morphometrics already
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
        
         %maximo=Nframe; % what is this for?
         %min=1;
         
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
    CCResults = fociAnalysis(dirname,stackname,Ncell,frame,limits,paramFit);

    %Save results 
    foci_name=[fociresfolder,'cc_res','_xy',num2str(num_xy),'_.mat']; 
    save(foci_name,'CCResults');

 
     
end       
end
