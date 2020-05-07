%Tristan Ursell
%cell Geneaology by overlap integral and connected component analysis
%Feb 2012

%Several bugs corrected - Andres Florez 04/20/20

function cells = family_cl(G_parent,pert_same,pert_child,frame_diff)

%calculate total number of regions in the stack
info1=imfinfo(G_parent);
n_frames=length(info1);
sx=info1.Width;
sy=info1.Height;
n_regs=zeros(1,n_frames);
G_im=cell([1,n_frames]);

%pre-load sparse images
for i=1:n_frames
    G_im{i}=sparse(double(imread(G_parent,i)));
    n_regs(i)=max(G_im{i}(:));
end

%index list
ind1=cumsum(n_regs);
ind2=[0,ind1];
Nregs=sum(n_regs);

%calculate the frame and object numbers for each individual element
fr=zeros(1,Nregs);
obj=zeros(1,Nregs);
for i=1:Nregs
    %calculate the frame for this cell
    fr(i)=find(i<=ind1,1,'first');
    %calculate the number of the cell in the frame
    obj(i)=i-ind2(fr(i));
end

%precalculate region areas
areas=zeros(Nregs,1);
Xbox1=zeros(Nregs,1);
Ybox1=zeros(Nregs,1);
Xbox2=zeros(Nregs,1);
Ybox2=zeros(Nregs,1);
maxD=zeros(Nregs,1);
for i=1:n_frames
    if n_regs(i)>0
        prop1=regionprops(uint16(full(G_im{i})),'Area','BoundingBox');
        
        %get indices for storing info
        if i==1
            curr_inds=1:ind1(1);
        else
            curr_inds=(ind1(i-1)+1):ind1(i);
        end
        
        %record areas
        areas(curr_inds)=[prop1.Area];
        
        %get bounding box
        btemp=cat(1,prop1.BoundingBox);
        Xbox1(curr_inds)=btemp(:,1);
        Xbox2(curr_inds)=Xbox1(curr_inds)+btemp(:,3);
        Ybox1(curr_inds)=btemp(:,2);
        Ybox2(curr_inds)=Ybox1(curr_inds)+btemp(:,4);
        
        %maximum dimension of segmented region
        maxD(curr_inds)=sqrt((Xbox1(curr_inds)-Xbox2(curr_inds)).^2+(Ybox1(curr_inds)-Ybox2(curr_inds)).^2);
    end
end
%round off bounding box values
Xbox1=round(Xbox1);
Xbox2=floor(Xbox2);
Ybox1=round(Ybox1);
Ybox2=floor(Ybox2);

%calculate region centroid weighting vectors
%x_weight=1:1:sx; %part of distance calc
%y_weight=1:1:sy; %part of distance calc

%{
%initialize grouping matrix
try
    G_groups=zeros(Nregs,Nregs);
    %D_groups=zeros(Nregs,Nregs); %part of distance calc
    set(handles.text_process,'String','Creating cell lineage...');
catch
    %calculate max number of non-zero sparse elements
    nzmax=Nregs*(2*frame_diff+1);
    
    clear G_groups
    G_groups=sparse([1,Nregs],[1,Nregs],[1,1],Nregs,Nregs,nzmax);
    set(handles.text_process,'String','Creating cell lineage (sparse)...');

    %D_groups=sparse([1,Nregs],[1,Nregs],[1,1],Nregs,Nregs,znmax);
end
%}

%initialize sparse grouping matrix
%calculate max number of non-zero sparse elements
nzmax=Nregs*(2*frame_diff+1);
clear G_groups
G_groups=sparse([1,Nregs],[1,Nregs],[1,1],Nregs,Nregs,nzmax);
disp('Creating cell lineage (sparse)...');
drawnow

for i=1:Nregs
    %get current single cell BW image
    sub1=G_im{fr(i)}==obj(i);
    
    %match frame_diff condition
    cond1_1=abs(fr(i)-fr)<=frame_diff;
    cond1_2=abs(fr(i)-fr)~=0;
    
    %find bounding boxes that overlap
    cond2_1=or(and(Xbox2>Xbox1(i),Xbox2(i)>Xbox1),and(Xbox1<Xbox2(i),Xbox1(i)<Xbox2));
    cond2_2=or(and(Ybox2>Ybox1(i),Ybox2(i)>Ybox1),and(Ybox1<Ybox2(i),Ybox1(i)<Ybox2));
    inds1=find(and(and(cond1_1,cond1_2),and(cond2_1,cond2_2)'));
    
    %get centroid for current contour
    %x1=sum(x_weight.*sum(sub1,1))/sum(sub1(:)); %part of distance calc
    %y1=sum(y_weight.*sum(sub1,2)')/sum(sub1(:)); %part of distance calc
    
    for j=1:length(inds1)
        %get current index
        ind0=inds1(j);
        
        %get current single cell BW image
        sub2=G_im{fr(ind0)}==obj(ind0);
        
        %get centroid for current contour
        %x2=sum(x_weight.*sum(sub2,1))/sum(sub2(:)); %part of distance calc
        %y2=sum(y_weight.*sum(sub2,2)')/sum(sub2(:)); %part of distance calc
        
        %calculate distance matrix
        %D_groups(i,ind0)=sqrt((x1-x2)^2+(y1-y2)^2); %part of distance calc
        
        %generate sub-index list
        xind0=round(min(Xbox1(i),Xbox1(ind0))):floor(max(Xbox2(i),Xbox2(ind0)));
        yind0=round(min(Ybox1(i),Ybox1(ind0))):floor(max(Ybox2(i),Ybox2(ind0)));
        
        %get overlap image
        sub3=sub1(yind0,xind0).*sub2(yind0,xind0);
        
        %calculate overlap, groups are formed only if the mutual overlap is
        %above the cutoff (pert_same)
        overlap_A=sum(sum(sub3))./[areas(i),areas(ind0)];
        G_groups(i,ind0)=and(overlap_A(1)>=pert_same,overlap_A(2)>=pert_same);
        G_groups(ind0,i)=and(overlap_A(1)>=pert_same,overlap_A(2)>=pert_same); %symmetrize
    end
    disp(['Cell lineage: ' num2str(round(100*i/Nregs)) '%'])
end

%{
%test code for graph vizualization
xy_grid=zeros(Nregs,2);
dt=2;
dx=2;
for i=1:Nregs
    xy_grid(i,1)=fr(i)*dt;
    %xy_grid(i,2)=dx*(obj(i)-n_regs(fr(i))/2);
    xy_grid(i,2)=dx*obj(i);
end
%}

%run connected component analysis
try %attempt regular group analysis (i.e. non-sparse)
    disp('Using sparse connected-component analysis.')
    [conn_groups,~]=graph_analysis_sparse_morphometrics_cl(G_groups,'min_conn',0,'max_link',n_frames); %calling right function fixed - Andres Florez
catch er1 %attempt group analysis with sparse representation
    if strcmp(er1.identifier,'MATLAB:nomem')
        disp('Using non-sparse connected-component analysis.')
        [conn_groups,~]=graph_analysis_morphometrics_cl(full(G_groups),'min_conn',0,'max_link',n_frames); %calling right function fixed - Andres Florez
    else
        error(['In graph_analysis*.m, ' er1.message])
    end
end
%create output structure
ngrps=length(conn_groups);

%perform check that no cells have more 'frame' entries than there are
%frames' -- if they do, this indicates two distinct cell lineages have been
%overlapped
fr_app=[conn_groups.num_els]>n_frames;
if sum(fr_app)>0
    warning(['Each of cells ' num2str(find(fr_app)) ' may be a combination of more than one lineage.'])
end

%{
%perform union check -- makes sure that the elements in each group are
%distinct from the elements in every other group (just a check)
%this version actually tells you which groups overlap with which
set_overlap=zeros(ngrps,ngrps);
for i=1:ngrps
    for j=1:ngrps
        if i~=j
            set_overlap(i,j)=~isempty(intersect([conn_groups(i).elements],[conn_groups(j).elements]));
        end
    end
end

if sum(set_overlap(:))>0
    warning('Cell lineage groups are overlapping -- something is up.')
end

%this version simply tells you that groups overlap (so no specificity), but requires less memory
for i=1:ngrps
    for j=1:ngrps
        if i~=j
            if ~isempty(intersect([conn_groups(i).elements],[conn_groups(j).elements]))
                warning('Cell lineage groups are overlapping -- something is up.')
                break
            end
        end
    end
end
%}

for i=1:ngrps
    %get frame and object list
    cells(i).frames=fr([conn_groups(i).elements]);
    cells(i).bw_label=obj([conn_groups(i).elements]);
end

%handle orphans
Fr_orphans=fr(sum(G_groups)==0);
obj_orphans=obj(sum(G_groups)==0);

end %This was missing - Andres Florez












