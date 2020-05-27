function  CCResults =  fociAnalysis(stackname,kymofolder,gc_fitfolder,kmeansfolder,fociresfolder,Ncell,frame,limits,paramFit,timeStep,Dparameter,exp_cut,noiseTh)
   

%loading .tif stack and setting xy number
tsStack = tiffread(stackname);
xy_pos = char(regexp(stackname,'(_xy\w*).','match'));
xy_pos = [erase(xy_pos,"."),'_'];
    
for N=1:Ncell

try    
             
     fit_data=zeros(1,size(tsStack,2));
     d_up=zeros(1,size(tsStack,2));
     d_bt=zeros(1,size(tsStack,2));
     x_pole_up=zeros(1,size(tsStack,2));
     x_pole_bt=zeros(1,size(tsStack,2));
     y_pole_up=zeros(1,size(tsStack,2));
     y_pole_bt=zeros(1,size(tsStack,2));   
     l_cell=zeros(1,size(tsStack,2));
     foci=zeros(1,size(tsStack,2));
     foci2=zeros(1,size(tsStack,2)); 
     foci3=zeros(1,size(tsStack,2)); 
     foci4=zeros(1,size(tsStack,2)); 
     tb=zeros(1,length(foci));
     tc=zeros(1,length(foci));
     a1=[];
     a2=[];
     a3=[];
     
     start=limits(N,1);
     finish=limits(N,2);
     cont=finish-start;

     allCN = [frame(finish).object.cellID]; %finding ids per cell 
     ind = find(allCN == N);

     if cont<paramFit  %consecutive points for a single cell to include the trajectory
          continue
     else  
         %%% Build up the mother trajectory for kymograph. It uses 
         %%% the centerline of each cell for the last frame 

         %%%start

%         try
            x_ztr=round(frame(finish).object(ind).centerline(:,2));
%         catch
%            warning('Problem');
%            continue
%         end
         y_ztr=round(frame(finish).object(ind).centerline(:,1));
         l_k=round(length(x_ztr)/2);
         kymo=zeros(3*cont,length(x_ztr));
         temp=1;
    end
        %%%ends

        %%%Generates and ROI using regionprops and BoundingBox
        %%%in the last frame for each cell. This is
        %%%used for wavelet foci counting. 
        %%%start


%    try
        xs=round(frame(finish).object(ind).Xcont);
        ys=round(frame(finish).object(ind).Ycont);
        label=num2str(N);
%    catch
%        warning(label)
%        continue
%    end


    AB=tsStack(start).data;
    Im_s=zeros(size(AB));

    for p=1:length(xs)
        Im_s(ys(p), xs(p))=1; 
    end

    se = strel('square',3);
    Im_s = imclose(Im_s,se);     
    Im_s=imfill(Im_s,'holes');
    stats = regionprops(Im_s,'BoundingBox');

    xMin = ceil(stats.BoundingBox(1));
    xMax = xMin + stats.BoundingBox(3) - 1;
    yMin = ceil(stats.BoundingBox(2));
    yMax = yMin + stats.BoundingBox(4) - 1;

    % ends for for each frame each cell. 
  
      %Foci counting
        for k=start:finish           
        
        try    
            
            
            
%               try
                    allCN = [frame(k).object.cellID]; % comma separated list expansion 
                    ind = find(allCN == N);
                    label=num2str(N);
%               catch
%                    warning(label)
%                    continue
%               end

%               try
                    xs=round(frame(k).object(ind).Xcont);
                    ys=round(frame(k).object(ind).Ycont);
                    label=num2str(N);
%               catch
%                     warning(label)
%                     continue
%               end
               
              % pole calculations
               pole_2=frame(k).object(ind).pole2; 

               if pole_2==1
                   pole_2=frame(k).object(ind).pole1;   
               end

%               try % added by Andres
                x_pole_up(k)=xs(1);    
                y_pole_up(k)=ys(1);
                x_pole_bt(k)=xs(pole_2);        
                y_pole_bt(k)=ys(pole_2);
                l_cell(k)=frame(k).object(ind).length;      
%               catch
%                   continue
%               end
  
                
                                
                A=tsStack(k).data;                                   
                AA=A;
                AB=A;                                                                                                  
                A = imgaussfilt(A,2); 
                I=A;

                %%Generates binary per cell per frame
                Im1=zeros(size(I));

                for p=1:length(xs)
                    Im1(ys(p), xs(p))=1; 
                end

                se = strel('square',3);
                Im1 = imclose(Im1,se);     
                Im2=imfill(Im1,'holes');
                loc=find(Im2==0);
                loc1=find(Im2==1);

                % Generate the 3 trayectories based on the centerline (and
                % 2 parallel trajectories). The signals is integrated for
                % thress trajectories. It is referenced to the main
                % trajectory (last frame) using centroid.

                x=round(frame(k).object(ind).centerline(:,2));
                y=round(frame(k).object(ind).centerline(:,1));
                x_m=round(frame(k).object(ind).pill_mesh(:,2));
                y_m=round(frame(k).object(ind).pill_mesh(:,1));

                %x_ms2=round(frame(k).object(ind).pill_mesh(:,4));
                %y_ms2=round(frame(k).object(ind).pill_mesh(:,3));

                delta_x= round(abs(x-x_m)/2.5);
                delta_y= round(abs(y-y_m)/2.5);
                %[x,y]=cleanup(x,y);
                %l_cell(k)=length(x);

                lin_ind = sub2ind(size(I),x,y);
                lin_ind2 = sub2ind(size(I),x+delta_x,y+delta_y);
                lin_ind3 = sub2ind(size(I),x-delta_x,y-delta_y);

                a1=double(I(lin_ind));
                a2=double(I(lin_ind2));
                a3=double(I(lin_ind3));
                a=a1+a2+a3;
                aa=A(lin_ind);
                
                %some cells are problematic, use this to prevent stopping
                %the calculation
%                 try
                    b=sgolayfilt(a,4,11); 
%                 catch
%                    continue
%                 end

                [pks1,locs1]=findpeaks(b,'MinPeakDistance',5);
                counter=0;

                % DIEGO'S FOCI COUNTING: Since singal2noise is low, local
                % variability of each peak was assessed by the derivate of the
                % signal over 6 neighboors pixels. Lower than Diego_parameter is
                % rejected.

                for q=1:length(pks1)
                      try
                           bbb=(b(locs1(q)-3:locs1(q)+3));
                      catch
                      continue
                      end
  
                      if sum(abs(diff(bbb)))>Dparameter
                         counter=counter+1;    
                      end
                end

               foci(k)=counter;
               %%referencing to mother trajectory
               halved=round(length(b)/2);
               r=1;

               if l_k<halved
                   continue
               end    

               for q=l_k-halved+1:1:l_k+halved-1
                   
                    kymo(3*temp,q)=aa(r);
                    kymo(3*temp-1,q)=aa(r);
                    kymo(3*temp-2,q)=aa(r);
                    r=r+1;
               end
               %%%% Wavelet counting
               AA(loc)=0;                                 
               AA=AA(yMin:yMax,xMin:xMax);
               AB=AB(yMin:yMax,xMin:xMax);                                
               foci3(k)=wavelet_foci(AA,AB,noiseTh);
               %%%% End wavelet 
               temp=temp+1;
               
        
        catch        
        continue
        end       
               
               
        end


    aux2=find(l_cell==0);
    l_cell(aux2)=NaN;         
    
    time=timeStep*(0:1:length(l_cell)-1); 
    [fit, gof, t_t]=createFit_exp(time, l_cell,paramFit,exp_cut);

    if t_t==1      
        continue                                
    end    
    
    % Saving growth curves and exponential fits
    gc_fitname  = [gc_fitfolder,'gcfit',xy_pos,'cell_',label,'_','.fig'];
    savefig(gc_fitname); %save plot                   
    close all

    d_up=[sqrt(diff(x_pole_up).^2+diff(y_pole_up).^2),NaN];
    d_bt=[sqrt(diff(x_pole_bt).^2+diff(y_pole_bt).^2),NaN];                         
    d_up(aux2)=NaN;  
    d_bt(aux2)=NaN;                            
    foci(aux2)=NaN;
    foci3(aux2)=NaN;
    foci4(aux2)=NaN;

    fit_data(1)=log(2)/fit.b;
    fit_data(2)=gof.rsquare;
    
    %Performing clustering on foci data
    kmeans_name = [kmeansfolder,'kmeansFoci',xy_pos,'cell_',label,'_','.fig'];
    [foci4, ~, tb, tc, tau]=cell_cycle_main2(foci,foci3,kmeans_name);
    output=[l_cell;fit_data;d_up;d_bt;foci3;foci;tb;tc;tau;foci4];
    CCResults(N).out=output;

    % Generating and saving kymographs
    kymo_name = [kymofolder,'kymograph',xy_pos,'cell_',label,'_','.tif'];
    I_a=adapthisteq(mat2gray(kymo));
    imwrite(I_a,kymo_name);

    
catch        
continue
end     
end

%Save results 
foci_name=[fociresfolder,'cc_res',xy_pos,'.mat']; 
save(foci_name,'CCResults');

end

    