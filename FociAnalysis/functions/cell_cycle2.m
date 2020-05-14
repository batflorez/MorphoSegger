function [foci4,foci2,tb,tc] = cell_cycle2(foci,foci3_N)


tc1=0;
tb1=0;

tc2=0;
tb2=0;

tc4=0;
tb4=0;

tc8=0;
tb8=0;

tc16=0;
tb16=0;

tc32=0;
tb32=0;


tb=zeros(1,length(foci));
tc=zeros(1,length(foci));

foci4=zeros(1,length(foci));

foci2=foci3_N;

foci3=foci2;


aux_g=0;


cont=1;            

            for i=2:1:length(foci)-2
                
                    if isnan(foci2(i)==1)
                        
                        cont=cont+1;   
                        
                        continue   

                    end
                
                
                
                    if aux_g==0

                            
                        
                        
                            if length(find(foci2(cont:cont+10)==1))<3

                                    aux_g=2;

                                    continue
                            
                           
                            end
                        
                            if foci2(i)==1

                               
                                
                                if foci2(i+1)==2
                                        
                                    
                                        
                                    aux_g=1;
                                        
                                    continue
                                    
                                    
                                end    
                                
                                
                            end    
                        

                            if foci(i)==0 

                                foci2(i)=0;
                                foci3_N(i)=0;
                                
                                
                                if (foci2==2)

                                    aux_g=1;
                                    
                                    continue

                                end    

                            end
                     
                           
                        
                        
                    

                
                    
                
                
                
                    elseif aux_g==1

                        
                        
                        
                        if foci2(i)==2

                               
                                
                                if (foci2(i+1)>3) && (foci2(i+1)<=5)
                                        
                                    
                                        
                                    aux_g=2;
                                        
                                    continue
                                    
                                    
                                end    
                                
                                
                        end    
                        

                        if foci(i)==0 

                                foci2(i)=0;
                                foci3_N(i)=0;
                                
                                
                                if (foci2(i+1)>3) && (foci2(i+1)<=5)

                                    aux_g=2;
                                    
                                    continue

                                end    

                        end
                     
                       if length(find(foci2>=1 & foci2<=2))<4

                            aux_g=2;
                            
                            continue
                            
                           
                       end
                        
                        
                    
%                  

                   




                    elseif aux_g==2

                        
                        if length(find(foci2>3 & foci2<=5))<4

                           aux_g=4;
                           
                           continue
                           
                        end    
                       
                    
                        if (foci2(i)>3) && (foci2(i)<=5)

                               
                                
                                if (foci2(i+1)>=6) && (foci2(i+1)<=10)
                                        
                                    aux_g=4;
                                    
                                    continue
                                        
                                end    
                                
                                
                        end    
                        

                        if foci(i)<1 

                                foci2(i)=0;
                                foci3_N(i)=0;
                                
                                
                                if (foci2(i+1)>=6) && (foci2(i+1)<=10)

                                    aux_g=4;
                                    
                                    continue

                                end    

                        end
                   
                       
                        
                   
                   
                   

                   
%                  
                    elseif aux_g==4

        
                    
                        if (foci2(i)>=6) && (foci2(i)<=10)

                               
                                
                                if (foci2(i+1)>=11) && (foci2(i+1)<=20)
                                        
                                    aux_g=8;
                                    
                                    continue
                                        
                                end    
                                
                                
                        end    
                        

                        if foci(i)<=2

                                foci2(i)=0;
                                foci3_N(i)=0;
                                
                                
                                if (foci2(i+1)>=11) && (foci2(i+1)<=20)

                                    aux_g=8;
                                    
                                    continue

                                end    

                        end
                   
                       if length(find(foci2>=6 & foci2<=10))<4

                                aux_g=8; 
                                
                                continue
                       end     
                        
                   
                   
                   
                   
                   
                    elseif aux_g==8

        
                    
                        if (foci2(i)>=11) && (foci2(i)<=20)

                               
                                
                                if (foci2(i+1)>=21) && (foci2(i+1)<=36)
                                        
                                    aux_g=16;
                                    
                                    continue
                                        
                                end    
                                
                                
                        end    
                        

                        if foci(i)<6

                                foci2(i)=0;
                                foci3_N(i)=0;
                                
                                
                                if (foci2(i+1)>=21) && (foci2(i+1)<=36)

                                    aux_g=16;
                                    
                                    continue

                                end    

                        end
                        
                        if length(find(foci2>=11 & foci2<=20))<4
    
                            aux_g=16;
                        end
                        
                        
                  
                  
                   
                   
                   
                   
                    elseif aux_g==16
                

                        if foci(i)<7

                                foci2(i)=0;
                                foci3_N(i)=0;
                                
                                
                                continue
      

                        end
                         
                   end 
           end
            
            
            aux_z=find(isnan(foci));
            
            foci3_N(aux_z)=NaN;
            
            %foci_nn=(foci+foci3_N)/2;
            
            %foci_nn=(foci+foci3_N)/2;
            
            %foci_nn=medfilt1(foci_nn,3);
            
            [fit_s, ~]=createFit_smoothingspline(foci,0.975);

            x=[1:length(foci)];
            foci_n1=fit_s(x)';
            
            [fit_s2, ~]=createFit_smoothingspline(foci3_N,0.975);

            %x=[1:length(foci3)];
            foci_n2=fit_s2(x)';
            
            %[fit_s3, ~]=createFit_smoothingspline(foci3_N,0.975);

            %x=[1:length(foci3_N)];
            %foci_n3=fit_s3(x)';
            
            foci_n=(foci_n1+foci_n2)/2;

            %foci_n=(foci+foci2)/2;
            
            foci_n(aux_z)=NaN;
            
            aux=find(foci3_N==0);
            foci_n(aux)=0;
            
            aux_g=0;
            
            cont=1;    
            
            for i=1:1:length(foci)-2
                
                
                if isnan(foci_n(i)==1)
                    
                    cont=cont+1;   
                 
                    continue   

                end
                
                
                if aux_g==0

        
                    if length(find(foci_n(cont:10+cont)>=0.5 & foci_n(cont:10+cont)<1.5))<3 && foci_n(i)~=0

                            aux_g=1;
                            
                            continue
                    end  
                    
                    
                    if (foci_n(i)>=0.5) && (foci_n(i)<1.5) 

                            tc1=tc1+1;
                                
                            if (foci_n(i+1)>=1.5) && (foci_n(i+1)<3)
                                
                                if (foci_n(i+2)>=1.5) && (foci_n(i+2)<3)
                                
                                aux_g=1;
                                
                                continue
                                
                                end
                                        
                            end    
                                
                                
                     end    

                     if foci_n(i)==0 

                            
                            tb1=tb1+1;
                            
                            
                                
                            if (foci_n(i+1)>=1.5) && (foci_n(i+1)<3)

                                aux_g=1;
                                
                                continue

                            end
                            
                            if (foci_n(i+1)>=0.5) && (foci_n(i+1)<1.5)
                                     
                                
                                     if (foci_n(i-1)>=0.5) && (foci_n(i-1)<1.5)
                                     
                                         tb1=tb1-1;
                                     
                                         tc1=tc1+1;
                                     
                                     end
                                        
                                     
                            end  
                            
                            if (foci_n(i+1)>=0.5) && (foci_n(i+1)<1.5)

                                if (foci_n(i+2)>=1.5) && (foci_n(i+2)<3)

                                    aux_g=1;
                                    
                                    %tc2=tc2+1;

                                    continue
                                
                                end

                            end
 
                            

                     end
                
                      
                         
                     
                
                
                
                
                elseif aux_g==1

                   
                    
                                if (foci_n(i)>=1.5) && (foci_n(i)<3) 

                                        tc2=tc2+1;

                                        if (foci_n(i+1)>=3) && (foci_n(i+1)<5.5)

                                            if (foci_n(i+2)>=3) && (foci_n(i+2)<5.5)


                                                    aux_g=2;

                                                    continue

                                            end    



                                        end    


                                 end    

                                 if foci_n(i)==0 

                                        tb2=tb2+1;



                                        if (foci_n(i+1)>=3) && (foci_n(i+1)<5.5)

                                            aux_g=2;

                                            continue

                                        end

                                        if (foci_n(i+1)>=1.5) && (foci_n(i+1)<3)


                                                 if (foci_n(i-1)>=1.5) && (foci_n(i-1)<3)

                                                     tb2=tb2-1;
                                                     tc2=tc2+1;

                                                 end


                                        end     





                                         if (foci_n(i+1)>=0.5) && (foci_n(i+1)<3)
            
                                            if (foci_n(i+2)>=3) && (foci_n(i+2)<5.5)
            
                                                aux_g=2;
                                                
                                                %tc4=tc4+1;
            
                                                continue
                                            
                                            end
            
                                        end




                                 end

                                 if length(find(foci_n>=1.5 & foci_n<3))<4 && foci_n(i)~=0

                                        aux_g=2;

                                        continue

                                 end  
    
                     
                
                
                
                
                
                
                elseif aux_g==2
                
                    
               
                    
                    
                    if (foci_n(i)>=3) && (foci_n(i)<5.5) 

                            tc4=tc4+1;
                            
                            if (foci_n(i+1)>=5.5) && (foci_n(i+1)<10)
                                     
                                
                                     if (foci_n(i+2)>=5.5) && (foci_n(i+2)<10)
                                        
                                         aux_g=4;
                                     
                                         continue
                                     
                                     end
                                        
                                     
                            end    
                                     
                    end    

                    if foci_n(i)==0 
             
                            
                            tb4=tb4+1;
                            
                            if (foci_n(i+1)>=5.5) && (foci_n(i+1)<10)
                                        
                                     aux_g=4;
                                     
                                     continue
                                        
                            end  
                            
                            if (foci_n(i+1)>=3) && (foci_n(i+1)<5.5)
                                     
                                
                                     if (foci_n(i-1)>=3) && (foci_n(i-1)<5.5)
                                     tb4=tb4-1;
                                     tc4=tc4+1;
                                     end
                                        
                                     
                            end   
                            
                            
%                             if (foci_n(i+1)>=0.5) && (foci_n(i+1)<5.5)
% 
%                                 if (foci_n(i+2)>=5.5) && (foci_n(i+2)<10)
% 
%                                     aux_g=1;
%                                     
%                                     tc8=tc8+1;
% 
%                                     continue
%                                 
%                                 end
% 
%                             end
                            

                    end
                
                    
                    if length(find(foci_n>=2.5 & foci_n<5.5))<4 && foci_n(i)~=0  

                            aux_g=4;
                            continue
                    end     
                        
                
                
                
                
                

                elseif aux_g==4    
                    
                    
                    
                    if (foci_n(i)>=5.5) && (foci_n(i)<10) 

                        tc8=tc8+1;
                        
                        if (foci_n(i+1)>=10) && (foci_n(i+1)<20)
                                     
                                     if (foci_n(i+2)>=10) && (foci_n(i+2)<20)   
                                        
                                         aux_g=8;
                                        
                                         continue
                                     
                                     end
                                        
                        end  
                       

                    end    

                    if foci_n(i)==0 

                        tb8=tb8+1;  
                        
                         if (foci_n(i+1)>=10) && (foci2(i+1)<20)
                                        
                                     aux_g=8;
                                     continue
                                        
                         end
                         


                    end

                    if length(find(foci_n>=5.5 & foci_n<10))<4 && foci_n(i)~=0

                            aux_g=8;
                    end 

                 
                
                
                
                
                
                elseif aux_g==8
                
                
                
                    if (foci_n(i)>=10) && (foci_n(i)<20) 
                        
                        tc16=tc16+1;  
                        
                        if (foci_n(i+1)>=20) && (foci_n(i+1)<36)
                                     
                                     if (foci_n(i+2)>=20) && (foci_n(i+2)<36)   
                                     
                                         aux_g=16;
                                         
                                         continue
                                     
                                     end
                                        
                        end  
 
                      
                    end    

                    if foci_n(i)==0 

                        tb16=tb16+1;  
                        
                        if (foci_n(i+1)>=20) && (foci_n(i+1)<36)
                                        
                                     aux_g=16;
                                     continue
                                        
                        end  

                    end
                
                 
                    if length(find(foci_n>=10 & foci_n<20))<4 && foci_n(i)~=0

                            aux_g=16;
                            
                            continue
                    end         
                    
                    
                
                
               
                
                elseif aux_g==16

                     if (foci_n(i)>=20) && (foci_n(i)<36) 

                        tc32=tc32+1;  

                    end    

                     if foci_n(i)==0

                        tb32=tb32+1;
                        continue


                     end

                end

            end
            
            
    tb(1)=tb1;
    tb(2)=tb2;
    tb(3)=tb4;
    tb(4)=tb8;
    tb(5)=tb16;
    tb(6)=tb32;
    
    tc(1)=tc1;
    tc(2)=tc2;
    tc(3)=tc4;
    tc(4)=tc8;
    tc(5)=tc16;
    tc(6)=tc32;
   
    foci4=foci_n;
    
    foci4(aux)=0;
    
        aux_fn=find(foci_n>=0, 1);
    
    if isempty(aux_fn)
       
        aux_fn=1;
        
    end    
    
       if tc1>0 && tb1<=1

        
            if tc2>0 && tc4>0 && tc8>8 
        
                   aux_stop=find(foci_n>=11,1);

                   foci4=foci4(aux_fn+tb1+tc1:aux_stop-1);
                   
                   foci4=[-3,foci4,NaN*ones(1,(length(foci)-length(foci4)-1))];
    
            end
   
  
        
       elseif tc1>0 && tb1>=1
    
            if tc2>0 && tc4>0 && tc8>0
                
                   aux_stop=find(foci_n>=11,1);
                   
                   foci4=foci4(aux_fn+tc1:aux_stop-1);
                   
                   foci4=[-2,foci4,NaN*ones(1,(length(foci)-length(foci4)-1))];
               

            end

           
       end

    foci4=round(foci4);
    
end