function [foci4,foci2,tb,tc] = cell_cycle_main(foci,foci3_N)

foci_t=foci;

        cont=1;
        for i=2:1:length(foci)-2
                
                    
                    if isnan(foci(i)==1)
                        
                        cont=cont+1;   
                        
                        continue                 
                        

                    end
                    
                    foci_t=foci_t(cont:cont+10);
                    foci_t=foci_t(foci_t>0);
                    foci_t=mean(foci_t);
                    
                    
                    
                    if foci_t<=1.5

                                [foci4,foci2,tb,tc] = cell_cycle2(foci,foci3_N);   
                                
                                break;
                    
                    else
                        
                               [foci4,foci2,tb,tc] = cell_cycle(foci,foci3_N);
                               
                               break
                            
                           
                    end
                        
                    
                    
                    
        end         
                    


    
end
