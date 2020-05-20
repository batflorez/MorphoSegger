function [foci4,foci2,tb,tc,tau_cluster] = cell_cycle_main2(foci,foci3_N,n_f)

foci_t=foci;
aux=find(foci_t>=0);  
                    
foci_t=foci_t(aux(1):aux(1)+10);
foci_t=foci_t(foci_t>0);
foci_t=mean(foci_t);

if foci_t<=1.5
    [foci4,foci2,foci3_N,tb,tc] = cell_cycle2(foci,foci3_N);  
else
    [foci4,foci2,foci3_N,tb,tc] = cell_cycle(foci,foci3_N);
end

tau_cluster=round(clustering_kmeans((foci+foci3_N)/2));
savefig(n_f); %save clustering figures

close all
                    
      
end                    
