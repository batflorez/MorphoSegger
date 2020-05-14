function [tau_cluster] = clustering_kmeans(foci)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

tau_cluster=zeros(1,length(foci));

foci_n=foci(foci>0);

x=[1:1:length(foci_n)];

aux=foci_n>=2 & foci_n<=10;

foci_n=foci_n(aux);

x=x(aux);

pos_a=diff(x)<=2;

x=x(pos_a);
    
foci_n=foci_n(pos_a);

if length(foci_n)<30
    
    return
    
end    


XXX=[x',foci_n'];

%aux=find(foci_n>0);


% est=round(log(foci_n(end))/log(2));
% 
% if est<3
%     
%     return
% 
% else

%E = evalclusters(XXX,'kmeans','silhouette','klist',[est:est+2]);

%n_cluster=E.OptimalK;

[idx,~] = kmeans(XXX,3,'MaxIter',100,'Replicates',10);%,'Start','sample');%,'Distance','cosine');'Distance','cosine'

%n_cluster+delta

for i=1:3%n_cluster
    
    a_vec=XXX(idx==i,2);
 
    %a_pos=XXX(idx==i,1);
    
    %pos_a=diff(a_pos)<=2;
    
    %a_v=a_vec(pos_a);

    
    %p_aux=a_pos<=2; 
    
    %a_vec(p_aux);
      
    
    tau=length(a_vec);
    
    m_tau=mean(a_vec);
    
    
    if m_tau>=2 && m_tau<=3
        
        tau_cluster(2)=tau;
        
    elseif m_tau>=3 && m_tau<=5.5
        
        tau_cluster(3)=tau;
    
    elseif m_tau>5.5 && m_tau<=10

        
        tau_cluster(4)=tau;
        
    end
    
    
end
        
figure;

plot(XXX(idx==1),XXX(idx==1,2),'r.','MarkerSize',20)

hold on

plot(XXX(idx==2),XXX(idx==2,2),'b.','MarkerSize',20)

hold on

plot(XXX(idx==3),XXX(idx==3,2),'c.','MarkerSize',20)







end

