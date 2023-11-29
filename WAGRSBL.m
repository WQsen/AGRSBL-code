function [source_power,Azimuth,source_id,iter,Noise_variance,Mu_x]=WAGRSBL(Y,Azimuth_grid,grid_interval,Aa,f0,D,c,K)
%% input:
%% Y: Array output data
%% Azimuth_grid: Space discontinue Azimuth angle grid
%% grid_interval: Interval of the adjacent angle grid
%% Aa: Compeleted array steering matrix
%% K: Source number,if K is not a prior,we let K=M-1.
%% output:
%% Source_power_re:signal space power
%% Azimuth:All space angle,include DOA
%% j1:iter number in pre-estimate
%% j2:iter number in secondly-estimate
%% source_id_: DOA index in Azimuth

[M,L]=size(Y(:,:,1));
F=size(Y,3);
b=10^(-7);
ro_=1;
es=10^(-12);
a_=1;
Ry=[];
for ff=1:F
Ry(:,:,ff)=(Y(:,:,ff)*Y(:,:,ff)')/L;     %样本协方差矩阵
end
tol1=0.001;    %迭代停止门限
Azimuth_grid_=Azimuth_grid;
N=length(Azimuth_grid);
noise_variance=10^(-3)*mean(var(mean(Y,3)));
alpha0=1/noise_variance;   %噪声方差倒数
source_power=zeros(N,1);
for ff=1:1:F
source_power =source_power+ sum(abs(Aa(:,:,ff)'*Y(:,:,ff)),2)/(L*M); % 最小二乘初始化
end
source_power=(source_power)/F;
Aa_=Aa;
jmax=500;
e=1;
j1=0;
Length_add=0;
gri_f=1.5;
I_=2*floor( log2(grid_interval/gri_f) ); % SBL过程最大迭代数
while (j1<=I_) %迭代次
    j1=j1+1;                  % number of iter
    source_power_old=source_power;
    gamma=diag(source_power);
    Mu_x_norm=zeros(length(source_power),1);
    YA_trace=0;
    Mu_x=[];sigma_x_ii=[];
    for ff=1:F
    B(:,:,ff)=Aa_(:,:,ff)*gamma*Aa_(:,:,ff)';
    sigma_y(:,:,ff)=1/alpha0*eye(M)+B(:,:,ff);
    sigma_y_inv(:,:,ff)=inv(sigma_y(:,:,ff));
    Aa_sigma_y_inv=Aa_(:,:,ff)'* sigma_y_inv(:,:,ff);
    Mu_x(:,:,ff)=gamma*Aa_sigma_y_inv*Y(:,:,ff);           % 
    Aa_sigma_y_Aa=[];
    for i=1:length(source_power)
        Aa_sigma_y_Aa(i)=Aa_sigma_y_inv(i,:)*Aa_(:,i,ff);
    end
    sigma_x_ii(:,ff)=diag(gamma)-diag(gamma).* Aa_sigma_y_Aa.'.*diag(gamma);  %
     YA_trace=YA_trace+(((norm(Y(:,:,ff)-Aa_(:,:,ff)*Mu_x(:,:,ff),'fro'))^2)+L*trace(B(:,:,ff)-B(:,:,ff)*sigma_y_inv(:,:,ff)*B(:,:,ff))); %
    %% M-SBL
    Mu_x_norm=Mu_x_norm+ sum(abs(Mu_x(:,:,ff)).^2, 2);  %对Mu_x逐行平方求和
    end
    source_power=(L*( -F+sum(real(sigma_x_ii),2)./source_power)+sqrt(L^2*( -F+sum(real(sigma_x_ii),2)./source_power).^(2)+4*b*Mu_x_norm))/(2*b)+es;
    
    %% noise variance update1:
    alpha0=(F*M*L)/YA_trace; %
    
    [~,spd_]=sort(source_power,'descend');
    spd=spd_(1:min((M-1),length(source_power)-1));
    source_id=spd; %
        source_power_=source_power;
        source_power_old_=source_power_old;
        gird_add=[];index_ad=[];
        for i=1:length(source_id)
            if source_id(i)>1 &&source_id(i)<length(Azimuth_grid_)
                gri_r=Azimuth_grid_(source_id(i)+1)-Azimuth_grid_(source_id(i)) ;
                gri_l=Azimuth_grid_(source_id(i))-Azimuth_grid_(source_id(i)-1) ;
                if gri_r>gri_f || gri_l>gri_f
                        if isempty(intersect(gird_add,(Azimuth_grid_(source_id(i)+1)+Azimuth_grid_(source_id(i)))/2) ) && isempty(intersect(gird_add,(Azimuth_grid_(source_id(i))+Azimuth_grid_(source_id(i)-1))/2) )    
                    for ff=1:F
                    if source_power(source_id(i)+1)>source_power(source_id(i)-1)
                        Sigyk_inv(:,:,ff)=(sigma_y(:,:,ff)-Aa_(:,source_id(i):source_id(i)+1,ff)*diag([source_power(source_id(i):source_id(i)+1)])*Aa_(:,source_id(i):source_id(i)+1,ff)')^(-1);
                    else  
                        Sigyk_inv(:,:,ff)=(sigma_y(:,:,ff)-Aa_(:,source_id(i)-1:source_id(i),ff)*diag([source_power(source_id(i)-1:source_id(i))])*Aa_(:,source_id(i)-1:source_id(i),ff)')^(-1);
                    end
                  SRS(:,:,ff)=L*Sigyk_inv(:,:,ff)*Ry(:,:,ff)* Sigyk_inv(:,:,ff);
                    end
                    %%  搜索最大值进行插值
                        ag=[];sp=[];L_L=[];
                       ag=Azimuth_grid_(source_id(i)-1): (Azimuth_grid_(source_id(i)+1)-Azimuth_grid_(source_id(i)-1))/5:Azimuth_grid_(source_id(i)+1);
                   if source_power(source_id(i)+1)>source_power(source_id(i)-1)
                           ag=setdiff(ag,[Azimuth_grid_(source_id(i)), Azimuth_grid_(source_id(i)+1)]); %不要包括已有网格点
                   else  
                           ag=setdiff(ag,[Azimuth_grid_(source_id(i)-1),Azimuth_grid_(source_id(i))]); %不要包括已有网格点
                   end
                   aa_s=[];Qi_s=[];Si=[];
                   for ff=1:F
                       aa_s(:,:,ff)=exp(1j*2*pi*f0(ff)*D'*sin(  ag*pi/180)/c);
                        Qi_s(:,ff)=diag((aa_s(:,:,ff)'* SRS(:,:,ff)*aa_s(:,:,ff))); Si(:,ff)=diag((aa_s(:,:,ff)'*Sigyk_inv(:,:,ff)*aa_s(:,:,ff)));
                  
                     sp(:,ff)= ( - 2*b -L*Si(:,ff)+sqrt( ( (b*2+L*Si(:,ff) )).^2-b*4*(b+L*Si(:,ff)-Qi_s(:,ff) ))) ./(b*2*Si(:,ff) );
                      end
                        L_L= sum(real( -L*log( ( 1+sp.*Si))+( Qi_s)./( sp.^(-1)+ Si )-b.*sp),2);    % 
                  
                        [agr_max]=find((L_L)==max(L_L));
                        theta_s=ag(agr_max);
                        
      % grid and power add
                     if  theta_s>=Azimuth_grid_(source_id(i)) && (Azimuth_grid_(source_id(i)+1)-Azimuth_grid_(source_id(i)))>gri_f
                        if isempty(intersect(gird_add,Azimuth_grid_(source_id(i))+(Azimuth_grid_(source_id(i)+1)-Azimuth_grid_(source_id(i)))/2 ))
                            index_ad=[index_ad,source_id(i)];
                            index_ad=sort(index_ad,'ascend');
                            ia=find(index_ad==source_id(i));
                            i1=ia-1;
                            source_power_=[source_power_(1:source_id(i)+i1-1);source_power_(source_id(i)+i1)/2; source_power_(source_id(i)+i1)/2;source_power_(source_id(i)+i1+1:end)];
                            source_power_old_=[source_power_old_(1:source_id(i)+i1-1);source_power_old_(source_id(i)+i1)/2; source_power_old_(source_id(i)+i1)/2;source_power_old_(source_id(i)+i1+1:end)];
                            gird_add=union(gird_add,Azimuth_grid_(source_id(i))+(Azimuth_grid_(source_id(i)+1)-Azimuth_grid_(source_id(i)))/2);
                        end
                        
                     elseif  theta_s<Azimuth_grid_(source_id(i)) && Azimuth_grid_(source_id(i))-Azimuth_grid_(source_id(i)-1)>gri_f
                        if isempty(intersect(gird_add,Azimuth_grid_(source_id(i))-(Azimuth_grid_(source_id(i))-Azimuth_grid_(source_id(i)-1))/2))
                            index_ad=[index_ad,source_id(i)];
                            index_ad=sort(index_ad,'ascend');
                            ia=find(index_ad==source_id(i));
                            i1=ia-1;
                            source_power_=[source_power_(1:source_id(i)+i1-1);source_power_(source_id(i)+i1)/2; source_power_(source_id(i)+i1)/2;source_power_(source_id(i)+i1+1:end)];
                            source_power_old_=[source_power_old_(1:source_id(i)+i1-1);source_power_old_(source_id(i)+i1)/2; source_power_old_(source_id(i)+i1)/2;source_power_old_(source_id(i)+i1+1:end)];
                            gird_add=union(gird_add,Azimuth_grid_(source_id(i))-(Azimuth_grid_(source_id(i))-Azimuth_grid_(source_id(i)-1))/2);
                        end
                     end           
                     end
                end
            end
        end
        
        if ~isempty(gird_add)
            Azimuth_grid_=union(Azimuth_grid_,gird_add);
        end
        Aa_=[];
         for ff=1:F
        Aa_(:,:,ff)=exp(1j*2*pi*f0(ff)*D'*sin(Azimuth_grid_*pi/180)/c);
         end

        source_power=[];
        source_power_old=[];
        source_power=source_power_;
        source_power_old=source_power_old_;
        
        Length_add=Length_add+length(gird_add);
end

j2=j1;

while  true
    j2=j2+1;                  % number of iter
    source_power_old=source_power;
    gamma=diag(source_power);
    Mu_x_norm=zeros(length(source_power),1);
    Mu_x=[];sigma_x_ii=[];
    for ff=1:F
    B(:,:,ff)=Aa_(:,:,ff)*gamma*Aa_(:,:,ff)';
    sigma_y(:,:,ff)=1/alpha0*eye(M)+B(:,:,ff);
    sigma_y_inv(:,:,ff)=inv(sigma_y(:,:,ff));
    Aa_sigma_y_inv=Aa_(:,:,ff)'* sigma_y_inv(:,:,ff);
    Mu_x(:,:,ff)=gamma*Aa_sigma_y_inv*Y(:,:,ff);           % eq.(11)
    Aa_sigma_y_Aa=[];
    for i=1:length(source_power)
        Aa_sigma_y_Aa(i)=Aa_sigma_y_inv(i,:)*Aa_(:,i,ff);
    end
    sigma_x_ii(:,ff)=diag(gamma)-diag(gamma).* Aa_sigma_y_Aa.'.*diag(gamma);  %
    %% M-SBL
    Mu_x_norm=Mu_x_norm+ sum(abs(Mu_x(:,:,ff)).^2, 2);  %对Mu_x逐行平方求和
    end
    source_power=(L*( -F+sum(real(sigma_x_ii),2)./source_power)+sqrt(L^2*( -F+sum(real(sigma_x_ii),2)./source_power).^(2)+4*b*Mu_x_norm))/(2*b)+es;
      source_power_=source_power;
     source_id=[]; 
   source_power1=[];source_power2=[];source_power3=[];
   source_power1 = source_power_(1:end-2);
   source_power2 = source_power_(2:end-1);
   source_power3 = source_power_(3:end);
   IX=[]; PeakLoc0=[];SortIdx=[];SigLoc=[];source_id=[];SigLoc_=[];source_id_=[];
   IX = find((source_power2 > source_power1) + (source_power2 >= source_power3) > 1);
   PeakLoc0 = IX+1;   %PeakLoc0 is index of spectrum
  [~,SortIdx] = sort(source_power_(PeakLoc0),'descend');%descend the spectrum，SortIdx 
   %% noise
   KP = min(length(PeakLoc0),K);
%     KP = min(length(PeakLoc0),M-1);
   SigLoc = PeakLoc0(SortIdx(1:KP));
   source_id_=sort(SigLoc,'ascend');
   signal_id=source_id_;
   

   KP_=length(signal_id);
   Ac=[];
   Ac=[Aa_(:,signal_id,:)]; 
   P=[]; 
   AP_RY=0;
   for ff=1:F
   P(:,:,ff)=Ac(:,:,ff)*(Ac(:,:,ff)'*Ac(:,:,ff))^(-1)*Ac(:,:,ff)';                           
   AP_RY=AP_RY+abs(trace((eye(M)-P(:,:,ff))*Ry(:,:,ff)));
   end
   alpha0=F*(M-K)/AP_RY;
    Noise_variance=alpha0^(-1);
    if  norm((source_power-source_power_old),2)/norm(source_power_old,2)<tol1  || j2>500   %
        break
    end
end
iter=j1+j2;
ds=1;
%%  refined DOA original:
source_id=[];
source_power1=[];source_power2=[];source_power3=[];
source_power1 = source_power(1:end-2);
source_power2 = source_power(2:end-1);
source_power3 = source_power(3:end);
IX=[]; PeakLoc0=[];SortIdx=[];SigLoc=[];source_id=[];SigLoc_=[];source_id_=[];
IX = find((source_power2 > source_power1) + (source_power2 >= source_power3) > 1);
PeakLoc0 = IX+1;   %PeakLoc0 is index of spectrum
[~,SortIdx] = sort(source_power(PeakLoc0),'descend');%descend the spectrum，SortIdx 
% KP = min(length(PeakLoc0),M-1);
  KP = min(length(PeakLoc0),K);
SigLoc = PeakLoc0(SortIdx(1:KP));
source_id=sort(SigLoc,'ascend');
Azimuth=Azimuth_grid_;
theta_r1=[];
Azimuth_refine1=[];
theta_refine1=[]; % 存放精确化估计的角度
gamma_k=[];
sigma_yk_inv=[];
r_step=0.1;
for i=1:length(source_id)
%     %%
    ar=0;
    %% 考虑谱泄露问题，去除两个能量幅值。
    Aa_k=[];qs=[];
for ff=1:F
    if real(source_power (source_id(i)-1)) < real(source_power(source_id(i)+1))
        Aa_k= [Aa_(:, source_id(i):source_id(i)+1,ff)];         % exclude the angle index
        gamma_k=source_power;
        gamma_k= [gamma_k(source_id(i):source_id(i)+1)];
        gamma_k=diag(gamma_k);
        sigma_yk=sigma_y(:,:,ff)-Aa_k* gamma_k*Aa_k';    %
        sigma_yk_inv(:,:,ff)=inv( sigma_yk);
        qs(:,:,ff)=sigma_yk_inv(:,:,ff)*Ry(:,:,ff)*sigma_yk_inv(:,:,ff)*L;
    else
        Aa_k= [Aa_(:, source_id(i)-1:source_id(i),ff)];         % exclude the angle index
        gamma_k=source_power;
        gamma_k= [gamma_k(source_id(i)-1:source_id(i))];
        gamma_k=diag(gamma_k);
        sigma_yk=sigma_y(:,:,ff)-Aa_k* gamma_k*Aa_k';    %
        sigma_yk_inv(:,:,ff)=inv( sigma_yk);
        qs(:,:,ff)=sigma_yk_inv(:,:,ff)*Ry(:,:,ff)*sigma_yk_inv(:,:,ff)*L;
    end
end
    source_power_rere=[];theta_r1=[];
    for agr= Azimuth(source_id(i)-1) :r_step:Azimuth(source_id(i)+1)  % 
        ar=ar+1;
        aa=[];qk=[];zk=[];
        for ff=1:F
         aa(:,ff)=exp(1j*2*pi*f0(ff)*D'*sin(agr*pi/180)/c); % 
        qk(:,ff)=aa(:,ff)'*qs(:,:,ff)*aa(:,ff);zk(:,ff)=aa(:,ff)'*sigma_yk_inv(:,:,ff)*aa(:,ff);
    
       source_power_rere(ar,ff)= ( -2*b-L*zk(:,ff)+sqrt( ((2*b+L*zk(:,ff)))^2-4*b*(b+L*zk(:,ff)-qk(:,ff) )) )/(2*b*zk(:,ff) );   % 
         end
        theta_r1(ar)= sum(real(-L*log( ( 1+( source_power_rere(ar,:)).*zk))+( qk)./( ( source_power_rere(ar,:)).^(-1)+ zk )-b*(source_power_rere(ar,:))) );    
    end
    [agr_max(i)]=find((theta_r1)==max(theta_r1));
    Azimuth_refine=[];
    Azimuth_refine= Azimuth(source_id(i)-1) :r_step:Azimuth(source_id(i)+1) ;
    theta_refine(i)=Azimuth_refine(agr_max(i));
        Azimuth(source_id(i))=theta_refine(i);
end

end