clear all
close all

M=10;  %阵元个数, element number
d=1;  %阵元间距半波长
c=1500;               % 声速, sound speed
f0=750;   %载频, frequency
snr=0;              % 信噪比
L=30; % 快拍数,number of snapshots
error=0;error1=0;
test=400; % 试验统计,number of independent test
timetol1=0;itertol=0;
noise_error=0;
i0=0;
Azimuth_G=[];
SPALL=[];
K_error=0;
for i=1:test
    i0=i0+1;
    grid_interval=12;  % 网格间距, grid interval
    Azimuth_grid=-90:grid_interval:90;         % all hypothetical angles
    % 随机角度集1：
      u=unifrnd(-1,1);
    AL=[ -7+u,3+u,10+u,30+u];
    % fixed DOAs：
%     AL=[ -7.63,3.25,10.32, 30.83];
    K=size(AL,2);         % 信源数量，signal number
    %% ULA array manifold
    A=[];
    for m=1:M
        for k=1:K
            A(m,k)=exp(1j*2*pi*d*(m-1)*f0/c*sin(AL(1,k)*pi/180));
        end
    end
    
    %% the ULA array completed manifold matrix :
    for m=1:M
        for ang= 1:length(Azimuth_grid)
            Az(ang)=Azimuth_grid(ang)*pi/180;
            Aa(m,ang)=exp(1j*2*pi*f0*(m-1)/c*sin(Az(ang)));
        end
    end
    
    %% 1.Uncorrelated signal : Guassian
    S=(randn(K,L)+1j*randn(K,L));     % 产生信号，基带处理
    Vj=diag(sqrt(  10^(snr/10)./diag(1/L*(S*S') ) ) );
    S=Vj*S;     % 产生信号, signal
    noise=sqrt(1/2)*(randn(M,L)+1j*randn(M,L));  % 加噪声, noise
    Y=A*S+noise;
    
    tstart1 = tic;
    %% AGR_SBL:
%     [Source_power_re1,Azimuth1,~,iter,Noise_variance]=AGRSBL(Y,Azimuth_grid,grid_interval,Aa,f0,c,K); % K is known
        [Source_power_re1,Azimuth1,~,iter,Noise_variance]=AGR_SBL_(Y,Azimuth_grid,grid_interval,Aa,f0,c); % K is unknown
    %
    Azimuth_G(i0,1:length(Source_power_re1))=Azimuth1;
    %     SPALL(i0,1:length(Source_power_re1))=Source_power_re1;
    time1=toc(tstart1);
    timetol1=timetol1+time1;  % total time
    [~, peakindex1] = findpeaks(real(Source_power_re1),'sortstr','descend');
    KP=min(length(peakindex1),K);
    source_id1=sort( peakindex1(1:KP),'ascend');
    DOA_estimate1=Azimuth1(source_id1);
    
    itertol=itertol+iter; %总迭代次数
    error_sum1=0;
    if length( DOA_estimate1)>=K
        for ii=1:length(AL)
            error_sum1=error_sum1+((DOA_estimate1(ii)-AL(ii)))^2;
        end
    else %极端情况，信号估计DOA小于真实DOA个数
        for ii=1:length(DOA_estimate1)
            error_sum1=error_sum1+((DOA_estimate1(ii)-AL(ii)))^2;
            
        end
        for ii=length(DOA_estimate1)+1:K
            error_sum1=error_sum1+((DOA_estimate1(length(DOA_estimate1))-AL(ii)))^2;
        end
    end
    error1=error1+error_sum1/length(AL);
    

end
    hold on
    power1=(real(Source_power_re1)/max(real(Source_power_re1)));
    plot(Azimuth1,(real(power1)),'g-','Linewidth',0.5);

    hold on
    scatter(AL(1),1);
    scatter(AL(2),1);
    scatter(AL(3),1);
    scatter(AL(4),1);
RMSE1=sqrt(error1/test); % root mean square error
itert_m=itertol/test;
ylim([0 1])
%% 归一化作图
ylim([0 1])
box on
xlim([-90 90])
legend('AGRSBL')
set(gca,'FontName','Helvetica','FontSize',16);
set(gca,'linewidth',1,'FontSize',20);
xlabel('DOA (degree)');ylabel('Unified spectrum','fontsize',20);

%% _AGR_SBL(K is unknown):
function [source_power,Azimuth,source_id,iterall,Noise_variance]=AGR_SBL_(Y,Azimuth_grid,grid_interval,Aa,f0,c)
%% input:
%% Y: Array output data
%% Azimuth_grid: Space discontinue Azimuth angle grid
%% grid_interval: Interval of the adjacent angle grid
%% Aa: Compeleted array steering matrix

%% output:
%% Source_power_re:signal space power
%% Azimuth:All space angle,include DOAgrid_interval
%% iterall:iter number in pre-estimate
%% source_id_: DOA index in Azimuth
[M,L]=size(Y);
b=10^(-7);
ro_=1;
es=10^(-12);
a_=1;
Ry=[];
Ry=(Y*Y')/L;     %样本协方差矩阵
tol1=0.001;    %迭代停止门限
Azimuth_grid_=Azimuth_grid;
N=length(Azimuth_grid);
% noise_variance1=10^(-2)*mean(var(Y)); %噪声初始化
noise_variance=10^(-2)*(norm(Y,'fro'))^2/(M*L);
alpha0=1/noise_variance;   %噪声方差倒数
source_power = sum(abs(Aa'*Y),2)/(L*M); % 最小二乘初始化
Aa_=Aa;
jmax=500;
e=1;
j1=0;
g=5;
Length_add=0;
gri_f=1.5; % 最小细分间距
I_=ceil( log2(grid_interval/gri_f) )+3; % AGR max iteration number
while (j1<=I_) %迭代次
    j1=j1+1;                  % number of iter
    source_power_old=source_power;
    gamma=diag(source_power);
    B=Aa_*gamma*Aa_';
    sigma_y=1/alpha0*eye(M)+B;
    sigma_y_inv=inv(sigma_y);
    Aa_sigma_y_inv=Aa_'* sigma_y_inv;
    Mu_x=gamma*Aa_sigma_y_inv*Y;           % 
    
    Aa_sigma_y_Aa=[];
    for i=1:length(source_power)
        Aa_sigma_y_Aa(i)=Aa_sigma_y_inv(i,:)*Aa_(:,i);
    end
    sigma_x_ii=diag(gamma)-diag(gamma).* Aa_sigma_y_Aa.'.*diag(gamma);  %
    %% M-SBL
    Mu_x_norm= sum(abs(Mu_x).^2, 2);  %对Mu_x逐行平方求和
    %% real Laplace
    source_power=(L*( -1+(real(sigma_x_ii))./source_power)+sqrt(L^2*( -1+(real(sigma_x_ii))./source_power).^(2)+4*b*Mu_x_norm))/(2*b)+es;
    % noise variance update1:
    alpha0=(M*L)/(((norm(Y-Aa_*Mu_x,'fro'))^2)+L*trace(B-B*sigma_y_inv*B)); % 如果此trace报错，那是因为和Bellhop的trace产生冲突，把Bellhop的相关路径删除后解决问题
    %% noise etimation
    [~,spd_]=sort(source_power,'descend');
    spd=spd_(1:min((M-1),length(source_power)));
    source_id=spd; % M-1个幅值，

    %% AGR,process
    source_power_=source_power;
    source_power_old_=source_power_old;
    gird_add=[];index_ad=[];
    for i=1:length(source_id)
        if source_id(i)>1 &&source_id(i)<length(Azimuth_grid_)
            gri_r=Azimuth_grid_(source_id(i)+1)-Azimuth_grid_(source_id(i)) ;
            gri_l=Azimuth_grid_(source_id(i))-Azimuth_grid_(source_id(i)-1) ;
            if gri_r>=2*gri_f || gri_l>=2*gri_f
               if isempty(intersect(gird_add,(Azimuth_grid_(source_id(i)+1)+Azimuth_grid_(source_id(i)))/2) ) && isempty(intersect(gird_add,(Azimuth_grid_(source_id(i))+Azimuth_grid_(source_id(i)-1))/2) )    
                    if source_power(source_id(i)+1)>source_power(source_id(i)-1)
                        Sigyk_inv(:,:,i)=(sigma_y-Aa_(:,source_id(i):source_id(i)+1)*diag([source_power(source_id(i):source_id(i)+1)])*Aa_(:,source_id(i):source_id(i)+1)')^(-1);
                    else
                        Sigyk_inv(:,:,i)=(sigma_y-Aa_(:,source_id(i)-1:source_id(i))*diag([source_power(source_id(i)-1:source_id(i))])*Aa_(:,source_id(i)-1:source_id(i))')^(-1);
                    end
                    
                    %%  搜索最大值进行插值
                    ag=[];sp=[];L_L=[];
                    ag=Azimuth_grid_(source_id(i)-1): (Azimuth_grid_(source_id(i)+1)-Azimuth_grid_(source_id(i)-1))/g :Azimuth_grid_(source_id(i)+1);
                    ag=setdiff(ag,[Azimuth_grid_(source_id(i)-1),Azimuth_grid_(source_id(i)),Azimuth_grid_(source_id(i)+1)]); %不要包括已有网格点
                    aa_s=exp(1j*2*pi*f0*[0:M-1]'*sin( ag*pi/180)/c);
                    Qi_s=sum( abs(( aa_s'*Sigyk_inv(:,:,i)*Y )).^2,2);
                    Si=real(diag((aa_s'*Sigyk_inv(:,:,i)*aa_s)));
                    sp= ( - 2*b -L*Si+sqrt( ( (L*Si)).^2+b*4*(Qi_s))) ./(b*2*Si);
                    L_L= ( -L*log( ( 1+sp.*Si))+( Qi_s)./( sp.^(-1)+ Si )...
                        -b.*sp);    % eq.(30)
                    [agr_max]=find((L_L)==max(L_L));
                    theta_s=ag(agr_max);
                    if  theta_s>=Azimuth_grid_(source_id(i)) && (Azimuth_grid_(source_id(i)+1)-Azimuth_grid_(source_id(i)))/2>=gri_f
                        if isempty(intersect(gird_add,Azimuth_grid_(source_id(i))+(Azimuth_grid_(source_id(i)+1)-Azimuth_grid_(source_id(i)))/2 ))
                            index_ad=[index_ad,source_id(i)];
                            index_ad=sort(index_ad,'ascend');
                            ia=find(index_ad==source_id(i));
                            i1=ia-1;
                            source_power_=[source_power_(1:source_id(i)+i1-1);source_power_(source_id(i)+i1)/2; source_power_(source_id(i)+i1)/2;source_power_(source_id(i)+i1+1:end)];
                            source_power_old_=[source_power_old_(1:source_id(i)+i1-1);source_power_old_(source_id(i)+i1)/2; source_power_old_(source_id(i)+i1)/2;source_power_old_(source_id(i)+i1+1:end)];
                            gird_add=union(gird_add,Azimuth_grid_(source_id(i))+(Azimuth_grid_(source_id(i)+1)-Azimuth_grid_(source_id(i)))/2);
                        end
                        
                    elseif  theta_s<Azimuth_grid_(source_id(i)) && (Azimuth_grid_(source_id(i))-Azimuth_grid_(source_id(i)-1))/2>=gri_f
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
    
    Aa_=exp(1j*pi*[0:M-1]'*sin(Azimuth_grid_*pi/180));
    source_power=[];
    source_power_old=[];
    source_power=source_power_;
    source_power_old=source_power_old_;
    
    Length_add=Length_add+length(gird_add);
    
    
end


%%
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


j2=j1;
j3=0;
Noise_variance=alpha0^(-1);
while  true
    j3=j3+1;
    j2=j2+1;                  % number of iter
    source_power_old=source_power;
    gamma=diag(source_power);
    B=Aa_*gamma*Aa_';
    sigma_y=1/alpha0*eye(M)+B;
    sigma_y_inv=inv(sigma_y);
    Aa_sigma_y_inv=Aa_'* sigma_y_inv;
    Mu_x=gamma*Aa_sigma_y_inv*Y;           %
    for i=1:length(source_power)
        Aa_sigma_y_Aa(i)=Aa_sigma_y_inv(i,:)*Aa_(:,i);
    end
    sigma_x_ii=diag(gamma)-diag(gamma).* Aa_sigma_y_Aa.'.*diag(gamma);  %
    %% M-SBL
    Mu_x_norm= sum(abs(Mu_x).^2, 2);  %对Mu_x逐行平方求和
    %% real Laplace
    source_power=(L*( -1+ (real(sigma_x_ii))./source_power)+sqrt(L^2*( -1+(real(sigma_x_ii))./source_power).^(2)+4*b.*Mu_x_norm))./(2*b)+es;
    %% noise etimation
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
    K_ = min(length(PeakLoc0),M-1);
    SigLoc = PeakLoc0(SortIdx(1:K_));
    
    source_id_=sort(SigLoc,'ascend');
    signal_id=source_id_;
    Ac=[Aa_(:,signal_id)];
    P=[];
    P=Ac*(Ac'*Ac)^(-1)*Ac';
    alpha0=1/abs(trace((eye(M)-P)*Ry)/(M-length(signal_id)));% eq.(18)
    Noise_variance=alpha0^(-1);
    if  norm((source_power-source_power_old),2)/norm(source_power_old,2)<tol1  || j2>500   % L2 norm
        break
    end
    
end
iterall=j2;
%% off-grid DOA estimtion,  refined DOA :
source_id=[];
source_power1=[];source_power2=[];source_power3=[];
source_power1 = source_power(1:end-2);
source_power2 = source_power(2:end-1);
source_power3 = source_power(3:end);
IX=[]; PeakLoc0=[];SortIdx=[];SigLoc=[];source_id=[];SigLoc_=[];source_id_=[];
IX = find((source_power2 > source_power1) + (source_power2 >= source_power3) > 1);
PeakLoc0 = IX+1;   %PeakLoc0 is index of spectrum
[~,SortIdx] = sort(source_power(PeakLoc0),'descend');%descend the spectrum，SortIdx 是谱峰降序顺序在这些谱峰里的下标
K_ = min(length(PeakLoc0),M-1);
SigLoc = PeakLoc0(SortIdx(1:K_));
source_id=sort(SigLoc,'ascend');
Azimuth=Azimuth_grid_;
theta_r1=[];
Azimuth_refine1=[];
theta_refine1=[]; % 存放精确化估计的角度
gamma_k=[];
sigma_yk_inv=[];
r_step=0.2;
start2=tic;
Aa_re_=Aa_;
gamma_re1=diag(source_power);
%
sigma_y_re1=sigma_y;
for i=1:length(source_id)
    %%
    Aa_k=[];
    if real(source_power (source_id(i)-1)) < real(source_power(source_id(i)+1))
        Aa_k= [Aa_(:, source_id(i):source_id(i)+1)];         %
        gamma_k=source_power;
        gamma_k= [gamma_k(source_id(i):source_id(i)+1)];
        gamma_k=diag(gamma_k);
    else
        Aa_k= [Aa_(:, source_id(i)-1:source_id(i))];         %
        gamma_k=source_power;
        gamma_k= [gamma_k(source_id(i)-1:source_id(i))];
        gamma_k=diag(gamma_k);
    end
    sigma_yk=sigma_y-Aa_k* gamma_k*Aa_k';    %
    sigma_yk_inv=inv( sigma_yk);
    sgY=sigma_yk_inv*Y;
    
    source_power_rere=[];theta_r1=[];
    ar=0;
    for agr= Azimuth(source_id(i)-1) :r_step:Azimuth(source_id(i)+1)  %
        ar=ar+1;
        aa=[];
        aa=exp(1j*2*pi*f0*[0:M-1]'*sin(agr*pi/180)/c); %
        qk=sum(abs(aa'*sgY).^2);zk=real(aa'*sigma_yk_inv*aa);
        source_power_rere(ar)= ( -2*b-L*zk+sqrt( ((2*b+L*zk))^2-4*b*(b+L*zk-qk)) )/(2*b*zk);    %
        theta_r1(ar)= ( -L*log( ( 1+( source_power_rere(ar))*zk))+( qk)/( ( source_power_rere(ar))^(-1)+ zk )...
            -b*(source_power_rere(ar)));    %
    end
    [agr_max(i)]=find((theta_r1)==max(theta_r1));
    Azimuth_refine=[];
    Azimuth_refine= Azimuth(source_id(i)-1) :r_step:Azimuth(source_id(i)+1) ;
    theta_refine(i)=Azimuth_refine(agr_max(i));
    %% 细分二次的精细化搜索:
    r_step2=0.05;
    ar=0;
    for agr2=  theta_refine(i)-r_step :r_step2: theta_refine(i)+r_step
        ar=ar+1;
        aa=[];
        aa=exp(1j*2*pi*f0*[0:M-1]'*sin(agr2*pi/180)/c);
        qk=sum(abs(aa'*sgY).^2);zk=real(aa'*sigma_yk_inv*aa);
        source_power_rere(ar)= ( -2*b-L*zk+sqrt( ((2*b+L*zk))^2-4*b*(b+L*zk-qk)) )/(2*b*zk);
        theta_r2(ar)= ( -L*log( ( 1+( source_power_rere(ar))*zk))+( qk)/( ( source_power_rere(ar))^(-1)+ zk )...
            -b*(source_power_rere(ar)));
        
    end
    [agr_max(i)]=find((theta_r2)==max(theta_r2));
    Azimuth_refine=[];
    Azimuth_refine= theta_refine(i)-r_step :r_step2: theta_refine(i)+r_step;
    theta_refine2(i)=Azimuth_refine(agr_max(i));
    Azimuth(source_id(i))=theta_refine2(i);
end

end

%% _AGR_SBL(K is known):
function [source_power,Azimuth,source_id,iterall,Noise_variance]=AGRSBL(Y,Azimuth_grid,grid_interval,Aa,f0,c,K)
%% input:
%% Y: Array output data
%% Azimuth_grid: Space discontinue Azimuth angle grid
%% grid_interval: Interval of the adjacent angle grid
%% Aa: Compeleted array steering matrix
%% K: Source number,if K is not a prior,we let K=M-1.

%% output:
%% Source_power_re:signal space power
%% Azimuth:All space angle,include DOAgrid_interval
%% iterall:iter number
%% source_id_: DOA index in Azimuth
[M,L]=size(Y);
b=10^(-7);
ro_=1;
es=10^(-10);
a_=1;
Ry=[];
Ry=(Y*Y')/L;     %样本协方差矩阵
tol1=0.001;    %迭代停止门限
Azimuth_grid_=Azimuth_grid;
N=length(Azimuth_grid);
% noise_variance=10^(-2)*mean(var(Y)); %噪声初始化
noise_variance=10^(-2)*(norm(Y,'fro'))^2/(M*L);
alpha0=1/noise_variance;   %噪声方差倒数
source_power = sum(abs(Aa'*Y),2)/(L*M); % 最小二乘初始化
Aa_=Aa;
jmax=500;
e=1;
j1=0;
g=5;
Length_add=0;
ai=1; %每个网格迭代次数
gri_f=grid_interval/2^3; %1.5
I_=ceil( log2(grid_interval/gri_f) )+3;  % AGR max iteration number
while (j1<=I_) %迭代次
    j1=j1+1;                  % number of iter
    source_power_old=source_power;
    gamma=diag(source_power);
    B=Aa_*gamma*Aa_';
    sigma_y=1/alpha0*eye(M)+B;
    sigma_y_inv=inv(sigma_y);
    Aa_sigma_y_inv=Aa_'* sigma_y_inv;
    Mu_x=gamma*Aa_sigma_y_inv*Y;           %
    
    Aa_sigma_y_Aa=[];
    for i=1:length(source_power)
        Aa_sigma_y_Aa(i)=Aa_sigma_y_inv(i,:)*Aa_(:,i);
    end
    sigma_x_ii=diag(gamma)-diag(gamma).* Aa_sigma_y_Aa.'.*diag(gamma);  %
    %% M-SBL
    Mu_x_norm= sum(abs(Mu_x).^2, 2);  %对Mu_x逐行平方求和
    %% real Laplace
    source_power=(L*( -1+(real(sigma_x_ii))./source_power)+sqrt(L^2*( -1+(real(sigma_x_ii))./source_power).^(2)+4*b*Mu_x_norm))/(2*b)+es;
    % noise variance update1:
    alpha0=(M*L)/(((norm(Y-Aa_*Mu_x,'fro'))^2)+L*trace(B-B*sigma_y_inv*B)); %
    %% AGR process:
    [~,spd_]=sort(source_power,'descend');
    spd=spd_(1:min((M-1),length(source_power)));
    source_id=spd; % M-1个幅值
    source_power_=source_power;
    source_power_old_=source_power_old;
    gird_add=[];index_ad=[];
    for i=1:length(source_id)
        if source_id(i)>1 &&source_id(i)<length(Azimuth_grid_)
            gri_r=Azimuth_grid_(source_id(i)+1)-Azimuth_grid_(source_id(i)) ;
            gri_l=Azimuth_grid_(source_id(i))-Azimuth_grid_(source_id(i)-1) ;
            if (gri_r>=2*gri_f || gri_l>=2*gri_f)
               if isempty(intersect(gird_add,(Azimuth_grid_(source_id(i)+1)+Azimuth_grid_(source_id(i)))/2) ) && isempty(intersect(gird_add,(Azimuth_grid_(source_id(i))+Azimuth_grid_(source_id(i)-1))/2) )    
                    %% 考虑谱泄露问题这里去除的是谱峰及其一侧次大值：
                    if source_power(source_id(i)+1)>source_power(source_id(i)-1)
                        Sigyk_inv(:,:,i)=(sigma_y-Aa_(:,source_id(i):source_id(i)+1)*diag([source_power(source_id(i):source_id(i)+1)])*Aa_(:,source_id(i):source_id(i)+1)')^(-1);
                    else
                        Sigyk_inv(:,:,i)=(sigma_y-Aa_(:,source_id(i)-1:source_id(i))*diag([source_power(source_id(i)-1:source_id(i))])*Aa_(:,source_id(i)-1:source_id(i))')^(-1);
                    end
                    
                    %%  搜索最大值进行插值
                    ag=[];sp=[];L_L=[];
                    ag=Azimuth_grid_(source_id(i)-1): (Azimuth_grid_(source_id(i)+1)-Azimuth_grid_(source_id(i)-1))/g :Azimuth_grid_(source_id(i)+1);
                    ag=setdiff(ag,[Azimuth_grid_(source_id(i)-1),Azimuth_grid_(source_id(i)),Azimuth_grid_(source_id(i)+1)]); %不要包括已有网格点
                    aa_s=exp(1j*2*pi*f0*[0:M-1]'*sin( ag*pi/180)/c);
                    Qi_s=sum( abs(( aa_s'*Sigyk_inv(:,:,i)*Y )).^2,2);
                    Si=real(diag((aa_s'*Sigyk_inv(:,:,i)*aa_s)));
                    sp= ( - 2*b -L*Si+sqrt( ( (L*Si)).^2+b*4*(Qi_s))) ./(b*2*Si);
                    L_L= ( -L*log( ( 1+sp.*Si))+( Qi_s)./( sp.^(-1)+ Si )...
                        -b.*sp);    % eq.(30)
                    [agr_max]=find((L_L)==max(L_L));
                    theta_s=ag(agr_max);
                    % grid and power add
                    if  theta_s>Azimuth_grid_(source_id(i)) && Azimuth_grid_(source_id(i)+1)-Azimuth_grid_(source_id(i))>gri_f
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
    
    Aa_=exp(1j*pi*[0:M-1]'*sin(Azimuth_grid_*pi/180));
    source_power=[];
    source_power_old=[];
    source_power=source_power_;
    source_power_old=source_power_old_;
    Azimuth_grid=Azimuth_grid_;
    Length_add=Length_add+length(gird_add);
    
    
end

j2=j1;
while  true
    j2=j2+1;                  % number of iter
    source_power_old=source_power;
    gamma=diag(source_power);
    B=Aa_*gamma*Aa_';
    sigma_y=1/alpha0*eye(M)+B;
    sigma_y_inv=inv(sigma_y);
    Aa_sigma_y_inv=Aa_'* sigma_y_inv;
    Mu_x=gamma*Aa_sigma_y_inv*Y;           %
    for i=1:length(source_power)
        Aa_sigma_y_Aa(i)=Aa_sigma_y_inv(i,:)*Aa_(:,i);
    end
    sigma_x_ii=diag(gamma)-diag(gamma).* Aa_sigma_y_Aa.'.*diag(gamma);  %
    %% M-SBL
    Mu_x_norm= sum(abs(Mu_x).^2, 2);  %对Mu_x逐行平方求和
    %% real Laplace
    source_power=(L*( -1+(real(sigma_x_ii))./source_power)+sqrt(L^2*( -1+(real(sigma_x_ii))./source_power).^(2)+4*b.*Mu_x_norm))./(2*b)+es;
    %% noise etimation
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
    % noise1
    KP = min(length(PeakLoc0),K);
    SigLoc = PeakLoc0(SortIdx(1:KP));
    source_id_=sort(SigLoc,'ascend');
    signal_id=source_id_;
    Ac=[Aa_(:,signal_id)];
    P=[];
    P=Ac*(Ac'*Ac)^(-1)*Ac';
    alpha0=1/abs(trace((eye(M)-P)*Ry)/(M-KP));
    Noise_variance=alpha0^(-1);
    if  norm((source_power-source_power_old),2)/norm(source_power_old,2)<tol1  || j2>500   % L2 norm
        break
    end
    
end
iterall=j2;
%% off-grid DOA estimtion:
source_id=[];
source_power1=[];source_power2=[];source_power3=[];
source_power1 = source_power(1:end-2);
source_power2 = source_power(2:end-1);
source_power3 = source_power(3:end);
IX=[]; PeakLoc0=[];SortIdx=[];SigLoc=[];source_id=[];SigLoc_=[];source_id_=[];
IX = find((source_power2 > source_power1) + (source_power2 >= source_power3) > 1);
PeakLoc0 = IX+1;   %PeakLoc0 is index of spectrum
[~,SortIdx] = sort(source_power(PeakLoc0),'descend');%descend the spectrum，SortIdx 是谱峰降序顺序在这些谱峰里的下标
KP = min(length(PeakLoc0),K);
SigLoc = PeakLoc0(SortIdx(1:KP));
source_id=sort(SigLoc,'ascend');
Azimuth=Azimuth_grid_;
theta_r1=[];
Azimuth_refine1=[];
theta_refine1=[]; % 存放精确化估计的角度
gamma_k=[];
sigma_yk_inv=[];
r_step=0.2;
start2=tic;
Aa_re_=Aa_;
gamma_re1=diag(source_power);
%
sigma_y_re1=sigma_y;
for i=1:length(source_id)
    %%
    Aa_k=[];
    if real(source_power (source_id(i)-1)) < real(source_power(source_id(i)+1))
        Aa_k= [Aa_(:, source_id(i):source_id(i)+1)];         %
        gamma_k=source_power;
        gamma_k= [gamma_k(source_id(i):source_id(i)+1)];
        gamma_k=diag(gamma_k);
    else
        Aa_k= [Aa_(:, source_id(i)-1:source_id(i))];         %
        gamma_k=source_power;
        gamma_k= [gamma_k(source_id(i)-1:source_id(i))];
        gamma_k=diag(gamma_k);
    end
    sigma_yk=sigma_y-Aa_k* gamma_k*Aa_k';    %
    sigma_yk_inv=inv( sigma_yk);
    sgY=sigma_yk_inv*Y;
    
    source_power_rere=[];theta_r1=[];
    ar=0;
    for agr= Azimuth(source_id(i)-1) :r_step:Azimuth(source_id(i)+1)  %
        ar=ar+1;
        aa=[];
        aa=exp(1j*2*pi*f0*[0:M-1]'*sin(agr*pi/180)/c); %
        qk=sum(abs(aa'*sgY).^2);zk=real(aa'*sigma_yk_inv*aa);
        source_power_rere(ar)= ( -2*b-L*zk+sqrt( ((2*b+L*zk))^2-4*b*(b+L*zk-qk)) )/(2*b*zk);    %
        theta_r1(ar)= ( -L*log( ( 1+( source_power_rere(ar))*zk))+( qk)/( ( source_power_rere(ar))^(-1)+ zk )...
            -b*(source_power_rere(ar)));    %
    end
    [agr_max(i)]=find((theta_r1)==max(theta_r1));
    Azimuth_refine=[];
    Azimuth_refine= Azimuth(source_id(i)-1) :r_step:Azimuth(source_id(i)+1) ;
    theta_refine(i)=Azimuth_refine(agr_max(i));
    %% 细分二次的精细化搜索:
    r_step2=0.05;
    ar=0;
    for agr2=  theta_refine(i)-r_step :r_step2: theta_refine(i)+r_step
        ar=ar+1;
        aa=[];
        aa=exp(1j*2*pi*f0*[0:M-1]'*sin(agr2*pi/180)/c);
        qk=sum(abs(aa'*sgY).^2);zk=real(aa'*sigma_yk_inv*aa);
        source_power_rere(ar)= ( -2*b-L*zk+sqrt( ((2*b+L*zk))^2-4*b*(b+L*zk-qk)) )/(2*b*zk);
        theta_r2(ar)= ( -L*log( ( 1+( source_power_rere(ar))*zk))+( qk)/( ( source_power_rere(ar))^(-1)+ zk )...
            -b*(source_power_rere(ar)));
        
    end
    [agr_max(i)]=find((theta_r2)==max(theta_r2));
    Azimuth_refine=[];
    Azimuth_refine= theta_refine(i)-r_step :r_step2: theta_refine(i)+r_step;
    theta_refine2(i)=Azimuth_refine(agr_max(i));
    Azimuth(source_id(i))=theta_refine2(i);
end

end
