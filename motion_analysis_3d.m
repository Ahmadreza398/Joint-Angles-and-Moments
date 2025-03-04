clc
clear
close all

%% load motion data
[file path]=uigetfile('*.xlsx');
data=xlsread([path file]);
frame=data(:,1)-data(1,1);
FR=120; % frame rate, Hz
t=frame/FR; %time, sec
n=length(t);
% pelvis markers, meters
LASIS=0.001*data(:,72:74);
RASIS=0.001*data(:,75:77);
LPSIS=0.001*data(:,78:80);
RPSIS=0.001*data(:,81:83);
% knee markers
RLKNE=0.001*data(:,36:38);
RMKNE=0.001*data(:,111:113);
% ankle markers
RLANK=0.001*data(:,42:44);
RMANK=0.001*data(:,114:116);
% foot markers
RMET1=0.001*data(:,126:128);
RMET5=0.001*data(:,78:80);
RTO=0.001*data(:,48:50);
RHE=0.001*data(:,45:47);

% tracking markers
% thigh markers
RTHI=0.001*data(:,33:35);
RTHIBU=0.001*data(:,51:53);
RTHIFU=0.001*data(:,54:56);
% shank markers
RTIB=0.001*data(:,39:41);
RTIBBU=0.001*data(:,63:65);
RTIBFU=0.001*data(:,66:68);

%% load force data
Fdata=xlsread([path file],2);
index=10:10:10*n;
Fdata=Fdata(index,:);
GRF=-Fdata(:,3:5); Fx=GRF(:,1); Fy=GRF(:,2); Fz=GRF(:,3);
Mz=0.001*Fdata(:,8); 
COP=0.001*Fdata(:,9:11); Cx=COP(:,1); Cy=COP(:,2); Cz=COP(:,3);
g=[0 0 -9.81];

%% segment LCS

Opelvis=0.5*(RASIS+LASIS); %pelvis center
% Opelvis=0.25*(RASIS+LASIS+RPSIS+LPSIS); %pelvis center
OkneeR=0.5*(RMKNE+RLKNE); %knee center
OankleR=0.5*(RMANK+RLANK); %ankle center
OheelR=RHE; %heel center

OthighR=(RTHI+RTHIBU+RTHIFU)/3;
OshankR=(RTIB+RTIBBU+RTIBFU)/3;

for i=1:n
    % pelvis LCS
    i_pel(i,:)=(RASIS(i,:)-Opelvis(i,:))/norm(RASIS(i,:)-Opelvis(i,:));
    v_pel(i,:)=(Opelvis(i,:)-0.5*(RPSIS(i,:)+LPSIS(i,:)))/norm(Opelvis(i,:)...
        -0.5*(RPSIS(i,:)+LPSIS(i,:)));
    k_pel(i,:)=cross(i_pel(i,:),v_pel(i,:));
    j_pel(i,:)=cross(k_pel(i,:),i_pel(i,:));
    Rpel(:,:,i)=[i_pel(i,:);j_pel(i,:);k_pel(i,:)];
    
    % hip joint center
    OhipR1(:,i)=[0.36*norm(RASIS(i,:)-LASIS(i,:));...
        -0.19*norm(RASIS(i,:)-LASIS(i,:));-0.3*norm(RASIS(i,:)-LASIS(i,:))];
    OhipL1(:,i)=[-0.36*norm(RASIS(i,:)-LASIS(i,:));...
        -0.19*norm(RASIS(i,:)-LASIS(i,:));-0.3*norm(RASIS(i,:)-LASIS(i,:))];
    OhipR(i,:)=Rpel(:,:,i)'*OhipR1(:,i)+Opelvis(i,:)';
    OhipL(i,:)=Rpel(:,:,i)'*OhipL1(:,i)+Opelvis(i,:)';
    
    %thigh LCS
    k_th(i,:)=(OhipR(i,:)-0.5*(RLKNE(i,:)+RMKNE(i,:)))/norm(OhipR(i,:)-0.5*(RLKNE(i,:)+RMKNE(i,:)));
    v_th(i,:)=(RLKNE(i,:)-RMKNE(i,:))/norm(RLKNE(i,:)-RMKNE(i,:));
    j_th(i,:)=cross(k_th(i,:),v_th(i,:));
    i_th(i,:)=cross(j_th(i,:),k_th(i,:));
    Rth(:,:,i)=[i_th(i,:);j_th(i,:);k_th(i,:)];
    
    % thigh LCS tracking
    [i_thT(i,:) j_thT(i,:) k_thT(i,:)]=local_axis(RTHIBU(i,:),RTHIFU(i,:),RTHI(i,:));
    RthT(:,:,i)=[i_thT(i,:);j_thT(i,:);k_thT(i,:)];
    RthT_th(:,:,i)=Rth(:,:,i)'*RthT(:,:,i);
    
    %shank LCS
    k_sh(i,:)=(OkneeR(i,:)-0.5*(RLANK(i,:)+RMANK(i,:)))/norm(OkneeR(i,:)-0.5*(RLANK(i,:)+RMANK(i,:)));
    v_sh(i,:)=(RLANK(i,:)-RMANK(i,:))/norm(RLANK(i,:)-RMANK(i,:));
    j_sh(i,:)=cross(k_sh(i,:),v_sh(i,:));
    i_sh(i,:)=cross(j_sh(i,:),k_sh(i,:)); 
    Rsh(:,:,i)=[i_sh(i,:);j_sh(i,:);k_sh(i,:)];
    
    % shank LCS tracking
    [i_shT(i,:) j_shT(i,:) k_shT(i,:)]=local_axis(RTIBBU(i,:),RTIBFU(i,:),RTIB(i,:));
    RshT(:,:,i)=[i_shT(i,:);j_shT(i,:);k_shT(i,:)];
    RshT_sh(:,:,i)=Rsh(:,:,i)'*RshT(:,:,i);
    
    %foot LCS
    k_fo(i,:)=(OankleR(i,:)-0.5*(RMET1(i,:)+RMET5(i,:)))/norm(OankleR(i,:)-0.5*(RMET1(i,:)+RMET5(i,:)));
    v_fo(i,:)=(RLANK(i,:)-RMANK(i,:))/norm(RLANK(i,:)-RMANK(i,:));
    j_fo(i,:)=cross(k_fo(i,:),v_fo(i,:));
    i_fo(i,:)=cross(j_fo(i,:),k_fo(i,:));
    Rfo(:,:,i)=[i_fo(i,:);j_fo(i,:);k_fo(i,:)];
    
    % foot LCS2 (HEEL)
    j_he(i,:)=(OheelR(i,:)-RTO(i,:))/norm(OheelR(i,:)-RTO(i,:));
    v_he(i,:)=(0.5*(RLANK(i,:)+RMANK(i,:))-OheelR(i,:))/norm(0.5*(RLANK(i,:)+RMANK(i,:))-OheelR(i,:));
    i_he(i,:)=cross(j_he(i,:),v_he(i,:));
    k_he(i,:)=cross(i_he(i,:),j_he(i,:));
    Rhe(:,:,i)=[i_he(i,:);j_he(i,:);k_he(i,:)];
end

%% segment and joint angular calculations
for i=1:n
    % pelvis angular calculations
    [a_pel(i,:) b_pel(i,:) g_pel(i,:)]=EC_angles(Rpel(:,:,i));
    [alpha_pel(i,:) beta_pel(i,:) gama_pel(i,:)]=angles(Rpel(:,:,i));
    % thigh angular calculations
    [a_th(i,:) b_th(i,:) g_th(i,:)]=EC_angles(Rth(:,:,i));
    [alpha_th(i,:) beta_th(i,:) gama_th(i,:)]=angles(Rth(:,:,i));
    % shank angular calculations
    [a_sh(i,:) b_sh(i,:) g_sh(i,:)]=EC_angles(Rsh(:,:,i));
    [alpha_sh(i,:) beta_sh(i,:) gama_sh(i,:)]=angles(Rsh(:,:,i));
    % foot angular calculations
    [a_fo(i,:) b_fo(i,:) g_fo(i,:)]=EC_angles(Rfo(:,:,i));
    [alpha_fo(i,:) beta_fo(i,:) gama_fo(i,:)]=angles(Rfo(:,:,i));

    % heel-toe angular calculations
    [a_he(i,:) b_he(i,:) g_he(i,:)]=EC_angles(Rhe(:,:,i));
    [alpha_he(i,:) beta_he(i,:) gama_he(i,:)]=angles(Rhe(:,:,i));

    % hip joint angular calculation
    [a_h(i,:) b_h(i,:) g_h(i,:)]=EC_angles(Rpel(:,:,i)'*Rth(:,:,i));
    [alpha_h(i,:) beta_h(i,:) gama_h(i,:)]=angles(Rpel(:,:,i)'*Rth(:,:,i));

    % knee joint angular calculation
    [a_k(i,:) b_k(i,:) g_k(i,:)]=EC_angles(RthT(:,:,i)'*RshT(:,:,i));
    [alpha_k(i,:) beta_k(i,:) gama_k(i,:)]=angles(RthT(:,:,i)'*RshT(:,:,i));
    
    %knee joint angle (tracking)
    [a_kT(i,:) b_kT(i,:) g_kT(i,:)]=EC_angles(RthT(:,:,i)'*Rsh(:,:,i));
    [alpha_kT(i,:) beta_kT(i,:) gama_kT(i,:)]=angles(Rth(:,:,i)'*Rsh(:,:,i));
        
    % ankle joint angular calculation
    [a_a(i,:) b_a(i,:) g_a(i,:)]=EC_angles(Rsh(:,:,i)'*Rfo(:,:,i));
    [alpha_a(i,:) beta_a(i,:) gama_a(i,:)]=angles(Rsh(:,:,i)'*Rfo(:,:,i));
    
    % pelvis helical angle
    teta_pel(i,:)=asind(0.5*sqrt((Rpel(2,3,i)-Rpel(3,2,i))^2+...
        (Rpel(1,3,i)-Rpel(3,1,i))^2+(Rpel(1,2,i)-Rpel(2,1,i))^2));
    % thigh helical angle
    teta_th(i,:)=asind(0.5*sqrt((Rth(2,3,i)-Rth(3,2,i))^2+...
        (Rth(1,3,i)-Rth(3,1,i))^2+(Rth(1,2,i)-Rth(2,1,i))^2));
end
ang_th=pi/180*[smooth(alpha_th) smooth(beta_th) smooth(gama_th)]; % radian
ang_sh=pi/180*[smooth(alpha_sh) smooth(beta_sh) smooth(gama_sh)];
ang_fo=pi/180*[smooth(alpha_fo) smooth(beta_fo) smooth(gama_fo)];

ang_h=pi/180*[smooth(alpha_h) smooth(beta_h) smooth(gama_h)];
ang_k=pi/180*[smooth(alpha_k) smooth(beta_k) smooth(gama_k)];
ang_a=pi/180*[smooth(alpha_a) smooth(beta_a) smooth(gama_a)];

% angular velocity and acceleration (GCS)
for i=1:3
    % segment
    [omg_th(:,i) alfa_th(:,i)]=va(t,ang_th(:,i));
    [omg_sh(:,i) alfa_sh(:,i)]=va(t,ang_sh(:,i));
    [omg_fo(:,i) alfa_fo(:,i)]=va(t,ang_fo(:,i));
    % joint
    [omg_h(:,i) alfa_h(:,i)]=va(t,ang_h(:,i));
    [omg_k(:,i) alfa_k(:,i)]=va(t,ang_k(:,i));
    [omg_a(:,i) alfa_a(:,i)]=va(t,ang_a(:,i));
end

%% segment mass and center of mass (COM)
% segment mass
M=75; %total body mass, Kg
m_fo=0.0145*M; %foot mass
m_sh=0.0465*M; %foot mass
m_th=0.1*M; %foot mass

% segment length
for i=1:n
    L_th(i,:)=norm(OkneeR(i,:)-OhipR(i,:));
    L_sh(i,:)=norm(OankleR(i,:)-OkneeR(i,:));
    L_fo(i,:)=norm(RTO(i,:)-OankleR(i,:));
end
% segment radii
for i=1:n
    Rp_th(i,:)=norm(OhipR(i,:)-RASIS(i,:));
    Rd_th(i,:)=norm(OkneeR(i,:)-RLKNE(i,:));
    x_th(i,:)=Rd_th(i,:)/Rp_th(i,:);
    Rp_sh(i,:)=Rd_th(i,:);
    Rd_sh(i,:)=norm(OankleR(i,:)-RLANK(i,:));
    x_sh(i,:)=Rd_sh(i,:)/Rp_sh(i,:);
    Rp_fo(i,:)=Rd_sh(i,:);
    Rd_fo(i,:)=0.5*norm(RMET1(i,:)-RMET5(i,:));
    x_fo(i,:)=Rd_fo(i,:)/Rp_fo(i,:);
end

% COM coefficient (c)
C=@(x) (1+2*x+3*x^2)/(4*(1+x+x^2));
c=@(x) C(x)*(x<=1)+(1-C(x))*(x>1);
% for i=1:n
%     c_th(i,:)=COM(x_th(i));
%     c_sh(i,:)=COM(x_sh(i));
%     c_fo(i,:)=COM(x_fo(i));
% end

for i=1:n
    c_th(i,:)=c(x_th(i));
    r_th(:,i)=[0;0;-c_th(i)*L_th(i)]; %LCS
    rcm_th(i,:)=Rth(:,:,i)'*r_th(:,i)+OhipR(i,:)';
    c_sh(i,:)=c(x_sh(i));
    r_sh(:,i)=[0;0;-c_sh(i)*L_sh(i)];
    rcm_sh(i,:)=Rsh(:,:,i)'*r_sh(:,i)+OkneeR(i,:)';    
    c_fo(i,:)=c(x_fo(i));
    r_fo(:,i)=[0;0;-c_fo(i)*L_fo(i)];
    rcm_fo(i,:)=Rfo(:,:,i)'*r_fo(:,i)+OankleR(i,:)';    
end

% segment COM velocity and acceleration
for i=1:3
   [vcm_th(:,i) acm_th(:,i)]=va(t,rcm_th(:,i));
   [vcm_sh(:,i) acm_sh(:,i)]=va(t,rcm_sh(:,i));
   [vcm_fo(:,i) acm_fo(:,i)]=va(t,rcm_fo(:,i));
end

%% segment moment of inertia
for i=1:n
    I_th(:,:,i)=MOI(m_th,L_th(i),Rp_th(i),Rd_th(i));
    I_sh(:,:,i)=MOI(m_sh,L_sh(i),Rp_sh(i),Rd_sh(i));
    I_fo(:,:,i)=MOI(m_fo,L_fo(i),Rp_fo(i),Rd_fo(i));
end    

%% net joint force
for i=1:n
    Fa(i,:)=m_fo*(acm_fo(i,:)-g)-GRF(i,:);
    Fk(i,:)=m_sh*(acm_sh(i,:)-g)+Fa(i,:);
    Fh(i,:)=m_th*(acm_th(i,:)-g)+Fk(i,:);
end

%% net joint moment
for i=1:n
    
    % GCS to LCS for angular velocity and acceleration
    omgfo(:,i)=Rfo(:,:,i)*omg_fo(i,:)';
    alfafo(:,i)=Rfo(:,:,i)*alfa_fo(i,:)';
    omgsh(:,i)=Rsh(:,:,i)*omg_sh(i,:)';
    alfash(:,i)=Rsh(:,:,i)*alfa_sh(i,:)';
    omgth(:,i)=Rth(:,:,i)*omg_th(i,:)';
    alfath(:,i)=Rth(:,:,i)*alfa_th(i,:)';    
    
    % inertial torque (moment)
    Tinfo(:,i)=I_fo(:,:,i)*alfafo(:,i)+cross(omgfo(:,i),I_fo(:,:,i)*omgfo(:,i));
    Tin_fo(:,i)=Rfo(:,:,i)'*Tinfo(:,i); %LCS to GCS
    Tinsh(:,i)=I_sh(:,:,i)*alfash(:,i)+cross(omgsh(:,i),I_sh(:,:,i)*omgsh(:,i));
    Tin_sh(:,i)=Rsh(:,:,i)'*Tinsh(:,i); %LCS to GCS
    Tinth(:,i)=I_th(:,:,i)*alfath(:,i)+cross(omgth(:,i),I_th(:,:,i)*omgth(:,i));
    Tin_th(:,i)=Rth(:,:,i)'*Tinth(:,i); %LCS to GCS 
    
    % GRF torque (force plate)
    Tz(i)=Mz(i)-Cx(i)*Fy(i)-Cy(i)*Fx(i);
    Tgrf(:,i)=[0;0;Tz(i)];
    
    % ankle torque (moment) distances
    ra_f(:,i)=rcm_fo(i,:)'-OankleR(i,:)';
    ra_grf(:,i)=COP(i,:)'-OankleR(i,:)';
    
    % ankle joint torque (moment)
    Ta(i,:)=Tin_fo(:,i)-Tgrf(:,i)-cross((ra_grf(:,i)-ra_f(:,i)),GRF(i,:)')...
        +cross(ra_f(:,i),Fa(i,:)');
    
    % knee torque (moment) distances
    rk_sh(:,i)=rcm_sh(i,:)'-OkneeR(i,:)';
    rk_a(:,i)=OankleR(i,:)'-OkneeR(i,:)';
    
    % knee joint torque (moment)
    Tk(i,:)=Tin_sh(:,i)+Ta(i,:)'+cross(rk_sh(:,i),Fk(i,:)')...
        +cross((rk_a(:,i)- rk_sh(:,i)),Fa(i,:)');
    
    % hip torque (moment) distances
    rh_th(:,i)=rcm_th(i,:)'-OhipR(i,:)';
    rh_k(:,i)=OkneeR(i,:)'-OhipR(i,:)';
    
    % hip joint torque (moment)
    Th(i,:)=Tin_th(:,i)+Tk(i,:)'+cross(rh_th(:,i),Fh(i,:)')...
        +cross((rh_k(:,i)-rh_th(:,i)),Fk(i,:)');
end   
    
%% joint power
for i=1:n
    Pa(i,:)=dot(Ta(i,:),omg_a(i,:)); % sum of ankle powers
    P_a(i,:)=Ta(i,:).*omg_a(i,:); % ankle powers (Px,Py,Pz)
    Pk(i,:)=dot(Tk(i,:),omg_k(i,:)); % sum of knee powers
    P_k(i,:)=Tk(i,:).*omg_k(i,:); % knee powers (Px,Py,Pz)
    Ph(i,:)=dot(Th(i,:),omg_h(i,:)); % sum of hip powers
    P_h(i,:)=Th(i,:).*omg_h(i,:); % hip powers (Px,Py,Pz)
end

%% showing motion
plotx=[RASIS(:,1) RPSIS(:,1) LASIS(:,1) LPSIS(:,1) RLKNE(:,1) RMKNE(:,1) RLANK(:,1) RMANK(:,1) RMET1(:,1) RMET5(:,1) RTO(:,1) RHE(:,1)];
ploty=[RASIS(:,2) RPSIS(:,2) LASIS(:,2) LPSIS(:,2) RLKNE(:,2) RMKNE(:,2) RLANK(:,2) RMANK(:,2) RMET1(:,2) RMET5(:,2) RTO(:,2) RHE(:,2)];
plotz=[RASIS(:,3) RPSIS(:,3) LASIS(:,3) LPSIS(:,3) RLKNE(:,3) RMKNE(:,3) RLANK(:,3) RMANK(:,3) RMET1(:,3) RMET5(:,3) RTO(:,3) RHE(:,3)];

figure('WindowState','Maximized');
for i=1:n
plot3(plotx(i,:),ploty(i,:),plotz(i,:),'*',...
    [OhipR(i,1) OhipL(i,1)],[OhipR(i,2) OhipL(i,2)],[OhipR(i,3) OhipL(i,3)],'o');
hold on;
% segment COM
plot3([rcm_th(i,1) rcm_sh(i,1) rcm_fo(i,1)],[rcm_th(i,2) rcm_sh(i,2) rcm_fo(i,2)]...
    ,[rcm_th(i,3) rcm_sh(i,3) rcm_fo(i,3)],'blacko')
% pelvis LCS
plot_axis(Opelvis(i,:),i_pel(i,:),j_pel(i,:),k_pel(i,:),RASIS(i,:),LASIS(i,:),0.5*(RPSIS(i,:)+LPSIS(i,:)));
% thigh LCS
plot_axis(OhipR(i,:),i_th(i,:),j_th(i,:),k_th(i,:),OhipR(i,:),RLKNE(i,:),RMKNE(i,:));
% thigh LCS tracking
plot_axis(OthighR(i,:),i_thT(i,:),j_thT(i,:),k_thT(i,:),RTHI(i,:),RTHIBU(i,:),RTHIFU(i,:));

% shank LCS
plot_axis(OkneeR(i,:),i_sh(i,:),j_sh(i,:),k_sh(i,:),OkneeR(i,:),RLANK(i,:),RMANK(i,:));
% shank LCS tracking
plot_axis(OshankR(i,:),i_shT(i,:),j_shT(i,:),k_shT(i,:),RTIB(i,:),RTIBBU(i,:),RTIBFU(i,:));

% foot LCS
plot_axis(OankleR(i,:),i_fo(i,:),j_fo(i,:),k_fo(i,:),OankleR(i,:),RMET1(i,:),RMET5(i,:));
% GRF
quiver3(Cx(i),Cy(i),Cz(i),0.001*Fx(i),0.001*Fy(i),0.001*Fz(i));
axis equal; view([1 0.5 1]);
pause(0.01);
hold off;
end

%% plot angular results
% Euler-Cardan angles
figure('WindowState','Maximized');
subplot(2,2,1)
plot(t,[a_pel b_pel g_pel]); 
title('pelvis Euler-Cardan angles, deg');
legend('alpha','beta','gama'); xlabel('time, sec');
subplot(2,2,2)
plot(t,[a_th b_th g_th]); 
title('thigh Euler-Cardan angles, deg');
legend('alpha','beta','gama'); xlabel('time, sec');
subplot(2,2,3)
plot(t,[a_sh b_sh g_sh]); 
title('shank Euler-Cardan angles, deg');
legend('alpha','beta','gama'); xlabel('time, sec');
subplot(2,2,4)
plot(t,[a_fo b_fo g_fo]); 
title('foot Euler-Cardan angles, deg');
legend('alpha','beta','gama'); xlabel('time, sec');

figure('WindowState','Maximized');
subplot(3,3,1)
plot(t,[alpha_h beta_h gama_h]); 
title('hip angles, deg');
legend('alpha','beta','gama'); xlabel('time, sec');
subplot(3,3,2)
plot(t,[alpha_k beta_k gama_k]); 
title('knee angles, deg');
legend('alpha','beta','gama'); xlabel('time, sec');
subplot(3,3,3)
plot(t,[alpha_a beta_a gama_a]); 
title('ankle angles, deg');
legend('alpha','beta','gama'); xlabel('time, sec');

subplot(3,3,4)
plot(t,omg_h); 
title('hip angluar velocities, rad/s');
legend('omgx','omgy','omgz'); xlabel('time, sec');
subplot(3,3,5)
plot(t,omg_k); 
title('knee angluar velocities, rad/s');
legend('omgx','omgy','omgz'); xlabel('time, sec');
subplot(3,3,6)
plot(t,omg_a); 
title('ankle angluar velocities, rad/s');
legend('omgx','omgy','omgz'); xlabel('time, sec');

subplot(3,3,4)
plot(t,omg_h); 
title('hip angluar velocities, rad/s');
legend('omgx','omgy','omgz'); xlabel('time, sec');
subplot(3,3,5)
plot(t,omg_k); 
title('knee angluar velocities, rad/s');
legend('omgx','omgy','omgz'); xlabel('time, sec');
subplot(3,3,6)
plot(t,omg_a); 
title('ankle angluar velocities, rad/s');
legend('omgx','omgy','omgz'); xlabel('time, sec');

subplot(3,3,7)
plot(t,alfa_h); 
title('hip angluar accelerations, rad/s^2');
legend('alfax','alfay','alfaz'); xlabel('time, sec');
subplot(3,3,8)
plot(t,alfa_k); 
title('knee angluar accelerations, rad/s^2');
legend('alfax','alfay','alfaz'); xlabel('time, sec');
subplot(3,3,9)
plot(t,alfa_a); 
title('ankle angluar accelerations, rad/s^2');
legend('alfax','alfay','alfaz'); xlabel('time, sec');

%% plot joint force results
figure('WindowState','Maximized');
subplot(1,3,1)
plot(t,Fh); title('Hip joint force, N');
legend('Fx','Fy','Fz'); xlabel('time, sec');
axis([0 1.5 -1200 200]);
subplot(1,3,2)
plot(t,Fk); title('Knee joint force, N');
legend('Fx','Fy','Fz'); xlabel('time, sec');
axis([0 1.5 -1200 200]);
subplot(1,3,3)
plot(t,Fa); title('Ankle joint force, N');
legend('Fx','Fy','Fz'); xlabel('time, sec');
axis([0 1.5 -1200 200]);

%% plot joint moments and powers

figure('WindowState','Maximized');
subplot(2,3,1)
plot(t,Th); title('hip joint moment, N.m');
legend('Mx','My','Mz'); xlabel('time, sec');
axis([0.2 1.2 -200 200]);
subplot(2,3,2)
plot(t,Tk); title('Knee joint moment, N.m');
legend('Mx','My','Mz'); xlabel('time, sec');
axis([0.2 1.2 -200 200]);
subplot(2,3,3)
plot(t,Ta); title('Ankle joint moment, N.m');
legend('Mx','My','Mz'); xlabel('time, sec');
axis([0.2 1.2 -200 200]);
subplot(2,3,4)
plot(t,P_h); title('hip joint power, Watts');
legend('Px','Py','Pz'); xlabel('time, sec');
axis([0.2 1.2 -200 200]);
subplot(2,3,5)
plot(t,P_k); title('Knee joint power, Watts');
legend('Px','Py','Pz'); xlabel('time, sec');
axis([0.2 1.2 -200 200]);
subplot(2,3,6)
plot(t,P_a); title('Ankle joint power, Watts');
legend('Px','Py','Pz'); xlabel('time, sec');
axis([0.2 1.2 -200 200]);