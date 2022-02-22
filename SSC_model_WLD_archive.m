clear all
load('anuga_data_multi_wind')
runid='model_test'

t_full=vertcat(t_full,inf);
timesteps=44000;
gridsize=10;

dt=10;

timestart=0;
time_counter=min(find(t_full>timestart));
threshold_h=.005; %meters, conc=0 if below this depth
[depth,xvel,yvel,elev,Xgrid2,Ygrid2]=get_gridded_data(u_full,v_full,h_full,x_full,y_full,elev_full,time_counter);

F_xvel = scatteredInterpolant(Xgrid2(:),Ygrid2(:),xvel(:),'linear','nearest'); %interpolate linearly, extrapolate nearest neighbor
F_yvel = scatteredInterpolant(Xgrid2(:),Ygrid2(:),yvel(:),'linear','nearest'); %interpolate linearly, extrapolate nearest neighbor
F_depth = scatteredInterpolant(Xgrid2(:),Ygrid2(:),depth(:),'linear','nearest'); %interpolate linearly, extrapolate nearest neighbor
F_eta= scatteredInterpolant(Xgrid2(:),Ygrid2(:),elev(:),'linear','nearest'); 
x=[min(Xgrid2(:))+gridsize:gridsize:max(Xgrid2(:))];
y=[min(Ygrid2(:))+gridsize:gridsize:max(Ygrid2(:))];
[Xgrid,Ygrid]=meshgrid(x,y);
u=F_xvel(Xgrid,Ygrid);
v=F_yvel(Xgrid,Ygrid);
h=F_depth(Xgrid,Ygrid);
eta=F_eta(Xgrid,Ygrid);
[cf,manning]=compute_friction(h,eta);
uedge=.5*u(:,1:end-1)+.5*u(:,2:end);
vedge=.5*v(1:end-1,:)+.5*v(2:end,:);
Nx=length(x);
Ny=length(y);
c=zeros(Ny,Nx);

r_x=zeros(size(uedge));
r_y=zeros(size(vedge));
flux_x=zeros(size(uedge));
flux_y=zeros(size(vedge));


cnew=zeros(Ny,Nx);
ctemp=zeros(Ny,Nx);
ctemp2=zeros(Ny,Nx);
E=zeros(Ny,Nx);
w_s=1*10^-4;


c(isnan(h))=NaN;
c(isnan(u))=NaN;
h(isnan(u))=NaN;

uedge(find(isnan(uedge)))=0;
vedge(find(isnan(vedge)))=0;


source_idx=find(Xgrid>6.524*10^5 &Xgrid<6.529*10^5&Ygrid>3.27084*10^6);
source_conc=1*10^-5;

save_interval=1000;
c_save=zeros(size(c,1),size(c,2),floor(timesteps/save_interval));
h_save=zeros(size(c,1),size(c,2),floor(timesteps/save_interval));
u_save=zeros(size(c,1),size(c,2),floor(timesteps/save_interval));
v_save=zeros(size(c,1),size(c,2),floor(timesteps/save_interval));
frict_save=zeros(size(c,1),size(c,2),floor(timesteps/save_interval));



for k=1:timesteps
    if dt*k+timestart>t_full(time_counter+1)
        time_counter=time_counter+1;
        h_old=h;
        [depth,xvel,yvel,elev,Xgrid2,Ygrid2]=get_gridded_data(u_full,v_full,h_full,x_full,y_full,elev_full,time_counter);
        F_xvel = scatteredInterpolant(Xgrid2(:),Ygrid2(:),xvel(:),'linear','nearest'); %interpolate linearly, extrapolate nearest neighbor
        F_yvel = scatteredInterpolant(Xgrid2(:),Ygrid2(:),yvel(:),'linear','nearest'); %interpolate linearly, extrapolate nearest neighbor
        F_depth = scatteredInterpolant(Xgrid2(:),Ygrid2(:),depth(:),'linear','nearest'); %interpolate linearly, extrapolate nearest neighbor

        u=F_xvel(Xgrid,Ygrid);
        v=F_yvel(Xgrid,Ygrid);
        h=F_depth(Xgrid,Ygrid);
        [cf,manning]=compute_friction(h,eta);
        uedge=.5*u(:,1:end-1)+.5*u(:,2:end);
        vedge=.5*v(1:end-1,:)+.5*v(2:end,:);
        k
        c=c.*h_old./h;
        c(h<=threshold_h)=0;
    end

    if rem(k,2)==0
    %x direction
    for j=3:Ny-2
        for i=2:Nx-2
            if uedge(j,i)>0
                r_x(j,i)=(c(j,i)-c(j,i-1))/(c(j,i+1)-c(j,i));
            else
                r_x(j,i)=(c(j,i+2)-c(j,i+1))/(c(j,i+1)-c(j,i));
            end
        end
    end
    %van leer flux limiter
    thi_x=(r_x+abs(r_x))./(1+abs(r_x));
    thi_x(find(isnan(thi_x)))=0;
    for j=3:Ny-2
        for i=2:Nx-2
            flux_x(j,i)=.5*uedge(j,i)*((1+sign(uedge(j,i)))*c(j,i)+(1-sign(uedge(j,i)))*c(j,i+1))+...
                .5*abs(uedge(j,i))*(1-abs(uedge(j,i)*dt/gridsize))*thi_x(j,i)*(c(j,i+1)-c(j,i));
        end
    end
    flux_x(find(isnan(flux_x)))=0;

    ctemp(:,2:end-1)=c(:,2:end-1)+dt*(flux_x(:,1:end-1)-flux_x(:,2:end))/gridsize...
        -dt*c(:,2:end-1).*(uedge(:,1:end-1)-uedge(:,2:end))/gridsize;
    ctemp(2,:)=ctemp(3,:); ctemp(1,:)=ctemp(3,:); ctemp(:,1)=ctemp(:,3); ctemp(:,2)=ctemp(:,3); ctemp(Ny,:)=ctemp(Ny-2,:); ctemp(Ny-1,:)=ctemp(Ny-2,:); ctemp(:,Nx)=ctemp(:,Nx-2); ctemp(:,Nx-1)=ctemp(:,Nx-2);
    ctemp(source_idx)=source_conc;
    %y direction
    for j=2:Ny-2
        for i=3:Nx-2
            if vedge(j,i)>0
                r_y(j,i)=(ctemp(j,i)-ctemp(j-1,i))/(ctemp(j+1,i)-ctemp(j,i));
            else
                r_y(j,i)=(ctemp(j+2,i)-ctemp(j+1,i))/(ctemp(j+1,i)-ctemp(j,i));
            end
        end
    end
    %van leer flux limiter
    thi_y=(r_y+abs(r_y))./(1+abs(r_y));
    thi_y(find(isnan(thi_y)))=0;

    for j=2:Ny-2
        for i=3:Nx-2
            flux_y(j,i)=.5*vedge(j,i)*((1+sign(vedge(j,i)))*ctemp(j,i)+(1-sign(vedge(j,i)))*ctemp(j+1,i))+...
                .5*abs(vedge(j,i))*(1-abs(vedge(j,i)*dt/gridsize))*thi_y(j,i)*(ctemp(j+1,i)-ctemp(j,i));
        end
    end
    flux_y(find(isnan(flux_y)))=0;

    ctemp2(2:end-1,:)=ctemp(2:end-1,:)+dt*(flux_y(1:end-1,:)-flux_y(2:end,:))/gridsize...
        -dt*ctemp(2:end-1,:).*(vedge(1:end-1,:)-vedge(2:end,:))/gridsize;
    else
    %y direction
    for j=2:Ny-2
        for i=3:Nx-2
            if vedge(j,i)>0
                r_y(j,i)=(c(j,i)-c(j-1,i))/(c(j+1,i)-c(j,i));
            else
                r_y(j,i)=(c(j+2,i)-c(j+1,i))/(c(j+1,i)-c(j,i));
            end
        end
    end
    %van leer flux limiter
    thi_y=(r_y+abs(r_y))./(1+abs(r_y));
    thi_y(find(isnan(thi_y)))=0;

    for j=2:Ny-2
        for i=3:Nx-2
            flux_y(j,i)=.5*vedge(j,i)*((1+sign(vedge(j,i)))*c(j,i)+(1-sign(vedge(j,i)))*c(j+1,i))+...
                .5*abs(vedge(j,i))*(1-abs(vedge(j,i)*dt/gridsize))*thi_y(j,i)*(c(j+1,i)-c(j,i));
        end
    end
    flux_y(find(isnan(flux_y)))=0;

    ctemp(2:end-1,:)=c(2:end-1,:)+dt*(flux_y(1:end-1,:)-flux_y(2:end,:))/gridsize...
        -dt*c(2:end-1,:).*(vedge(1:end-1,:)-vedge(2:end,:))/gridsize;
    ctemp(2,:)=ctemp(3,:); ctemp(1,:)=ctemp(3,:); ctemp(:,1)=ctemp(:,3); ctemp(:,2)=ctemp(:,3); ctemp(Ny,:)=ctemp(Ny-2,:); ctemp(Ny-1,:)=ctemp(Ny-2,:); ctemp(:,Nx)=ctemp(:,Nx-2); ctemp(:,Nx-1)=ctemp(:,Nx-2);
    ctemp(source_idx)=source_conc;
     %x direction
    for j=3:Ny-2
        for i=2:Nx-2
            if uedge(j,i)>0
                r_x(j,i)=(ctemp(j,i)-ctemp(j,i-1))/(ctemp(j,i+1)-ctemp(j,i));
            else
                r_x(j,i)=(ctemp(j,i+2)-ctemp(j,i+1))/(ctemp(j,i+1)-ctemp(j,i));
            end
        end
    end
    %van leer flux limiter
    thi_x=(r_x+abs(r_x))./(1+abs(r_x));
    thi_x(find(isnan(thi_x)))=0;

    for j=3:Ny-2
        for i=2:Nx-2
            flux_x(j,i)=.5*uedge(j,i)*((1+sign(uedge(j,i)))*ctemp(j,i)+(1-sign(uedge(j,i)))*ctemp(j,i+1))+...
                .5*abs(uedge(j,i))*(1-abs(uedge(j,i)*dt/gridsize))*thi_x(j,i)*(ctemp(j,i+1)-ctemp(j,i));
        end
    end
    flux_x(find(isnan(flux_x)))=0;
    ctemp2(:,2:end-1)=ctemp(:,2:end-1)+dt*(flux_x(:,1:end-1)-flux_x(:,2:end))/gridsize...
        -dt*ctemp(:,2:end-1).*(uedge(:,1:end-1)-uedge(:,2:end))/gridsize;
    end
    ctemp2(source_idx)=source_conc;
    ctemp2(2,:)=ctemp2(3,:); ctemp2(1,:)=ctemp2(3,:); ctemp2(:,1)=ctemp2(:,3); ctemp2(:,2)=ctemp2(:,3); ctemp2(Ny,:)=ctemp2(Ny-2,:); ctemp2(Ny-1,:)=ctemp2(Ny-2,:); ctemp2(:,Nx)=ctemp2(:,Nx-2); ctemp2(:,Nx-1)=ctemp2(:,Nx-2);


    cnew(2:end-1,2:end-1)=ctemp2(2:end-1,2:end-1)+(dt./(h(2:end-1,2:end-1))).*(-ctemp2(2:end-1,2:end-1))*w_s;
    cnew(find(cnew<0))=0;
    cnew(h<=threshold_h)=0;
    c=cnew;
    c(source_idx)=source_conc;
    c(2,:)=c(3,:); c(1,:)=c(3,:); c(:,1)=c(:,3); c(:,2)=c(:,3); c(Ny,:)=c(Ny-2,:); c(Ny-1,:)=c(Ny-2,:); c(:,Nx)=c(:,Nx-2); c(:,Nx-1)=c(:,Nx-2);
   
    if rem(k,save_interval)==0;

        c_save(:,:,k/save_interval)=c;
        h_save(:,:,k/save_interval)=h;
        u_save(:,:,k/save_interval)=u;
        v_save(:,:,k/save_interval)=v;

        k
    end
    if rem(k,2000)==0
        save([runid], '-v7.3')
    end
end

t_save=[timestart+dt*save_interval:dt*save_interval:timestart+dt*save_interval*length(c_save(1,1,:))];
save([runid], '-v7.3')
