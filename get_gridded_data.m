function [depth,xvel,yvel,elev,Xgrid,Ygrid]=get_gridded_data(u,v,h,x_full,y_full,elev,row)
%% grid cell averaging
xmin=6.38*10^5-20;
xmax=6.62*10^5+20;
ymin=3.252*10^6-20;
ymax=3.271*10^6+20;
gridsize=20;
x=x_full;
y=y_full;
xedge=[xmin:gridsize:xmax];
yedge=[ymin:gridsize:ymax];

depth_avg=nan(length(yedge)-1,length(xedge)-1);
u_avg=nan(length(yedge)-1,length(xedge)-1);
v_avg=nan(length(yedge)-1,length(xedge)-1);
elev_avg=nan(length(yedge)-1,length(xedge)-1);
[Xgrid,Ygrid]=meshgrid(xedge(1:end-1)+gridsize/2,yedge(1:end-1)+gridsize/2);

h=h(row,:);
u=u(row,:);
v=v(row,:);

for i=1:length(xedge)-1 %average data points that lie within the same grid cell
    for j=1:length(yedge)-1
        f=find(x>=xedge(i)&x<xedge(i+1)&y>=yedge(j)&y<yedge(j+1));
        if ~isempty(f)
            depth_avg(j,i)=nanmean(h(f)); 
            u_avg(j,i)=nanmean(u(f)./h(f));
            v_avg(j,i)=nanmean(v(f)./h(f));
            elev_avg(j,i)=nanmean(elev(f));

        end
    end
    
end

%% interpolate

f=find(~isnan(u_avg));
F_xvel = scatteredInterpolant(Xgrid(f),Ygrid(f),u_avg(f),'linear','nearest'); %interpolate linearly, extrapolate nearest neighbor
F_yvel = scatteredInterpolant(Xgrid(f),Ygrid(f),v_avg(f),'linear','nearest'); %interpolate linearly, extrapolate nearest neighbor
F_elev = scatteredInterpolant(Xgrid(f),Ygrid(f),elev_avg(f),'linear','nearest'); %interpolate linearly, extrapolate nearest neighbor
F_depth = scatteredInterpolant(Xgrid(f),Ygrid(f),depth_avg(f),'linear','nearest'); %interpolate linearly, extrapolate nearest neighbor

depth=F_depth(Xgrid,Ygrid);
xvel=F_xvel(Xgrid,Ygrid);
yvel=F_yvel(Xgrid,Ygrid);
elev=F_elev(Xgrid,Ygrid);

end