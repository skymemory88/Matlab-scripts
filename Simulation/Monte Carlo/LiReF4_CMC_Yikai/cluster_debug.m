ipos = params.pos;
isp = eSpin;
ipos(cluster,:) = [];
isp(cluster,:) = [];
figure
quiver3(ipos(:,1), ipos(:,2), ipos(:,3), isp(:,1), isp(:,2), isp(:,3),'k');
pbaspect([1 1 int16(params.dims(3)/params.dims(1))])
hold on
quiver3(params.pos(cluster,1), params.pos(cluster,2), params.pos(cluster,3),...
    eSpin(cluster,1), eSpin(cluster,2), eSpin(cluster,3),'r');
pbaspect([1 1 int16(params.dims(3)/params.dims(1))])