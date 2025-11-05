function u=z_nodes(Coord,t)
for j=1:size(Coord,1),
    curcoords=[Coord(j,:),t];
    u(j)=ze(curcoords);
end