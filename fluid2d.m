function fluid2d()

GRAVITY=0.01;
DT=0.1;
TIME_STEP_TOTAL=1000;
PARTICLE_NUM=10000;
GRID_H=100;
GRID_W=100;

function outbound = is_out_of_bound(pos)
    outbound=(pos(1)<0||pos(1)>GRID_W-1)||(pos(2)<0||pos(2)>GRID_H-1);
end

function pos = force_push_inbound_pos(pos)
    if pos(1)<0
        pos(1)=0.001;
    end
    if pos(1)>GRID_W-1
        pos(1)=GRID_W-1.001;
    end
    if pos(2)<0
        pos(2)=0.001;
    end
    if pos(2)>GRID_H-1
        pos(2)=GRID_H-1.001;
    end
end

function vel = force_push_inbound_vel(pos,vel)
    if pos(1)<0
        vel(1)=-vel(1)*0.5;
    end
    if pos(1)>GRID_W-1
        vel(1)=-vel(1)*0.5;
    end
    if pos(2)<0
        vel(2)=-vel(2)*0.5;
    end
    if pos(2)>GRID_H-1
        vel(2)=-vel(2)*0.5;
    end
end

%function p_v = interpolate_v(p_pos)
%    p_x=ceil(p_pos(1));
%    p_y=ceil(p_pos(2));
%end

clear global;
close all;

grid=zeros(GRID_W,GRID_H);
grid_v=zeros(GRID_W,GRID_H,2);
grid_v_copy=zeros(GRID_W,GRID_H,4);

particles=zeros(PARTICLE_NUM,2);
for i=1:PARTICLE_NUM
    particles(i,1)=rand*GRID_W;
    particles(i,2)=rand*GRID_H;
end
particle_v=zeros(PARTICLE_NUM,2);

for i=1:TIME_STEP_TOTAL
    % update particle velocity
    for j=1:PARTICLE_NUM
        %disp(particle_v(j,2))
        % add gravity
        particle_v(j,2)=particle_v(j,2)+GRAVITY*DT;
        % add grid velocity
        particle_v(j,1)=particle_v(j,1)+grid_v_copy(ceil(particles(j,2)),ceil(particles(j,1)),1)+grid_v_copy(ceil(particles(j,2)),ceil(particles(j,1)),2);
        particle_v(j,2)=particle_v(j,2)+grid_v_copy(ceil(particles(j,2)),ceil(particles(j,1)),3)+grid_v_copy(ceil(particles(j,2)),ceil(particles(j,1)),4);
    end
    % update position
    particles=particles+particle_v*DT;
    % clear grid
    grid=zeros(GRID_W,GRID_H);
    grid_v=zeros(GRID_W,GRID_H,2);
    grid_v_copy=zeros(GRID_W,GRID_H,4);
    % update particles in grid
    for j=1:PARTICLE_NUM
        if is_out_of_bound(particles(j,:))
            % if particle is out of bound, push the particle back
            particle_v(j,:)=force_push_inbound_vel(particles(j,:),particle_v(j,:));
            particles(j,:)=force_push_inbound_pos(particles(j,:));
        end
        % push particles to grid
        grid(ceil(particles(j,2)),ceil(particles(j,1)))=grid(ceil(particles(j,2)),ceil(particles(j,1)))+1;
        % push particle velocity to corresponding grid velocity
        grid_v(ceil(particles(j,2)),ceil(particles(j,1)),1)=grid_v(ceil(particles(j,2)),ceil(particles(j,1)),1)+particle_v(j,1);
        grid_v(ceil(particles(j,2)),ceil(particles(j,1)),2)=grid_v(ceil(particles(j,2)),ceil(particles(j,1)),2)+particle_v(j,2);
    end
    % make grid incompressible
    for j=1:GRID_H
        for k=1:GRID_W
            d=0;
            cnt=0;
            if j~=1
                grid_v_copy(j-1,k,3)=grid_v(j-1,k,2);
                d=d-grid_v(j-1,k,2);
                cnt=cnt+1;
            end
            if j~=GRID_H
                grid_v_copy(j+1,k,4)=grid_v(j+1,k,2);
                d=d+grid_v(j+1,k,2);
                cnt=cnt+1;
            end
            if k~=1
                grid_v_copy(j,k-1,1)=grid_v(j,k-1,1);
                d=d-grid_v(j,k-1,1);
                cnt=cnt+1;
            end
            if k~=GRID_W
                grid_v_copy(j,k+1,2)=grid_v(j,k+1,1);
                d=d+grid_v(j,k+1,1);
                cnt=cnt+1;
            end
            d=d/cnt;
            disp(d)
            if j~=1
                grid_v_copy(j-1,k,3)=grid_v_copy(j-1,k,3)+d;
            end
            if j~=GRID_H
                grid_v_copy(j+1,k,4)=grid_v_copy(j+1,k,4)-d;
            end
            if k~=1
                grid_v_copy(j,k-1,1)=grid_v_copy(j,k-1,1)+d;
            end
            if k~=GRID_W
                grid_v_copy(j,k+1,2)=grid_v_copy(j,k+1,2)-d;
            end
        end
    end
    
    imshow(grid/20);
    %h.GridVisible='off';
    %disp(max(grid_v_copy,[],'all'))
    drawnow

end

end
