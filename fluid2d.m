function fluid2d()

GRAVITY=-10;
DT=0.1;
CO=0.8;
NUM_ITER=50;
TIME_STEP_TOTAL=1000;
PARTICLE_NUM=1000;
GRID_H=100;
GRID_W=150;
PARTICLE_PER_GRID=40;

grid=zeros(GRID_H,GRID_W);
grid_type=ones(GRID_H,GRID_W);
grid_v_x=zeros(GRID_H,GRID_W+1);
grid_v_y=zeros(GRID_H+1,GRID_W);
grid_v_x_w=zeros(GRID_H,GRID_W+1);
grid_v_y_w=zeros(GRID_H+1,GRID_W);

particles=zeros(PARTICLE_NUM,2);
for i=1:PARTICLE_NUM
    particles(i,1)=rand*GRID_W;
    particles(i,2)=rand*GRID_H;
end
particle_v=zeros(PARTICLE_NUM,2);
particle_v_copy=zeros(PARTICLE_NUM,2);

function addGravity()
    particle_v(:,2)=particle_v(:,2)+GRAVITY*DT;
end

function outbound = is_out_of_bound(pos)
    outbound=(pos(1)<=0||pos(1)>=GRID_W)||(pos(2)<=0||pos(2)>=GRID_H);
end

function pos = force_push_inbound_pos(pos)
    if pos(1)<0
        pos(1)=0.001;
    end
    if pos(1)>GRID_W
        pos(1)=GRID_W-0.001;
    end
    if pos(2)<0
        pos(2)=0.001;
    end
    if pos(2)>GRID_H
        pos(2)=GRID_H-0.001;
    end
end

function vel = force_push_inbound_vel(pos,vel)
    if pos(1)<0
        vel(1)=-vel(1)*0.5;
    end
    if pos(1)>GRID_W
        vel(1)=-vel(1)*0.5;
    end
    if pos(2)<0
        vel(2)=-vel(2)*0.5;
    end
    if pos(2)>GRID_H
        vel(2)=-vel(2)*0.5;
    end
end

function valid = is_valid_grid(x,y)
    valid = ~is_out_of_bound([x,y]);
    if valid
        valid = valid && grid_type(y,x);
    end
end

function particle_vx_to_grid(vx,pos)
    x=pos(1);
    y=pos(2);
    dx=x-floor(x);
    dy=y-floor(y);
    x=ceil(x);
    y=ceil(y);
    i1=zeros(2);
    i2=zeros(2);
    i3=zeros(2);
    i4=zeros(2);

    if dx<0.5
        dx=dx+0.5;
        i1=[x-1,y+1];
        i2=[x,y+1];
        i3=[x,y];
        i4=[x-1,y];
    else
        dx=dx-0.5;
        i1=[x,y+1];
        i2=[x+1,y+1];
        i3=[x+1,y];
        i4=[x,y];
    end

    w1=(1-dx)*(1-dy);
    w2=dx*(1-dy);
    w3=dx*dy;
    w4=(1-dx)*dy;

    s1=0;
    s2=0;
    s3=0;
    s4=0;

    if is_valid_grid(i1(1),i1(2))
        s1=1;
    end
    if is_valid_grid(i2(1),i2(2))
        s2=1;
    end
    if is_valid_grid(i3(1),i3(2))
        s3=1;
    end
    if is_valid_grid(i4(1),i4(2))
        s4=1;
    end

    w=w1*s1+w2*s2+w3*s3+w4*s4;

    if is_valid_grid(i1(1),i1(2))
        grid_v_x(i1(2),i1(1))=grid_v_x(i1(2),i1(1))+vx*w1/w;
        grid_v_x_w(i1(2),i1(1))=grid_v_x_w(i1(2),i1(1))+w1/w;
    end
    if is_valid_grid(i2(1),i2(2))
        grid_v_x(i2(2),i2(1))=grid_v_x(i2(2),i2(1))+vx*w2/w;
        grid_v_x_w(i2(2),i2(1))=grid_v_x_w(i2(2),i2(1))+w2/w;
    end
    if is_valid_grid(i3(1),i3(2))
        grid_v_x(i3(2),i3(1))=grid_v_x(i3(2),i3(1))+vx*w3/w;
        grid_v_x_w(i3(2),i3(1))=grid_v_x_w(i3(2),i3(1))+w3/w;
    end
    if is_valid_grid(i4(1),i4(2))
        grid_v_x(i4(2),i4(1))=grid_v_x(i4(2),i4(1))+vx*w4/w;
        grid_v_x_w(i4(2),i4(1))=grid_v_x_w(i4(2),i4(1))+w4/w;
    end
end

function particle_vy_to_grid(vy,pos)
    x=pos(1);
    y=pos(2);
    dx=x-floor(x);
    dy=y-floor(y);
    x=ceil(x);
    y=ceil(y);
    i1=zeros(2);
    i2=zeros(2);
    i3=zeros(2);
    i4=zeros(2);

    if dy<0.5
        dy=dy+0.5;
        i1=[x,y];
        i2=[x+1,y];
        i3=[x+1,y-1];
        i4=[x,y-1];
    else
        dy=dy-0.5;
        i1=[x,y+1];
        i2=[x+1,y+1];
        i3=[x+1,y];
        i4=[x,y];
    end

    w1=(1-dx)*(1-dy);
    w2=dx*(1-dy);
    w3=dx*dy;
    w4=(1-dx)*dy;

    s1=0;
    s2=0;
    s3=0;
    s4=0;

    if is_valid_grid(i1(1),i1(2))
        s1=1;
    end
    if is_valid_grid(i2(1),i2(2))
        s2=1;
    end
    if is_valid_grid(i3(1),i3(2))
        s3=1;
    end
    if is_valid_grid(i4(1),i4(2))
        s4=1;
    end

    w=w1*s1+w2*s2+w3*s3+w4*s4;

    if is_valid_grid(i1(1),i1(2))
        grid_v_y(i1(2),i1(1))=grid_v_y(i1(2),i1(1))+vy*w1/w;
        grid_v_y_w(i1(2),i1(1))=grid_v_y_w(i1(2),i1(1))+w1/w;
    end
    if is_valid_grid(i2(1),i2(2))
        grid_v_y(i2(2),i2(1))=grid_v_y(i2(2),i2(1))+vy*w2/w;
        grid_v_y_w(i2(2),i2(1))=grid_v_y_w(i2(2),i2(1))+w2/w;
    end
    if is_valid_grid(i3(1),i3(2))
        grid_v_y(i3(2),i3(1))=grid_v_y(i3(2),i3(1))+vy*w3/w;
        grid_v_y_w(i3(2),i3(1))=grid_v_y_w(i3(2),i3(1))+w3/w;
    end
    if is_valid_grid(i4(1),i4(2))
        grid_v_y(i4(2),i4(1))=grid_v_y(i4(2),i4(1))+vy*w4/w;
        grid_v_y_w(i4(2),i4(1))=grid_v_y_w(i4(2),i4(1))+w4/w;
    end
end

function grid_vx_to_particle(pos,index)
    x=pos(1);
    y=pos(2);
    dx=x-floor(x);
    dy=y-floor(y);
    x=ceil(x);
    y=ceil(y);
    i1=zeros(2);
    i2=zeros(2);
    i3=zeros(2);
    i4=zeros(2);

    if dx<0.5
        dx=dx+0.5;
        i1=[x-1,y+1];
        i2=[x,y+1];
        i3=[x,y];
        i4=[x-1,y];
    else
        dx=dx-0.5;
        i1=[x,y+1];
        i2=[x+1,y+1];
        i3=[x+1,y];
        i4=[x,y];
    end

    w1=(1-dx)*(1-dy);
    w2=dx*(1-dy);
    w3=dx*dy;
    w4=(1-dx)*dy;

    s1=0;
    s2=0;
    s3=0;
    s4=0;

    if is_valid_grid(i1(1),i1(2))
        s1=1;
    end
    if is_valid_grid(i2(1),i2(2))
        s2=1;
    end
    if is_valid_grid(i3(1),i3(2))
        s3=1;
    end
    if is_valid_grid(i4(1),i4(2))
        s4=1;
    end

    w=w1*s1+w2*s2+w3*s3+w4*s4;

    if is_valid_grid(i1(1),i1(2))
        particle_v(index,1)=particle_v(index,1)+grid_v_x(i1(2),i1(1))*w1/w/grid_v_x_w(i1(2),i1(1));
    end
    if is_valid_grid(i2(1),i2(2))
        particle_v(index,1)=particle_v(index,1)+grid_v_x(i2(2),i2(1))*w2/w/grid_v_x_w(i2(2),i2(1));
    end
    if is_valid_grid(i3(1),i3(2))
        particle_v(index,1)=particle_v(index,1)+grid_v_x(i3(2),i3(1))*w3/w/grid_v_x_w(i3(2),i3(1));
    end
    if is_valid_grid(i4(1),i4(2))
        particle_v(index,1)=particle_v(index,1)+grid_v_x(i4(2),i4(1))*w4/w/grid_v_x_w(i4(2),i4(1));
    end
end

function grid_vy_to_particle(pos,index)
    x=pos(1);
    y=pos(2);
    dx=x-floor(x);
    dy=y-floor(y);
    x=ceil(x);
    y=ceil(y);
    i1=zeros(2);
    i2=zeros(2);
    i3=zeros(2);
    i4=zeros(2);

    if dy<0.5
        dy=dy+0.5;
        i1=[x,y];
        i2=[x+1,y];
        i3=[x+1,y-1];
        i4=[x,y-1];
    else
        dy=dy-0.5;
        i1=[x,y+1];
        i2=[x+1,y+1];
        i3=[x+1,y];
        i4=[x,y];
    end

    w1=(1-dx)*(1-dy);
    w2=dx*(1-dy);
    w3=dx*dy;
    w4=(1-dx)*dy;

    s1=0;
    s2=0;
    s3=0;
    s4=0;

    if is_valid_grid(i1(1),i1(2))
        s1=1;
    end
    if is_valid_grid(i2(1),i2(2))
        s2=1;
    end
    if is_valid_grid(i3(1),i3(2))
        s3=1;
    end
    if is_valid_grid(i4(1),i4(2))
        s4=1;
    end

    w=w1*s1+w2*s2+w3*s3+w4*s4;

    if is_valid_grid(i1(1),i1(2))
        particle_v(index,2)=particle_v(index,2)+grid_v_y(i1(2),i1(1))*w1/w/grid_v_y_w(i1(2),i1(1));
    end
    if is_valid_grid(i2(1),i2(2))
        particle_v(index,2)=particle_v(index,2)+grid_v_y(i2(2),i2(1))*w2/w/grid_v_y_w(i2(2),i2(1));
    end
    if is_valid_grid(i3(1),i3(2))
        particle_v(index,2)=particle_v(index,2)+grid_v_y(i3(2),i3(1))*w3/w/grid_v_y_w(i3(2),i3(1));
    end
    if is_valid_grid(i4(1),i4(2))
        particle_v(index,2)=particle_v(index,2)+grid_v_y(i4(2),i4(1))*w4/w/grid_v_y_w(i4(2),i4(1));
    end
end

function solve_incompressible()
    for ite=1:NUM_ITER
    end
end

clear global;
close all;

for i=1:TIME_STEP_TOTAL
    % add gravity to particle
    addGravity();

    % clear grid
    grid=zeros(GRID_H,GRID_W);
    grid_v_x=zeros(GRID_H,GRID_W+1);
    grid_v_y=zeros(GRID_H+1,GRID_W);
    grid_v_x_w=zeros(GRID_H,GRID_W+1);
    grid_v_y_w=zeros(GRID_H+1,GRID_W);

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
        particle_vx_to_grid(particle_v(j,1),particles(j,:));
        particle_vy_to_grid(particle_v(j,2),particles(j,:));
    end
    % clear particle velocity after interpolation
    particle_v_copy=particle_v;
    particle_v=zeros(PARTICLE_NUM,2);
    
    % make grid incompressible

    % add particle velocity back
    for j=1:PARTICLE_NUM
        grid_vx_to_particle(particles(j,:),j);
        grid_vy_to_particle(particles(j,:),j);
    end

    particle_v=CO*particle_v_copy+(1-CO)*particle_v;

    % update position
    particles=particles+particle_v*DT;
    
    imshow(flipud(grid));
    %h.GridVisible='off';
    %disp(max(grid_v_y,[],'all'))
    drawnow

end

end
