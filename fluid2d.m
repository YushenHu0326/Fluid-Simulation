function fluid2d()

GRAVITY=1;
DT=0.1;
NUM_ITER=50;
TIME_STEP_TOTAL=1000;
PARTICLE_NUM=10000;
GRID_H=100;
GRID_W=100;
PARTICLE_PER_GRID=10;

grid=zeros(GRID_W,GRID_H);
grid_type=zeros(GRID_W,GRID_H);
grid_v=zeros(GRID_W,GRID_H,4);
grid_v_copy=zeros(GRID_W,GRID_H,4);

particles=zeros(PARTICLE_NUM,2);
for i=1:PARTICLE_NUM
    particles(i,1)=rand*GRID_W;
    particles(i,2)=rand*GRID_H;
end
particle_v=zeros(PARTICLE_NUM,2);

function addGravity()
    particle_v(:,2)=particle_v(:,2)+GRAVITY*DT;
end

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

function particle_vel_to_grid(vel,pos)
    x=pos(1);
    y=pos(2);
    dx=x-floor(x);
    dy=y-floor(y);
    x=ceil(x);
    y=ceil(y);
    w1=(1-dx)*(1-dy);
    w2=dx*(1-dy);
    w3=dx*dy;
    w4=(1-dx)*dy;
    if dx<0.5
        if x>1
            grid_v(x-1,y,3)=grid_v(x-1,y,3)+vel(2)*w4/(w1+w2+w3+w4);
            grid_v(x,y,3)=grid_v(x,y,3)+vel(2)*w3/(w1+w2+w3+w4);
            grid_v(x-1,y,4)=grid_v(x-1,y,4)+vel(2)*w1/(w1+w2+w3+w4);
            grid_v(x,y,4)=grid_v(x,y,4)+vel(2)*w2/(w1+w2+w3+w4);
        else
            grid_v(x,y,3)=grid_v(x,y,3)+vel(2)*w3/(w2+w3);
            grid_v(x,y,4)=grid_v(x,y,4)+vel(2)*w2/(w2+w3);
        end
    else
        if x<GRID_W
            grid_v(x,y,3)=grid_v(x,y,3)+vel(2)*w4/(w1+w2+w3+w4);
            grid_v(x+1,y,3)=grid_v(x+1,y,3)+vel(2)*w3/(w1+w2+w3+w4);
            grid_v(x,y,4)=grid_v(x,y,4)+vel(2)*w1/(w1+w2+w3+w4);
            grid_v(x+1,y,4)=grid_v(x+1,y,4)+vel(2)*w2/(w1+w2+w3+w4);
        else
            grid_v(x,y,3)=grid_v(x,y,3)+vel(2)*w4/(w1+w4);
            grid_v(x,y,4)=grid_v(x,y,4)+vel(2)*w1/(w1+w4);
        end
    end

    if dy<0.5
        if y>1
            grid_v(x,y-1,1)=grid_v(x,y-1,1)+vel(1)*w4/(w1+w2+w3+w4);
            grid_v(x,y,1)=grid_v(x,y,1)+vel(1)*w1/(w1+w2+w3+w4);
            grid_v(x,y-1,2)=grid_v(x,y-1,2)+vel(1)*w3/(w1+w2+w3+w4);
            grid_v(x,y,2)=grid_v(x,y,2)+vel(1)*w2/(w1+w2+w3+w4);
        else
            grid_v(x,y,1)=grid_v(x,y,1)+vel(1)*w1/(w1+w2);
            grid_v(x,y,2)=grid_v(x,y,2)+vel(1)*w2/(w1+w2);
        end
    else
        if y<GRID_H
            grid_v(x,y,1)=grid_v(x,y,1)+vel(1)*w4/(w1+w2+w3+w4);
            grid_v(x,y+1,1)=grid_v(x,y+1,1)+vel(1)*w1/(w1+w2+w3+w4);
            grid_v(x,y,2)=grid_v(x,y,2)+vel(1)*w3/(w1+w2+w3+w4);
            grid_v(x,y+1,2)=grid_v(x,y+1,2)+vel(1)*w2/(w1+w2+w3+w4);
        else
            grid_v(x,y,1)=grid_v(x,y,1)+vel(1)*w4/(w1+w4);
            grid_v(x,y,2)=grid_v(x,y,2)+vel(1)*w3/(w1+w4);
        end
    end
end

function grid_vel_to_particle(index,pos)
    x=pos(1);
    y=pos(2);
    dx=x-floor(x);
    dy=y-floor(y);
    x=ceil(x);
    y=ceil(y);
    w1=(1-dx)*(1-dy);
    w2=dx*(1-dy);
    w3=dx*dy;
    w4=(1-dx)*dy;
    if dx<0.5
        if x>1
            particle_v(index,2)=grid_v(x-1,y,3)*w4;
            particle_v(index,2)=particle_v(index,2)+grid_v(x,y,3)*w3;
            particle_v(index,2)=particle_v(index,2)+grid_v(x-1,y,4)*w1;
            particle_v(index,2)=particle_v(index,2)+grid_v(x,y,4)*w2;
            particle_v(index,2)=particle_v(index,2)/(w1+w2+w3+w4);
        else
            particle_v(index,2)=grid_v(x,y,3)*w3;
            particle_v(index,2)=particle_v(index,2)+grid_v(x,y,4)*w2;
            particle_v(index,2)=particle_v(index,2)/(w2+w3);
        end
    else
        if x<GRID_W
            particle_v(index,2)=grid_v(x,y,3)*w4;
            particle_v(index,2)=particle_v(index,2)+grid_v(x+1,y,3)*w3;
            particle_v(index,2)=particle_v(index,2)+grid_v(x,y,4)*w1;
            particle_v(index,2)=particle_v(index,2)+grid_v(x+1,y,4)*w2;
            particle_v(index,2)=particle_v(index,2)/(w1+w2+w3+w4);
        else
            particle_v(index,2)=grid_v(x,y,3)*w4;
            particle_v(index,2)=particle_v(index,2)+grid_v(x,y,4)*w1;
            particle_v(index,2)=particle_v(index,2)/(w1+w4);
        end
    end

    if dy<0.5
        if y>1
            particle_v(index,1)=grid_v(x,y-1,1)*w4;
            particle_v(index,1)=particle_v(index,1)+grid_v(x,y,1)*w1;
            particle_v(index,1)=particle_v(index,1)+grid_v(x,y-1,2)*w3;
            particle_v(index,1)=particle_v(index,1)+grid_v(x,y,2)*w2;
            particle_v(index,1)=particle_v(index,1)/(w1+w2+w3+w4);
        else
            particle_v(index,1)=grid_v(x,y,1)*w1;
            particle_v(index,1)=particle_v(index,1)+grid_v(x,y,2)*w2;
            particle_v(index,1)=particle_v(index,1)/(w1+w2);
        end
    else
        if y<GRID_H
            particle_v(index,1)=grid_v(x,y,1)*w4;
            particle_v(index,1)=particle_v(index,1)+grid_v(x,y+1,1)*w1;
            particle_v(index,1)=particle_v(index,1)+grid_v(x,y,2)*w3;
            particle_v(index,1)=particle_v(index,1)+grid_v(x,y+1,2)*w2;
            particle_v(index,1)=particle_v(index,1)/(w1+w2+w3+w4);
        else
            particle_v(index,1)=grid_v(x,y,1)*w4;
            particle_v(index,1)=particle_v(index,1)+grid_v(x,y,2)*w3;
            particle_v(index,1)=particle_v(index,1)/(w3+w4);
        end
    end
end

clear global;
close all;

for i=1:TIME_STEP_TOTAL
    % add gravity to particle
    addGravity();

    % clear grid
    grid=zeros(GRID_W,GRID_H);
    grid_v=zeros(GRID_W,GRID_H,4);

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
        particle_vel_to_grid(particle_v(j,:),particles(j,:));
    end
    % make grid incompressible
    for n=1:NUM_ITER
        for j=1:GRID_H
            for k=1:GRID_W
                if grid(j,k)>PARTICLE_PER_GRID
                    d=0;
                    cnt=0;
                    if j~=1
                        d=d-grid_v(j-1,k,3);
                        cnt=cnt+1;
                    end
                    if j~=GRID_H
                        d=d+grid_v(j+1,k,4);
                        cnt=cnt+1;
                    end
                    if k~=1
                        d=d-grid_v(j,k-1,1);
                        cnt=cnt+1;
                    end
                    if k~=GRID_W
                        d=d+grid_v(j,k+1,2);
                        cnt=cnt+1;
                    end
                    d=d/cnt;
                    %disp(d)
                    if j~=1
                        grid_v(j-1,k,3)=grid_v(j-1,k,3)+d;
                    end
                    if j~=GRID_H
                        grid_v(j+1,k,4)=grid_v(j+1,k,4)-d;
                    end
                    if k~=1
                        grid_v(j,k-1,1)=grid_v(j,k-1,1)+d;
                    end
                    if k~=GRID_W
                        grid_v(j,k+1,2)=grid_v(j,k+1,2)-d;
                    end
                end
            end
        end
    end

    for j=1:PARTICLE_NUM
        % update particle velocity based on corresponding grid velocity
        grid_vel_to_particle(j,particles(j,:));
    end

    % update position
    particles=particles+particle_v*DT;
    
    imshow(grid/20);
    %h.GridVisible='off';
    %disp(max(grid_v_copy,[],'all'))
    drawnow

end

end
