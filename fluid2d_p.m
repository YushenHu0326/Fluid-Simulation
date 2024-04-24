function fluid2d_p()

GRAVITY=-10;
DT=0.1;
O=1.9;
K=1;
CO=0.5;
NUM_ITER=10;
TIME_STEP_TOTAL=1000;
PARTICLE_NUM=500;
GRID_H=10;
GRID_W=10;
PARTICLE_PER_GRID=10;
PARTICLE_R=4;

GRID_HASH_SIZE=100;

img=zeros(GRID_H*5,GRID_W*5);

grid=zeros(GRID_H,GRID_W);
grid_type=ones(GRID_H,GRID_W);
grid_v_x=zeros(GRID_H,GRID_W+1);
grid_v_y=zeros(GRID_H+1,GRID_W);
grid_v_x_w=zeros(GRID_H,GRID_W+1);
grid_v_y_w=zeros(GRID_H+1,GRID_W);
grid_hash=zeros(GRID_H,GRID_W,GRID_HASH_SIZE);
unhashed_particles=zeros(PARTICLE_NUM+1);

particles=zeros(PARTICLE_NUM,2);
for i=1:PARTICLE_NUM
    particles(i,1)=rand*GRID_W/2;
    particles(i,2)=rand*GRID_H/2;
end
particle_v=zeros(PARTICLE_NUM,2);
particle_v_copy=zeros(PARTICLE_NUM,2);

function add_gravity()
    particle_v(:,2)=particle_v(:,2)+GRAVITY*DT;
end

function outbound = is_out_of_bound(pos)
    outbound=(pos(1)<=0||pos(1)>=GRID_W)||(pos(2)<=0||pos(2)>=GRID_H);
end

function solve_bound_collision()
    for i=1:PARTICLE_NUM
        x=particles(i,1);
        y=particles(i,2);
        if x<0
            particles(i,1)=0;
            particle_v(i,1)=-particle_v(i,1)/10;
        end
        if x>GRID_W
            particles(i,1)=GRID_W;
            particle_v(i,1)=-particle_v(i,1)/10;
        end
        if y<0
            particles(i,2)=0;
            particle_v(i,2)=-particle_v(i,2)/10;
        end
        if y>GRID_H
            particles(i,2)=GRID_H;
            particle_v(i,2)=-particle_v(i,2)/10;
        end
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

    w01=1;
    w02=1;
    w03=1;  
    w04=1;

    if is_valid_grid(i1(1),i1(2))
        if grid_v_x_w(i1(2),i1(1))~=0
            w01=grid_v_x_w(i1(2),i1(1));
        end
        particle_v(index,1)=particle_v(index,1)+grid_v_x(i1(2),i1(1))*w1/w/w01;
    end
    if is_valid_grid(i2(1),i2(2))
        if grid_v_x_w(i2(2),i2(1))~=0
            w02=grid_v_x_w(i2(2),i2(1));
        end
        particle_v(index,1)=particle_v(index,1)+grid_v_x(i2(2),i2(1))*w2/w/w02;
    end
    if is_valid_grid(i3(1),i3(2))
        if grid_v_x_w(i3(2),i3(1))~=0
            w03=grid_v_x_w(i3(2),i3(1));
        end
        particle_v(index,1)=particle_v(index,1)+grid_v_x(i3(2),i3(1))*w3/w/w03;
    end
    if is_valid_grid(i4(1),i4(2))
        if grid_v_x_w(i4(2),i4(1))~=0
            w04=grid_v_x_w(i4(2),i4(1));
        end
        particle_v(index,1)=particle_v(index,1)+grid_v_x(i4(2),i4(1))*w4/w/w04;
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

    w01=1;
    w02=1;
    w03=1;
    w04=1;

    if is_valid_grid(i1(1),i1(2))
        if grid_v_y_w(i1(2),i1(1))~=0
            w01=grid_v_y_w(i1(2),i1(1));
        end
        particle_v(index,2)=particle_v(index,2)+grid_v_y(i1(2),i1(1))*w1/w/w01;
    end
    if is_valid_grid(i2(1),i2(2))
        if grid_v_y_w(i2(2),i2(1))~=0
            w02=grid_v_y_w(i2(2),i2(1));
        end
        particle_v(index,2)=particle_v(index,2)+grid_v_y(i2(2),i2(1))*w2/w/w02;
    end
    if is_valid_grid(i3(1),i3(2))
        if grid_v_y_w(i3(2),i3(1))~=0
            w03=grid_v_y_w(i3(2),i3(1));
        end
        particle_v(index,2)=particle_v(index,2)+grid_v_y(i3(2),i3(1))*w3/w/w03;
    end
    if is_valid_grid(i4(1),i4(2))
        if grid_v_y_w(i4(2),i4(1))~=0
            w04=grid_v_y_w(i4(2),i4(1));
        end
        particle_v(index,2)=particle_v(index,2)+grid_v_y(i4(2),i4(1))*w4/w/w04;
    end
end

function solve_incompressible()
    for ite=1:NUM_ITER
        for y=1:GRID_H
            for x=1:GRID_W
                if grid(y,x)>PARTICLE_PER_GRID
                    d=-grid_v_x(y,x);
                    d=d+grid_v_x(y,x+1);
                    d=d-grid_v_y(y,x);
                    d=d+grid_v_y(y+1,x);
                    d=d*O;
                    d=d-K*(grid(y,x)-PARTICLE_PER_GRID);
                    s=0;
                    if is_valid_grid(x-1,y)
                        s=s+1;
                    end
                    if is_valid_grid(x+1,y)
                        s=s+1;
                    end
                    if is_valid_grid(x,y-1)
                        s=s+1;
                    end
                    if is_valid_grid(x,y+1)
                        s=s+1;
                    end

                    %disp(d)

                    if is_valid_grid(x-1,y)
                        grid_v_x(y,x)=grid_v_x(y,x)+d/s;
                    end
                    if is_valid_grid(x+1,y)
                        grid_v_x(y,x+1)=grid_v_x(y,x+1)-d/s;
                    end
                    if is_valid_grid(x,y-1)
                        grid_v_y(y,x)=grid_v_y(y,x)+d/s;
                    end
                    if is_valid_grid(x,y+1)
                        grid_v_y(y+1,x)=grid_v_y(y+1,x)-d/s;
                    end
                end
            end
        end
    end
end

function hash_particles()
    grid_hash=zeros(GRID_H,GRID_W,GRID_HASH_SIZE);
    unhashed_particles=zeros(PARTICLE_NUM+1);
    for i=1:PARTICLE_NUM
        x=ceil(particles(i,1));
        y=ceil(particles(i,2));
        if is_out_of_bound([x,y])
            continue;
        end
        if grid_hash(y,x,1)<GRID_HASH_SIZE-1
            grid_hash(y,x,1)=grid_hash(y,x,1)+1;
            grid_hash(y,x,grid_hash(y,x,1)+1)=i;
        else
            unhashed_particles(1)=unhashed_particles(1)+1;
            unhashed_particles(unhashed_particles(1)+1)=i;
        end
    end
end

function particle_collision()
    
        for i=1:GRID_H
            for j=1:GRID_H
                k=grid_hash(i,j,1);
                if k>1
                    for m=2:k
                        if m<k
                            for n=m+1:k
                                d=particles(m,:)-particles(n,:);
                                %disp(norm(d))
                                if norm(d)<2*PARTICLE_R
                                    if norm(d)==0
                                        d=[1,1];
                                    end
                                    dp=d/norm(d)*(PARTICLE_R-norm(d)/2)/50;
                                    %if isnan(dp(1)) || isnan(dp(2))
                                        %disp(normalize(d))
                                    %end
                                    particle_v(m,1)=particle_v(m,1)+dp(1);
                                    particle_v(m,2)=particle_v(m,2)+dp(2);
                                    particle_v(n,1)=particle_v(n,1)-dp(1);
                                    particle_v(n,2)=particle_v(n,2)-dp(2);
                                end
                            end
                        end
                    end
                end
            end
        end

end

clear global;
close all;

for t=1:TIME_STEP_TOTAL
    % clear grid
    grid=zeros(GRID_H,GRID_W);
    grid_v_x=zeros(GRID_H,GRID_W+1);
    grid_v_y=zeros(GRID_H+1,GRID_W);
    grid_v_x_w=zeros(GRID_H,GRID_W+1);
    grid_v_y_w=zeros(GRID_H+1,GRID_W);
    
    % add gravity to particle
    add_gravity();
    % update particles in grid
    solve_bound_collision();
    % calculate particle collision
    hash_particles();
    particle_collision();

    for i=1:PARTICLE_NUM
        % push particles to grid
        if is_out_of_bound(particles(i,:))
            continue;
        end
        grid(ceil(particles(i,2)),ceil(particles(i,1)))=grid(ceil(particles(i,2)),ceil(particles(i,1)))+1;
        % push particle velocity to corresponding grid velocity
        particle_vx_to_grid(particle_v(i,1),particles(i,:));
        particle_vy_to_grid(particle_v(i,2),particles(i,:));
    end

    % clear particle velocity after interpolation
    %particle_v_copy=particle_v;
    %particle_v=zeros(PARTICLE_NUM,2);
    
    % make grid incompressible
    solve_incompressible();

    % add particle velocity back
    for i=1:PARTICLE_NUM
        %grid_vx_to_particle(particles(i,:),i);
        %grid_vy_to_particle(particles(i,:),i);
    end

    %particle_v=CO*particle_v_copy+(1-CO)*particle_v;

    % update position
    particles=particles+particle_v*DT;

    img=zeros(GRID_H*5,GRID_W*5);
    for i=1:PARTICLE_NUM
        if ~is_out_of_bound(particles(i,:))
            img(ceil(particles(i,2)*5),ceil(particles(i,1)*5))=img(ceil(particles(i,2)*5),ceil(particles(i,1)*5))+1;
        end
    end
    
    imshow(flipud(img));
    %imshow(flipud(grid));
    %h.GridVisible='off';
    %disp(min(particle_v(:),[],'all'))
    %set(gca, 'XLim', [0,10], 'YLim', [0,10])
    %plot(particles(:,1),particles(:,2),'ro', 'MarkerSize', 3)
    %set(gca, 'XLim', [0,10], 'YLim', [0,10])
    drawnow

end

end