function fluid2d_p()

GRAVITY=-10;
DT=0.05;
TIME_STEP_TOTAL=1000;
PARTICLE_NUM=500;
GRID_H=20;
GRID_W=20;
GRID_HASH_SIZE=100;

% particle properties
M=0.1;
M2=M*M;
R=1;
R2=R*R;
R3=R*R2;
R4=R*R3;
R5=R*R4;
GAS_CONSTANT=2;
RESTDENSITY=2;
VISCOSITY=0;

img=zeros(GRID_H*5,GRID_W*5,3);

grid_hash=zeros(GRID_H,GRID_W,GRID_HASH_SIZE);
unhashed_particles=zeros(PARTICLE_NUM+1);

particles=zeros(PARTICLE_NUM,2);
for i=1:PARTICLE_NUM
    particles(i,1)=rand*GRID_W/2+GRID_W/4;
    particles(i,2)=rand*GRID_H/2+GRID_W/4;
end
particle_v=zeros(PARTICLE_NUM,2);
particle_a=zeros(PARTICLE_NUM,2);
particle_d=zeros(PARTICLE_NUM,2);
particle_p=zeros(PARTICLE_NUM,2);

function outbound = is_out_of_bound(pos)
    outbound=(pos(1)<=0||pos(1)>=GRID_W)||(pos(2)<=0||pos(2)>=GRID_H);
end

function solve_bound_collision()
    for i=1:PARTICLE_NUM
        x=particles(i,1);
        y=particles(i,2);
        if x<0
            particles(i,1)=0.5;
            particle_v(i,1)=-particle_v(i,1)/5;
        end
        if x>GRID_W
            particles(i,1)=GRID_W-0.5;
            particle_v(i,1)=-particle_v(i,1)/5;
        end
        if y<0
            particles(i,2)=0.5;
            particle_v(i,2)=-particle_v(i,2)/5;
        end
        if y>GRID_H
            particles(i,2)=GRID_H-0.5;
            particle_v(i,2)=-particle_v(i,2)/5;
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
            if grid_hash(y,x,1)==0
                grid_hash(y,x,1)=2;
            else
                grid_hash(y,x,1)=grid_hash(y,x,1)+1;
            end
            grid_hash(y,x,grid_hash(y,x,1))=i;
        else
            unhashed_particles(1)=unhashed_particles(1)+1;
            unhashed_particles(unhashed_particles(1)+1)=i;
        end
    end
end

function k = poly6kernel(distance_sqr)
    x=1-distance_sqr/R2;
    k=315/(64*pi*R3)*x*x*x;
end

function k = spikykerneld1(distance)
    x=1-distance/R;
    k=-45/(pi*R4)*x*x;
end

function k = spikykerneld2(distance)
    x=1-distance/R;
    k=90/(pi*R5)*x;
end

function k = spikykernelg(distance,direction)
    k=spikykerneld1(distance)*direction;
end

function compute_density()
    for i=1:PARTICLE_NUM
        x=ceil(particles(i,1));
        y=ceil(particles(i,2));

        sum=0;
        index=grid_hash(y,x,1);
        if index>1
            ii=2;
            while ii<=index
                j=grid_hash(y,x,ii);
                d=norm(particles(i,:)-particles(j,:));
                if d<R
                    sum=sum+poly6kernel(d*d*0.004);
                end
            ii=ii+1;
            end
        end

        particle_d(i)=sum*M+0.0001;
        particle_p(i)=GAS_CONSTANT*(particle_d(i)-RESTDENSITY);
    end
end

function compute_acc()
    for i=1:PARTICLE_NUM
        x=ceil(particles(i,1));
        y=ceil(particles(i,2));

        pressure=zeros(2);
        viscosity=zeros(2);

        D2=particle_d(i)*particle_d(i);

        index=grid_hash(y,x,1);
        if index>1
            ii=2;
            while ii<=index
                j=grid_hash(y,x,ii);
                d=norm(particles(i,:)-particles(j,:));
                if d~=0
                    if d<R && particle_d(j)>0.0001
                        pd=(particles(i,:)-particles(j,:))/d;
                        pc=M2*spikykernelg(d,pd);
                        pc=pc*(particle_p(i)/D2+particle_p(j)/(particle_d(j)*particle_d(j)));

                        vc=VISCOSITY*M2*(particle_v(j,:)-particle_v(i,:))/particle_d(j);
                        vc=vc*spikykerneld2(d);

                        pressure=pressure+pc;
                        viscosity=viscosity+vc;
                    end
                end
                ii=ii+1;
            end
        end

        particle_a(i,1)=-pressure(1)/M+viscosity(1)/M;
        particle_a(i,2)=GRAVITY-pressure(2)/M+viscosity(2)/M;
        %disp(pressure(2))
    end
end

clear global;
close all;

for t=1:TIME_STEP_TOTAL
    % clear grid
    grid=zeros(GRID_H,GRID_W);
    
    % update particles in grid
    solve_bound_collision();

    hash_particles();

    compute_density();

    compute_acc();

    particle_v=particle_v+particle_a*DT;

    % update position
    particles=particles+particle_v*DT;

    img=zeros(GRID_H*5,GRID_W*5,3);
    for i=1:PARTICLE_NUM
        if ~is_out_of_bound(particles(i,:))
            %img(ceil(particles(i,2)*5),ceil(particles(i,1)*5))=img(ceil(particles(i,2)*5),ceil(particles(i,1)*5))+1;
            img(ceil(particles(i,2)*5),ceil(particles(i,1)*5),1)=1;
            img(ceil(particles(i,2)*5),ceil(particles(i,1)*5),2)=img(ceil(particles(i,2)*5),ceil(particles(i,1)*5),2)+norm(particle_p(i))/20;
        end
    end
    
    imshow(flipud(img));
    %imshow(flipud(grid));
    %h.GridVisible='off';
    %disp(max(particle_a(:,2),[],'all'))
    %set(gca, 'XLim', [0,10], 'YLim', [0,10])
    %plot(particles(:,1),particles(:,2),'ro', 'MarkerSize', 3)
    %set(gca, 'XLim', [0,10], 'YLim', [0,10])
    drawnow

end

end
