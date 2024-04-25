function fluid2d_p()

GRAVITY=-10;
DT=0.01;
TIME_STEP_TOTAL=5000;
MAX_PARTICLE_NUM=1000;
PARTICLE_NUM=500;
GRID_H=1;
GRID_W=1;

% particle properties
M=0.1;
DAMP=1;
SMOOTHING_LEN=0.02;
DENS_TO_P=0.05;

img=zeros(500,500,3);

generate_method=1;
present_heatmap=0;
has_obstacle=0;

particle_v=zeros(PARTICLE_NUM,2);

if generate_method==0
    particles=zeros(PARTICLE_NUM,2);
    px=SMOOTHING_LEN;
    py=SMOOTHING_LEN;
    for i=1:PARTICLE_NUM
        particles(i,1)=px+SMOOTHING_LEN;
        particles(i,2)=py;
        px=px+SMOOTHING_LEN;
        if px>0.2
            py=py+SMOOTHING_LEN;
            px=SMOOTHING_LEN;
        end
    end
else
    PARTICLE_NUM=1;
    particles(1,1)=0;
    particles(1,2)=GRID_H/2;
    particle_v=zeros(PARTICLE_NUM,2);
    particle_v(1,1)=2+rand;
end

particle_a=zeros(PARTICLE_NUM,2);
particle_d=zeros(PARTICLE_NUM,2);
particle_p=zeros(PARTICLE_NUM,2);

obstacle=[0.5,0.2,0.3];

function outbound = is_out_of_bound(pos)
    outbound=(pos(1)<=0||pos(1)>=GRID_W)||(pos(2)<=0||pos(2)>=GRID_H);
end

function solve_bound_collision(obstacle_on)
    for i=1:PARTICLE_NUM
        x=particles(i,1);
        y=particles(i,2);
        if x<0
            particles(i,1)=0;
            particle_a(i,1)=100;
        end
        if x>GRID_W
            particles(i,1)=GRID_W;
            particle_a(i,1)=-100;
        end
        if y<0
            particles(i,2)=0;
            particle_a(i,2)=100;
        end
        if y>GRID_H
            particles(i,2)=GRID_H;
            particle_a(i,2)=-100;
        end

        if obstacle_on
            if norm(particles(i,:)-[obstacle(1),obstacle(2)])<obstacle(3)
                d=particles(i,:)-[obstacle(1),obstacle(2)];
                particles(i,:)=d/norm(d)*obstacle(3)+[obstacle(1),obstacle(2)];
                particle_a(i,:)=particle_a(i,:)+d*100;
            end
        end
    end
end

function k = kernel(r,h)
    k=1/(h^2*pi)*exp(norm(r)^2/h^2);	
end

function k = kernelg(r,h)
	n=-2*exp(-norm(r)^2/h^2)/(h^4*pi) ;
	k=n.*r;	
end

function compute_density()
    for i=1:PARTICLE_NUM
        particle_d(i)=M*kernel(0,SMOOTHING_LEN);
        for j=i+1:PARTICLE_NUM
            d=particles(i,:)-particles(j,:);
            if norm(d)<SMOOTHING_LEN
                dens=M*kernel(d,SMOOTHING_LEN);
                particle_d(i)=particle_d(i)+dens;
                %particle_d(j)=particle_d(j)+dens;
            end
        end

        particle_p(i)=DENS_TO_P*particle_d(i);
    end
end

function compute_acc()
    for i=1:PARTICLE_NUM
        for j=i+1:PARTICLE_NUM
            d=particles(i,:)-particles(j,:);
            if norm(d)<SMOOTHING_LEN
                a=-M*(particle_p(i)/particle_d(i)^2+particle_p(j)/particle_d(j)^2)*kernelg(d,SMOOTHING_LEN);
                particle_a(i,1)=particle_a(i,1)+a(1);
                particle_a(i,2)=particle_a(i,2)+a(1);
                particle_a(j,1)=particle_a(j,1)-a(2);
                particle_a(j,2)=particle_a(j,2)-a(2);
            end
        end
        particle_a(i,1)=particle_a(i,1)-DAMP*particle_v(i,1);
        particle_a(i,2)=particle_a(i,2)-DAMP*particle_v(i,2);
    end
end

clear global;
close all;

for t=1:TIME_STEP_TOTAL
    if generate_method==1
        if PARTICLE_NUM<MAX_PARTICLE_NUM
            PARTICLE_NUM=PARTICLE_NUM+1;
            particles(t+1,1)=0;
            particles(t+1,2)=GRID_H/2;
            particle_v(t+1,1)=2+rand;
        end
    end
    particle_a=zeros(PARTICLE_NUM,2);
    particle_d=zeros(PARTICLE_NUM,2);
    particle_p=zeros(PARTICLE_NUM,2);

    particle_a(:,2)=GRAVITY;
    % update particles in grid
    solve_bound_collision(has_obstacle);
    compute_density();
    compute_acc();
    particle_v=particle_v+particle_a*DT;

    % update position
    particles=particles+particle_v*DT;
    
    if present_heatmap
        img=zeros(500,500,3);
        for i=1:PARTICLE_NUM
            if ~is_out_of_bound(particles(i,:))
                img(ceil(particles(i,2)/GRID_W*500),ceil(particles(i,1)/GRID_H*500))=img(ceil(particles(i,2)/GRID_W*500),ceil(particles(i,1)/GRID_H*500))+1;
                img(ceil(particles(i,2)/GRID_W*500),ceil(particles(i,1)/GRID_H*500),1)=1;
                img(ceil(particles(i,2)/GRID_W*500),ceil(particles(i,1)/GRID_H*500),2)=img(ceil(particles(i,2)/GRID_W*500),ceil(particles(i,1)/GRID_H*500),2)+particle_d(i)/1000;
            end
        end
        imshow(flipud(img));
        drawnow
    else
        plot(particles(:,1),particles(:,2),'b.');
        if has_obstacle
            viscircles([obstacle(1),obstacle(2)],obstacle(3));
        end
        axis equal
        xlim([0,GRID_W]);
        ylim([0,GRID_H]);
        drawnow
        obstacle(2)=obstacle(2)+0.001;
    end
    %imshow(flipud(grid));
    %h.GridVisible='off';
    %disp(max(particle_a(:,2),[],'all'))
    %set(gca, 'XLim', [0,10], 'YLim', [0,10])
    %plot(particles(:,1),particles(:,2),'ro', 'MarkerSize', 3)
    %set(gca, 'XLim', [0,10], 'YLim', [0,10])


end

end
