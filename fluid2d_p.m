function fluid2d_p()

GRAVITY=-10;
DT=0.01;
TIME_STEP_TOTAL=5000;
PARTICLE_NUM=500;
GRID_H=1;
GRID_W=1;

% particle properties
M=0.1;
DAMP=1;
SMOOTHING_LEN=0.01;
DENS_TO_P=0.04;

%img=zeros(GRID_H*5,GRID_W*5,3);

particles=zeros(PARTICLE_NUM,2);
px=SMOOTHING_LEN;
py=SMOOTHING_LEN+0.2;
for i=1:PARTICLE_NUM
    particles(i,1)=px+SMOOTHING_LEN*2;
    particles(i,2)=py;
    px=px+SMOOTHING_LEN*2;
    if px>GRID_W/4
        py=py+SMOOTHING_LEN*2;
        px=SMOOTHING_LEN;
    end
end
particle_v=zeros(PARTICLE_NUM,2);
particle_a=zeros(PARTICLE_NUM,2);
particle_d=zeros(PARTICLE_NUM,2);
particle_p=zeros(PARTICLE_NUM,2);

function solve_bound_collision()
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
                particle_d(j)=particle_d(j)+dens;
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
    particle_a=zeros(PARTICLE_NUM,2);
    particle_a(:,2)=GRAVITY;
    % update particles in grid
    solve_bound_collision();
    compute_density();
    compute_acc();
    particle_v=particle_v+particle_a*DT;

    % update position
    particles=particles+particle_v*DT;

    %img=zeros(GRID_H*5,GRID_W*5,3);
    %for i=1:PARTICLE_NUM
        %disp(particles(i,:))
        %if ~is_out_of_bound(particles(i,:))
            %disp(particles(i,:))
            %img(ceil(particles(i,2)*5),ceil(particles(i,1)*5))=img(ceil(particles(i,2)*5),ceil(particles(i,1)*5))+1;
            %img(ceil(particles(i,2)*5),ceil(particles(i,1)*5),1)=1;
            %img(ceil(particles(i,2)*5),ceil(particles(i,1)*5),2)=img(ceil(particles(i,2)*5),ceil(particles(i,1)*5),2)+particle_d(i)/50;
        %end
    %end
    
    %imshow(flipud(img));
    %imshow(flipud(grid));
    %h.GridVisible='off';
    %disp(max(particle_a(:,2),[],'all'))
    %set(gca, 'XLim', [0,10], 'YLim', [0,10])
    %plot(particles(:,1),particles(:,2),'ro', 'MarkerSize', 3)
    %set(gca, 'XLim', [0,10], 'YLim', [0,10])

    plot(particles(:,1),particles(:,2),'b.');
    axis equal
    xlim([0,GRID_W]);
    ylim([0,GRID_H]);
    drawnow

end

end
