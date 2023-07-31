close all;
clear all;
clc;

rng('shuffle');

f1 = fopen('internprobdata.out','w');
fclose(f1);

xmin = 0;
xmax = 2.8; %200 nm wide in both x and y directions, goes from 0 to 2.8um

discretespots = 2.8/0.01; %number of sites in one direction
Delta_x = (xmax - xmin)/discretespots; %finding lattice spacing
D = 0.2; %EGFR diffusivity
Delta_t = Delta_x^2/(4*D); %Delta t given the lattice spacing and diffusivity

%barrier heights for entering/exiting domains, in kBT
energy_enterclath_yeslig = 4.4;
energy_leaveclath_yeslig = 10;%
clathintern_rate = 0.05;%1/20; %rate at which EGFR in clathrin domain is internalized
ligunbind_rate = 0.05; %rate of ligand unbinding from egfr
clath_radius = 0.050; %radius of clathrin domains

numruns = 100;

internalizedcount = 0;
for ri=1:numruns
    %place clathrin domains
    clear clath_x clath_y
    overlap = 1;
    clath_x(1) = rand()*(xmax - xmin - 2*clath_radius) + xmin + clath_radius;
    clath_y(1) = rand()*(xmax - xmin - 2*clath_radius) + xmin + clath_radius;
    
    for i=2:10
        overlap = 1;
        while(overlap == 1)
            clath_x(i) = rand()*(xmax - xmin - 2*clath_radius) + xmin + clath_radius;
            clath_y(i) = rand()*(xmax - xmin - 2*clath_radius) + xmin + clath_radius;
    
            overlap = 0;
            for j=1:(i-1)
                distance = sqrt((clath_x(i) - clath_x(j))^2 + (clath_y(i) - clath_y(j))^2);
                if(distance <= (2*clath_radius))
                    overlap = 1;
                end
            end
        end
    end
    %this finishes the placement of clathrin domains
    
    %next bunch of lines starts the EGFR at a site that is outside of
    %clathrin domains
    overlap = 1;
    while(overlap == 1)
        xpos = Delta_x*floor((xmin + (xmax - xmin)*rand())/Delta_x);
        ypos = Delta_x*floor((xmin + (xmax - xmin)*rand())/Delta_x);
    
        overlap = 0;
        for j=1:10
            distance = sqrt((xpos - clath_x(j))^2 + (ypos - clath_y(j))^2);
            if(distance <= (clath_radius))
                overlap = 1;
            end
        end
    end
    %ends EGFR placement
    
    %initial conditions
    t = 0; %initial time
    clathrinintern = 0; %whether EGFR has been internalized
    ligunbind = 0;
    clathfirst = 0;
    count = 0; %counting number of timesteps
    
    clath_before = 0; %whether EGFR was in clathrin after last timestep (0 means no, 1 means yes)
    
    t_clath_yeslig = 0; %initializing time that EGFR is in a clathrin domain with a ligand bound
    lastclathenter = -1; %initializing tracking of the last time at which an EGFR entered a clathrin domain
    clathentercount = 0;
    clathleavecount = 0;
    
    %simulation runs until EGFR is internalized in clathrin or the
    %ligand unbinds
    while((clathrinintern == 0) && (ligunbind == 0))
        xpos_old = xpos; %recording last x position
        ypos_old = ypos; %recording last y position
    
        clath_before = 0; %recording whether in clathrin before next step
        for i=1:10
            distance = sqrt((xpos - clath_x(i))^2 + (ypos - clath_y(i))^2);
            if(distance < clath_radius)
                clath_before = 1;
            end
        end
    
        rannum = rand(); %deciding which direction to step
        if(rannum < 0.25)
            xpos = xpos - Delta_x;
        elseif(rannum < 0.5)
            xpos = xpos + Delta_x;
        elseif(rannum < 0.75)
            ypos = ypos - Delta_x;
        else
            ypos = ypos + Delta_x;
        end
    
        if(xpos > xmax) %enforcing periodic boundary conditions
            xpos = xpos - xmax;
        elseif(xpos < xmin)
            xpos = xmax + (xpos - xmin);
        end
    
        if(ypos > xmax)
            ypos = ypos - xmax;
        elseif(ypos < xmin)
            ypos = xmax + (ypos - xmin);
        end
    
        clath_after = 0;  %determining if step has finished in a clathrin domain
        for i=1:10
            distance = sqrt((xpos - clath_x(i))^2 + (ypos - clath_y(i))^2);
            if(distance < clath_radius)
                clath_after = 1;
                lastclathin1 = i;
            end
        end
   
        %enforcing energy barrier if exited a clathrin domain with bound
        %ligand
        if((clath_before == 1) && (clath_after == 0))
            if(rand() > exp(-energy_leaveclath_yeslig))
                xpos = xpos_old;
                ypos = ypos_old;
                clath_after = 1;
            end
        end
    
        %enforcing energy barrier if entered a clathrin domain with bound
        %ligand
        if((clath_before == 0) && (clath_after == 1))
            if(rand() > exp(-energy_enterclath_yeslig))
                xpos = xpos_old;
                ypos = ypos_old;
                clath_after = 0;
            else
                lastclathenter = t;
                lastclathin = lastclathin1;
            end
        end

        %record clathrin enter data
        if((clath_before == 0) && (clath_after == 1))
            clathentercount = clathentercount + 1;
            clathentertime(clathentercount) = t;
            if(clathentercount == 1)
                clathenterwait(clathentercount) = t;
                firstenterwait(ri) = t;
            else
                clathenterwait(clathentercount) = t - clathentertime(clathentercount - 1);
            end
        end

        %record clathrin exit data and move clathrin domain
        if((clath_before == 1) && (clath_after == 0))
            clathleavecount = clathleavecount + 1;
            clathleavetime(clathleavecount) = t;
            clathleavewait(clathleavecount) = t - clathentertime(clathentercount);

            %move clathrin domain
            overlap = 1;
            while(overlap == 1)
                clath_x(lastclathin) = rand()*(xmax - xmin - 2*clath_radius) + xmin + clath_radius;
                clath_y(lastclathin) = rand()*(xmax - xmin - 2*clath_radius) + xmin + clath_radius;
        
                overlap = 0;
                for j=1:10
                    if(j~=lastclathin)
                        distance = sqrt((clath_x(lastclathin) - clath_x(j))^2 + (clath_y(lastclathin) - clath_y(j))^2);
                        if(distance <= (2*clath_radius))
                            overlap = 1;
                        end
                    end
                end

                distance = sqrt((clath_x(lastclathin) - xpos)^2 + (clath_y(lastclathin) - ypos)^2);
                if(distance <= (clath_radius))
                    overlap = 1;
                end
            end
        end

        %record data if first clathrin entry
        if(clathfirst == 0)
            if(clath_after == 1)
                clathtime = t;
                clathfirst = 1;
            end
        end
    
        %if ligand is bound, determine if ligand unbinds
        if(rand() < (ligunbind_rate*Delta_t))
            ligunbind = 1;
        end
    
        %if ligand-bound in clathrin domain, determine if internalized
        if(clath_after == 1)
            if(rand() < (clathintern_rate*Delta_t))
                clathrinintern = 1;
            end
        end

        %if both ligand unbound and internalized, pick one proportional
        %to rates
        if((ligunbind == 1) && (clathrinintern == 1))
            if(rand() < (ligunbindrate/(ligunbindrate+clathintern_rate)))
                ligunbind = 1;
                clathrinintern = 0;
            else
                ligunbind = 0;
                clathrinintern = 1;
            end
        end
    
        %record any time spent in clathrin domains
        if(clath_after == 1)
            t_clath_yeslig = t_clath_yeslig + Delta_t;
        end
    
        %increment time forward
        t = t + Delta_t;
    
        count = count + 1;
    end
    
    t_rec(ri) = t; %record time until internalized

    if(clathfirst == 1)
        clathtime_rec(ri) = clathtime;
    else
        clathtime_rec(ri) = -1;
    end
    
    if(clathrinintern == 1)
        internalizedcount = internalizedcount + 1;
    end
    
    t_clath_yeslig_frac(ri) = t_clath_yeslig/t; %record fraction of total time that was spend in clathrin with ligand bound
    
    t_clath_yeslig_rec(ri) = t_clath_yeslig;

    lastclathenter_rec(ri) = lastclathenter; %record last time to enter a clathrin domain
    lastclathlifetime_rec(ri) = t - lastclathenter; %record time from last entry into clathrin until internalization

    clathentercount_rec(ri) = clathentercount;
    clathleavecount_rec(ri) = clathleavecount;
    if(clathentercount == 0)
        clathenterwait_avg(ri) = -1;
    else
        clathenterwait_avg(ri) = mean(clathenterwait(1:clathentercount));
    end
    if(clathleavecount == 0)
        clathleavewait_avg(ri) = -1;
    else
        clathleavewait_avg(ri) = mean(clathleavewait(1:clathleavecount));
    end
end

%calculate clathrin averages
totalclathenter = 0;
totalclathleave = 0;
clathenterwait_avgall = 0;
clathleavewait_avgall = 0;
for ri=1:numruns
    totalclathenter = totalclathenter + clathentercount_rec(ri);
    clathenterwait_avgall = clathenterwait_avgall + clathentercount_rec(ri)*clathenterwait_avg(ri);
    totalclathleave = totalclathleave + clathleavecount_rec(ri);
    clathleavewait_avgall = clathleavewait_avgall + clathleavecount_rec(ri)*clathleavewait_avg(ri);
end
clathenterwait_avgall = clathenterwait_avgall/totalclathenter;
clathleavewait_avgall = clathleavewait_avgall/totalclathenter;

internalized_count_avg = mean(internalizedcount)/numruns;
t_clath_yeslig_avg = mean(t_clath_yeslig_rec);

f1 = fopen('internprobdata.out','a');
fprintf(f1,'%E %E %E %E\n', ligunbind_rate,internalized_count_avg,...
    t_clath_yeslig_avg,clathenterwait_avgall);
fclose(f1);