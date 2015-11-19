function [velocities,atoms_old] = initialize_velocities(atoms,latvec,kb_T,dt);
%establish random velocities at temperature T
%Evan Reed
%February 12, 2010
%MatSci 331 HW #3

[natoms,temp]=size(atoms);
velocities=random('Uniform',-1,1,[natoms,3]);

%make center of mass velocity zero
vel_cm=sum(velocities)/natoms;
for k=1:natoms
    velocities(k,:)=velocities(k,:)-vel_cm;
end

kinetic_energy=0.5*sum(sum(velocities.*velocities));
target_kinetic_energy=1.5*kb_T*natoms;
scale_factor=sqrt(target_kinetic_energy/kinetic_energy);
velocities=velocities*scale_factor;

%set previous atom positions for use with original Verlet 
atoms_old=atoms-velocities*dt;






