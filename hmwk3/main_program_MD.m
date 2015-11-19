%Calculate LJ vibration frequencies
%Evan Reed
%January 25, 2010
%MatSci 331 HW #3

%use LJ units
%atom positions are LJ units (not scaled)

clear;
kb_T=0.2;
dt=2e-2;
nsteps=1e3;

%set size of computational cell
L=3;
M=3;
N=3;


%potential minimum is 2^(1/6)
%set lattice constant (cubic primitive cell)
lattice=sqrt(2)*2^(1/6);
%lattice=lattice*0.95;  %can scale lattice constant



rcut=1.3;

%set lattice vectors
latvec=[L*lattice 0 0; 0 M*lattice 0; 0 0 N*lattice];

%set up computational cell for perfect xtal
atoms=setup_cell(L,M,N,latvec);
[natoms,temp]=size(atoms);
   
%initialize velocities
[velocities,atoms_old]=initialize_velocities(atoms,latvec,kb_T,dt);


for time=1:nsteps
    if (rcut*2<latvec(1,1) && rcut*2<latvec(2,2) && rcut*2<latvec(3,3))
      [pot_e(time),forces]=calc_energy_faster(atoms,latvec,rcut,1);
    else
        [pot_e(time),forces]=calc_energy(atoms,latvec,rcut,1);
    end
    [atoms_new,velocities]=integrate(atoms,atoms_old,velocities,forces,latvec,dt);
    [instantaneous_kb_T(time),kin_e(time)]=calc_ke(velocities,lattice);
    total_energy(time)=kin_e(time)+pot_e(time);
    saved_velocities(time,:,:)=velocities;
    instantaneous_kb_T(time)
    atoms_old=atoms;
    atoms=atoms_new;
end

    





