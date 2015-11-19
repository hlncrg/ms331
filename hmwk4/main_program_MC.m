%Calculate LJ vibration frequencies
%Evan Reed
%January 29, 2011
%MatSci 331 HW #3

%use LJ units
%atom positions are LJ units (not scaled)

clear;
rejected_configs=0; %number of configurations rejected through Boltzmann criterion
kb_T=2;
beta=1.0/kb_T;
nsteps=1e6;
mc_max_move=0.2; %maximum spatial move in x, y, and z directions for single atom random displacement (in LJ units)
delta_kb_T=-kb_T/nsteps; %do simulated annealing
%delta_kb_T=0;

%initialize random number generator to produce same sequence every time the
%code is run.  This enables reproduction of results.
s=RandStream('mt19937ar');
RandStream.setDefaultStream(s);
reset(s,1);

%set size of computational cell
L=3;
M=3;
N=3;


%potential minimum is 2^(1/6)
%set lattice constant (cubic primitive cell)
lattice=sqrt(2)*2^(1/6);
%lattice=sqrt(2)*2^(1/6)*2.0;
%lattice=lattice*0.95;  %can scale lattice constant



rcut=1.3;

%set lattice vectors
latvec=[L*lattice 0 0; 0 M*lattice 0; 0 0 N*lattice];

%set up computational cell for perfect xtal
atoms=setup_cell(L,M,N,latvec);
[natoms,temp]=size(atoms);
atoms_start=atoms; %for use in calculating diffusion coefficient

%calculate the energy of the perfect crystal
if (rcut*2<latvec(1,1) && rcut*2<latvec(2,2) && rcut*2<latvec(3,3))
     [etot_perfect_xtal,forces]=calc_energy_faster(atoms,latvec,rcut,1);
     etot_perfect_xtal=etot_perfect_xtal/natoms;
     etot=etot_perfect_xtal;
else
    disp('Cell too small for calc_energy_faster');
    quit;
end

for step=1:nsteps
    
    step
    [etot_sequence(step),atoms,rejected_configs] = MC_move(atoms,latvec,natoms,mc_max_move,rcut,beta,etot,rejected_configs);
    etot=etot_sequence(step)
    rejected_configs/step
    r_squared(step)=calc_r_squared(atoms,atoms_start,latvec);

    kb_T=kb_T+delta_kb_T;
    beta=1/kb_T;
end
    


figure;
plot(etot_sequence,'k');


figure;
plot(r_squared,'k');





