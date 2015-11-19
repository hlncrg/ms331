%Calculate relaxed vacancy energy
%Evan Reed
%January 25, 2010
%MatSci 331 HW #2

%First, calculte the perfect crystal energy

clear;

%set size of computational cell
L=1;
M=L;
N=L;


%potential minimum is 2^(1/6)
%set lattice constant (cubic primitive cell)
lattice=sqrt(2)*2^(1/6);
lattice=lattice*0.95;  %can scale lattice constant



rcut=3.0;
force_tol=5e-3;
%force_tol=1e5;
alpha=2.5e-4;

%set lattice vectors
latvec=[L*lattice 0 0; 0 M*lattice 0; 0 0 N*lattice];

%set up computational cell for perfect xtal
atoms=setup_cell(L,M,N);
[natoms,temp]=size(atoms);
   
%calculate the energy of the perfect crystal
[etot_perfect_xtal,forces]=calc_energy(atoms,latvec,rcut,1);
etot_perfect_xtal=etot_perfect_xtal/natoms;


%loop over computational cell sizes
for L=3:3
    M=L;
    N=L;
    
    %set lattice vectors and create computational cell
    latvec=[L*lattice 0 0; 0 M*lattice 0; 0 0 N*lattice];
    clear atoms;
    atoms=setup_cell(L,M,N);

    %add an interstitial
    [natoms,temp]=size(atoms);
    atoms(natoms+1,:)=[0.5/L 0.5/M 0.5/N];

    %create a vacancy by removing an atom
    %[natoms,n]=size(atoms);
    %atoms=atoms(2:natoms,:);
   
    %minimize eneregy with respect to atom positions
    [etot,atoms_minimized]=minimize_energy(atoms,latvec,rcut,force_tol,alpha);

    %save some values for output
    vacancy_energy(L)=etot_perfect_xtal*(natoms-1)-etot;
    max_displacement(L)=max(max((atoms-atoms_minimized)*latvec));
end


etot_perfect_xtal
vacancy_energy
max_displacement




