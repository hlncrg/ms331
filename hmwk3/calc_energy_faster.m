function [etot,forces]=calc_energy_faster(atoms,latvec,rcut,force_flag)
%Calculate LJ energy of an arrangment of atoms with PBC
%This routine requires 2*rcut < computational cell lengths 
%Evan Reed
%Feb. 12, 2010
%MatSci 331 HW #3


[natoms,dummy]=size(atoms);


rcut_sqr=rcut*rcut;
rcut_6=rcut_sqr*rcut_sqr*rcut_sqr;
rcut_12=rcut_6*rcut_6;
ecut=(-1/rcut_6+1/rcut_12);
            

%define number of periodic copies of computational cell to consider
lmax=ceil(rcut/latvec(1,1));
mmax=ceil(rcut/latvec(2,2));
nmax=ceil(rcut/latvec(3,3));


etot=0;
forces=zeros(natoms,3);
%loop over atoms
for i=1:(natoms-1)
for j=(i+1):natoms
        %calculate displacement vector

        disp=atoms(i,:)-atoms(j,:);
        disp=disp-[round(disp(1)/latvec(1,1)) round(disp(2)/latvec(2,2)) round(disp(3)/latvec(3,3))]*latvec;

        %square of distance between atoms
        d_sqr=disp*disp';

        %only calculate energy for atoms within cutoff distance
        %don't calculate interactions between the same atom
        if (d_sqr<rcut_sqr && d_sqr > 0)
            inv_r_6=1/(d_sqr*d_sqr*d_sqr);
            inv_r_12=inv_r_6*inv_r_6;
            etot=etot-inv_r_6+inv_r_12-ecut;           
       
            %calcualte forces
            if (force_flag)
                fac=24*(2*inv_r_12-inv_r_6)/d_sqr;
                forces(i,:)=forces(i,:)+fac*disp;
                forces(j,:)=forces(j,:)-fac*disp;
                
            end
        end
            

end
end

etot=etot*4;
