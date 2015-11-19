function atoms = setup_cell(L,M,N);
%setup fcc computational cell
%Evan Reed
%January 25, 2010
%MatSci 331 HW #2


%define primitive cell
natoms=0;
basis(1,:)=[0 0 0];
basis(2,:)=[0.5 0.5 0];
basis(3,:)=[0 0.5 0.5];
basis(4,:)=[0.5 0 0.5];
nbasis=4;

%make periodic copies of the primitive cell
for l=0:(L-1)
for m=0:(M-1)
for n=0:(N-1)
    for k=1:nbasis
        atoms(natoms+k,:)=basis(k,:)+[l m n];
    end
    natoms=natoms+nbasis;
end
end
end
%make scaled atom coordinates, i.e. all atom positions range between 0 and 1
atoms(:,1)=atoms(:,1)/L;
atoms(:,2)=atoms(:,2)/M;
atoms(:,3)=atoms(:,3)/N;