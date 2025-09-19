c23456789012345678901234567890123456789012345678901234567890123456789012
      program readco

C*****************************************************
C*    CODE by Ozge Kurkcuoglu,                       *
C*    Istanbul Technical University, 2025	     *					  *
C*						     *
C*    YOU NEED TO SPECIFY:                           *
C*    -DIMENSION npar1, npar2(=3*npar1), npar3=20    *
C*    npar1 should be larger than number of residues *
C*    -CUTOFF DISTANCE rcutt=10 A is sufficient      *
C*    -WHICH RESIDUES ARE TO BE COARSE-GRAINED cg-cgx*
C*     cg and cgx are res indices in the network     *
C*    -EXACT TOTAL NUMBER OF AMINO ACIDS  rno        *
C*    -EXACT TOTAL NUMBER OF ATOMS res		     *
C*                                                   *
C*    If the script cannot read the pdb file, please *
C*    check the format of the example input file     *						     *
C*****************************************************

      integer npar1

      PARAMETER (npar1=4000)
      parameter (npar2=12000)
      parameter (npar3=20)


      character*6 dummy
      character*3 resnam,rest(npar1),att,ratom(npar1,55)
      character*3 restypp(npar1),restype(npar1),rrestyp(npar1)
      character*3 atnam,resatom(npar1,55),atom(npar1),c
      character*2 chain,chainold,reschain(npar1),rchain(npar1)
      character*2 chainn(npar1),cchain(npar1)
      character*2 at,resat(npar1),rcg(npar1)

      integer ires,ichain,irestot,m,n,cg,skip,cgx
      integer natoms(npar1),res,rno
      integer atnum,resnum,resold
      integer resnumo(npar1),mass(npar1,55),resn(npar1)
      integer dim1,dim2,nat(npar1)


      real x(npar1,55),y(npar1,55),z(npar1,55)
      real xx(npar1),yy(npar1),zz(npar1)
      real xx2(npar1),yy2(npar1),zz2(npar1)
      real xdist1(npar1),ydist1(npar1),zdist1(npar1)
      real xdist2(npar1),ydist2(npar1),zdist2(npar1)
      real xc(npar1,55),yc(npar1,55),zc(npar1,55)

      real bfactor(npar1,55),bfact(npar1),mw2(npar1)
      real t(3,3),mw(npar1),bfac(npar1),rcut(npar1)
	real centx(npar1),centy(npar1),centz(npar1),centbf(npar1)
	real contact(npar1,npar1),r(npar1,npar1)
	real molwt(npar1),ma(npar1),mm(npar2)
	real w(3),v(3,3),ez(npar1,3),a(npar1,3)
	real prod(3,3),prod2(3,npar1)
	real xxnew(npar1),yynew(npar1),zznew(npar1)
	real cm(npar2),ds(npar1),smo(npar1)
        real cmod(npar1,npar3),csum(npar1,npar3)
        real det,ezz(3,npar1),aaa(3,npar1),vv(3,3)
        	
      real dr2(npar1,npar2),sumdr2(npar2),alpha(npar2)
      real sum(npar2),collect(npar2)

      common /test1/ w1(npar2),v1(npar2,npar2),dn(npar2,npar2)

c--------cutoff distance-----------------------------------

	rcutt=10

c--------exact total number of amino acids is rno-----------
c--------total number of atoms is res-----------------------

	rno=611	
	rrn=rno
        res=4728
        
c--------from 1 to cg and from cgx+1 to res is low-resolution--------
c--------from cg+1 to cgx is high resolution-------------------------
c--------if you type same number for cg and cgx, it will not mcg-----
	cg=610
	cgx=610
c-----------------------------------------------------------
	skip=1

        OPEN (10,file='7aku-drugA.pdb',status='OLD')

        OPEN (11,file='check.coor4')
        open (49,file='coordinates_com.out')
        open (12,file='principal_axes1.out')
	open (13,file='principal_axes2.out')
	open (14,file='min_coordinates.out')

        open (6000,file='node_mass.out')
	open(21,file='eigenvalues.out')
	open(22,file='eigen_cumulative.out')

 	open(31,file='forceconst.out')
	open(32,file='bfactors.out')
	open(33,file='slowest.out')
	open(34,file='slowest_cumulative.out')
	open(35,file='fastest_cumulative.out')

	open (61,file='mod1a.pdb')
        open (62,file='mod1b.pdb')
      	open (63,file='mod2a.pdb')
        open (64,file='mod2b.pdb')
      	open (65,file='mod3a.pdb')
      	open (66,file='mod3b.pdb')
      	open (67,file='mod4a.pdb')
      	open (68,file='mod4b.pdb')
      	
        open (110,file='first100modes.out')
        
c*********************************************************************


      do j=1,npar1
         xx(j)=0.
         yy(j)=0.
         zz(j)=0.
         xx2(j)=0.
         yy2(j)=0.
         zz2(j)=0.
	rcut(j)=0.
      enddo
  
      do i=1,npar1
      do j=1,55
         x(i,j)=0.
         y(i,j)=0.
         z(i,j)=0.
      enddo
      enddo

 50   Format (A6) 
 52   Format(i6,3x,A3,2x,4f8.3)    
 55      format(a4,9x,a2,2x,a3,2x,i4,f12.3,1x,f7.3,1x,f7.3,
     :1x,f5.2,f6.2)
 56   Format(A3)
 57   Format(A6,I5,2X,a2,2x,A3,2X,I4,4X,3F8.3,1x,f4.2,2x,f5.1)
 77   format (A6,I5,2X,a2,2x,A3,2X,I4,4X,3F8.3,6X,f6.2)
 78   format (A6,I5,2X,a1,3x,A3,A2,I4,4X,3F8.3,6X,f6.2,1x,i3)
 79   format (A6,I5,2X,A3)

      i=0
	icount=0


 5    read(10,50)dummy

      if(dummy.ne.'ATOM') goto 5
c	if(dummy.ne.'HETATM') goto 5

      backspace (10)

	do k=1,res
 12   i=i+1


 31   read(10,50)dummy
c      if(dummy.eq.'HETATM') goto 31
  	 if(dummy.eq.'CONECT') goto 33
	if(dummy.eq.'MASTER') goto 33
         if(dummy.eq.'ANISOU') goto 31
      if(dummy.eq.'TER') then
         iter=iter+1
         goto 31
      endif 

      if(dummy.eq.'SIGATM') goto 31
      backspace (10)


	read(10,77)dummy,atnum,at,resnam,resnum,xa,ya,za,bf

	backspace(10)

	read(10,79)dummy,atnum,att
        write(11,77)dummy,atnum,at,resnam,resnum,xa,ya,za,bf


	icount=icount+1

	if(at.eq.'N '.or.at.eq.'NE'.or.at.eq.'NH'.or.at.eq.
     :    'ND'.or.at.eq.
     :       'NZ'.or.at.eq.'N1'.or.at.eq.'N2'.or.at.eq.'N3'
     :    .or.at.eq.'N4'
     :       .or.at.eq.'N5'.or.at.eq.'N6'.or.at.eq.'N7'.or.at.
     :       eq.'N8'.or.at.eq.
     :   'N9') then
	m=14
        endif
  
 	
	if(at.eq.'C '.or.at.eq.'CA'.or.at.eq.'CB'.or.at.eq.
     :   'CD'.or.at.eq.
     :     'CE'.or.at.eq.'CG'.or.at.eq.'CH'.or.at.eq.'CZ'.or.at.eq.'C1'
     :    .or.at.eq.'C2'.or.at.eq.'C3'.or.at.eq.'C4'.or.at.
     :     eq.'C5'.or.at.eq.
     :     'C6'.or.at.eq.'C7'.or.at.eq.'C8'.or.at.eq.'C9'.or.at.eq.'CC'
     :    .or.at.eq.'CF'.or.at.eq.'CZ') then
	m=12
        endif

	if(at.eq.'O '.or.at.eq.'OD'.or.at.eq.'OE'.or.at.
     :    eq.'OG'.or.at.eq.'OH'
     :    .or.at.eq.'OX'.or.at.eq.'O1'.or.at.eq.'O2'.or.
     :    at.eq.'O3'.or.at.eq.
     :    'O4'.or.at.eq.'O5'.or.at.eq.'O6'.or.at.eq.'O7'.or.at.eq.
     :   'O8'.or.at.eq
     :    .'O9'.or.at.eq.'OT') then
	m=16
        endif
	

	if(at.eq.'S '.or.at.eq.'SD'.or.at.eq.'SG') then
	m=32
        endif
	
	if(at.eq.'P '.or.at.eq.'P1') then
	m=31
        endif
	

	if(atnum.eq.1) then
         resold=resnum
         chainold=chain
         ires=1
         ichain=1
      endif
 
      if(resnum.ne.resold)then
         natoms(ires)=i-1
         ires=ires+1
         i=1
         resold=resnum
      endif

      if(chain.ne.chainold)then
         ichain=ichain+1
         chainold=chain
      endif


      x(ires,i)=xa
      y(ires,i)=ya
      z(ires,i)=za

	
	mass(ires,i)=m
        restypp(ires)=resnam
	rest(ires)=resnam
        resnumo(ires)=resnum
	resatom(ires,i)=at
        ratom(ires,i)=att
        bfactor(ires,i)=bf
        reschain(ires)=chain
	reschain(ires)=chainn(i)

	xx(icount)=xa
	yy(icount)=ya
	zz(icount)=za

c      go to 12
	enddo

 33   maxat=atnum
      natoms(ires)=i
	
      

c ------------low-resolution--------------------------------------------
c In the low resolution, one-node-per-residue (CA) is taken for proteins
c three-nodes-per-nucleotide (P, C4*, C2') are taken for RNA and DNA
        	
	do k=1,cg,skip
	centx(i)=0
	centy(i)=0
	centz(i)=0
	say=say+1
	enddo
	
	
	icount=0
	nn=0
	do k=1,cg

	if(restypp(k).ne.'A  '.and.restypp(k).ne.'T  '
     :     .and.restypp(k).ne.'U  '.and.restypp(k).ne.'G  '
     :     .and.restypp(k).ne.'C  '.and.restypp(k).ne.'GUA'
     :     .and.restypp(k).ne.'CYT'.and.restypp(k).ne.'THY'
     :     .and.restypp(k).ne.'ADE'.and.restypp(k).ne.'  T'
     :     .and.restypp(k).ne.'  C'.and.restypp(k).ne.'  G'
     :     .and.restypp(k).ne.'  A'.and.restypp(k).ne.'  U')
     :     then

	icount=icount+1
	nn=say
	do j=1,natoms(k)
	mw(icount)=mw(icount)+mass(k,j)
	mw2(icount)=mw2(icount)+mass(k,j)*mass(k,j)
	bfac(icount)=bfac(icount)+mass(k,j)*mass(k,j)*bfactor(k,j)
	centx(icount)=centx(icount)+x(k,j)*mass(k,j)
	centy(icount)=centy(icount)+y(k,j)*mass(k,j)
	centz(icount)=centz(icount)+z(k,j)*mass(k,j)

	resat(icount)=resatom(k,j)

	xc(icount,j)=x(k,j)
	yc(icount,j)=y(k,j)
	zc(icount,j)=z(k,j)


	enddo

	bfac(icount)=bfac(icount)/mw(icount)/mw(icount)
        centx(icount)=centx(icount)/mw(icount)
	centy(icount)=centy(icount)/mw(icount)
	centz(icount)=centz(icount)/mw(icount)

	rrestyp(icount)=restypp(k)
	resn(icount)=resnumo(k)
       reschain(icount)=reschain(k)
	rcg(icount)='CA'
	nat(icount)=natoms(k)
	molwt(icount)=mw(icount)

	endif

	if(restypp(k).eq.'A  '.or.restypp(k).eq.'T  '
     :     .or.restypp(k).eq.'U  '.or.restypp(k).eq.'G  '
     :     .or.restypp(k).eq.'C  '.or.restypp(k).eq.'GUA'
     :     .or.restypp(k).eq.'CYT'.or.restypp(k).eq.'THY'
     :     .or.restypp(k).eq.'ADE'.or.restypp(k).eq.'  G'
     :     .or.restypp(k).eq.'  C'.or.restypp(k).eq.'  A'
     :     .or.restypp(k).eq.'  T'.or.restypp(k).eq.'  U')
     :       then

	nn=nn+3
c --- phosphate group----------------------------------
        icount=icount+1
	do j=1,3

	mw(icount)=mw(icount)+mass(k,j)
	mw2(icount)=mw2(icount)+mass(k,j)*mass(k,j)
	bfac(icount)=bfac(icount)+mass(k,j)*mass(k,j)*bfactor(k,j)
	centx(icount)=centx(icount)+x(k,j)*mass(k,j)
	centy(icount)=centy(icount)+y(k,j)*mass(k,j)
	centz(icount)=centz(icount)+z(k,j)*mass(k,j)

	resat(icount)=resatom(k,j)

	xc(icount,j)=x(k,j)
	yc(icount,j)=y(k,j)
	zc(icount,j)=z(k,j)
	enddo
	bfac(icount)=bfac(icount)/mw(icount)/mw(icount)
        centx(icount)=centx(icount)/mw(icount)
	centy(icount)=centy(icount)/mw(icount)
	centz(icount)=centz(icount)/mw(icount)
	rcg(icount)='P '
	rrestyp(icount)=restypp(k)
	resn(icount)=resnumo(k)
       reschain(icount)=reschain(k)
	nat(icount)=3
	molwt(icount)=mw(icount)

c --- sugar group--------------------------------------	
	icount=icount+1

	do j=4,10

	mw(icount)=mw(icount)+mass(k,j)
	mw2(icount)=mw2(icount)+mass(k,j)*mass(k,j)
	bfac(icount)=bfac(icount)+mass(k,j)*mass(k,j)*bfactor(k,j)
	centx(icount)=centx(icount)+x(k,j)*mass(k,j)
	centy(icount)=centy(icount)+y(k,j)*mass(k,j)
	centz(icount)=centz(icount)+z(k,j)*mass(k,j)

	resat(icount)=resatom(k,j)

	xc(icount,j)=x(k,j)
	yc(icount,j)=y(k,j)
	zc(icount,j)=z(k,j)
	enddo
	bfac(icount)=bfac(icount)/mw(icount)/mw(icount)
        centx(icount)=centx(icount)/mw(icount)
	centy(icount)=centy(icount)/mw(icount)
	centz(icount)=centz(icount)/mw(icount)
	rcg(icount)='C4'
	rrestyp(icount)=restypp(k)
	resn(icount)=resnumo(k)
       reschain(icount)=reschain(k)
	nat(icount)=8
	molwt(icount)=mw(icount)

c --- base---------------------------------------------

	icount=icount+1

	do j=11,natoms(k)

	mw(icount)=mw(icount)+mass(k,j)
	mw2(icount)=mw2(icount)+mass(k,j)*mass(k,j)
	bfac(icount)=bfac(icount)+mass(k,j)*mass(k,j)*bfactor(k,j)
	centx(icount)=centx(icount)+x(k,j)*mass(k,j)
	centy(icount)=centy(icount)+y(k,j)*mass(k,j)
	centz(icount)=centz(icount)+z(k,j)*mass(k,j)

	resat(icount)=resatom(k,j)

	xc(icount,j)=x(k,j)
	yc(icount,j)=y(k,j)
	zc(icount,j)=z(k,j)
	enddo
	bfac(icount)=bfac(icount)/mw(icount)/mw(icount)
        centx(icount)=centx(icount)/mw(icount)
	centy(icount)=centy(icount)/mw(icount)
	centz(icount)=centz(icount)/mw(icount)
	rcg(icount)='C2'
	rrestyp(icount)=restypp(k)
	resn(icount)=resnumo(k)
       reschain(icount)=reschain(k)
	nat(icount)=natoms(k)-10
	molwt(icount)=mw(icount)

	endif


        enddo


c-------------------------------------------------------------------

c--------high resolution---------------------------------
	

	
	do k=cg+1,cgx
	do j=1,natoms(k)
	
	icount=icount+1
	ncount=ncount+1

	centx(icount)=x(k,j)
	centy(icount)=y(k,j)
	centz(icount)=z(k,j)

	bfac(icount)=bfactor(k,j)
	reschain(icount)=reschain(k)
	rcg(icount)=resatom(k,j)
	
	rrestyp(icount)=restypp(k)
	resn(icount)=resnumo(k)

	xc(icount,1)=x(k,j)
	yc(icount,1)=y(k,j)
	zc(icount,1)=z(k,j)

	molwt(icount)=mass(k,j)
	enddo
	enddo
c-------------------------------------------------------------


c ---------low resolution-------------------------------------

	do k=cgx+1,rno

      	if(restypp(k).ne.'A  '.and.restypp(k).ne.'T  '
     :     .and.restypp(k).ne.'U  '.and.restypp(k).ne.'G  '
     :     .and.restypp(k).ne.'C  '.and.restypp(k).ne.'GUA'
     :     .and.restypp(k).ne.'CYT'.and.restypp(k).ne.'THY'
     :     .and.restypp(k).ne.'ADE'.and.restypp(k).ne.'  U'
     :    .and.restypp(k).ne.'  T'.and.restypp(k).ne.'  A'
     :    .and.restypp(k).ne.'  C'.and.restypp(k).ne.'  G')
     :     then



	icount=icount+1

	do j=1,natoms(k)
	mw(icount)=mw(icount)+mass(k,j)
	mw2(icount)=mw2(icount)+mass(k,j)*mass(k,j)
	bfac(icount)=bfac(icount)+mass(k,j)*mass(k,j)*bfactor(k,j)
	centx(icount)=centx(icount)+x(k,j)*mass(k,j)
	centy(icount)=centy(icount)+y(k,j)*mass(k,j)
	centz(icount)=centz(icount)+z(k,j)*mass(k,j)

	resat(icount)=resatom(k,j)

	xc(icount,j)=x(k,j)
	yc(icount,j)=y(k,j)
	zc(icount,j)=z(k,j)


	enddo

	bfac(icount)=bfac(icount)/mw(icount)/mw(icount)
        centx(icount)=centx(icount)/mw(icount)
	centy(icount)=centy(icount)/mw(icount)
	centz(icount)=centz(icount)/mw(icount)

	rrestyp(icount)=restypp(k)
	resn(icount)=resnumo(k)
       reschain(icount)=reschain(k)
	rcg(icount)='CA'
	nat(icount)=natoms(k)
	molwt(icount)=mw(icount)

	endif
	
 	if(restypp(k).eq.'A  '.or.restypp(k).eq.'T  '
     :     .or.restypp(k).eq.'U  '.or.restypp(k).eq.'G  '
     :     .or.restypp(k).eq.'C  '.or.restypp(k).eq.'GUA'
     :     .or.restypp(k).eq.'CYT'.or.restypp(k).eq.'THY'
     :     .or.restypp(k).eq.'ADE'.or.restypp(k).eq.'  U'
     :   .or.restypp(k).eq.'  A'.or.restypp(k).eq.'  C'
     :   .or.restypp(k).eq.'  G'.or.restypp(k).eq.'  T')
     :    then

c --- phosphate group----------------------------------
        icount=icount+1
	do j=1,3

	mw(icount)=mw(icount)+mass(k,j)
	mw2(icount)=mw2(icount)+mass(k,j)*mass(k,j)
	bfac(icount)=bfac(icount)+mass(k,j)*mass(k,j)*bfactor(k,j)
	centx(icount)=centx(icount)+x(k,j)*mass(k,j)
	centy(icount)=centy(icount)+y(k,j)*mass(k,j)
	centz(icount)=centz(icount)+z(k,j)*mass(k,j)

	resat(icount)=resatom(k,j)

	xc(icount,j)=x(k,j)
	yc(icount,j)=y(k,j)
	zc(icount,j)=z(k,j)
	enddo
	bfac(icount)=bfac(icount)/mw(icount)/mw(icount)
        centx(icount)=centx(icount)/mw(icount)
	centy(icount)=centy(icount)/mw(icount)
	centz(icount)=centz(icount)/mw(icount)
	rcg(icount)='P '
	rrestyp(icount)=restypp(k)
	resn(icount)=resnumo(k)
       reschain(icount)=reschain(k)
	nat(icount)=3
	molwt(icount)=mw(icount)

c --- sugar group--------------------------------------	
	icount=icount+1

	do j=4,10

	mw(icount)=mw(icount)+mass(k,j)
	mw2(icount)=mw2(icount)+mass(k,j)*mass(k,j)
	bfac(icount)=bfac(icount)+mass(k,j)*mass(k,j)*bfactor(k,j)
	centx(icount)=centx(icount)+x(k,j)*mass(k,j)
	centy(icount)=centy(icount)+y(k,j)*mass(k,j)
	centz(icount)=centz(icount)+z(k,j)*mass(k,j)

	resat(icount)=resatom(k,j)

	xc(icount,j)=x(k,j)
	yc(icount,j)=y(k,j)
	zc(icount,j)=z(k,j)
	enddo
	bfac(icount)=bfac(icount)/mw(icount)/mw(icount)
        centx(icount)=centx(icount)/mw(icount)
	centy(icount)=centy(icount)/mw(icount)
	centz(icount)=centz(icount)/mw(icount)
	rcg(icount)='C4'
	rrestyp(icount)=restypp(k)
	resn(icount)=resnumo(k)
       reschain(icount)=reschain(k)
	nat(icount)=7
	molwt(icount)=mw(icount)

c --- base---------------------------------------------

	icount=icount+1

	do j=11,natoms(k)

	mw(icount)=mw(icount)+mass(k,j)
	mw2(icount)=mw2(icount)+mass(k,j)*mass(k,j)
	bfac(icount)=bfac(icount)+mass(k,j)*mass(k,j)*bfactor(k,j)
	centx(icount)=centx(icount)+x(k,j)*mass(k,j)
	centy(icount)=centy(icount)+y(k,j)*mass(k,j)
	centz(icount)=centz(icount)+z(k,j)*mass(k,j)

	resat(icount)=resatom(k,j)

	xc(icount,j)=x(k,j)
	yc(icount,j)=y(k,j)
	zc(icount,j)=z(k,j)
	enddo
	bfac(icount)=bfac(icount)/mw(icount)/mw(icount)
        centx(icount)=centx(icount)/mw(icount)
	centy(icount)=centy(icount)/mw(icount)
	centz(icount)=centz(icount)/mw(icount)
	rcg(icount)='C2'
	rrestyp(icount)=restypp(k)
	resn(icount)=resnumo(k)
       reschain(icount)=reschain(k)
	nat(icount)=natoms(k)-10
	molwt(icount)=mw(icount)

	endif



        enddo
c-------------------------------------------------------------

         
c--------new number of nodes is cgx+icount---------------------

	rno=icount
	write(6,*)'total number of nodes',rno
c-------------------------------------------------------------

	do k=nn+1,nn+ncount
	nat(k)=1
	enddo

c-------------------------------------------------------------
	
c     contact number between two different nodes 
c     (high-high resolution), (low-low resolution), (high-low resolution) 
c     is count

	do i=1,rno
	do j=1,rno

	do k=1,nat(i)
	do l=1,nat(j)


	dx=xc(i,k)-xc(j,l)
	dy=yc(i,k)-yc(j,l)
	dz=zc(i,k)-zc(j,l)

	r2=sqrt(dx*dx+dy*dy+dz*dz)
	if(i.ne.j.and.r2.le.rcutt) then
	contact(i,j)=contact(i,j)+1
	endif

	enddo
	enddo

	enddo
	enddo


	do i=1,rno
	write(6000,*)rcg(i),molwt(i)
      	mwt=mwt+molwt(i)

        enddo


c --------------------------------------------------------------------

	do i=1,rno
        
	   cx=cx+centx(i)*molwt(i)
           cy=cy+centy(i)*molwt(i)
           cz=cz+centz(i)*molwt(i)
       
        enddo

	cx=cx/mwt
	cy=cy/mwt
	cz=cz/mwt

	do i=1,rno
           centx(i)=-cx+centx(i)
           centy(i)=-cy+centy(i)
           centz(i)=-cz+centz(i)
        enddo

	write(6,*)'c.o.m.',cx,cy,cz
	write(49,*)'Coordinates with respect to mass center'

	do i=1,rno
           write(49,52)i,rrestyp(i),centx(i),centy(i),centz(i),molwt(i)
	
        enddo


c ------Generation of the (3nx3n) mass matrix------------------------
	
	do i=1,rno
	ma(i)=1/sqrt(molwt(i))

	enddo

	do i=1,3*rno
	mm(i)=0.
	enddo

 	do j=1,rno

c diagonals (for j)
              mm(3*j-2)=ma(j)
              mm(3*j-1)=ma(j)
              mm(3*j)=ma(j)

	enddo



c----------------------------------------------------------------------
c     Calculation of the principal axes of the molecule

c     Coordinates will be saved in matrix a (of size resnumx3)
c     which will be copied to matrix ez, 
c     before submission to the subroutine

c     dim2 is the number of rows in a or ez (=resnum)
c     dim1 is the number of columns (=3)

        do i=1,rno
           a(i,1)=centx(i)
           a(i,2)=centy(i)
           a(i,3)=centz(i)
        enddo

	do i=1,rno
	do j=1,3
           ez(i,j)=a(i,j)
        enddo
        enddo
	

	dim1=3
	dim2=rno

	
C********************************************************************
        call svdcmp(ez,dim2,dim1,dim2,dim1,w,v)
C********************************************************************


C     THE MATRICES THAT RETURN FROM THE SUBROUTINE ARE
C     EZ[dim2,dim1] IS THE V MATRIX
C     V[dim1,dim1]  IS THE U MATRIX
C     W[dim1] vector HAS THE DIAGONAL ELEMENTS OF SIGMA (S) MATRIX,
C     if we adopt the notation A = U S VT  

      write(12,*)'Performing checks on principal axes calculations'

c     checking if VVT is identity matrix

        write(12,*)
        write(12,*)'if VVT is identity matrix'

	do i=1,dim1
	do j=1,dim1
           summ=0.
           do k=1,dim1
              summ=summ+v(i,k)*v(j,k)
           enddo
           write(12,*)'i,j,summ=  ',i,j,summ
	enddo
        enddo
	write(12,*)      

	write(12,*) 'singular values and elements of u(i,j)'
	write(12,*) 
        do i=1,dim1
           write(12,113)i,w(i),(v(i,j),j=1,dim1)
        enddo
	write(12,*)

  113	format(i4,f12.6,5x,3f9.3)

c---------------------------------------------- 

c       checking whether A =  U S VT

 	do i=1,dim1
 	do j=1,dim1
           prod(i,j)=v(i,j)*w(j)
       enddo
        enddo

 	do i=1,dim1
 	do j=1,dim2
           prod2(i,j)=0.
           do k=1,dim1
              prod2(i,j)=prod2(i,j)+prod(i,k)*ez(j,k)
           enddo   
        enddo
        enddo

	do i=1,dim1
	do j=1,dim2
           dif=prod2(i,j)-a(j,i)
           if(dif.gt.0.01) write(6,*)i,j,'dif=',a(j,i),dif
        enddo
        enddo

c---------------------------------------------- 


c	PRINCIPAL AXES

        write(13,*)'Principal axes'
	do i=1,3
           write(13,*)(w(i)/4.*v(j,i),j=1,3)
        enddo

         det=v(1,1)*v(2,2)*v(3,3)-v(1,1)*v(1,3)*v(3,2)
     :      -v(1,2)*v(2,1)*v(3,3)
     :  +v(1,2)*v(2,3)*v(3,1)+v(1,3)*v(2,1)*v(3,2)
     :      -v(1,3)*v(2,2)*v(3,1)

          if(det.lt.0) then

        w(3)=-w(3)
      do i=1,3
      
         v(i,3)=-v(i,3)
      enddo

         write(13,*)'Principal axes 2'
        do i=1,3
         write(13,*)(w(i)/4.*v(j,i),j=1,3)
        enddo
      
      endif
      

	write(13,*) 'Singular values and elements of u(i,j)'
	write(13,*) 
        do i=1,dim1
           write(13,113)i,w(i),(v(i,j),j=1,dim1)
        enddo
	write(13,*)

        close(13)

c---------------------------------------------- 
c	TRANSFORMATION OF THE COORDINATES TO THEIR REPRESENTATION
c	IN THE FRAME SPANNED BY THE PRINCIPAL AXES


c	transformation matrix is UT, and after operating on A it
c	simply yields S VT. Therefore the transformed coordinates are
c	the elements of VT, multiplied by eigenvalues.


c       take the transpose of v(3,3)

	do i=1,3
	do j=1,3
	vv(i,j)=v(j,i)
	enddo
	enddo

	do i=1,dim2
	do j=1,3
	aaa(j,i)=a(i,j)
	enddo
	xxnew(i)=0.
	yynew(i)=0.
	zznew(i)=0.
	enddo
	
 	do j=1,dim1
           do i=1,dim2
        ezz(j,i)=vv(j,1)*aaa(1,i)+vv(j,2)*aaa(2,i)+vv(j,3)*aaa(3,i)
           enddo
	enddo

	do i=1,dim2
           xxnew(i)=ezz(1,i)
           yynew(i)=ezz(2,i)
           zznew(i)=ezz(3,i)
        
	enddo

	dummy='ATOM  '

	do  i=1,rno
           write(14,57)dummy,i,rcg(i),rrestyp(i),resn(i),xxnew(i)
     :     ,yynew(i),zznew(i),1.00,bfac(i)
        enddo

	write(14,56)'TER'

c-----------------------------------------------------------------------

c	BEGINNING OF ANM CALCULATIONS

	res1=3*rno
	
c-----------------------------------------------------------------------

c	GENERATION OF THE HESSIAN

	icontact=0

        do j=1,rno
        do k=1,rno
           ga=1.
	
           bx=xxnew(j)-xxnew(k)
           by=yynew(j)-yynew(k)
           bz=zznew(j)-zznew(k)
           r(j,k)=sqrt(bx*bx+by*by+bz*bz)

           if(j.ne.k) then


c diagonals (for j)
              dn(3*j-2,3*j-2)=dn(3*j-2,3*j-2)+
     :           contact(j,k)*ga*bx*bx/r(j,k)**2.
              dn(3*j-1,3*j-1)=dn(3*j-1,3*j-1)+
     :           contact(j,k)*ga*by*by/r(j,k)**2.
              dn(3*j,3*j)=dn(3*j,3*j)+
     :           contact(j,k)*ga*bz*bz/r(j,k)**2.
c off-diagonals of diagonal superelements (for j)
              dn(3*j-2,3*j-1)=dn(3*j-2,3*j-1)+
     :           contact(j,k)*ga*bx*by/r(j,k)**2.
              dn(3*j-2,3*j)=dn(3*j-2,3*j)+
     :           contact(j,k)*ga*bx*bz/r(j,k)**2.
              dn(3*j-1,3*j-2)=dn(3*j-1,3*j-2)+
     :           contact(j,k)*ga*by*bx/r(j,k)**2.
              dn(3*j-1,3*j)=dn(3*j-1,3*j)+
     :           contact(j,k)*ga*by*bz/r(j,k)**2.
              dn(3*j,3*j-2)=dn(3*j,3*j-2)+
     :           contact(j,k)*ga*bx*bz/r(j,k)**2.
              dn(3*j,3*j-1)=dn(3*j,3*j-1)+
     :           contact(j,k)*ga*by*bz/r(j,k)**2.
           endif
        enddo
        enddo


        do j=1,rno
        do k=1,rno
           ga=1.0
           bx=xxnew(j)-xxnew(k)
           by=yynew(j)-yynew(k)
           bz=zznew(j)-zznew(k)
           r(j,k)=sqrt(bx*bx+by*by+bz*bz)

           if(j.ne.k) then

c diagonals (for j&k)
              dn(3*j-2,3*k-2)=-contact(j,k)*ga*bx*bx/r(j,k)**2.
              dn(3*j-1,3*k-1)=-contact(j,k)*ga*by*by/r(j,k)**2.
              dn(3*j,3*k)=-contact(j,k)*ga*bz*bz/r(j,k)**2.
c off-diagonals (for j&k)
              dn(3*j-2,3*k-1)=-contact(j,k)*ga*bx*by/r(j,k)**2.
              dn(3*j-2,3*k)=-contact(j,k)*ga*bx*bz/r(j,k)**2.
              dn(3*j-1,3*k-2)=-contact(j,k)*ga*by*bx/r(j,k)**2.
              dn(3*j-1,3*k)=-contact(j,k)*ga*by*bz/r(j,k)**2.
              dn(3*j,3*k-2)=-contact(j,k)*ga*bx*bz/r(j,k)**2.
              dn(3*j,3*k-1)=-contact(j,k)*ga*by*bz/r(j,k)**2.
           endif
       enddo
       enddo
	write(6,*)'multiplication starts here'	

c Here the hessian matrix is multiplied by mass matrix as
c F=m**-1/2 * H * m**-1/2


	do i=1,3*rno
	do j=1,3*rno

 	dn(i,j)=mm(i)*dn(i,j)*mm(j) 	
   	
	enddo
	enddo



c-----------------------------------------------------------------------

C     INVERSION OF THE HESSIAN AND CORRELATIONS

      jres=res1

C********************************************************************
      call svdcmh(jres,jres,jres,jres)
C********************************************************************

c     dn(3nx3n) is  sent to the subroutine via the COMMON statement.

c     The Hessian (dn) is now transformed into H = U Lambda VT.
 
c     Here, the matrices that return from the subroutine via the 
c     COMMON statement are:
c     v1 is our matrix U(3nx3n). 
c     We should have dn transpose = v1.
c     w1(3n) gives the eigenvalues of the diagonal matrix Lambda.

c-----------------------------------------------------------------------
c     EIGENVALUES


      do i=1,jres
      do j=i,jres
         if(w1(i).lt.w1(j)) then
            top1=w1(i)
            w1(i)=w1(j)
            w1(j)=top1
            do k=1,jres     
               top1=v1(k,i)
               v1(k,i)=v1(k,j)
               v1(k,j)=top1
               top2=dn(k,i)
               dn(k,i)=dn(k,j)
               dn(k,j)=top2
            enddo
         endif
      enddo
      enddo

c     The eigenvalues are organized in descending order,
c     i.e. the last eigenvalues refer to the slowest modes.
c     The last 6 (= kc) eigenvalues should be zero.


c------------------------------------------------------------------
c     Counting the number of zero eigenvalues as k
        write(21,*)'mode no    eigenvalue   Sum(1/eigenvalue)'
      kc=0
      eigsum=0.
      do i=1,jres
         if(w1(i).le.1.0e-05) then
            kc=kc+1
         else
            eigsum=eigsum+1./w1(i)
         endif

         write(21,*)i,w1(i),eigsum
      enddo


      write(31,*)'kc',kc
      if(kc.ne.6) then
         write(6,*)'PROBLEM WITH ZERO EIGENVALUES!'
         stop
      endif

      eigcum=0.
      do i=1,jres-kc
         j=jres-kc-i+1
         aa=1./eigsum/w1(j)
         eigcum=eigcum+aa
         write(22,*)i,aa,eigcum,j
      enddo

      close(21)
      close(22)

c-------------------------------------------------------


c-------------------------------------------------------
c     B-FACTORS
c     Summing up over all modes to calculate B-factors

      do i=1,jres
         cm(i)=0.
         do k=1,jres-kc
           cm(i)= cm(i)+v1(i,k)*dn(i,k)/w1(k)
         enddo
      enddo

c     top1: summation of experimental B-factors
c     top2: summation of cm(k) from ANM calculations

      top1=0.
      top2=0.
      do k=1,jres
         top2=top2+cm(k)
      enddo

      do k=1,rno
         top1=top1+bfac(k)
      enddo

c     Rescaled ms fluctuations according experimental B-factors
c     Calculation of gamma as a result of scaling

      fact=top1/top2
      do k=1,jres
         cm(k)=cm(k)*fact
      enddo

      c1=8.*3.1415927**2.
      gamma=0.6*c1/(top1/top2)/3

c     Spring constant gamma changes for each residue pair. The following
c     calculation gives the average gamma in the elastic network.

      write(31,*)
      write(31,*)'top1,top2',top1,top2
      write(31,*)'top1/top2',(top1/top2)
      write(31,*)
      write(31,*)'gamma in (kcal/mol.A2)',gamma

      close(31)


c     Array of rescaled ms fluctuations for all residues
c     ds is r2, calculated as x2+y2+z2

      do i=1,rno
         k=(i-1)*3+1
         ds(i)=cm(k)+cm(k+1)+cm(k+2)
      enddo


c     smo(i): ms fluctuation smoothed over first neighbors of i

      smo(1)=(ds(1)+ds(2))/2.
      smo(rno)=(ds(rno)+ds(rno-1))/2.

      do i=2,rno-1
         smo(i)=(ds(i-1)+ds(i)+ds(i+1))/3.
      enddo
         write(32,*)'     index     Bf_exp    Bf_theo    Bf_theo_smooth'
	iskip=1
      do i=1,rno
         indexi=(i-1)*iskip+1
         write(32,*)indexi,bfac(i),ds(i),smo(i)
      enddo

      close(32)
      


c---------------------------------------------------------------------
c     Calculating ms fluctuations of slowest 20(=npar3) modes

c     cmod(i,j): ms fluctuation of residue i in slow mode j (=jres-kc-j+1)
c     csum(i,j): cumulative ms fluct of residue i in the slowest j modes

      write(49,*)'TOTAL NUMBER OF MODES= ',jres
      write(49,*)'Number of zero eigenvalues= ',kc
      write(49,*)
      write(49,*)'slowest_mode no: ',npar3

      ka1=jres-kc
      do j=1,npar3
         do i=1,rno
            k=(i-1)*3+1
            cmod(i,j)=0.
            csum(i,j)=0.
            cmod(i,j)= v1(k,ka1)*dn(k,ka1)/w1(ka1)+
     :            v1(k+1,ka1)*dn(k+1,ka1)/w1(ka1)+
     :            v1(k+2,ka1)*dn(k+2,ka1)/w1(ka1)
            if(i.eq.1)then
               csum(i,j)=cmod(i,j)
            else
               csum(i,j)=csum(i,j-1)+cmod(i,j)
            endif
         enddo
         write(49,*)'ka1 ',ka1
         ka1=ka1-1
      enddo


      do i=1,rno
         indexi=(i-1)*iskip+1
         write(33,58)indexi,(cmod(i,k),k=1,npar3)
         write(34,58)indexi,(csum(i,k),k=1,npar3)
      enddo

 58   Format(I4,(20f10.7,1x))

c---------------------------------------------------------------------
c---------------------------------------------------------------------
c     Calculating ms fluctuations of fastest 20(=npar3) modes

c     cmod(i,j): ms fluctuation of residue i in fast mode j (=jres-kc-j+1)
c     csum(i,j): cumulative ms fluct of residue i in the fast j modes

      write(49,*)
      write(49,*)'fastest_mode no: ',npar3
      do j=1,npar3
         ka1=j
         do i=1,rno
            k=(i-1)*3+1
            cmod(i,j)=0.
            csum(i,j)=0.
            cmod(i,j)= v1(k,ka1)*dn(k,ka1)/w1(ka1)+
     :            v1(k+1,ka1)*dn(k+1,ka1)/w1(ka1)+
     :            v1(k+2,ka1)*dn(k+2,ka1)/w1(ka1)
            if(i.eq.1)then
               csum(i,j)=cmod(i,j)
            else
               csum(i,j)=csum(i,j-1)+cmod(i,j)
            endif
         enddo
         write(49,*)'ka1 ',ka1
      enddo


      do i=1,rno
         indexi=(i-1)*iskip+1
         write(35,58)indexi,(csum(i,k),k=1,npar3)
      enddo


c----------------------------------------------------------------------
c     DISTORTED COORDINATES (BY FLUCTUATIONS) in PDB format
c     for 3-D representation by MIDAS

c     CALCULATION FOR INDIVIDUAL SLOWEST MODES 
c     starting with the slowest mode ka1=jres-kc


      c1=3./c1

      ka1=jres-kc
      ifile1=61


      do i=1,4

c        the components of the eigenvectors are rescaled by:
         rescale=(c1*top1/top2/w1(ka1))**0.5 

         write(49,*)
         write(49,*)'dist_coord: ka1,i',ka1,i
         write(49,*)'file',ifile1,ifile2
         write(49,*)'rescale',rescale

c        x component's index is 3j-2, y's is 3j-1, z's is 3j

         do j=1,rno
            ix=3*j-2
            iy=3*j-1
            iz=3*j
            xdist1(j)=xxnew(j)+v1(ix,ka1)*rescale
            ydist1(j)=yynew(j)+v1(iy,ka1)*rescale
            zdist1(j)=zznew(j)+v1(iz,ka1)*rescale
            xdist2(j)=xxnew(j)-v1(ix,ka1)*rescale
            ydist2(j)=yynew(j)-v1(iy,ka1)*rescale
            zdist2(j)=zznew(j)-v1(iz,ka1)*rescale

         enddo

         do j=1,rno
            write(ifile1,55)dummy,rcg(j),rrestyp(j),
     :          j,xdist1(j),ydist1(j),zdist1(j),1.00,smo(j)

         enddo

         close(ifile1)
     
         ifile1=ifile1+1
    

     
         do j=1,rno
            write(ifile1,55)dummy,rcg(j),rrestyp(j),
     :          j,xdist2(j),ydist2(j),zdist2(j),1.00,smo(j)
   
         enddo

         close(ifile1)
      
         ifile1=ifile1+1
      

         ka1=ka1-1

      enddo

c     First 100 modes with dispacement vectors

90     FORMAT(I6,1x,100(f8.5,1x))

	do j=1,rno
	ix=3*j-2
	iy=3*j-1
	iz=3*j
	write(110,90)j,(v1(ix,ka1),ka1=jres-kc,jres-kc-99,-1)
	write(110,90)j,(v1(iy,ka1),ka1=jres-kc,jres-kc-99,-1)
	write(110,90)j,(v1(iz,ka1),ka1=jres-kc,jres-kc-99,-1)

	enddo	



      END


c*************************************************************************
c     SUBROUTINE FOR FINDING THE PRINCIPAL AXES 
c*************************************************************************


      SUBROUTINE svdcmp(ag,m,n,mp,np,w,v)

C     This subroutine uses function: pythag

      INTEGER m,mp,n,np,NMAX,npar1
      PARAMETER (npar1=4000)
      PARAMETER (NMAX=1000)

      INTEGER i,its,j,jj,k,l,nm
      REAL ag(npar1,3),v(3,3),w(3)
      REAL anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX),pythag

c     here dim2=m,mp, and dim1=n,np
      g=0.0
      scale=0.0
      anorm=0.0

      do 25 i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0
        s=0.0
        scale=0.0
        if(i.le.m)then
          do 11 k=i,m
            xx7=ag(k,i)
            scale=scale+abs(xx7)
11        continue
          if(scale.ne.0.0)then
            do 12 k=i,m
              ag(k,i)=ag(k,i)/scale
              s=s+ag(k,i)*ag(k,i)
12          continue
            f=ag(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            ag(i,i)=f-g
            do 15 j=l,n
              s=0.0
              do 13 k=i,m
                s=s+ag(k,i)*ag(k,j)
13            continue
              f=s/h
              do 14 k=i,m
                ag(k,j)=ag(k,j)+f*ag(k,i)
14            continue
15          continue
            do 16 k=i,m
              ag(k,i)=scale*ag(k,i)
16          continue
          endif
        endif
        w(i)=scale *g
        g=0.0
        s=0.0
        scale=0.0
        if((i.le.m).and.(i.ne.n))then
          do 17 k=l,n
            scale=scale+abs(ag(i,k))
17        continue
          if(scale.ne.0.0)then
            do 18 k=l,n
              ag(i,k)=ag(i,k)/scale
              s=s+ag(i,k)*ag(i,k)
18          continue
            f=ag(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            ag(i,l)=f-g
            do 19 k=l,n
              rv1(k)=ag(i,k)/h
19          continue
            do 23 j=l,m
              s=0.0
              do 21 k=l,n
                s=s+ag(j,k)*ag(i,k)
21            continue
              do 22 k=l,n
                ag(j,k)=ag(j,k)+s*rv1(k)
22            continue
23          continue
            do 24 k=l,n
              ag(i,k)=scale*ag(i,k)
24          continue
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
25    continue

      do 32 i=n,1,-1
        if(i.lt.n)then
          if(g.ne.0.0)then
            do 26 j=l,n
              v(j,i)=(ag(i,j)/ag(i,l))/g
26          continue
            do 29 j=l,n
              s=0.0
              do 27 k=l,n
                s=s+ag(i,k)*v(k,j)
27            continue
              do 28 k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
28            continue
29          continue
          endif
          do 31 j=l,n
            v(i,j)=0.0
            v(j,i)=0.0
31        continue
        endif
        v(i,i)=1.0
        g=rv1(i)
        l=i
32    continue

      do 39 i=min(m,n),1,-1
        l=i+1
        g=w(i)
        do 33 j=l,n
          ag(i,j)=0.0
33      continue
        if(g.ne.0.0)then
          g=1.0/g
          do 36 j=l,n
            s=0.0
            do 34 k=l,m
              s=s+ag(k,i)*ag(k,j)
34          continue
            f=(s/ag(i,i))*g
            do 35 k=i,m
              ag(k,j)=ag(k,j)+f*ag(k,i)
35          continue
36        continue
          do 37 j=i,m
            ag(j,i)=ag(j,i)*g
37        continue
        else
          do 38 j= i,m
            ag(j,i)=0.0
38        continue
        endif
        ag(i,i)=ag(i,i)+1.0
39    continue

      do 49 k=n,1,-1
        do 48 its=1,30
          do 41 l=k,1,-1
            nm=l-1
            if((abs(rv1(l))+anorm).eq.anorm)  goto 2
            if((abs(w(nm))+anorm).eq.anorm)  goto 1
41        continue
1         c=0.0
          s=1.0
          do 43 i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm).eq.anorm) goto 2
            g=w(i)
            h=pythag(f,g)
            w(i)=h
            h=1.0/h
            c= (g*h)
            s=-(f*h)
            do 42 j=1,m
              y=ag(j,nm)
              z=ag(j,i)
              ag(j,nm)=(y*c)+(z*s)
              ag(j,i)=-(y*s)+(z*c)
42          continue
43        continue
2         z=w(k)
          if(l.eq.k)then
            if(z.lt.0.0)then
              w(k)=-z
              do 44 j=1,n
                v(j,k)=-v(j,k)
44            continue
            endif
            goto 3
          endif
          if(its.eq.30) pause 'no convergence in svdcmp'
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
          g=pythag(f,1.0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0
          s=1.0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=pythag(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
45          continue
            z=pythag(f,h)
            w(j)=z
            if(z.ne.0.0)then
              z=1.0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do 46 jj=1,m
              y=ag(jj,j)
              z=ag(jj,i)
              ag(jj,j)= (y*c)+(z*s)
              ag(jj,i)=-(y*s)+(z*c)
46          continue
47        continue
          rv1(l)=0.0
          rv1(k)=f
          w(k)=x
48      continue
3       continue
49    continue

      return
      END

c***********************************************************************

c***********************************************************************
 
      FUNCTION pythag(a,b)

      REAL a,b,pythag
      REAL absa,absb

      absa=abs(a)
      absb=abs(b)

      if(absa.gt.absb)then
        pythag=absa*sqrt(1.+(absb/absa)**2)
      else
        if(absb.eq.0.)then
          pythag=0.
        else
          pythag=absb*sqrt(1.+(absa/absb)**2)
        endif
      endif

      return
      END

c***********************************************************************
c     SUBROUTINE FOR INVERTING THE HESSIAN MATRIX
c***********************************************************************


      SUBROUTINE svdcmh(m,n,mp,np)

C     This subroutine uses function: pythag

      INTEGER m,mp,n,np,npar2
      PARAMETER (npar2=12000)

      INTEGER i,its,j,jj,k,l,nm
      REAL anorm,c,f,g,h,s,scale,x,y,z,rv1(npar2),pythag
     
      COMMON /test1/ w1(npar2),v1(npar2,npar2),dn(npar2,npar2)

      g=0.0
      scale=0.0
      anorm=0.0

      do 25 i=1,n
	write(6,*)'loop 25 i=',i
        l=i+1
        rv1(i)=scale*g
        g=0.0
        s=0.0
        scale=0.0
        if(i.le.m)then
          do 11 k=i,m
c	write(6,*)'loop 11 k=',k
	xx7=dn(k,i)
            scale=scale+abs(xx7)
11        continue
          if(scale.ne.0.0)then
            do 12 k=i,m
              dn(k,i)=dn(k,i)/scale
              s=s+dn(k,i)*dn(k,i)
12          continue
            f=dn(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            dn(i,i)=f-g
            do 15 j=l,n
              s=0.0
              do 13 k=i,m
                s=s+dn(k,i)*dn(k,j)
13            continue
              f=s/h
              do 14 k=i,m
                dn(k,j)=dn(k,j)+f*dn(k,i)
14            continue
15          continue
            do 16 k=i,m
              dn(k,i)=scale*dn(k,i)
16          continue
          endif
        endif
        w1(i)=scale *g
        g=0.0
        s=0.0
        scale=0.0
        if((i.le.m).and.(i.ne.n))then
          do 17 k=l,n
            scale=scale+abs(dn(i,k))
17        continue
          if(scale.ne.0.0)then
            do 18 k=l,n
              dn(i,k)=dn(i,k)/scale
              s=s+dn(i,k)*dn(i,k)
18          continue
            f=dn(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            dn(i,l)=f-g
            do 19 k=l,n
              rv1(k)=dn(i,k)/h
19          continue
            do 23 j=l,m
              s=0.0
              do 21 k=l,n
                s=s+dn(j,k)*dn(i,k)
21            continue
              do 22 k=l,n
                dn(j,k)=dn(j,k)+s*rv1(k)
22            continue
23          continue
            do 24 k=l,n
              dn(i,k)=scale*dn(i,k)
24          continue
          endif
        endif
        anorm=max(anorm,(abs(w1(i))+abs(rv1(i))))
25    continue

      do 32 i=n,1,-1
	write(6,*)'loop 32 i=',i
        if(i.lt.n)then
          if(g.ne.0.0)then
            do 26 j=l,n
              v1(j,i)=(dn(i,j)/dn(i,l))/g
26          continue
            do 29 j=l,n
              s=0.0
              do 27 k=l,n
                s=s+dn(i,k)*v1(k,j)
27            continue
              do 28 k=l,n
                v1(k,j)=v1(k,j)+s*v1(k,i)
28            continue
29          continue
          endif
          do 31 j=l,n
            v1(i,j)=0.0
            v1(j,i)=0.0
31        continue
        endif
        v1(i,i)=1.0
        g=rv1(i)
        l=i
32    continue

      do 39 i=min(m,n),1,-1
        write(6,*)'loop 39  i=',i
        l=i+1
        g=w1(i)
        do 33 j=l,n
          dn(i,j)=0.0
33      continue
        if(g.ne.0.0)then
          g=1.0/g
          do 36 j=l,n
            s=0.0
            do 34 k=l,m
              s=s+dn(k,i)*dn(k,j)
34          continue
            f=(s/dn(i,i))*g
            do 35 k=i,m
              dn(k,j)=dn(k,j)+f*dn(k,i)
35          continue
36        continue
          do 37 j=i,m
            dn(j,i)=dn(j,i)*g
37        continue
        else
          do 38 j= i,m
            dn(j,i)=0.0
38        continue
        endif
        dn(i,i)=dn(i,i)+1.0
39    continue

      do 49 k=n,1,-1
	write(6,*)'loop 49 k=',k
        do 48 its=1,30
          do 41 l=k,1,-1
            nm=l-1
            if((abs(rv1(l))+anorm).eq.anorm)  goto 2
            if((abs(w1(nm))+anorm).eq.anorm)  goto 1
41        continue
1         c=0.0
          s=1.0
          do 43 i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm).eq.anorm) goto 2
            g=w1(i)
            h=pythag(f,g)
            w1(i)=h
            h=1.0/h
            c= (g*h)
            s=-(f*h)
            do 42 j=1,m
              y=dn(j,nm)
              z=dn(j,i)
              dn(j,nm)=(y*c)+(z*s)
              dn(j,i)=-(y*s)+(z*c)
42          continue
43        continue
2         z=w1(k)
          if(l.eq.k)then
            if(z.lt.0.0)then
              w1(k)=-z
              do 44 j=1,n
                v1(j,k)=-v1(j,k)
44            continue
            endif
            goto 3
          endif
          if(its.eq.30) pause 'no convergence in svdcmp'
          x=w1(l)
          nm=k-1
          y=w1(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
          g=pythag(f,1.0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0
          s=1.0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w1(i)
            h=s*g
            g=c*g
            z=pythag(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 jj=1,n
              x=v1(jj,j)
              z=v1(jj,i)
              v1(jj,j)= (x*c)+(z*s)
              v1(jj,i)=-(x*s)+(z*c)
45          continue
            z=pythag(f,h)
            w1(j)=z
            if(z.ne.0.0)then
              z=1.0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do 46 jj=1,m
              y=dn(jj,j)
              z=dn(jj,i)
              dn(jj,j)= (y*c)+(z*s)
              dn(jj,i)=-(y*s)+(z*c)
46          continue
47        continue
          rv1(l)=0.0
          rv1(k)=f
          w1(k)=x
48      continue
3       continue
49    continue

      return
      END

