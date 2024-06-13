
!                        /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
!!!!!!!!!!!!!!!!!!!!!    COMPUTATIONAL MODEL FOR GENERATING TROPHIC NETWORKS     !!!!!!!!!!!!!!!!!!!!!!!
!                        \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

!    Accompanying the paper entitled "AUTOMATED GENERATION OF REAL-WORLD ECOLOGICAL NETWORKS AND ITS POTENTIAL APPLICATION IN CONSERVATION BIOLOGY" 

!    By: 
!    5 Authors
!    
!    This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or any later version.

!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details: 
!    For more details, see <http://www.gnu.org/licenses/>.

!    This program is written in Fortran language. This is the simplest version of the code. 
!    INPUT: completed species list listed in the filelist.txt fle
!    OUTPUTS: binary trophic networks (one per replicate), node lists (one per ecosystem), topological descriptors (one in total)

!    To compile: Download and install a fortran compiler: In Linux; type in terminal: "sudo apt-get install gfortran"  
!    (give permissions to files and executables if necessary).
!    Go to the folder where this program is contained (*.f90 file + filelist.txt + *.csv files)
!    Then type "f95 Networker.f90 -fbounds-check -o networker". To launch it, type: "./Networker". 
!    Output files will be created automatically in the same folder.
!    Output NODEK*.csv and EDGEK*.csv datafiles can be plotted using the external graphical software Gephi: <https://www.gephi.org>

program networker

implicit none

integer :: i,j,k,ii,jj,kk,n,ios,ich1,ich2,ret,replica,replicas          ! counters
integer :: nspp,come,come2,iii,jjj,kkk,sp1,sp2,im,edg,nod,force,forcemax,whois(300),dishr,sizhr
integer,allocatable :: NWS(:,:),NWSniche(:,:)                           ! network of species and genera & diets(types 1-5) of gen and species
real*4,allocatable  :: dis(:,:)                                         ! dietary vectors
real*4  :: a,b,c,w,x,y,z,xx,yy,zz,zzz,logx,logy,oki,toki,niche,noniche         ! topological descriptors
real*4  :: sizmaxs,sizmins,avgx,avgn,R,Rniche,sptot,sueltosf,sueltosc,sueltost ! topological descriptors
real*4  :: L1,L2,CC,CCC,CCCC,L1niche,L2niche,Cniche,CCniche,CCCniche,CCCCniche ! topological descriptors
real,parameter     :: pi=3.1415926535, e=2.7182818284, delta=1d-6       ! some useful constants
real*4,allocatable :: siz(:),ELES(:,:)                                  ! sizes of gen and species
character(len=50)  :: attr(6),node(2000)                                ! taxonomic attributes
character(len=150) :: diet,habitat,dietS(2000)                          ! diet and habitat
character(len=40)  :: arx1                                              ! character variables for file managing 
character(len=60)  :: arx2,arx3                                         ! character variables for file managing  
character(len=16)  :: arxaux                                            ! character variables for file managing 
character(len=1000):: arxi,arxu                                         
character(len=1000),allocatable :: arximax(:)
integer            :: matrixdiet(7,7),pred(2000,8),prey(2000),problematico(3000)          ! cateforias generales y especificas
integer,allocatable :: matrix(:,:),matrixniche(:,:)
real*4 :: AVL,CON,IND,OUTD,ASYM                      ! average links per species, connectance, average indegree and outdegree
real*4 :: TOP,INTS,BAS,CAN                           ! fraction of spp that are top, intermediate, basal (without resources) and cannibals.
real*4 :: MINL,MAXL,AVGL                             ! minimum, maximum and average (Cohen)chain lenght  Average (Cohen)chain lenght 
integer :: NCHAIN,NLOOPS,morpho                      ! minimum and maximum (Cohen)chain lenght, number of trophic chains and loops
real*4 :: AVLniche,CONniche,INDniche,OUTDniche,ASYMniche ! average links per species, connectance, average indegree and outdegree
real*4 :: TOPniche,INTSniche,BASniche,CANniche       ! fraction of spp that are top, intermediate, basal (without resources) and cannibals.
real*4 :: MINLniche,MAXLniche,AVGLniche              ! minimum, maximum and average (Cohen)chain lenght  Average (Cohen)chain lenght 
real*4 :: DD1,DD2,DD3                                ! degree distribution
real*4 :: degree(1000,3),degreeniche(1000,3)         ! degree distribution
integer :: NCHAINniche,NLOOPSniche                   ! minimum and maximum (Cohen)chain lenght, number of trophic chains and loops

 integer, allocatable :: seed(:)                     ! setting the seed for random number generator (not accessible)
 integer sizer                                       ! this way all replicates give the same result
 call random_seed(size=sizer)                        !
 allocate(seed(sizer))                               !
 call random_seed(put=seed)                          !

                   ! DIETARY MATRIX (WHO EATS WHOM) ......
                   !A H C D S P P    ! 0 cannot eat  
                   !U E A E C A A    ! 1 eats, size does not matter
                   !T R R T A 1 2    ! 2 eats, size matters, 3 inverse niche model applies (for animal-animal aprasites)
matrixdiet(1,1:7)=(/0,1,0,1,0,1,0/)  ! eats AUTotroph
matrixdiet(2,1:7)=(/0,0,2,0,1,0,1/)  ! eats HERbivore
matrixdiet(3,1:7)=(/0,0,2,0,1,0,1/)  ! eats CARnivore
matrixdiet(4,1:7)=(/0,0,2,0,1,0,1/)  ! eats DETritivore
matrixdiet(5,1:7)=(/0,0,2,0,1,0,1/)  ! eats SCAvenger
matrixdiet(6,1:7)=(/0,0,0,0,0,0,0/)  ! PA1-Rasites plants      
matrixdiet(7,1:7)=(/0,0,3,0,1,0,1/)  ! PA2-Rasites animals

forcemax=20                                                             ! maximum trials before brute-force to link all 
open(90,file='filelist.txt'   , action='read'  ,iostat=ios)             ! list of files to open
open(77,file='morphospace.dat', action='write' ,iostat=ios)             ! output file wih topological descriptors
morpho=1 ; niche=0.0 ; noniche=0.0 ; sptot=0.0

do i=1,3                                                                ! number of datasets in filelist.txt MANUAL CHECK
  read(90,*)arx1 
  open(88,file=arx1, action='read'  ,iostat=ios)                        ! file to open C1
  nspp=0                                                                ! setting the counter
  do while(ios.eq.0)                                                    ! counting the number of SPECIES
    read(88,*,iostat=ios)attr(1:2) ; nspp=nspp+1                        ! counting the number of SPECIES
  end do ; nspp=nspp-1 ; nod=1 ; arxi=' ' ; arxu=' '
  rewind(88)
  write(*,*)'    '
  write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,*)'INPUTFILE',i,arx1,nspp
  arx2=' ';arx3=' ';arx2(1:40)=arx1(1:40);arx3(1:40)=arx1(1:40)! NETWORK FILES  
  do j=1,40-3                                                  ! NETWORK FILES
    if(arx1(j:j+3).eq.'.csv')then                              ! NETWORK FILES
      arx2(j:j+9)='_NODEK.csv'; arx3(j:j+15)='_EDGEK_R____.csv'! NETWORK FILES
    end if                                    !123456789012345 ! NETWORK FILES
  end do                      ; arxaux(1:16)='_EDGEK_R____.csv'! NETWORK FILES   
  open(101,file=arx2 , action='write',iostat=ios)              ! NETWORK FILES: NODES
  write(101,"(A52)")'Id,Label,Kingdom,Phylum,Class,Order,Family,Diet,Size'
                    !123456789012345678901234567890123456789012345678901234567890123456789012345678901234
            
  !!!!!!!!!!!!!!!!!!!!!!!!! READ THE DATA !!!!!!!!!!!!!!!!!!!!!!
  allocate(siz(nspp)) ; allocate(dis(nspp,7))                           ! 6 dietary categories (coarse graining)
  allocate(arximax(nspp))
  allocate(eles(nspp,2)) ; eles=0.0
  dietS=' ' ; siz=0.0 ; sizmaxs=0.0 ; sizmins=1e10 ; dis=0.0 ; node=' '
  pred=0 ; prey=0
  do j=1,nspp
    x=0.0 ; y=0.0 ; attr=' ' ; diet=' ' ; habitat=' ' ; ich1=0 ; ich2=0 
    read(88,*,iostat=ios)attr(1:6),x,y,diet,habitat                     ! 123456789012345678901234567890
    ich1=iachar(diet(1:1)) ; ich2=iachar(habitat(1:1)) ; toki=toki+1.0  ! handling missing fields
      if(habitat(1:3).eq.'MOR')then          ; goto 15 ; end if                   
      if((ich1.ne.32).and.(ich1.ge.97).and.(ich2.le.90))then
      habitat=' '            ; backspace(88) ; goto 15 ; end if      
      if((ich1.eq.32).and.(ich2.le.90))then
      habitat=' ' ; diet=' ' ; backspace(88) ; goto 15 ; end if
    !end if   
    15 siz(j)=log(x) 
    if((x.eq.0.0).and.(attr(2)(1:6).eq.'Animal'))then 
      call random_number(x) ; siz(j)=log(x+1d-4)                        ! if some organisms still have missing data fter ALGORITHM 1
      write(*,*)'asignosdize',siz(j),attr(1)(1:30)
    end if                                                              ! because Log(0)=Nan, but animalas must have some size   
    if(diet(1:len(diet)).eq.' ')then                                    ! random diet assignation for missing data
      if(attr(2)(1:6).eq.'Animal')then   ; call getdiet(attr(3:6),diet(1:len(diet))) 
      else ; diet(1:14)= 'photoautotroph'; end if 
    end if  
    dietS(j)=diet     
    if((siz(j).gt.sizmaxs).and.(attr(2)(1:6).eq.'Animal'))then          ! for the normalization [0,1]
    sizmaxs=siz(j); end if                                              ! for the normalization [0,1]
    if((siz(j).lt.sizmins).and.(x.gt.0.0).and.(attr(2)(1:6).ne.'Planta')&
    .and.(attr(2)(1:5).ne.'Fungi').and.(attr(2)(1:6).ne.'Bacter'))then 
    sizmins=siz(j); end if                                              ! for the normalization [0,1]
    if((attr(2)(1:4).eq.'Plan').or.(attr(2)(1:4).eq.'Fung').or.(attr(2)(1:4).eq.'Bact').or.(attr(2)(1:4).eq.'Chro'))then 
    siz(j)=0.0 ; end if    
    node(j)(1:50)=attr(1)(1:50)
   
    dis(j,1:7)=0.0
    do k=1,len(diet)-8   !!! CATEGORIC COARSE GRAINING !!! CATEGORIC COARSE GRAINING !!! CATEGORIC COARSE GRAINING !!! CATEGORIC COARSE GRAINING 
      if(diet(k:k+8).eq.'autotroph')then ;  dis(j,1)=dis(j,1)+1.0 ; end if   ! CATEGORY 1: AUTOTROPH (QUIMIO OR PHOTO)
      if(diet(k:k+7).eq.'symbiotr' )then ;  dis(j,1)=dis(j,1)+1.0 ; end if   ! CATEGORY 1: AUTOTROPH (QUIMIO OR PHOTO)
                        !1234567890                      1234567890                    1234567890
      if(diet(k:k+8).eq.'herbivore')then ;  dis(j,2)=dis(j,2)+1.0 ; end if   ! CATEGORY 2: HERBIVORE
      if(diet(k:k+4).eq.'fruit'    )then ;  dis(j,2)=dis(j,2)+1.0 ; end if   ! CATEGORY 2: HERBIVORE
      if(diet(k:k+5).eq.'leaves'   )then ;  dis(j,2)=dis(j,2)+1.0 ; end if   ! CATEGORY 2: HERBIVORE     
      if(diet(k:k+4).eq.'plant'    )then ;  dis(j,2)=dis(j,2)+1.0 ; end if   ! CATEGORY 2: HERBIVORE
      if(diet(k:k+3).eq.'seed'     )then ;  dis(j,2)=dis(j,2)+1.0 ; end if   ! CATEGORY 2: HERBIVORE
      if(diet(k:k+3).eq.'alga'     )then ;  dis(j,2)=dis(j,2)+1.0 ; end if   ! CATEGORY 2: HERBIVORE
      if(diet(k:k+3).eq.'leaf'     )then ;  dis(j,2)=dis(j,2)+1.0 ; end if   ! CATEGORY 2: HERBIVORE
      if(diet(k:k+3).eq.'wood'     )then ;  dis(j,2)=dis(j,2)+1.0 ; end if   ! CATEGORY 2: HERBIVORE   
      if(diet(k:k+3).eq.'root'     )then ;  dis(j,2)=dis(j,2)+1.0 ; end if   ! CATEGORY 2: HERBIVORE
      if(diet(k:k+5).eq.'flower'   )then ;  dis(j,2)=dis(j,2)+1.0 ; end if   ! CATEGORY 2: HERBIVORE   
      if(diet(k:k+5).eq.'pollen'   )then ;  dis(j,2)=dis(j,2)+1.0 ; end if   ! CATEGORY 2: HERBIVORE 
      if(diet(k:k+5).eq.'nectar'   )then ;  dis(j,2)=dis(j,2)+1.0 ; end if   ! CATEGORY 2: HERBIVORE
      if(diet(k:k+3).eq.'weed'     )then ;  dis(j,2)=dis(j,2)+1.0 ; end if   ! CATEGORY 2: HERBIVORE
      if(diet(k:k+5).eq.'lichen'   )then ;  dis(j,2)=dis(j,2)+1.0 ; end if   ! CATEGORY 2: HERBIVORE
      if(diet(k:k+5).eq.'nectar'   )then ;  dis(j,2)=dis(j,2)+1.0 ; end if   ! CATEGORY 2: HERBIVORE
      if(diet(k:k+3).eq.'fung'     )then ;  dis(j,2)=dis(j,2)+1.0 ; end if   ! CATEGORY 2: HERBIVORE
      if(diet(k:k+4).eq.'omniv'    )then ;  dis(j,2)=dis(j,2)+1.0 ; end if   ! CATEGORY 2: HERBIVORE
      if(diet(k:k+4).eq.'grass'    )then ;  dis(j,2)=dis(j,2)+1.0 ; end if   ! CATEGORY 2: HERBIVORE
      if(diet(k:k+4).eq.'phyto'    )then ;  dis(j,2)=dis(j,2)+1.0 ; end if   ! CATEGORY 2: HERBIVORE
      if(diet(k:k+5).eq.'microo'   )then ;  dis(j,2)=dis(j,2)+1.0 ; end if   ! CATEGORY 2: HERBIVORE
      if(diet(k:k+4).eq.'grass'    )then ;  dis(j,2)=dis(j,2)+1.0 ; end if   ! CATEGORY 2: HERBIVORE
                         !1234567890                      1234567890                    1234567890
      if(diet(k:k+8).eq.'carnivore')then ;  dis(j,3)=dis(j,3)+1.0 ; 
        if(attr(2)(1:6).eq.'Planta')then  ; dis(j,1)=1            ; end if   ! correction for insectivore plants  ...
                                                                    end if   ! CATEGORY 3: CARNIVORE
      if(diet(k:k+2).eq.'egg'      )then ;  dis(j,3)=dis(j,3)+1.0 ; end if   ! CATEGORY 3: CARNIVORE
      if(diet(k:k+6).eq.'vertebr'  )then ;  dis(j,3)=dis(j,3)+1.0 ; end if   ! CATEGORY 3: CARNIVORE
      if(diet(k:k+3).eq.'worm'     )then ;  dis(j,3)=dis(j,3)+1.0 ; end if   ! CATEGORY 3: CARNIVORE      
      if(diet(k:k+5).eq.'insect'   )then ;  dis(j,3)=dis(j,3)+1.0 ; end if   ! CATEGORY 3: CARNIVORE  
      if(diet(k:k+3).eq.'fish'     )then ;  dis(j,3)=dis(j,3)+1.0 ; end if   ! CATEGORY 3: CARNIVORE
      if(diet(k:k+8).eq.'crustacea')then ;  dis(j,3)=dis(j,3)+1.0 ; end if   ! CATEGORY 3: CARNIVORE 
      if(diet(k:k+5).eq.'echino'   )then ;  dis(j,3)=dis(j,3)+1.0 ; end if   ! CATEGORY 3: CARNIVORE  
      if(diet(k:k+6).eq.'amphibi'  )then ;  dis(j,3)=dis(j,3)+1.0 ; end if   ! CATEGORY 3: CARNIVORE
      if(diet(k:k+6).eq.'mollusc'  )then ;  dis(j,3)=dis(j,3)+1.0 ; end if   ! CATEGORY 3: CARNIVORE   
      if(diet(k:k+5).eq.'mammal'   )then ;  dis(j,3)=dis(j,3)+1.0 ; end if   ! CATEGORY 3: CARNIVORE  
      if(diet(k:k+3).eq.'bird'     )then ;  dis(j,3)=dis(j,3)+1.0 ; end if   ! CATEGORY 3: CARNIVORE
      if(diet(k:k+6).eq.'reptile'  )then ;  dis(j,3)=dis(j,3)+1.0 ; end if   ! CATEGORY 3: CARNIVORE   
      if(diet(k:k+6).eq.'arthrop'  )then ;  dis(j,3)=dis(j,3)+1.0 ; end if   ! CATEGORY 3: CARNIVORE 
      if(diet(k:k+2).eq.'zoo'      )then ;  dis(j,3)=dis(j,3)+1.0 ; end if   ! CATEGORY 3: CARNIVORE
      if(diet(k:k+4).eq.'omniv'    )then ;  dis(j,3)=dis(j,3)+1.0 ; end if   ! CATEGORY 3: CARNIVORE     
      if(diet(k:k+6).eq.'predato'  )then ;  dis(j,3)=dis(j,3)+1.0 ; end if   ! CATEGORY 3: CARNIVORE 
      if(diet(k:k+5).eq.'microo'   )then ;  dis(j,3)=dis(j,3)+1.0 ; end if   ! CATEGORY 3: CARNIVORE
                         !1234567890                      1234567890                    1234567890
      if(diet(k:k+5).eq.'detrit'   )then ;  dis(j,4)=dis(j,4)+1.0  ; end if  ! CATEGORY 4: DETRITIVORE    
      if(diet(k:k+6).eq.'sedimen'  )then ;  dis(j,4)=dis(j,4)+1.0  ; end if  ! CATEGORY 4: DETRITIVORE     
      if(diet(k:k+2).eq.'mor'      )then ;  dis(j,4)=dis(j,4)+1.0  ; end if  ! CATEGORY 4: DETRITIVORE        
      if(diet(k:k+5).eq.'litter'   )then ;  dis(j,4)=dis(j,4)+1.0  ; end if  ! CATEGORY 4: DETRITIVORE        
      if(diet(k:k+3).eq.'clay'     )then ;  dis(j,4)=dis(j,4)+1.0  ; end if  ! CATEGORY 4: DETRITIVORE        
      if(diet(k:k+6).eq.'suspend'  )then ;  dis(j,4)=dis(j,4)+1.0  ; end if  ! CATEGORY 4: DETRITIVORE    
      if (diet(k:k+3).eq.'filt'    )then ;  dis(j,4)=dis(j,4)+1.0  ; end if  ! CATEGORY 4: DETRITIVORE                                                               
                         !1234567890                      1234567890                    1234567890            
      if(diet(k:k+6).eq.'scaveng'  )then ;  dis(j,5)=dis(j,5)+1.0  ; end if  ! CATEGORY 5: SCAVENGER  
      if(diet(k:k+6 ).eq.'carrion' )then ;  dis(j,5)=dis(j,5)+1.0  ; end if  ! CATEGORY 5: SCAVENGER
      if((diet(k:k+6).eq.'parasit').and.(attr(2)(1:6).eq.'Animal'))then
      dis(j,6)=dis(j,6)+1.0 ; end if ! CATEGORY 6: PARASITE OF PLANTS (USUALLY OTHER PLANTS) 
      if((diet(k:k+6).eq.'parasit').and.(attr(2)(1:6).eq.'Planta'))then
      dis(j,7)=dis(j,7)+1.0 ; end if ! CATEGORY 7: PARASITE OF ANIMALS (ANIMALS, SOME FUNGI)    
      
      !! SPECIFIC DIETS CARNIVORES !!!!!!! SPECIFIC DIETS CARNIVORES !!!!!!! SPECIFIC DIETS CARNIVORES !!!!!!! SPECIFIC DIETS CARNIVORES !!!!!
      if((diet(k:k+5).eq.'insect') .or.(diet(k:k+6).eq.'arthrop').or.(diet(k:k+5).eq.'invert' ))then; pred(j,1)=1 ; end if ! INSECTIVORE
      if((diet(k:k+3).eq.'worm')   .or.(diet(k:k+6).eq.'mollusc').or.(diet(k:k+5).eq.'invert' ))then; pred(j,2)=1 ; end if ! MOLLUSCIVORE
      if((diet(k:k+6).eq.'crustac').or.(diet(k:k+5).eq.'echino') .or.(diet(k:k+5).eq.'invert' ))then; pred(j,3)=1 ; end if ! DUROPHAGOUS
      if( diet(k:k+3).eq.'fish')                                                                then; pred(j,4)=1 ; end if ! PISCIVORE
      if((diet(k:k+6).eq.'amphibi').or.(diet(k:k+7).eq.'|vertebr')                              )then;pred(j,5)=1 ; end if ! "FROG"VORE
      if((diet(k:k+6).eq.'reptile').or.(diet(k:k+7).eq.'|vertebr')                              )then;pred(j,6)=1 ; end if ! HERPETIVORE
      if((diet(k:k+2).eq.'egg')    .or.(diet(k:k+3).eq.'bird')   .or.(diet(k:k+7).eq.'|vertebr'))then;pred(j,7)=1 ; end if ! BIRD-VORE-OOPHAGOUS
      if((diet(k:k+5).eq.'mammal' ).or.(diet(k:k+7).eq.'|vertebr')                              )then;pred(j,8)=1 ; end if ! "MAMMALI"-VORE
      if((diet(k:k+6).eq.'carnivo').or.(diet(k:k+6).eq.'omnivor'))then                 ! unnecessary if there is already something more specific                                
      if(sum(pred(j,1:8)).eq.0)then                  ; pred(j,:)=1 ; end if ; end if   ! unnecessary if there is already something more specific  
      if((sum(pred(j,:)).eq.0).and.(dis(j,3).eq.1))then; pred(j,:)=1 ; end if          ! sometimes happens. to prevent spare nodes...
    end do ! of categorization "k"    
       
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!   TAXONOMICAL COARSE-GRAINING  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    if((attr(3)(1:7).eq.'Arthrop').and.(attr(4)(1:9).ne.'Crustacea'))then; prey(j)=1 ;end if ! type of prey: insect, spyders, and relatives
    if((attr(3)(1:7).eq.'Mollusc').or. (attr(3)(1:8).eq.'Annelida' ))then; prey(j)=2 ;end if ! type of prey: snail, mollusc, worm 
    if((attr(4)(1:7).eq.'Crustac' ).or.(attr(3)(1:5).eq.'Echino'   ))then; prey(j)=3 ;end if ! type of prey: crustacea, echinoderms, etc   
    
    if((attr(4)(1:8).eq.'Actinopt').or.(attr(4)(1:9).eq.'Elasmobra') &
   .or.(attr(4)(1:8).eq.'Sarcopte').or.(attr(4)(1:9).eq.'Cephalasp'))then; prey(j)=4 ;end if ! type of prey: fish
   
    if (attr(4)(1:8).eq.'Amphibia')                                  then; prey(j)=5 ;end if ! type of prey: frogs, tadpoles
    if (attr(4)(1:8).eq.'Reptilia')                                  then; prey(j)=6 ;end if ! type of prey: lizards, snakes, reptiles
    if (attr(4)(1:4).eq.'Aves')                                      then; prey(j)=7 ;end if ! type of prey: birds and nestlings
    if (attr(4)(1:8).eq.'Mammalia')                                  then; prey(j)=8 ;end if ! type of prey: mammals 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    y=maxval(dis(j,:))!sum(dis(j,:))                                    ! normalization and relativization of diets 
    dis(j,:)=dis(j,:)/y                                                 ! normalization and relativization of diets 
    dis(j,:)=(dis(j,:))**2                                              ! optional. It increases specialization 
    if((dis(j,1).eq.1).and.(dis(j,3).eq.1))then                                    ! correction for insectivore plants ... insects and frogs only
    dis(j,3)=1 ; siz(j)=sizmins ; pred(j,:)=0 ; pred(j,1)=1 ; pred(j,5)=1 ; end if ! correction for insectivore plants ...
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dishr=0 
    do k=1,5 ; if(dis(j,k).eq.1.0)then ; dishr=k ; end if ; end do      ! dishr=human readable diet (optional). Only considers preferential diet ...
    write(arxi,"(I4,A2,6(A50,A2),I1)")nod,',',attr(1),',',attr(2),&
    ',',attr(3),',',attr(4),',',attr(5),',',attr(6),',',dishr           !,',',sizhr
    write(*,*)'&& ',j,nspp,nod,',',attr(1)(1:15),',',attr(2)(1:15),',',attr(3)(1:15),'SIZZ',siz(j),'DIS',dis(j,:),'HR=',dishr
    call removespaces(arxi,arxu)                                        ! write the attributes of the secies in a *csv, Gephi-readable format
    arxu=trim(arxu) ; arximax(j)(1:1000)=arxu(1:1000)                   ! this will be used later ....
    nod=nod+1                                                           ! goes to the next node/species in the network  
  end do                                                                ! end of species loop (j)
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!! GO FOR NETWORKS !!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!! GO FOR NETWORKS !!!!!!!!!!!!!!!!!!!  
  !!!!!!!!!!!!!!!!!!!!!!!!! GO FOR NETWORKS !!!!!!!!!!!!!!!!!!!
  sizmins=sizmins*0.6666666                                             ! to avoid zero-sizes (it does not matter as it will be made relative)
  do j=1,size(siz) ; if(siz(j).ne.0.0)then ; siz(j)=siz(j)-sizmins           ; end if ; end do ! for the normalization [0,1] 
  do j=1,size(siz) ; if(siz(j).ne.0.0)then ; siz(j)=siz(j)/(sizmaxs-sizmins) ; end if ; end do ! for the normalization [0,1] 

  allocate(NWS(nspp,nspp)) ; allocate(NWSniche(nspp,nspp))              ! allocating networks ...
  16 CCC=0.12                                                           ! Directed/desired connectance (CC) ! Average from empirical networks (see SI)
  call random_number(xx) ; call random_number(yy)                       ! Box-Muller algorithm N(0,sdev)  
  zz=0.05*sqrt(-2*log(xx))*cos(2*pi*yy)                                 ! Box-Muller algorithm N(0,sdev)
  CCC=CCC+zz 
  if((CCC.lt.0.05).or.(CCC.gt.0.25))then ; goto 16 ; end if  
  !!!! write NODEK files   !!!! write NODEK files   !!!! write NODEK files   !!!! write NODEK files   !!!! write NODEK files 
  !!!! write NODEK files   !!!! write NODEK files   !!!! write NODEK files   !!!! write NODEK files   !!!! write NODEK files 
  !nod=1
  do j=1,nspp
    sizhr=0 ;     sizhr=int(siz(j)*9.0)       
    write(*,*)'lasizesh ',siz(j),sizhr
    do k=1,1000
      arxu=' '
      if(arximax(j)(k:k).eq.' ')then ; write(arximax(j)(k:k+1),"(A1,I1)")',',sizhr 
      arxu(1:k+1)=arximax(j)(1:k+1) 
      write(101,"(A1000)")arxu ; flush(101) 
      exit ; end if      
    end do
    write(*,"(A11,A20)")'arxu2 to10=',arxu(1:20)
  end do
 close(101)                                                             ! NODES files
  
  !!!!!  REPLICATES VARY (APPROXIMATELY) AS R=2exp ||Log2(N/2)||. 
  if        (nspp.le.32  )                    then ; replicas=32        ! number of network replicates 
    else if((nspp.gt.32  ).and.(nspp.le.128 ))then ; replicas=16        ! number of network replicates 
    else if((nspp.gt.128 ).and.(nspp.le.512 ))then ; replicas=8         ! number of network replicates 
    else if((nspp.gt.512 ).and.(nspp.le.1024))then ; replicas=4         ! number of network replicates 
    else if (nspp.gt.1024)                    then ; replicas=2         ! number of network replicates 
  end if
  
  !!!!
  deallocate(arximax)
  do replica=1,replicas                                                 ! NETWORK FILES: EDGES
    sptot=sptot+real(nspp)  ; force=0 ; problematico=0 ; eles=0.0
    allocate(matrix(nspp,nspp)) ; allocate(matrixniche(nspp,nspp)) 
     do j=1,60-3 ;       if(arx3(j:j+3).eq.'_EDG')then                  ! NETWORK FILES: EDGES
       write(arxaux,"(A9,I3,A4)")arxaux(1:9),replica,arxaux(13:16)      ! max 999 replicates
       do im=1,16 ; if (arxaux(im:im)==" ") arxaux(im:im)="0";end do    ! composing filename
       arx3(j:j+15)=arxaux(1:16) ; goto 77                              ! NETWORK FILES: EDGES
     end if                                                             ! NETWORK FILES: EDGES
     end do														        ! NETWORK FILES: EDGES
    77 open(102,file=arx3 , action='write',iostat=ios)                  ! NETWORK FILES: EDGES
    write(102,"(A42)")'Source Target Type Id Label Interval Weight'     ! first output line (GEPHI-readable format)
    NWS=0 ; NWSniche=0 ; edg=1 ; sueltosf=0.0 ; sueltosc=0.0 ; sueltost=0.0  ; matrix=0 ; matrixniche=0
    
    do j=1,nspp       ! J come a K (VERIFY)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!          ! correction for insectivore plants ...
        if((dis(j,1).eq.1).and.(dis(j,3).eq.1))then                                    ! correction for insectivore plants ...
        dis(j,3)=2 ; siz(j)=sizmins ; pred(j,:)=0 ; pred(j,1)=1 ; pred(j,5)=1 ; end if ! correction for insectivore plants ...
      !!!!!!!!!!!!!!!! calculating the trophic range of the species J according to the niche model. 
      36 CC=CCC/(real(nspp))**2                                         ! The 2nd parameter of the NM (1st one is nspp)
      CC=2.0*CC                                                         ! Expected value of the beta distribution
      a=1.0                                                             ! alfa parameter of the beta distribution is 1 (as in Rivas & Martinez)
      b=(a*(1.0-CC))/CC                                                 ! beta parameter of the beta distribution to match the mean of 2CC. 
      63 call random_number(x)                                          ! to get a random number X [0,1] from the beta distribution, this is the input
      y=1.0-(1.0-x)**b                                                  ! the output from the beta distribution (checked)
      R=abs(y*siz(j))                                                   ! trophic range of the niche model (must be, at most, equal to sizes(i), hence the normalization     
      call random_number(xx) ; Rniche=abs(y*xx*0.75) 
      R=R*0.75                                                          ! balancing     
      call random_number(z)                                             ! to get the centroid Ci
      
      x=R/2.0                                                           ! left limit for the interval used to calculate centroid U(x,y)
      if(siz(j).lt.1.0-x)then ; y=siz(j)                                ! right limit for the interval used to calculate centroid U(x,y)
      else ; y=1.0-x; end if                                            ! right limit for the interval used to calculate centroid U(x,y)      
      C=x+z*abs(y-x)                                                    ! trophic optimum
      L1=C-(R/2.0) ; L2=C+(R/2.0)                                       ! upper and lower limits to the trophic range
      eles(j,1)=L1 ; eles(j,2)=L2
      CCCC=(CCC/sqrt(real(nspp)))*10.0                                  ! just to adjust the non-allometric conections
      call random_number(zz)
      zz=1.0/(1.0+zz)                                        
      CCCC=cccc*zz 
      
      x=Rniche/2.0                                                      ! left limit for the interval used to calculate centroid U(x,y)
      if(xx.lt.1.0-x)then ; y=xx                                        ! right limit for the interval used to calculate centroid U(x,y)
      else ; y=1.0-x; end if                                            ! right limit for the interval used to calculate centroid U(x,y)      
      Cniche=x+z*abs(y-x)                                               ! trophic optimum
      L1niche=Cniche-(Rniche/2.0) ; L2niche=Cniche+(Rniche/2.0)         ! upper and lower limits to the trophic range
      do k=1,nspp                                                       ! PURE NICHE MODEL ! PURE NICHE MODEL ! PURE NICHE MODEL        
        call random_number(zz)                                          ! PURE NICHE MODEL ! PURE NICHE MODEL ! PURE NICHE MODEL 
        if((zz.ge.L1).and.(zz.le.L2))then                               ! PURE NICHE MODEL ! PURE NICHE MODEL ! PURE NICHE MODEL  
          matrixniche(j,k)=1                                            ! PURE NICHE MODEL ! PURE NICHE MODEL ! PURE NICHE MODEL 
        end if
      end do
      
      write(*,*)j,node(j)(1:20),SIZ(J),nspp,' OPTIMUM ',C,' RANGE ',R,' LIMITS ',L1,L2  
      write(*,*)j,node(j)(1:20),CCCC,       ' niche parameters, diets          ',dis(j,:)
      if(dis(j,4).eq.1.0)then ; dis(j,1)=1.0 ; end if                       ! detritivores as photosintetizers      
      if((node(j)(1:6).eq.'Animal').and.(R.le.0.0))then                     ! a range cannot be negative
      goto 63 ; end if                                                      ! a range cannot be negative
      if((L1.gt.-1.0).and.(L1.le.1.0).and.(L2.gt.0.0).and.(L2.le.2.0))then  ! double check. Everything is OK 
      else 
        if (dis(j,3).eq.1)then 
          write(*,*)J,L1,L2,'ERROR: carnivore with no trophic range',sizmins! ERROR. carnivore with no trophic range assigned
          READ(*,*) 
        end if
      end if
      come2=0
      do k=1,nspp                                                       ! J eats K (VERIFY)             
        come=0     
        do jj=1,7                                                       ! coarse-graining matching
          do kk=1,7                                                     ! coarse-graining matching
            if(dis(j,jj).eq.0.0)then ; jjj=0 ; else ; jjj=1 ; end if    ! coarse-graining matching
            if(dis(k,kk).eq.0.0)then ; kkk=0 ; else ; kkk=1 ; end if    ! coarse-graining matching
            if(dis(k,kk).eq.0.0)then ; kkk=0 ; else ; kkk=1 ; end if    ! coarse-graining matching
            if((jjj.gt.0).and.(kkk.gt.0))then                           ! the two species have an assigned trophic category (non-zero)       
              if(matrixdiet(kk,jj).eq.1)then                
                call random_number(x)                                   ! no size restriction. Interaction with probability P(x)                
                zzz=CCCC*dis(j,jj)                                      ! keeping proportionality 
                if(x.le.zzz)then ; come=1 ; end if                      ! keeping the same connectance !
              end if
              if(matrixdiet(kk,jj).eq.2)then                            ! ALLOMETRIC NICHE MODEL APPLIES
                if((siz(k).ge.L1).and.(siz(k).le.L2))then               ! ALLOMETRIC NICHE MODEL APPLIES           
                  do sp1=1,8                                            ! adding specificity 
                  if((pred(j,sp1).eq.1).and.(prey(k).eq.sp1))then       ! predator-prey matching
                    call random_number(x)                               ! proportional readjustment  
                    if(x.lt.dis(j,jj))then                              ! proportional readjustment  
                      come=2        ; niche=niche+1.0 
                    end if                                              ! proportional readjustment  
                  end if                                                ! adding specificity
                  end do                                                ! adding specificity 
                end if
              end if            
            end if
          end do
        end do
        if(come.gt.0)then
          if(j.ne.k) then                                               ! non-cannibal interactions
            write(*,*)replica,nspp,j,k,node(j)(1:20),' se come a ',node(k)(1:20),siz(k),come,force
            come2=1
            matrix(j,k)=1                                               ! j eats k
            edg=edg+1                   ; noniche=noniche+1.0
          else                                                          ! cannibalism allowed but less likely
            call random_number(x) ; if(x.lt.0.5)then
            write(*,*)replica,nspp,j,k,node(j)(1:20),' se come b ',node(k)(1:20),siz(k),come,force           
            matrix(j,k)=1                                               ! j eats k
            edg=edg+1                   ; noniche=noniche+1.0  ; end if
          end if
        end if          
      end do                                                            ! end loop of the "k" especies      
      !!!!!!!!!!!!!!!!    FAST CONTROL  (double check)    !!!!!!!!!!!!!!!      
      if((come2.eq.0).and.(dis(j,1).eq.0.0).and.(dis(j,4).eq.0.0))then  ! non-producer and non-detritivore, but feeding on nothing
        siz(j)=siz(j)*1.1 
        if(dis(j,3).eq.1)then 
          17 call random_number(x) ; if((int(x*8+1).gt.8).or.(int(x*8+1).lt.1))then ; goto 17 ; end if
          pred(j,int(x*8+1))=1                                          ! relaxing trophic preferences
        end if  
        goto 36       
      end if      
    end do                                                              ! end loop of the "j" especies      
    
    do j=1,nspp                                                         ! SPARE NODES PURE NICHE-MODEL  ! SPARE NODES PURE NICHE-MODEL 
      if(sum(matrixniche(1:nspp,j)).eq.0)then                           ! SPARE NODES PURE NICHE-MODEL  ! SPARE NODES PURE NICHE-MODEL 
        do k=1,nspp                                                     ! SPARE NODES PURE NICHE-MODEL  ! SPARE NODES PURE NICHE-MODEL  
           if(sum(matrixniche(k,1:nspp)).eq.0)then                      ! SPARE NODES PURE NICHE-MODEL  ! SPARE NODES PURE NICHE-MODEL 
             call random_number(x) ; if(x.lt.0.5)then                   ! connecting the spare node as a consumer
               187 call random_number(x) ; if((int(x*nspp+1).gt.nspp).or.(int(x*nspp+1).lt.1))then ; goto 187 ; end if
               matrixniche(1:nspp,j)=1  
             else                                                       ! connecting the spare node as a producer
               188 call random_number(x) ; if((int(x*nspp+1).gt.nspp).or.(int(x*nspp+1).lt.1))then ; goto 188 ; end if
               matrixniche(k,1:nspp)=1  
             end if   
           end if                                                       ! SPARE NODES PURE NICHE-MODEL  ! SPARE NODES PURE NICHE-MODEL  
        end do    											            ! SPARE NODES PURE NICHE-MODEL  ! SPARE NODES PURE NICHE-MODEL 
      end if 												            ! SPARE NODES PURE NICHE-MODEL  ! SPARE NODES PURE NICHE-MODEL 
    end do 													            ! SPARE NODES PURE NICHE-MODEL  ! SPARE NODES PURE NICHE-MODEL 
    
    !!!!!!!!!!!!!!! checking spare nodes (things not eaten by anyone)
      
    whois=0 ; edg=1                                                     ! edg is just a counter here ...
    do k=1,nspp                                                         ! (:,K)=0 K is not eaten by anything, K is apex predator 
      if(sum(matrix(1:nspp,k)).eq.0)then                                ! K is not eaten by anything
        sueltosf=sueltosf+1.0
        if(dis(k,1).eq.1.0)then ; problematico(k)=1 ; end if            ! a plant not eaten by anything    
      end if
      if(sum(matrix(k,1:nspp)).eq.0)then                                ! K does not eat anything
        if(dis(k,1).eq.0.0)then ; problematico(k)=2 ; end if            ! an animal that does not eat anything
        sueltosc=sueltosc+1.0    
      else if(sum(matrix(k,1:nspp)).eq.1)then                           ! to fix obligate cannibalism            
        if(matrix(k,k).eq.1)then                                        ! to fix obligate cannibalism                    
          problematico(k)=2 ; sueltosc=sueltosc+1.0    
        end if
      end if
    end do     
    sueltosf=sueltosf/real(nspp) ; sueltosc=sueltosc/real(nspp) ; sueltost=sueltost/real(nspp)
    
    !! N=FORCEMAX LOOP TO ENSURE THE WHOLE CONNECTEDNESS OF THE NETWORK
    !! N=FORCEMAX LOOP TO ENSURE THE WHOLE CONNECTEDNESS OF THE NETWORK      
    do force=1,forcemax
      if(ccc.lt.0.25)then ; CCC=1.1*CCC ; end if                        ! increasing niche width
      do j=1,nspp                                                       ! J eats K (VERIFY)
        if(problematico(j).eq.2)then                                    ! an animal that does not eat nothing 
          177 call random_number(x) ; if((int(x*8+1).gt.8).or.(int(x*8+1).lt.1))then ; goto 177 ; end if
          pred(j,int(x*8+1))=1                                          ! relaxing trophic preferences
          if(force.eq.forcemax/2)then 
            178 call random_number(x) ; if((int(x*7+1).gt.7).or.(int(x*7+1).le.1))then ; goto 178 ; end if
            call random_number(z)     ; dis(j,int(x*7+1))=z             ! relaxing trophic preferences
          end if  
          CC=CCC/(real(nspp))**2 ; CC=2.0*CC ; a=1.0 ; b=(a*(1.0-CC))/CC   
          64 call random_number(x) ; y=1.0-(1.0-x)**b ; R=y*siz(j) ; R=abs(R*0.75) 
          if((node(j)(1:6).eq.'Animal').and.(R.le.0.0))then 
          goto 64 ; end if                                              ! a range cannot be negative
          call random_number(z) ; x=R/2.0                               ! left limit for the interval used to calculate centroid U(x,y)
          if(siz(j).lt.1.0-x)then ; y=siz(j)                            ! right limit for the interval used to calculate centroid U(x,y)
          else ; y=1.0-x; end if                                        ! right limit for the interval used to calculate centroid U(x,y)
          C=x+z*abs(y-x) ; C=0.5*C ; L1=C-(R/2.0) ; L2=C+(R/2.0)        ! upper and lower limits to the trophic range
          CCCC=(CCC/sqrt(real(nspp)))*10.0                              ! just to adjust the non-allometric conections   
          call random_number(zz)                                        ! new
          zz=1.0/(1.0+zz)                                               ! new
          CCCC=cccc*zz                                                  ! new
            CCCC=1.1*CCCC
            eles(j,1)=eles(j,1)*0.8 ; eles(j,2)=eles(j,2)*1.2
            L1=eles(j,1)            ; L2=eles(j,2)
            
          do k=1,nspp                                                   ! J eats K (VERIFY)             
            come=0     
            do jj=1,7
              do kk=1,7
                if(dis(j,jj).eq.0.0)then ; jjj=0 ; else ; jjj=1 ; end if
                if(dis(k,kk).eq.0.0)then ; kkk=0 ; else ; kkk=1 ; end if
                
                if((jjj.gt.0).and.(kkk.gt.0))then                       ! the two species have an assigned trophic category (non-zero)                 
                  if(matrixdiet(kk,jj).eq.1)then                
                    call random_number(x)                               ! no size restriction. Interaction with probability P(x)   
                    zzz=CCCC*dis(j,jj)                                  ! keeping proportionality 
                    if(x.le.zzz)then ; come=1 ; end if                  ! keeping the same connectance !
                  end if
                  if(matrixdiet(kk,jj).eq.2)then                        ! ALLOMETRIC NICHE MODEL APPLIES
                    if((siz(k).ge.L1).and.(siz(k).le.L2))then           ! ALLOMETRIC NICHE MODEL APPLIES           
                      do sp1=1,8                                        ! adding specificity 
                      if((pred(j,sp1).eq.1).and.(prey(k).eq.sp1))then   ! predator-prey matching
                        call random_number(x)                           ! proportional readjustment  
                        if(x.lt.dis(j,jj))then                          ! proportional readjustmen 
                          come=2        ; niche=niche+1.0
                        end if                                          ! proportional readjustmen 
                      end if                                            ! adding specificity
                      end do                                            ! adding specificity 
                    end if
                  end if            
                end if
              end do
            end do
                        
            if(come.gt.0)then
              if(j.ne.k) then                                           ! non-cannibal interactions
                problematico(j)=0 ; matrix(j,k)=1 ; edg=edg+1 ; noniche=noniche+1.0
              else                                                      ! cannibalism allowed but less likely
                call random_number(x) ; if(x.lt.0.5)then
                write(*,*)replica,nspp,j,k,node(j)(1:20),' eats G ',node(k)(1:20),siz(k),force,COME
                problematico(j)=0 ; matrix(j,k)=1 ; edg=edg+1 ; noniche=noniche+1.0  ; end if
              end if
            end if          
          end do                                                        ! end loop of the k species        
        else if(problematico(j).eq.1)then                               ! plant not eaten by anything
          179 call random_number(x) ; kk=int(x*nspp+1)
          if((kk.gt.nspp).or.(kk.lt.1))then ; goto 179 ; end if
          do k=kk,nspp                                                  ! conencted to a random herbivore if necessary
            if(dis(k,2).eq.1)then                                       ! herbivore 
               matrix(k,j)=1 ; edg=edg+1 ; noniche=noniche+1.0          ! k to j (inverted link)
                write(*,*)replica,nspp,j,k,node(k)(1:20),' eats H ',node(j)(1:20),siz(j),force,kk
               problematico(j)=0  ; exit
            end if
            if(k.eq.nspp)then ; goto 179 ; end if
          end do      
        end if         ! end loop of "problematic" nodes  
      end do           ! end loop of the j species      
    end do             ! end loop of the brut-force protocol to reconnect spare nodes
       
    !!!!!! WRITES THE "FINAL" MARIX after solving spare nodes and just before closing the output file ...                 
    edg=1               
    do j=1,nspp
      do k=1,nspp
        if(matrix(j,k).eq.1)then
          call writeline(edg,k,j,3)                                     ! 3 is for trophic edges (necessary for multi-layer extensions) 
          edg=edg+1          !     matrix(j,k)=1                        ! j eats k
        end if
      end do    
    end do
    close(102)                                                          ! EDGES files 
    
    write(*,*)replica,'SPARES=',sum(problematico),sueltost*real(nspp),sueltosf,sueltosc,sueltost,CCC
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!! deg distrib included
    degree=0.0 ; degreeniche=0.0 ; DD1=0.0 ; DD2=0.0 ; DD3=0.0
    do k=1,nspp
      jj=sum(matrix(k,1:nspp))      ; if(jj.gt.0)then  ; degree(jj,1)      =degree(jj,1) +1.0     ; end if
      kk=sum(matrix(1:nspp,k))      ; if(kk.gt.0)then  ; degree(kk,2)      =degree(kk,2) +1.0     ; end if
      jj=jj+kk                      ; if(jj.gt.0)then  ; degree(jj,3)      =degree(jj,3) +1.0     ; end if
      jjj=sum(matrixniche(k,1:nspp)); if(jjj.gt.0)then ; degreeniche(jjj,1)=degreeniche(jjj,1)+1.0; end if
      kkk=sum(matrixniche(1:nspp,k)); if(kkk.gt.0)then ; degreeniche(kkk,2)=degreeniche(kkk,2)+1.0; end if
      jjj=jjj+kkk                   ; if(jjj.gt.0)then ; degreeniche(jjj,3)=degreeniche(jjj,3)+1.0; end if
    end do   
    do j=1,1000
      x=degree(j,1)/maxval(degree(:,1)) ; xx=degreeniche(j,1)/maxval(degreeniche(:,1))
      y=degree(j,2)/maxval(degree(:,2)) ; yy=degreeniche(j,2)/maxval(degreeniche(:,2))
      z=degree(j,3)/maxval(degree(:,3)) ; zz=degreeniche(j,3)/maxval(degreeniche(:,3))
      DD1=DD1+(x-xx)**2  ; DD2=DD2+(y-yy)**2 ; DD3=DD3+(z-zz)**2
    end do
    DD1=sqrt(DD1) ; DD2=sqrt(DD2) ; DD3=sqrt(DD3)  
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! pure niche model
    AVLniche=0.0 ; CONniche=0.0 ; INDniche=0.0 ; OUTDniche=0.0 ; ASYMniche=0.0
    call basic1(nspp,matrixniche,AVLniche,CONniche,INDniche,OUTDniche,ASYMniche)
        
    TOPniche=0.0 ; INTSniche=0.0 ; BASniche=0.0 ; CANniche=0.0
    call basic2(nspp,matrixniche,TOPniche,INTSniche,BASniche,CANniche)
    
    AVGLniche=0.0 ; MINLniche=0 ; MAXLniche=0  ; NCHAINniche=0 ; NLOOPSniche=0
    !call basic4(i,replica,nspp,matrixniche,AVGLniche,MINLniche,MAXLniche,NCHAINniche,NLOOPSniche)
    ! topological assessment of networks (niche model + functional traits)    
    AVL=0.0 ; CON=0.0 ; IND=0.0 ; OUTD=0.0 ; ASYM=0.0
    call basic1(nspp,matrix,AVL,CON,IND,OUTD,ASYM)
        
    TOP=0.0 ; INTS=0.0 ; BAS=0.0 ; CAN=0.0
    call basic2(nspp,matrix,TOP,INTS,BAS,CAN)
    
    AVGL=0.0 ; MINL=0 ; MAXL=0  ; NCHAIN=0 ; NLOOPS=0
    !call basic4(i,replica,nspp,matrix,AVGL,MINL,MAXL,NCHAIN,NLOOPS)    ! COMPLEX TOPOLOGICAL ANALYSIS CARRIED OUT WITH EXTERNAL SOFTWARE
    
    write(*,*)i,replica,'properties',AVL,CON,IND,OUTD,ASYM,TOP,'intsbas',INTS,BAS,CAN
          
    !call random_number(x)
    write(77,*)arx2,morpho,i,nspp,replica,AVL,CON,IND,OUTD,ASYM,TOP,INTS,BAS,CAN,int(MINL),int(MAXL),AVGL,NCHAIN,NLOOPS,DD1,DD2,DD3&
    ,AVL-AVLniche,CON-CONniche,IND-INDniche,OUTD-OUTDniche,TOP-TOPniche,INTS-INTSniche,BAS-BASniche,CAN-CANniche
    flush(77)                             ! 5   6   7   8    9    10  11  12  13  14   15   16   17     18           
                                                                                                                         
    morpho=morpho+1
    deallocate(matrix) ; deallocate(matrixniche)
    
  end do           ! end loop replicates
  deallocate(siz) ; deallocate(NWS) ; deallocate(eles) ; deallocate(NWSniche)	; deallocate(dis) ; close(88)  
  
end do             ! end loop files /ecosystems

close(90) 

ret=SYSTEM('mv *NODES*.csv redes/') 
ret=SYSTEM('mv *EDGES*.csv redes/') 
end program

!!!!!!!!!!!!!!!!!!!            SUBROUTINES            SUBROUTINES            SUBROUTINES            SUBROUTINES            SUBROUTINES            SUBROUTINES
!!!!!!!!!!!!!!!!!!! SUBROUTINES            SUBROUTINES            SUBROUTINES            SUBROUTINES            SUBROUTINES            SUBROUTINES
!!!!!!!!!!!!!!!!!!!            SUBROUTINES            SUBROUTINES            SUBROUTINES            SUBROUTINES            SUBROUTINES            SUBROUTINES

subroutine removespaces(tres,atres)
  character(len=1000)::tres,atres                                 ! format issues ....  
  integer :: c1,c2,c3
  c2=1    ; atres=' '       
  do c1=1,len(tres)
    if(tres(c1:c1).ne.' ')then ; atres(c2:c2)=tres(c1:c1) ; c2=c2+1 ; end if
  end do
  atres=trim(atres)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BASIC 4 ! COMPLEX TOPOLOGICAL ANALYSIS CARRIED OUT WITH EXTERNAL SOFTWARE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BASIC 1
subroutine basic1(ns,matris,AVLs,CONs,INDs,OUTDs,ASYMs)
integer :: ns,i1,j1
real*4  :: avls,cons,inds,outds,asyms                                   ! same variables but "subroutinized"
integer :: matris(ns,ns)                                                ! same variables but "subroutinized"
integer :: idg1,odg1,idg2,odg2,nl

  avls=0.0 ; cons=0.0 ; inds=0.0 ; outds=0.0
  nl=sum(matris)
  AVLs=real(nl)/real(ns)                                                ! AVL: average number of links per species 
  CONs=real(nl)/real(ns)**2                                             ! CON: connectance =i/l**2
  idg1=0 ; odg1=0 ;  idg2=0 ; odg2=0
  do j1=1,ns
    if(sum(matris(1:ns,j1)).ne.0)then ; odg1=odg1+sum(matris(1:ns,j1)) ; odg2=odg2+1 ; end if ! by columns
    if(sum(matris(j1,1:ns)).ne.0)then ; idg1=idg1+sum(matris(j1,1:ns)) ; idg2=idg2+1 ; end if ! by rows
  end do
  INDs =real(idg1)/real(idg2)                                           ! IN-DEGREE.  By colums. Average inputs
  OUTDs=real(odg1)/real(odg2)                                           ! OUT-DEGREE. By colums. Average inputs
  ASYMS=INDS/OUTDS                                                      ! ASYMMETRY (Bascompte)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BASIC 2
subroutine basic2(ns,matris,TOPs,INTSs,BASs,CANs)
integer :: ns,i1,j1
real*4  :: tops,intss,bass,cans                                         ! same variables but "subroutinized"
integer :: matris(ns,ns)                                                ! same variables but "subroutinized"
integer :: idg1,odg1,idg2,odg2,nl

  nl=sum(matris) !; rows=0.0 ; cols=0.0
  tops=0.0 ; ints=0.0 ; bass=0.0 ; cans=0.0
  do j1=1,ns
    if(matris(j1,j1).eq.1)then ; cans=cans+1.0 ; end if                 ! self-interaction is cannibalism
    if(sum(matris(1:ns,j1)).eq.0)then ; tops=tops+1.0 ; end if          !by columns ! TOP PREDATORS 
    if(sum(matris(j1,1:ns)).eq.0)then ; bass=bass+1.0 ; end if          !by rows    ! BASAL, AUTOTROPH
  end do  
  CANS=CANS/real(ns) ; TOPS=TOPS/real(ns) ; BASS=BASS/real(ns) ; INTSS=1.0-(TOPS+BASS)  
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine writeline(edgi,sk,skk,slabel)                         ! it transforms the data into Gephi readable files 
integer :: edgi,sk,skk,slabel

if(sk.lt.10)then
  if(skk.lt.10)then
    if(edgi.lt.10)then                           ; write(102,"(I1,I2,A9,I2,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if((edgi.ge.10 ).and.(edgi.lt.100 ))then; write(102,"(I1,I2,A9,I3,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if((edgi.ge.100).and.(edgi.lt.1000))then; write(102,"(I1,I2,A9,I4,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if( edgi.gt.1000)then                   ; write(102,"(I1,I2,A9,I5,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    end if
  else if((skk.ge.10 ).and.(skk.lt.100))then
    if(edgi.lt.10)then                           ; write(102,"(I1,I3,A9,I2,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if((edgi.ge.10 ).and.(edgi.lt.100 ))then; write(102,"(I1,I3,A9,I3,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if((edgi.ge.100).and.(edgi.lt.1000))then; write(102,"(I1,I3,A9,I4,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if( edgi.gt.1000)then                   ; write(102,"(I1,I3,A9,I5,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    end if
  else if((skk.ge.100).and.(skk.lt.1000))then
    if(edgi.lt.10)then                           ; write(102,"(I1,I4,A9,I2,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if((edgi.ge.10 ).and.(edgi.lt.100 ))then; write(102,"(I1,I4,A9,I3,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if((edgi.ge.100).and.(edgi.lt.1000))then; write(102,"(I1,I4,A9,I4,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if( edgi.gt.1000)then                   ; write(102,"(I1,I4,A9,I5,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    end if
  end if
else if((sk.ge.10 ).and.(sk.lt.100))then
  if(skk.lt.10)then
    if(edgi.lt.10)then                           ; write(102,"(I2,I2,A9,I2,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if((edgi.ge.10 ).and.(edgi.lt.100 ))then; write(102,"(I2,I2,A9,I3,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if((edgi.ge.100).and.(edgi.lt.1000))then; write(102,"(I2,I2,A9,I4,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if( edgi.gt.1000)then                   ; write(102,"(I2,I2,A9,I5,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    end if
  else if((skk.ge.10 ).and.(skk.lt.100))then
    if(edgi.lt.10)then                           ; write(102,"(I2,I3,A9,I2,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if((edgi.ge.10 ).and.(edgi.lt.100 ))then; write(102,"(I2,I3,A9,I3,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if((edgi.ge.100).and.(edgi.lt.1000))then; write(102,"(I2,I3,A9,I4,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if( edgi.gt.1000)then                   ; write(102,"(I2,I3,A9,I5,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    end if
  else if((skk.ge.100).and.(skk.lt.1000))then
    if(edgi.lt.10)then                           ; write(102,"(I2,I4,A9,I2,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if((edgi.ge.10 ).and.(edgi.lt.100 ))then; write(102,"(I2,I4,A9,I3,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if((edgi.ge.100).and.(edgi.lt.1000))then; write(102,"(I2,I4,A9,I4,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if( edgi.gt.1000)then                   ; write(102,"(I2,I4,A9,I5,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    end if
  end if
else if((sk.ge.100).and.(sk.lt.1000))then
  if(skk.lt.10)then
    if(edgi.lt.10)then                           ; write(102,"(I3,I2,A9,I2,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if((edgi.ge.10 ).and.(edgi.lt.100 ))then; write(102,"(I3,I2,A9,I3,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if((edgi.ge.100).and.(edgi.lt.1000))then; write(102,"(I3,I2,A9,I4,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if( edgi.gt.1000)then                   ; write(102,"(I3,I2,A9,I5,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    end if
  else if((skk.ge.10 ).and.(skk.lt.100))then
    if(edgi.lt.10)then                           ; write(102,"(I3,I3,A9,I2,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if((edgi.ge.10 ).and.(edgi.lt.100 ))then; write(102,"(I3,I3,A9,I3,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if((edgi.ge.100).and.(edgi.lt.1000))then; write(102,"(I3,I3,A9,I4,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if( edgi.gt.1000)then                   ; write(102,"(I3,I3,A9,I5,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    end if
  else if((skk.ge.100).and.(skk.lt.1000))then
    if(edgi.lt.10)then                           ; write(102,"(I3,I4,A9,I2,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if((edgi.ge.10 ).and.(edgi.lt.100 ))then; write(102,"(I3,I4,A9,I3,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if((edgi.ge.100).and.(edgi.lt.1000))then; write(102,"(I3,I4,A9,I4,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    else if( edgi.gt.1000)then                   ; write(102,"(I3,I4,A9,I5,I2,A1,I2)")sk,skk,' directed',edgi,slabel,' ',1
    end if
  end if
end if

end subroutine


subroutine getdiet(attri,dieti)
character(len=50)  :: attri(4)
character(len=150) :: dieti
real*4             :: dx
   call random_number(dx)
   ! attri1=phylum ; attri2=class ; attri3=order ; attri4=family 
                         !1234567890                      1234567890123456789012345678901234567890         
   if(attri(3)(1:10).eq.'Proboscide')then ; dieti(1:40)='herbivore|fruit|leaves                  ' ; return ; end if  
   if(attri(3)(1:10).eq.'Primates  ')then ; dieti(1:40)='omnivore|fruit|leaves|insect            ' ; return ; end if              
   if(attri(3)(1:10).eq.'Lagomorpha')then ; dieti(1:40)='herbivore|leaves|root                   ' ; return ; end if   
   if(attri(3)(1:10).eq.'Rodentia  ')then ; dieti(1:40)='herbivore|seed|leaves|omnivore          ' ; return ; end if              
   if(attri(3)(1:10).eq.'Chiroptera')then    
     if(px.gt.0.5)then                    ; dieti(1:40)='insect|carnivore                        ' ; else
                                          ; dieti(1:40)='herbivore|fruit                         ' ; end if ; return ; end if   
   if(attri(3)(1:10).eq.'Carnivora ')then ; dieti(1:40)='carnivore|vertebrate|carrion            ' ; return ; end if              
   if(attri(3)(1:10).eq.'Perissodac')then ; dieti(1:40)='herbivore|leaves                        ' ; return ; end if   
   if(attri(3)(1:10).eq.'Artiodacty')then ; dieti(1:40)='herbivore|leaves                        ' ; return ; end if   
   if(attri(3)(1:10).eq.'Cetacea   ')then ; dieti(1:40)='fish|crustacean|plakto                  ' ; return ; end if                            
   if(attri(3)(1:10).eq.'Testudines')then ; dieti(1:40)='herbivore|fruit|leaves                  ' ; return ; end if            
   if(attri(2)(1:10).eq.'Amphibia  ')then ; dieti(1:40)='carnivore|worm|insect|arthropod|detritus' ; return ; end if        
   if(attri(2)(1:10).eq.'Reptilia  ')then ; dieti(1:40)='carnivore                               ' ; return ; end if
   if(attri(2)(1:10).eq.'Aves      ')then 
     if(px.gt.0.5)then                     ; dieti(1:40)='carnivore                               ' ; else
                                           ; dieti(1:40)='herbivore|seed|fruit                    ' ; end if ; return ; end if
   if(attri(1)(1:10).eq.'Chordata  ')then                                                         ! any other vertebrate order and all fishes ...
     if(px.gt.0.25)then                    ; dieti(1:40)='carnivore|omnivore                      ' ; else
                                           ; dieti(1:40)='herbivore                               ' ; end if ; return ; end if                                 
   if(attri(2)(1:10).eq.'Arachnida ')then  ; dieti(1:40)='carnivore|insect|arthropod              ' ; return ; end if                 
   if(attri(4)(1:10).eq.'Aphididae ')then  ; dieti(1:40)='plant|leaves                            ' ; return ; end if                  
   if(attri(3)(1:10).eq.'Lepidopter')then  ; dieti(1:40)='pollen|plant|leaves                     ' ; return ; end if          
   if(attri(3)(1:10).eq.'Orthoptera')then  ; dieti(1:40)='plant|leaves                            ' ; return ; end if                 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!! for all the remaining organisms     
                                             dieti(1:40)='omnivore                         '                                   
                                           
end subroutine
