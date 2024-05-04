subroutine readifiepieda(inname, frag_i, frag_j, pfrag_i, pfrag_j, dist, hfifie, mp2ifie, &
&prtype1, grimme, jung, hill, es, ex, ct, di, erest, qval, fdimesint, ifpair, pipair)

! 2020/05/17
! Author: Koji Okuwaki
implicit none
integer i, j, skip
character*1, dimension(100000000):: fdimes
character*40 head
integer, dimension(100000000):: frag_i, frag_j, pfrag_i, pfrag_j, fdimesint
double precision, dimension(100000000):: dist, hfifie,mp2ifie,prtype1,grimme,jung,hill
double precision, dimension(100000000):: es,ex,ct,di,erest,qval

integer ifpair, pipair
logical disp
character*200 inname
character*3 :: method = 'Non'

!! Main
! read(*,'(A)') inname
write(*,'(A)') trim(adjustl(inname))

! initialize
ifpair=0
pipair=0
frag_i(:) = 0
frag_j(:) = 0
dist(:) = 0.0
fdimesint(:) = 0
hfifie(:) = 0.0
mp2ifie(:) = 0.0
prtype1(:) = 0.0
grimme(:) = 0.0
jung(:) = 0.0
hill(:) = 0.0

pfrag_i(:) = 0
pfrag_j(:) = 0
es(:) = 0.0
ex(:) = 0.0
ct(:) = 0.0
di(:) = 0.0
erest(:) = 0.0
qval(:) = 0.0
disp=.false.

! Read File
open(17,file=trim(adjustl(inname)), status='old')

! IFIE - Section
do
    read(17,'(a)', end=999) head
    if (trim(adjustl(head)) == "Disp                    = ON") then
        disp=.true.
        print *, '-- DISP mode --'
    end if

    if (trim(adjustl(head)) == "## MP2-IFIE" .or. trim(adjustl(head))=="## MP3-IFIE" &
        .or. trim(adjustl(head))=="## MP4-IFIE") then
        method='MP'
        ! write(*, '(a)') 'Read IFIE Section'
        do skip=1,4  !skip 4lines
            read(17,'()')
        end do
        goto 110

    else if (trim(adjustl(head))=="## HF-IFIE") then
        ! write(*, '(a)') 'Read IFIE Section'
        method='HF'
        do skip=1,4  !skip 4lines
            read(17,'()')
        end do
        goto 110

    else if (trim(adjustl(head))=="## LRD-IFIE") then
        method='LRD'
        ! write(*, '(a)') 'Read IFIE Section'
        do skip=1,4  !skip 4lines
            read(17,'()')
        end do
        goto 110
    end if
end do
110 continue

! write (*, '(a)') trim(adjustl(method))
! Read IFIE val
i = 1
do
    ! (MP2:) J-PAIR    DIST     DIMER-ES   HF-IFIE    MP2-IFIE   PR-TYPE1   GRIMME     JUNG       HILL
    ! (MP3:) IJ-PAIR    DIST     DIMER-ES   HF-IFIE    MP2-IFIE   USER-MP2   MP3-IFIE   USER-MP3   PADE[2/1]
    if (trim(adjustl(method)) == "MP") then
        read(17, *,  err=200)frag_i(i), frag_j(i), dist(i), fdimes(i), &
            &hfifie(i), mp2ifie(i), prtype1(i), grimme(i), jung(i), hill(i)
    else if (trim(adjustl(method)) == "HF") then
        read(17, *,  err=200)frag_i(i), frag_j(i), dist(i), fdimes(i), hfifie(i)
    else if (trim(adjustl(method)) == "LRD") then
        read(17, *,  err=200)frag_i(i), frag_j(i), dist(i), fdimes(i), hfifie(i), mp2ifie(i)
    ! if(i == 100) write(*, '(2i7, f10.3, a3, 6f10.3)') frag_i(i),frag_j(i),dist(i),fdimes(i), &
        ! & hfifie(i),mp2ifie(i),prtype1(i),grimme(i),jung(i),hill(i)
    end if
    i = i + 1
end do

200 ifpair = i - 1
! write(*, '(a, i8)') 'ifpair = ', ifpair
backspace(17) ! back 1 line
read(17,'(a)') head

!! PIEDA -SECTION
if (trim(adjustl(head))=="## PIEDA") then   !この下にIFIE情報
    ! write(*, '(a)') 'Read PIEDA Section'
    do skip=1,4                        !4行とばす
        read(17,'()')
    end do
  goto 120
else if(trim(adjustl(head))=="## Mulliken") then
    write(*, '(a)') 'log do not have pieda data: end'
    goto 998
else
    goto 998
end if

!! Read PIEDA val
120 continue
i = 1
do
        ! IJ-PAIR       ES             EX             CT+mix         DI(MP2)        q(I=>J).
    if (disp .eqv. .false.) then
        read(17,*, err=210) pfrag_i(i),pfrag_j(i),es(i),ex(i),ct(i),di(i),qval(i)
    else
        read(17,*, err=210) pfrag_i(i),pfrag_j(i),es(i),ex(i),ct(i),di(i),erest(i),qval(i)
    end if

    ! if(i == 100) write(*, '(2i7, 5f10.3)') pfrag_i(i),pfrag_j(i),es(i),ex(i),ct(i),di(i),qval(i)
    i = i + 1
end do

210 pipair = i - 1
! write(*, '(a, i8)') 'pipair = ', pipair
! write(*, '(a)') 'Read end'

998 continue
! post process for ifie
fdimesint(:) = 0
do i = 1, ifpair
    if (hfifie(i) < -2)  then
        hfifie(i)=0.
        mp2ifie(i)=0.
        prtype1(i)=0.
        grimme(i)=0.
        jung(i)=0.
        hill(i)=0.
    end if

    if ((trim(adjustl(fdimes(i)))) == "F") then
        fdimesint(i) = 0
    else if ((trim(adjustl(fdimes(i)))) == "T") then
        fdimesint(i) = 1
    end if

    hfifie(i) =  hfifie(i) * 627.5095
    mp2ifie(i)=  mp2ifie(i) * 627.5095
    prtype1(i)=  prtype1(i) * 627.5095
    grimme(i) =  grimme(i) * 627.5095
    jung(i)   =  jung(i) * 627.5095
    hill(i)   =  hill(i) * 627.5095
end do

999 continue
end subroutine

!      ## HF-IFIE
! 
!            IJ-PAIR    DIST     DIMER-ES   HF-IFIE        HF-IFIE        HF-IFIE
!                       / A      APPROX.   / Hartree      / kcal/mol     / kJ/mol
!         --------------------------------------------------------------------------
!             2    1    0.000000   F      -15.064065   -9452.843160  -39550.695783

!      ## PIEDA
! 
!            IJ-PAIR       ES             EX             CT+mix         DI             q(I=>J)  
!                          / kcal/mol     / kcal/mol     / kcal/mol     / kcal/mol     / e      
!         -----------------------------------------------------------------------------------------------------
!             3    1       1.152159       0.015701      -0.718282       0.000000       0.001781
!             4    1       0.869703       0.025875      -0.141362       0.000000      -0.000038
!             4    2       0.672861       0.942428      -1.383603       0.000000       0.006733
!             5    1       2.328359       0.000275      -0.104993       0.000000      -0.000337
!             5    2     -10.950065       3.413130      -1.410642       0.000000       0.021846
!             5    3       2.594522       1.238201      -1.481940       0.000000       0.005841
! 
!      ## LRD-IFIE
! 
!            IJ-PAIR    DIST     DIMER-ES   HF-IFIE    LRD-IFIE
!                       / A      APPROX.   / Hartree  / Hartree
!         ----------------------------------------------------------------------------------------------------
!             2    1    0.000000   F      -15.064065  -0.001832
!             3    1    3.405247   F        0.000716  -0.000527
!             3    2    0.000000   F      -15.064590  -0.003540
!             4    1    3.399003   F        0.001202  -0.000583
!             4    2    2.756835   F        0.000369  -0.001841
!             4    3    0.000000   F      -15.071663  -0.003561
!             5    1    3.542065   F        0.003544  -0.000280
!             5    2    2.105612   F       -0.014259  -0.001511
!             5    3    2.756819   F        0.003746  -0.002613
!             5    4    0.000000   F      -15.064453  -0.004435
! 
!      ## PIEDA
! 
!            IJ-PAIR       ES             EX             CT+mix         DI(LRD)        q(I=>J)  
!                          / kcal/mol     / kcal/mol     / kcal/mol     / kcal/mol     / e      
!         -----------------------------------------------------------------------------------------------------
!             3    1       1.152159       0.015701      -0.718282      -0.330770       0.001781
!             4    1       0.869703       0.025875      -0.141362      -0.366012      -0.000038
!             4    2       0.672861       0.942428      -1.383603      -1.155489       0.006733
!             5    1       2.328359       0.000275      -0.104993      -0.175521      -0.000337
!             5    2     -10.950065       3.413130      -1.410642      -0.948301       0.021846
!             5    3       2.594522       1.238201      -1.481940      -1.639756       0.005841

