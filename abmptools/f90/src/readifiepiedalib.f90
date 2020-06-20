subroutine readifiepieda(inname, frag_i, frag_j, pfrag_i, pfrag_j, dist, hfifie, mp2ifie, &
&prtype1, grimme, jung, hill, es, ex, ct, di, qval, fdimesint, ifpair, pipair)

! 2020/05/17
! Author: Koji Okuwaki
implicit none
integer i, j, skip
character*1, dimension(10000000):: fdimes
character*25 head
integer, dimension(10000000):: frag_i, frag_j, pfrag_i, pfrag_j, fdimesint
double precision, dimension(10000000):: dist, hfifie,mp2ifie,prtype1,grimme,jung,hill
double precision, dimension(10000000):: es,ex,ct,di,qval

integer ifpair, pipair
character*200 inname

!! Main
! read(*,'(A)') inname
! write(*,'(A)') trim(adjustl(inname))

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
qval(:) = 0.0
! Read File
open(17,file=trim(adjustl(inname)), status='old')

! IFIE - Section
do
    read(17,'(a)', end=999) head
    if (trim(adjustl(head))=="## MP2-IFIE" .or. trim(adjustl(head))=="## MP3-IFIE" &
        .or. trim(adjustl(head))=="## MP4-IFIE") then
        ! write(*, '(a)') 'Read IFIE Section'
        do skip=1,4  !skip 4lines
            read(17,'()')
        end do
      goto 110
    end if
end do
110 continue

! Read IFIE val
i = 1
do
    ! (MP2:) J-PAIR    DIST     DIMER-ES   HF-IFIE    MP2-IFIE   PR-TYPE1   GRIMME     JUNG       HILL
    ! (MP3:) IJ-PAIR    DIST     DIMER-ES   HF-IFIE    MP2-IFIE   USER-MP2   MP3-IFIE   USER-MP3   PADE[2/1]
    read(17, *,  err=200)frag_i(i), frag_j(i), dist(i), fdimes(i), &
        &hfifie(i), mp2ifie(i), prtype1(i), grimme(i), jung(i), hill(i)
    ! if(i == 100) write(*, '(2i7, f10.3, a3, 6f10.3)') frag_i(i),frag_j(i),dist(i),fdimes(i), &
        ! & hfifie(i),mp2ifie(i),prtype1(i),grimme(i),jung(i),hill(i)
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
    write(*, '(a)') 'log dont have pieda data: end'
    goto 998
else
    goto 998
end if

!! Read PIEDA val
120 continue
i = 1
do
        ! IJ-PAIR       ES             EX             CT+mix         DI(MP2)        q(I=>J).
    read(17,*, err=210) pfrag_i(i),pfrag_j(i),es(i),ex(i),ct(i),di(i),qval(i)
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
