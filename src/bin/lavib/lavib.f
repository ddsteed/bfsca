C * JAN 01 2015 - RDS - REMOVE SYMMETRY NUMBER DEPENDENCY
C *
C * JUN 24 2006 - RDS - ADD SUBROUTINE GAULEG TO GENERATE GAUSS-LEGENDRE ABSCISSAS AND
C *                   - WEIGHTS WITH ANY DIMENSION
C *
C * JUN 22 2006 - RDS - MODIFY LAVIB.F SO THAT 
C *                   - INPUT POTENTIAL, EXCHANGE AND LA MESH COULD BE DIFFERENT
C *                   - TRAPEZOIDAL INTEGRALS GIVE THE SAME RESULTS AS THE ORIGINAL CODE
C *
C * JUN 13 2006 - RDS - MODIFIED ALL SUBROUTINES OF LA CAL. WITH NEW INTEGRAL FORMULA
C *                   - (Gauss-Legendre or trapezoidal points) TRAPEZOIDAL INTEGRAL
C *
C * JUN 06 2006 - RDS - FIX THE BUG OF DEFINING NEXDIM (IN SUBROUTINE OF INDATA)
C *
C * MAY 27 2006 - RDS - REWRITE CLEBEXCH AND KERNCALC TO COPE WITH E-N2
C *                   - ADD A JUDGEMENT TO ENSURE THE CONVERGENCE OF EXCHANGE KERNEL- nexdim >= nlproj
C *                   - FIX THE BUG OF VNUC [ADD (-1)**LAM TO LEFT ATOM]
C *
C * MAY 10 2006 - RDS - READ VIBRATIONAL SPECTROSCOPIC CONSTANTS FROM OUTER FILE
C *                   - DEFINED MULTI-DIMENSIONAL ARRAYS AS FIXED BOUNDARIES 
C *                   - IF THEY MAY BE PASSED TO MAIN PROGRAM (DYNAMIC ALLOCATED 
C *                   - ARRAYS ARE EASILY PASSED WRONG ADDRESSES.)
C *
C * MAY 02 1991 - WKT - THE CLOSED CHANNEL BUSINESS WORKS (BELIEVE IT OR NOT)
C *
C * APR 19 1991 - WKT - PRODUCES CORRECT VIBRATIONAL EXCITATION CROSS SECTIONS WITH EXACT EXCHANGE
C *
C * APR 12 1991 - WKT - WORKS FOR VIBRATIONAL EXCITATION FOR A LOCAL POTENTIAL
C *
C * Mar 26 1991 - WKT - THE R-MATRIX PROPAGATOR IS IN WORKING ORDER
C *
C * Mar 19 1991 - WKT - NON-LOCAL POTENTIAL WORKS TO 10 BOHR
C *                   - NOW WE ARE INSTALLING TH R-MATRIX PROPAGATOR FOR THE ASYMPTOTIC REGION
C *
C * Mar 01 1991 - WKT - THE NON-LOCAL POTENTIAL IS IN THE WORKS HERE
C *                   - THE NON-LOCAL EXCHANGE POTENTIAL NOW WORKS
C *
C * PROGRAM LAVIB
C *
C * This program uses the Direct iterative procedure in the
C * Linear Algebraic method to solve the Schrodinger
C * equation for e-Mol scattering
C *
C * Important notes:
C *   1) This program is written largely without common statments.  
C *      Common statements often make a program difficult to decypher
C *      so all variables are passed through the argument list of the 
C *      call to the subroutine.
C *
C *   2) This program is ONLY for v = 0 ---> v' since kw(i) is cal. 
C *                               *****
C *      by E_v(v') - E_v(0) in subroutine chener (Hao Feng, May 20, 2006)
C *
C * Note:  
C *   1) The R-matrix propagator routine will produce a symmetric 
C *      K-matrix almost no matter what it starts with, however, the LA 
C *      part of the program will not.  So if the K-matrix that is 
C *      produced at the end of the LA region is not symmetric, then 
C *      something is wrong, for example,
C *      NITER may be too small so there are not enough basis functions 
C *      to produce an accurate solution.
C *
C *   2) The LA max bound is at the exchange bound
C *
C *-/

      program lavib
      implicit none

      integer chanmax,itermax,nrgmax,lammax
      integer enermax,setmax,boundmax,vibmax,pwavemax
      integer nlmomax,exlmmax
      parameter(chanmax=165,itermax=10,nrgmax=100,enermax=120,
     $     lammax=21,setmax=1,boundmax=6,vibmax=15,pwavemax=11,nlmomax
     $     =16,exlmmax=60)
      integer potptsmax,ptsmax
      parameter(potptsmax=330,ptsmax=331)

C *   NCH: the total number of channels in the scattering problem
C * NCHSR: the number of channels in the Short Range region
C *        as the program goes to large r, nch is usually reduced
C *        so nchsr retains the number of channels that were used
C *        in the Short-Range region
C *
C * NITER: the number of iterations to be done
C *        i.e. the number of wavefunction basis vectors
C *
C * NBOUND: the number of bound target orbitals
C * NVOPEN: the number of OPEN vibrational channels... 
C *         it is calculated in the routine CHENER and is used in BESSL
      integer nch,nchsr,niter,nbound,nvopen

C * below is the mesh data... needs to be modified for
C * varying the mesh with channel index
C *
C * NREGEX: the number of regions used in calculation of the exchange
C *         potential and exchange integrals...  
C * NPTSEX: the number of radial points in the exchange potential
C *
C *  NREGPOT: the number of integration regions in the input potential
C *  NPTSPOT: the number of points in the input potential
C *
C * NREG: the number of integration regions
C * NPTS: the number of radial mesh points in the LA region
C *
C * For trapezoidal grid, NREG = NREGPOT >= NREGEX
C * For Gauss-Legendre grid, these regions make no sense
C * However, for both, NPTS = NPTSEXCH
      double precision step(nrgmax),stepex(nrgmax),steppot(nrgmax)
      integer nstep(nrgmax),nstepex(nrgmax),nsteppot(nrgmax)
      integer nreg,nregex,nregpot,nptsex,nptspot,npts

C * IRMAT = 1 means we'll be using the R-Matrix Propagator
      integer irmat

C * RMATSTEP: the step size to be used in the R-Matrix propagator 
C *           (if it is being used)
C *    RASYM: R at which we quit integrating
C *     REND: RMAX if no R-Matrix propagation, RASYM otherwise
C *           and is used when doing Born r-closure
      double precision rmatstep,rasym,rend

C * ALPHA0: the Spherical Polarizability
C * ALPHA2: the Non-SPherical Polarizability
C *      Q: the Quadrupole Moment
      double precision alpha0(vibmax,vibmax),alpha2(vibmax,vibmax)
      double precision q(vibmax,vibmax) 

C * NLAMS: number of lambdas in the long range region
      integer nlams

C * b: the solution to the Lin Alg eqns ... Ax=b
C *    for the LA method without the R-Matrix Propagator
C *    b = g1.  
C *    If, however, the R-Matrix propagator
C *    is used then b=(1/k)(n+Cj)  C=n'(a)/j'(a)
C *    you get the idea.
C *    See Collins & Schneider PRA 1981 (p2387) but be
C *    warned, there are lots of errors in the R-Matrix part
C *    of that paper.
      double precision b(chanmax,ptsmax)

C * g1, g2: the regular and irregular functions used in the solution
C *         of the integral equations 
      double precision g1(chanmax,ptsmax),g2(chanmax,ptsmax)

C * I1, I2: the I matrices from the Integral Equation Prop method
C *         Note that the I2 matrices are propagated in the reverse
C *         direction (that is toward the origin) to remove stability 
C *         problems
      double precision i1(chanmax,chanmax,ptsmax)
      double precision i2(chanmax,chanmax,ptsmax)
      double precision q1(chanmax,chanmax,ptsmax)
      double precision q2(chanmax,chanmax,ptsmax)

C * XHI: the unnormalized wavevector
C * PSI: normalized
      double precision xhi(chanmax,chanmax,ptsmax)
      double precision psi(chanmax,chanmax,ptsmax)
      double precision kpsi(chanmax,chanmax)

C * PSIFINAL: the solution set of wavevectors
C * XHIFINAL: the unnormalized set
      double precision psifin(itermax,chanmax,chanmax,ptsmax)
      double precision xhifin(itermax,chanmax,chanmax,ptsmax)

C * V: the interaction potential
      double precision v(chanmax,chanmax,ptsmax)

C * KMAT: the K-Matrix
      double precision kmat(chanmax,chanmax)
      double precision kmattemp(chanmax,chanmax)
      double precision eigshft(chanmax), eigsum
      double precision eigvec(chanmax,chanmax)

C * TMATI, TMATR: the real and imag parts of the T-Matrix
      double precision tmati(chanmax,chanmax), tmatr(chanmax,chanmax)

C *   CROSS: contains the cross sections temporarily having the
C *          same dimension as the T-Matrix
C * CROSSLA: has the cross sections calculated at the end of the LA region
C *   KW(i): the wavenumber of the ith channel
      double precision cross(enermax,vibmax), kw(chanmax)
      double precision crossla(enermax,vibmax)

C * ASYMPSI(n,n): the wavefunction matrix at the end of the Linear 
C *               Algebraic region
C *         RMAT: the R-Matrix 
C *      ASYMPSI: truncated and put into RMAT when the R-Matrix option is on
      double precision asympsi(chanmax,chanmax)
      double precision rmat(chanmax,chanmax)
      integer i,j,k,numv,iex

      integer l0,lmax,symlam

C * VLAM: contains the Legendre projections of the local potential
      double precision vlam(lammax,vibmax,vibmax,potptsmax)

C *   SYMLAM: the symmetry label capital Lambda Sigma ==> 0, Pi ==> 1....
C * NLAMDAIN: the number of lamdas read in for the local potential
C *           currently this number is used everywhere
C * NLAMLOC: the total number of legendre (electronic + nuclear)
C *          to be included in the local part of the input potential
C *          Currently these two numbers have to be the same for 
C *          Vibrational excitation calculations, but for RR calcs 
C *          NLAMLOC can be larger and the code will pad out using
C *          the analytic nuclear terms
      integer nlamdain
      integer nlamloc

C * FMAT: contains the Clebsch Gordan for the exchange terms.. (called 
C *       h_lam in Morrison and Collins PRA '78)
      double precision fmat(pwavemax,pwavemax,lammax)

C * rpot: r corresponding to input potential
C *  ity: = 0 Gauss-Legendre quadrature
C *       = 1 Trapezoidal quadrature
C *  rgs: Gauss-Legendre nodes
C *  wtt: Gauss-Legendre weights
      integer ity
      double precision rpot(potptsmax),rgs(ptsmax),wtt(ptsmax)

C * choice: the initial guess
C * ipot: the potential choice
      integer choice,ipot

C * NUMENER: the number of energies read in
C *   IENER: an index for counting the energies
      integer iener,numener

C * ENERGY(): the list of energies to be used
      double precision energy(enermax)
      double precision acoefs(chanmax*itermax,chanmax)

C * channel: indexer and Vibrational data
C *    NVIB: the number of vibrational states in the scattering CALCULATION
C *          this number has to be less than NVIBIN.  The valu of NVIB as
C *          readin from unit 5 is the number of vibrational channels to
C *          be used in the short range (LA) part of the scattering 
C *          calculation. the value of NVIB is then set to NVTRUNC for
C *          the R-Matrix Propagation.
C * NVIBSYM: the number of vibrational levels to include in the
C *          scattering calculation for each symmetry
C *  NVIBIN: the number of vibrational levels in the input potential
C *          AND in the long range moments (maybe that should be changed)
C *    NOFR: not currently used (set it to 1)
C *     NPW: carries the number of partial waves to be included in the
C *          scattering calculation (the vlaue of NPW is set to NLTRUNC
C *          for the R-Matrix propagation part of the program)
C *   NPWIN: stores the value of NPW that was used in the short range
C *          part of the calculation because the long range part needs
C *          it when getting angular matrix elements from the matrix FMAT

C * IRR = 1: if this is a Rigid Rotator Calculation or a VIBAV Calculation
C *          it means that the local part of the input potential is read
C *          in from a formatted file (unstead of the unformatted file)
C *   INFIL: The name of the input file (as usual) 
C * RCLOSURE = 1: if we are doing Born r-closure
      integer chin(chanmax,2),nvib,nvibin,nofr,npw,npwin
      integer nvibsym
      character*8 infil
      integer irr,ind,rclosure

C * the exchange stuff
C *
C * KERNLFIL: names of the files containing the kernel (the '.ker' will be appended internally here) 
C *  NAMELEN: containing the length in characters of the above
C *           filenames.  NAMELEN is calculated in the routine INDATA
      character*8 kernlfil
      integer namelen

C * EXUNIT: the unit associated with the exchang kernel. I set it to 55
C *  NCHEX: the number of partial wave l's in the sphercal projections
C *         that we want to include in the scattering calculation
C * NEXDIM: this is the number of partial waves used in the CREATION of 
C *         the Exchange kernel.  It is important to get this number 
C *         right so that the unformatted reads of the exchange kernel
C *         work out ok.
      integer exunit,nchex,nexdim

C *    NVIBX: the number of vib levels in the exchange kernel that are 
C *           retained in this calculation
C *   NPWAVX: the number of partial wave channels in the exchange kernel
C *           that are retained in the calculation
C *  NVIBXIN: the number of vibrational levels in the input potential
C * NPWAVXIN: the number of partial wave channels in the input potential
      integer nvibx,npwavx
      integer nvibxin,npwavxin

C * SPHRJ(): the spherical projections of the molecular orbitals of the
C *          target.
      double precision sphrj(boundmax,chanmax,ptsmax)

C * NLPROJ: the number of projections for the orbital
C * MLPROJ: the Ml for that orbital
C * L0PROJ: the L0 for the orbital
      integer nlproj(boundmax),mlproj(boundmax),l0proj(boundmax)
      double precision clebx1(boundmax,chanmax,exlmmax,nlmomax)
      double precision clebx2(boundmax,chanmax,exlmmax,nlmomax)
      double precision clebx3(boundmax,chanmax,exlmmax,nlmomax)
      double precision clebx4(boundmax,chanmax,exlmmax,nlmomax)
      double precision clebx5(boundmax,chanmax,exlmmax,nlmomax)

C * nbset: the number of different iteration numbers with which
C *        to run the program
C * numiters(20): contains the above numbers
      integer nbset,numiter(setmax),ibset

C * rmax: where the integration stops
      double precision rmax

C * truncation parameters for truncating the r-matrix to a small
C * number of channels in the asymptotic region
C * NVTRUNCI: the requested number of vibrational channels to keep
C *           in the outer region.  This number is the put into NVTRUNC 
C *           and used where appropriate in the program (subroutine 
C *           RTRUNC).  If there are fewer than NVTRUNCI channels 
C *           actually open then NVOPEN is put into NVTRUNC.
C *  NLTRUNC: the number of partial waves to keep in the asymptotic
C *           region.
C * NCHTRUNC = NVTRUNC*NLTRUNC = the total number of channels in the 
C *           outer region
      integer nvtrunci,nvtrunc,nltrunc,nchtrunc

      double precision faclog(500)
      double precision rac,u9
      integer ia,ib,ic,id,ie,if,l9(9),kselc

C * iswit: switch for units of vibrational constants
C *    we: harmonic vibrational spectroscopic constant
C *  wexe: non-harmonic vibrational spectroscopic constant
      integer iswit
      double precision we, wexe

C * variables for C-G coeff.
      common/faclg/faclog,rac,u9,ia,ib,ic,id,ie,if,l9,kselc

C * dummy variables for print
      integer ii, jj, nf, ipts

C * set up the factorial table
      kselc = 1
      call facset

C * start the main program

      call indata(nbset,numiter,l0,lmax,symlam,nlamdain,numener
     $     ,energy,iex,choice,ipot,nreg,nstep,step,nregex,nstepex,stepex
     $     ,nregpot,nsteppot,steppot,npts,nptsex,nptspot,nbound,nlproj
     $     ,mlproj,l0proj,irmat,rmax,rmatstep,rasym,alpha0,alpha2,q
     $     ,nlams,rpot,ity,rgs,wtt,nchex,kernlfil,namelen,nexdim,nofr
     $     ,nvibsym,nvibin,infil,nvtrunci,nltrunc,irr,nlamloc,nvibx
     $     ,npwavx,nvibxin,npwavxin,rclosure,iswit,we,wexe)

C * this routine gets the potential matrix elements into the Matrix V
C * the routine reads unformatted files when vibration is allowed
C * and the standard formatted "VLAM" output files otherwise
      ind = index(infil,' ')
      if(ind.eq.0)ind = 9
      if(irr.eq.1) then
C        write(6,*)'input file is ',infil(1:ind-1)//'.dat'
         write(6,*)'input file is ',infil(1:ind-1)
         call potread(vlam,steppot,alpha0(1,1),alpha2(1,1),q(1,1)
     $        ,nlamdain,nlamloc,nptspot,infil,nregpot,nsteppot)
      else
C        write(6,*)'input file is ',infil(1:ind-1)//'.vib'
         write(6,*)'input local potential file is ',infil(1:ind-1)
         call vpotread(vlam,infil,nlamdain,nvibin,nofr,nptspot)
      endif
      nlamdain = nlamloc

C * might want to consider rearranging the order of energy and symmetry
C * loops to take advantage of the calculation of the exchange kernel
C * which is energy independent, but depends on the number of 
C * channels (unless we set up variables for exchange channels)(!)
      do iener = 1,numener
         exunit = 55
         write(0,'(a,i2,a,i2,a,f6.3,a)') "   Now cal. No.",iener," of ",
     $        numener,"; Energy = ", energy(iener)*27.2116d0, " eV"

         npw = (lmax-l0)/2 + 1
         nvib = nvibsym
         nvtrunc = nvtrunci
         npwin = npw
         nch = nvib*npw
         write(6,*) " nch = ", nch
         nchsr = nch

C * set up the channel indexing matrix
         write(6,*)
         write(6,*)' The channel indexing is as follows'
         call chindex(chin,npw,nch)
         write(6,*)

C * Exchange stuff
C * if IEX=2 then we already have the exchange kernel
C * and just need to read it out of FORT.55 
         if (iex.eq.2) then
C * Open as unit EXUNIT the file containing the appropriate exchange
C * kernel.  Notice that since iex=2 the file is opened as 'old'
            open(unit=exunit,status='old',file =kernlfil(1:namelen),form
     $           ='unformatted')

C * Otherwise we have to calculate the exchange kernel
         else if (iex.eq.1) then
C * Judge if the number of partial waves is enough
            call judgenpwav(nbound,nlproj,nexdim)

C * Open as unit EXUNIT the file containing the appropriate exchange
C * kernel.  Notice that since iex=1 the file is opened as 'new'
            open(unit=exunit,status='new',file =kernlfil(1:namelen),form
     $           ='unformatted')

C * If we are just beginning, then read in the spherical projections
C * Of course, they only need be read in once for all syms and eners
            if (iener .eq. 1) then
               call sphjread(sphrj,nlproj,nbound,nptspot,nptsex,rgs
     $              ,nexdim)
            endif

            call clebexch(clebx1,clebx2,clebx3,clebx4,clebx5,symlam,l0
     $           ,nlproj,mlproj,l0proj,nbound ,nexdim)

C * Calculate the exchange Kernel for this symmetry
C * of course, you should only do this once for all energies
            call kerncalc(sphrj,xhi,clebx1,clebx2,clebx3,clebx4 ,clebx5
     $           ,nlproj,l0proj,exunit,l0,nexdim ,nbound,rgs,nptsex)

C * Once we have calculated the exchange kernels for each symmetry then
C * we don't need to do it anymore (since they are energy independent)
            iex = 2
         endif

C * the subroutine fmatrix calculates the angular coupling
C * coefficients for the local potential.. (Clebsch Gordons)
         call fmatrix(fmat,chin,l0,nlamdain,symlam,npw ,nch)

         do ibset = 1,nbset
            niter = numiter(ibset)
            write(6,*) 'Energy =  ',energy(iener)
            write(6,*) 'Number of Iterations: ',numiter(ibset)

C * kw(nch) is the energy of the particular channel
C * this routine sets up the channel energies
            call chener(energy(iener),kw,chin,nvib,nch,nvopen,iswit,we
     $           ,wexe)

            write(6,*)'there are',nvopen,' open channels'
            write(6,*)
            write(6,*)' The channel energies are as follows'
            do i = 1,nch
               if (chin(i,1).le.nvopen)then
                  write(6,*)i,kw(i),13.6058*kw(i)*kw(i)
               else
                  write(6,*)i,kw(i),-13.6058*kw(i)*kw(i)
               endif
            end do

C * we now set up the greens function arrays from the routine bessel
            call bessl(irmat,g1,g2,b,kw,l0,rgs,chin,nch,npts ,nvopen)

C * call pot to fill the v matrix for all channels at all mesh points
            call pot(v,fmat,vlam,chin,nlamdain,nlamloc,nvibin,nch
     $           ,nptspot,npts,rpot,rgs,npw)

C * now we submit the first guess to the solution... it will be 
C * unnormalized and we will call it xhi^0
            call frstgues(g1,xhi,choice,nch,npts)

C * now we start the crunch part of the code.  once we have the trial
C * initial function above, we generate the rest of the through the`
C * direct iterative variational procedure.  this is a short cut
C * (at least in the memory sense) to the linear algebraic matrices.
            do numv = 1,niter

C * First we rewind the file containing the exchange kernel
C * so that we can get the exchange potential
C * at every iteration
               rewind(exunit)

C * then we gram-schmidt orthonormalize these vectors with the already
C * orthonormal set currently in psi and we add them to psi
               call gramschm(xhi,xhifin,psifin,psi,numv,step,nstep ,nreg
     $              ,nch,npts,wtt,niter,ity)

C * first we calculate the q1 and q2 matrices from collins and schneider 
C * cpc paper

C * if iex=1 then we call the exchange version and we pass XHI for it to
C * use as a dummy variable for storing the exchange kernel
C * also, we send ibn the array i1 to be used as a dummy array
               if (iex .ne. 0) then
                  call qget3x(q1,q2,g1,g2,v,xhi,i1,exunit,psi,kpsi,
     $                 nstep,step,nreg,nregex,nptsex,nch,npts,nchex
     $                 ,nexdim,nvib,nvibx ,nvibxin,npw,npwavx,npwavxin
     $                 ,wtt,ity)
               else
                  call qget3(q1,q2,g1,g2,v,psi,nstep,step,nreg,nch,npts
     $                 ,wtt,ity)
               endif
                  
C * then from the q matrices we get the i matrices
               call imats(i1,i2,q1,q2,nch,npts)

C * then from the i matrices and the greens functions we get 
C * the new wavefunction basis set iterates (xhi^1)
               call newxhi(xhi,g1,g2,i1,i2,nch,npts,chin)
            end do

C * now we have a set of vectors to span the space of the solution
C * to the linear algebraic equations
C * so we now calculate the projection coefficients of the solution
C * on the orthonormal basis we created.
            call expcoefs(psifin,xhifin,acoefs,b,step,nstep,nreg,nch
     $           ,npts,niter,wtt,ity)

C * this routine calculates the final solution from the basis
C * of orthonormal functions stored in psifin... the solution
C * gets stored in psi
            call psifinal(psifin,psi,acoefs,asympsi,nch,npts,niter)

            k = 6
            if(niter.lt.k) k = niter

C * get the K-matrix
C * if we plan to propagate the r-matrix into the Asymptotic
C * region then that is done here before K-matrix extraction
            if (irmat .eq. 1) then
               write(6,*)
               write(6,*)
               write(6,*) ' The K-matrix before R-Mat prop:'
               write(20,*)
               write(20,*) ' The K-matrix before R-Mat prop:'
               call rkmat(kmat,kw,rmax,l0,asympsi,chin,nch)
               write(20,*)
               write(20,*) 'R-Matrix before truncation'
               do i = 1, nch
                  write(20,504) (asympsi(i,j), j = 1,nch)
               end do
               write(20,*)

C * -----------------------------------------------------------------
C * This little section of code truncates the R-Matrix before
C * propagation into the asymptotic region
               write(6,*)'number of open channels',nvopen
               if (nvtrunc .gt. nvopen) nvtrunc = nvopen
               nchtrunc = nvtrunc*nltrunc
               call rtrunc(asympsi,nvtrunc,nltrunc,nchtrunc,rmat,chin
     $              ,nch)
               nch = nchtrunc
               nvib = nvtrunc
               npw = nltrunc
               call chindex(chin,nltrunc,nch)
               call chener(energy(iener),kw,chin,nvib,nch,nvopen,iswit
     $              ,we,wexe)
C * ----------------------------------------------------------------

               write(20,*) 'the energies'
               write(20,504) (kw(i)*kw(i)*2.d0*13.6058d0, i=1,nch)
               write(20,*)
               write(20,*) 'after truncation'
               call rkmat(kmat,kw,rmax,l0,rmat,chin,nch)

C * Here we are going to get the T-Matrix and the Cross Sections
C * at the end of the LA region.

C * calculate the t matrix from the k matrix
               call tmatrix(kmat,tmati,tmatr,kmattemp,nch)

C * calculate the cross sections from the t matrix
               call crossec(tmati,tmatr,kw,crossla,symlam,nvib,npw,nch
     $              ,numener,iener)
               write(20,*)
               write(6,*)
               write(6,*) 'Now carry out the Propagation of the R-Mat'
               write(20,*)'Now carry out the Propagation of the R-Mat'
               write(6,*)

C * Call the R-Matrix propagator for the 
C * asymptotic propagation

C * 7-30-04 A. Feldt fix - change to pass vibmax instead of nvibin so
C *         that further usage gets the appropriate long range elements
C *         This fixes the dimensioning error that had obtained if
C *         nvibin < vibmax
               call rprop(rmat,rmax,rasym,rmatstep,kw,l0,alpha0,alpha2
     $              ,q,fmat,nlamdain,nlams,chin,vibmax ,nch,npwin)

C * convert the r-matrix to the K-matrix
               call rkmat(kmat,kw,rasym,l0,rmat,chin,nch)
               rend = rasym
            else

C * else calculate the K-matrix directly from the linear Algebraic
C * solution wavefunction
               call kmatrix(asympsi,kmat,kw,g1,g2,nch,npts)
               rend = rmax
            endif

C * now perform Born r-closure, if we so choose - added by A. Feldt
C * 8-11-04 (this required adding the variable 'rend' and setting
C * it appropriately above, too).  The choice of doing r-closure is
C * made by a variable in the very first line of input to the code.
            if (rclosure .eq. 1) then
               call addbornk(kmat,kw,rend,l0,chin,nch,alpha0,alpha2,q
     $              ,vibmax,fmat,npwin,nlamdain)
            endif

C * calculate the eigenphase sums
            call eigfas(kmat,kmattemp,eigshft,eigsum,eigvec,nch)

C * calculate the t matrix from the k matrix
            call tmatrix(kmat,tmati,tmatr,kmattemp,nch)

C * calculate the cross sections from the t matrix
            call crossec(tmati,tmatr,kw,cross,symlam,nvib,npw,nch
     $           ,numener,iener)

C * the routine OUTDATA takes care of most of the output
C * -- an effort to localize the write statements
C * -- probably a failure
            call outdata(kmat,tmati,tmatr,cross,crossla,kw,asympsi
     $           ,eigshft,eigvec,v,energy(iener),step,nstep,nreg,nch
     $           ,nchsr,npts,l0,symlam,nvib,npw,iener,numener)
         end do

C * We are now done with this symmetry so close the kernel file and we 
C * can go on to the next one
         close(exunit)
      end do

C * print wavefunction
      nf = 100
      do ii = 1, nvib
         do jj = 1, nvib
            nf = nf + 1
            do ipts = 1, nptsex
               write(nf,'(1x,f10.5,2x,1pe12.4)') rgs(ipts),
     &               psi(ii,jj,ipts)
            enddo
         enddo
      enddo

C * here we want to create some tables from the results of all the runs
C * this will eventually be turned into a separate routine
      write(50,*)'Cross Sections at the LA/R-Mat Boundary'
      call crossprt(crossla,energy,numener,nvib,50)
      write(51,*)'Cross Sections at the asymptotic boundary'
      call crossprt(cross,energy,numener,nvib,51)

 504  format(6(1pe12.4))
      stop
      end

C *-
C * this subroutine prints the output
      subroutine outdata(kmat,tmati,tmatr,cross,crossla,k,asympsi
     $     ,eigshft,eigvec,v,energy,step,nstep,nreg,nch,nchsr,npts,l0
     $     ,symlam,nvib,npw,iener,numener)

      integer npts,nch,nreg,nchsr
      integer l0,symlam,npw,nvib
CARR
      integer chanmax,enermax,symmax,vibmax,nrgmax,potptsmax
      parameter(chanmax=165,enermax=120,symmax=4,vibmax=15,nrgmax=100
     $     ,potptsmax=330)
      double precision eigshft(chanmax),eigvec(chanmax,chanmax)
      double precision kmat(chanmax,chanmax),tmati(chanmax,chanmax)
     $     ,tmatr(chanmax,chanmax)
      double precision cross(enermax,vibmax)
      double precision crossla(enermax,vibmax)
      double precision k(chanmax),energy
      double precision step(nrgmax)
      integer nstep(nrgmax)
      double precision asympsi(chanmax,chanmax),v(chanmax,chanmax
     $     ,potptsmax)

      double precision htoev
      integer i,j

      htoev = 2.d0*13.6058d0
      write(6,*)
      write(6,*) 'the following is the mesh'
      do i = 1,nreg
         write(6,*) nstep(i),step(i)
      end do
      write(6,*)
      write(6,*) 'the wavefn at rmax:'
      do i = 1,nch
         write(6,500) (asympsi(i,j),j=1,nch)
      end do
      write(6,*)
      write(6,*) 'k =  ',(k(i),i=1,nch)
      write(6,*)
      write(6,*) 'the potential at the origin'
      do i = 1,nchsr
         write(6,500) (v(i,j,1),j=1,nchsr)
      end do
      write(6,*)
      write(6,*) 'the potential at R-mat boundary'
      do i = 1,nchsr
         write(6,500) (v(i,j,npts),j=1,nchsr)
      end do
      write(6,*)
      write(6,*) ' the k-matrix:'
      write(70,503)symlam,l0,l0+(npw-1)*2,nch,energy*2.d0,nvib
 503  format(4i5,d15.7,i5)
      do i = 1,nch
         write(6,500) (kmat(i,j),j=1,nch)
         write(70,501)(kmat(i,j),j=1,nch)
      end do
 501  format(4d23.16)
C501  format(4d24.16)
      write(6,*)
      write(6,*) ' the eigenphase shifts are'
      write(6,500)(eigshft(i),i=1,nch)
      write(6,*)
      write(6,*) ' the eigenvectors are'
      do i = 1,nch
         write(6,500)(eigvec(i,j),j=1,nch)
      end do
      write(6,*)
      write(6,*)
      write(6,*) ' the imaginary part of the t-matrix:'
      do i = 1,nch
         write(6,500) (tmati(i,j),j=1,nch)
      end do
      write(6,*)
      write(6,*) ' the double precision part of the t-matrix:'
      do i = 1,nch
         write(6,500) (tmatr(i,j),j=1,nch)
      end do
      write(6,*)
      write(6,*) ' the cross section is at the LA boundary is :'
      write(6,511) energy*htoev, (crossla(iener,j),j=1,nvib)
      write(6,*)
      write(6,*) ' the cross section at the asymptotic region is :'
      write(6,511) energy*htoev, (cross(iener,j),j=1,nvib)
      write(19,*) ' the cross section is:'
      write(19,511) energy*htoev, (cross(iener,j),j=1,nvib)
 500  format(8(1pe10.2))
      write(6,*)
      write(20,*) ' the K-Matrix is:'
      write(20,509)energy
      do i = 1,nch
         write(20,510)(kmat(i,j),j=1,nch)
      end do
 509  format(f9.4)
 510  format(5(1pe14.6))
 511  format(f7.3,5f10.5)
      return
      end

C *-
C * This is the routine that reads in the Unit 5 crap
      Subroutine indata(nbset,numiter,l0,lmax,symlam,nlamdain,
     $     numener,energy,iex,choice,ipot,nreg,nstep,step,nregex,nstepex
     $     ,stepex,nregpot,nsteppot,steppot,npts,nptsex,nptspot,nbound
     $     ,nlproj ,mlproj,l0proj,irmat,rmax,rmatstep,rasym,alpha0
     $     ,alpha2,q ,nlams,rpot,ity,rgs,wtt,nchex,kernlfil,namelen
     $     ,nexdim,nofr,nvibsym,nvibin,infil,nvtrunc,nltrunc,irr
     $     ,nlamloc,nvibx,npwavx,nvibxin,npwavxin,rclosure,iswit,we
     $     ,wexe)

      implicit none

      integer chanmax,itermax,nrgmax,lammax,symmax,enermax
      integer setmax,boundmax,vibmax

      parameter(chanmax=165,itermax=10,nrgmax=100,lammax=21,symmax=4
     $     ,enermax=120,setmax=1,boundmax=6,vibmax=15)
      
      integer potptsmax,ptsmax
      parameter(potptsmax=330,ptsmax=331)

      integer nreg,nregex,nregpot,nptspot,nptsex,npts
      integer nvibin,nofr
      integer nvibxin,npwavxin
      integer npwavx,nvibx
      integer nvibsym
      character*8 infil
      double precision rmax
      double precision rmatstep,rasym
      double precision alpha0(vibmax,vibmax),alpha2(vibmax,vibmax)
      double precision q(vibmax,vibmax) 

      integer ity,nps,ic,jz,ipts,istep
      double precision rs,rfn,rfs,rttp,wttp,dell,del(nrgmax),r1,wx
      double precision rpot(potptsmax),rgs(ptsmax),wtt(ptsmax)
      double precision rnode(ptsmax),wnode(ptsmax)

      integer iswit
      double precision we, wexe

C * NLAMS number of lambdas in the long range region
      integer nlams
      integer irmat

C *      IRR = 1 ==> RR or VIBAV calculation
C * RCLOSURE = 1 ==> do Born r-closure
      integer irr,rclosure
      integer nbset,numiter(setmax)
      integer l0,lmax,symlam

C * EUNITS -- Units for the input energies
C *        = 0 ==> Hartree
C *        = 1 ==> eV
      integer nlamdain,nlamloc,numener,eunits
      double precision energy(enermax) 
      integer iex,choice,ipot

C * Exchange stuff
C * NEXDIM: the number of partial waves used in the CREATION
C *         of the exchange kernel
      integer nbound,nchex,nexdim
      integer nlproj(boundmax),mlproj(boundmax),l0proj(boundmax)

C * KERNLFIL: an array containing the names of the kernel files
C *  NAMELEN: an array containing the lengths of the above filenames
      character*8 kernlfil
      integer namelen
      integer nstep(nrgmax),nstepex(nrgmax),nsteppot(nrgmax)
      integer i,j,ind

C * the junk below has to be the same as on the potential mesh
C * so maybe we shouldn't read it in...
      double precision step(nrgmax),stepex(nrgmax),steppot(nrgmax)

C * truncation parameters for truncating the r-matrix to a small
C * number of channels in the asymptotic region
      integer nvtrunc,nltrunc

C * start the main program
      write(6,*)
      read(5,*) nbset,irmat,rclosure
      write(6,*) 'The Number of sets of iterations is: ',nbset
      write(6,*)
      if(rclosure.eq.1) then
         write(6,*)' Born r-closure will be performed on the K matrix'
      endif
      read(5,*) rmatstep,rasym,nlams
      if(irmat.eq.1) then
         write(6,*)' the R-matrix propagator in the asymptotic region'
         write(6,*)
         write(6,*) 'The R-Matrix is propagated to',rasym
         write(6,*) 'The R-matrix step size is',rmatstep
         write(6,*) 'The number of lamdas in the outer region is',nlams
      endif

C *  NVIBIN: the number of vibrational levels included in the input
C *          potential and the long range moments
C * NVIBSYM: the number of vibrational levels to include in the
C *          scattering calculation for each symmetry

C * NOFR: a formatting option in the wlam code set it to 1 always
C * IRR = 1 ==> RR or VIBAV calculation
      read(5,*) nvibin,nofr,irr
      write(6,*)
      write(6,*)'nvibin :',nvibin
      write(6,*)
      if(irr.eq.1) then
         write(6,*)' This is a RR or VIBAV calc. So the local part'
         write(6,*)' of the potential is from a formatted file'
         write(6,*)
      endif
      read(5,*) nvtrunc,nltrunc
      write(6,*)' Parameters for truncating the R-matrix in the'
      write(6,*) '   asymptotic region'
      write(6,*)'Number of Vib states:',nvtrunc
      write(6,*)'Number of partial waves:',nltrunc
      write(6,*)

C * INFIL: the name of the file containing the (unformatted) input
C *        local potential
      read(5,1500)infil
 1500 format(10a8)
      ind = index(infil,' ')
      if(ind.eq.0) ind = 9
      if(irr.eq.1) then
         write(6,*) ' The name of the (formatted) input local potential'
         write(6,*) ' is: ', infil(1:ind-1)//'.dat'
      else
         write(6,*) ' The name of the (unformatted) input local pot'
C        write(6,*) ' is: ', infil(1:ind-1)//'.vib'
         write(6,*) ' is: ', infil(1:ind-1)
      endif
      write(6,*) 'The Polarizabilities of the target:'
      write(6,*)
      write(6,*) ' Alpha0'
      write(6,*)
      do i = 1,nvibin
         read(5,*) (alpha0(i,j),j=1,nvibin)
         write(6,502) (alpha0(i,j),j=1,nvibin)
      end do
      write(6,*)
      write(6,*) ' Alpha2'
      write(6,*)
      do i = 1,nvibin
         read(5,*) (alpha2(i,j),j=1,nvibin)
         write(6,502) (alpha2(i,j),j=1,nvibin)
      end do
      write(6,*)
      write(6,*) ' Quadrupole Moment'
      write(6,*)
      do i = 1,nvibin
         read(5,*) (q(i,j),j=1,nvibin)
         write(6,502) (q(i,j),j=1,nvibin)
      end do
 502  format(5(1pd12.3))
      write(6,*)
      write(6,*)

      read(5,*) (numiter(i),i=1,nbset)
      write(6,*) 'The number of iterations in each trial:'
      write(6,*) (numiter(i),i=1,nbset)
      write(6,*)
      write(6,*) 'l0, lmax, symlam, vib chans:'

      read(5,*)  l0,lmax,symlam,nvibsym
      write(6,*) l0,lmax,symlam,nvibsym

      write(6,*)
      read(5,*) nlamdain,nlamloc
      write(6,*)'Number of Lamdas in input potential file:',nlamdain
      read(5,*) numener,eunits
      write(6,*)
      write(6,*) numener,' Energies'
      read(5,*) (energy(i),i=1,numener)
      if(eunits.eq.0) then
         write(6,*)' Input energies are in Hartrees'
         write(6,*) (energy(i),i=1,numener)
         write(6,*)
      else if(eunits .eq. 1) then
         write(6,*)' Input Energies are in eV'
         write(6,*) (energy(i),i=1,numener)
         write(6,*)
         do i = 1,numener
            energy(i) = energy(i)/27.2116d0
         end do
      endif

C * the exchange part... iex turns non-local exchange on
      read(5,*) iex
      write(6,*) 'Exchange: (1,2=exch, 0=none)',iex

C * exchange stuff
      write(6,*)
      read(5,*) nbound
      write(6,*)'Number of bound target MOs = ', nbound
      read(5,*) npwavxin
      read(5,*) nvibxin
      write(6,*)
      write(6,*)'Total number of vib chans in input Exch Kernel file:'
      write(6,*) nvibxin
      write(6,*)
      write(6,*)'Total number of Prtl Waves in input Exch Kernel file:'
      write(6,*) npwavxin
      write(6,*)
      read(5,*) npwavx
      read(5,*) nvibx
      write(6,*)'Number of vib chans in Exch Kernel:'
      write(6,*)' for this calculation'
      write(6,*) nvibx
      write(6,*)
      write(6,*)'Number of Prtl Waves in Exch Kernel:'
      write(6,*)' for this calculation'
      write(6,*) npwavx
      write(6,*)

C * NCHEX: the number of exchange channels included in the scattering
C *        calculation
      nexdim = nvibxin * npwavxin
      nchex  =   nvibx * npwavx

      write(6,500)
 500  format('   Orbital    nlproj    mlproj    l0proj')
 501  format(4i10)
      do i = 1, nbound
         read(5,*)      nlproj(i),mlproj(i),l0proj(i)
         write(6,501) i,nlproj(i),mlproj(i),l0proj(i)
      end do

C * now we read in the names of the files containing the exchange kernels
C * clearly, one will be read in for each scattering symmetry (right?)
      read(5,1500) kernlfil
      write(6,*) 'The names of the input kernel files are:'

      namelen = index(kernlfil,' ')-1
      if(namelen .eq. -1) namelen = 8
      write(6,*) kernlfil(1:namelen)

C * the following lines set up the interaction potential
C * what kind of first guess?
      write(6,*)
      read(5,*) choice
      write(6,*) 'Choice for initial guess: ',choice

C * which potential
      read(5,*) ipot
      write(6,*)
      write(6,*) 'Choice for system: ',ipot

C * Integration Mesh:
      read(5,*) nreg,nregpot,nregex
      write(6,*)
      write(6,*) 'There are ',nreg,' LA Regions'
      write(6,*) 'There are ',nregpot,' input potential regions'
      rmax = 0.d0
      r1 = 0.d0
      ipts = 0
      do i = 1,nregpot
         read(5,*) nstep(i),step(i)
         write(6,*) nstep(i),step(i)
         nsteppot(i) = nstep(i)
         steppot(i) = step(i)
         rmax = rmax+nsteppot(i)*steppot(i)

         do istep = 1, nsteppot(i)
            ipts = ipts + 1
            r1 = r1 + steppot(i)
            rpot(ipts) = r1
         enddo

      end do
      write(6,*)
      write(6,*) 'The local potential extends to ',rmax,' bohr'
      nptspot = 0
      do i = 1,nregpot
         nptspot = nptspot + nsteppot(i)
      end do
      if (nptspot .gt. potptsmax) stop 'potptsmax is too small!'
      write(6,*) 'There are ',nptspot,' input potential points'

C * ity = 0 - Trapezoidal pts and wts
C *       1 - Gauss-Legendre pts and wts (copied from Saha)
C *       2 - Gauss-Legendre pts and wts (cal. by gauleg)
C *       
C * nps = No. of pts. in a region
C *  rs = starting radius of region
C * rfn = ending radius of region

      write(6,*)
      WRITE (6,*) 'There are',nregex,' exhchange mesh regions'
      WRITE (6,*)
      WRITE (6,*)'ity nps    rs     rfn '
      ic = 0
      do i = 1, nregex
         READ (5,*) ity, nps, rs, rfn
         WRITE (6,'(i2,i5,f8.3,f8.3)')ity, nps, rs, rfn
         rfs = rfn - rs
         nstepex(i) = nps
         stepex(i) = rfs/nps
         if (ity .eq. 0) then
            del(i) = rfs/nps
            r1 = rs - del(i)
            if (i .eq. 1) then
               r1 = rs
               nps = nps - 1
            endif
            do jz = 1, nps
               dell = del(i)
               ic = ic + 1
               r1 = r1 + dell
               wx = 1.0d0
               if (ic .eq. 1) wx = 0.5d0
               if (jz .eq. 1 .and. i .ne. 1) then
                  wx = 0.5d0
                  dell = del(i) + del(i-1)
               endif
               wtt(ic) = wx*dell
               rgs(ic) = r1
            enddo
         else if (ity .eq. 1) then
            do jz = 1, nps
               ic = ic + 1
               if (nps .gt. 3) then
                  call lgndrx(nps,jz,wttp,rttp)
               else if (nps .le. 3) then
                  call lgndx2(nps,jz,wttp,rttp)
               endif
               wtt(ic) = wttp*rfs
               rgs(ic) = rttp*rfs + rs
            enddo
         else if (ity .eq. 2) then
            call gauleg(rs,rfn,rnode,wnode,nps)
            do jz = 1, nps
               ic = ic + 1
               rgs(ic) = rnode(jz)
               wtt(ic) = wnode(jz)
            enddo
         endif
      end do

      nptsex = ic + 1

C * LA mesh points are the same as the exchange points
      npts = nptsex
      if (nptsex .gt. ptsmax) stop 'ptsmax is too small'

      rgs(nptsex) = rfn
      wtt(nptsex) = 0.0d0
      write(6,*)
      write(6,*) 'Now the LA mesh points are:',npts
      write(6,*) 'Now the abscissas and weights are: (',nptsex
     $     ,' points)'
      do jz = 1, nptsex
         write(6,'(i5,f9.5,f14.10)') jz,rgs(jz),wtt(jz)
      enddo

      write(6,*)
      write(6,*) 'There are ',nptsex,' exchange potential points'

C * the LA max bound is at the exchange bound
      if (rgs(nptsex) .lt. rmax) rmax = rgs(nptsex)

      if (ity .eq. 0 .and. nregex .gt. nregpot) then
         write(6,*)
         write(6,*) "CAUSTION: now you input trapezoidal grids,"
         write(6,*) "but your exchange regions are greater than"
         write(6,*) "the input potential regions!"
         write(6,*)
      endif
      if (abs(rgs(nptsex)-rpot(nptspot)) .gt. 1.d-3) then
         write(6,*)
         write(6,*) "CAUSTION: now your rmax of exchange is greater"
         write(6,*) "than rmax of input potential!"
         write(6,*) " rmax of exchange = ",rgs(nptsex)
         write(6,*) "rmax of potential = ",rpot(nptspot)
         write(6,*)
      endif

      read(5,*) iswit, we, wexe
      write(6,*)
      write(6,*) 'The vibrational spectroscopic constants are: '
      if (iswit .eq. 1) then
         write(6,*) "  (cm-1)"
      else if (iswit .eq. 0) then
         write(6,*) "  (Hartree)"
      endif
      write(6,*) ' we = ', we, ' wexe = ', wexe
      write(6,*)
      return
      end

C *-
C * This is a subroutine to print out the cross sections in
C * a neat table. (HAHAHAHAHAHAHAHAHAHHAHAHA)
      subroutine crossprt(cross,energy,nener,nvib,ounit)
      implicit none
      integer nener,nvib

      integer enermax, vibmax
      parameter(enermax=120, vibmax=15)
      double precision cross(enermax, vibmax)
      double precision energy(enermax)

      integer ounit
      integer i,j,k
      double precision htoev

C * first get the sum of the symmetries for each vib exit channel
      htoev = 2.d0*13.6058d0
      do i = 1,nvib
         write(ounit,*)
         do j = 1,nener
            write(ounit,500) htoev*energy(j), cross(j,i)
         end do
      end do
 500  format(f7.3, f13.8)
      return
      end

C *-
C * This routine calculates the cross sections from the T-Matrices
C * It currently gives for each channel
      subroutine crossec(tmati,tmatr,k,cross,symlam,nvib,npw,nch,nener
     $     ,iener)

      implicit none
      integer nch,symlam,nvib,npw
      integer nener, iener
CARR
      integer chanmax,enermax,vibmax
      parameter(chanmax=165,enermax=120,vibmax=15)
      double precision tmati(chanmax,chanmax),tmatr(chanmax,chanmax)
      double precision k(chanmax),pi
      double precision cross(enermax,vibmax)

      integer i,j,l

      pi = 4.d0*datan(1.d0)
      do i = 1,nvib
         do j = 1,npw
            do l = 1,npw
               cross(iener,i) = cross(iener,i) +
     $              (tmati(j,(i-1)*npw+l)*tmati(j,(i-1)*npw+l) +
     $               tmatr(j,(i-1)*npw+l)*tmatr(j,(i-1)*npw+l))
            end do
         end do

         if(symlam.ne.0) cross(iener,i) = cross(iener,i)*2.d0
         cross(iener,i) = pi*cross(iener,i)/(k(1)*k(1))
         write(6,*)' in the routine Crossec'
         write(6,*) cross(iener,i)

      end do
      return
      end

C *-
      subroutine expcoefs(psifin,xhifin,acoefs,g1,step,nstep,nreg,nch
     $     ,npts,niter,wtt,ity)

      implicit none

      integer nch,npts,niter
      integer nreg,ireg,istep,ipts
      integer nstep(nreg)
      double precision step(nreg)
      double precision wtt(npts)
      integer ity
CARR
      integer chanmax,ptsmax,itermax
      parameter (chanmax=165,ptsmax=331,itermax=10)
c$$$      double precision psifin(niter,nch,nch,npts)
c$$$      double precision xhifin(niter,nch,nch,npts)
c$$$      double precision acoefs(nch*niter,nch)
c$$$      double precision g1(nch,npts)
c$$$      double precision psixhim(nch*niter,nch*niter)
c$$$      double precision psipsim(nch*niter,nch*niter)
      double precision psifin(itermax,chanmax,chanmax,ptsmax)
      double precision xhifin(itermax,chanmax,chanmax,ptsmax)
      double precision acoefs(chanmax*itermax,chanmax)
      double precision g1(chanmax,ptsmax)
      double precision psixhim(chanmax*itermax,chanmax*itermax)
      double precision psipsim(chanmax*itermax,chanmax*itermax)

      double precision psixhi, psib
      double precision psipsi


C * to carry out LU decomposition and create an inverse
CARR
c$$$      integer indx(nch*niter)
c$$$      double precision temp(nch*niter,nch*niter)
c$$$      double precision pxinv(nch*niter,nch*niter)
c$$$      double precision psibm(nch*niter,nch)
c$$$      double precision unitary(nch*niter,nch*niter)
      integer indx(chanmax*itermax)
      double precision temp(chanmax*itermax,chanmax*itermax)
      double precision pxinv(chanmax*itermax,chanmax*itermax)
      double precision psibm(chanmax*itermax,chanmax)
      double precision unitary(chanmax*itermax,chanmax*itermax)

      integer i,j,k,l
      integer ich,m

C * if the matrices get too big then you have to allow work space for 
C * the IMSL routine... the error message looks something like this:
C **   TERMINAL ERROR from DLINRG.  Insufficient workspace for current
C **            allocation(s). Correct by calling IWKIN from main 
C **            program with the three following statements: (REGARDLESS
C **            OF PRECISION)
C **   TERMINAL ERROR from DLINRG.  The workspace is based on N, where 
C **   N = 72.
C * the huge number of do loops below is used to create
C * a matrix of products of the normalized z vectors and the 
C * unnormalized x vectors... so we get a set of coeffs
C * nch*niter x nch*niter corresponding to all possible
C * projections of a z onto an x

C * surely it is possible to get these projections
C * in the orthonormalization procedure
      do i = 1,niter
         do j = 1,nch
            do k = 1,niter
               do l = 1,nch
                  psixhi = 0.0d0
                  psipsi = 0.0d0
                  do m = 1,nch
                     if (ity .eq. 0) then
                        ipts = 0
                        do ireg = 1,nreg
                           do istep = 1,nstep(ireg)
                              ipts = ipts+1
                              psipsi = psipsi+psifin(i,j,m,ipts)
     $                             *psifin(k,l,m,ipts)*step(ireg)
                              psixhi = psixhi+psifin(i,j,m,ipts)
     $                             *xhifin(k,l,m,ipts)*step(ireg)
                           end do
                        end do
                     else if (ity .eq. 1 .or. ity .eq. 2) then
                        do ipts = 1, npts
                           psipsi = psipsi+psifin(i,j,m,ipts)*psifin(k,l
     $                          ,m,ipts)*wtt(ipts)
                           psixhi = psixhi+psifin(i,j,m,ipts)*xhifin(k,l
     $                          ,m,ipts)*wtt(ipts)
                        enddo
                     endif
                  end do
                  psipsim(nch*(i-1)+j,nch*(k-1)+l) = psipsi
                  psixhim(nch*(i-1)+j,nch*(k-1)+l) = psixhi
               end do
            end do
         end do
      end do

c$$$      write(23,*)
c$$$      write(23,*) ' the overlap matrix:'
c$$$      write(23,*)
c$$$      write(23,*)
      write(6,*)
      write(6,*)

C * in the do loops below we project the b vectors (in ax=b) onto
C * the orthonormal set {z} and store the result in psibm
C * so psibm is nch*niter x nch   in size
      do ich = 1,nch
         do i = 1,niter
            do j = 1,nch
               psib = 0.0d0
               if (ity .eq. 0) then
                  ipts = 0
                  do ireg = 1,nreg
                     do istep = 1,nstep(ireg)
                        ipts = ipts+1
                        psib = psib+psifin(i,j,ich,ipts)*g1(ich,ipts)*
     $                       step(ireg)
                     end do
                  end do
               else if (ity .eq. 1 .or. ity .eq. 2) then
                  do ipts = 1, npts
                     psib = psib+psifin(i,j,ich,ipts)*g1(ich,ipts)
     $                    *wtt(ipts)
                  end do
               endif
               psibm(nch*(i-1)+j,ich) = psib
            end do
         end do
      end do

C * here we solve the small (matrix) problem ax=b
C * for the expansion coeffs of the solution (to the big AX=B)
C * in the basis of iterates

C * store the projections with the right coefficients
      do i = 1,(niter-1)*nch
         do j = 1,(niter-1)*nch
            psixhim(i,j) = -psixhim(i,j+nch)
         end do
         psixhim(i,i) = psixhim(i,i)+1.0d0
      end do

      call nrinverse(psixhim,pxinv,temp,nch*(niter-1),chanmax*itermax
     $     ,indx)

      call mymatmul(psixhim,pxinv,unitary,nch*(niter-1),chanmax
     $     *itermax)

C * ...from above and matrix multiplying that result by 
C * b (the projections <zi|bi>) as follows:
      do ich = 1,nch
         do i = 1,nch*(niter-1)
            acoefs(i,ich) = 0.d0
            do j = 1,nch*(niter-1)
               acoefs(i,ich) = acoefs(i,ich)+pxinv(i,j)*psibm(j,ich)
            end do
         end do
      end do
      return
      end

C *-
      subroutine frstgues(a,psi,choice,nch,npts)
      implicit none

      integer nch,npts
CARR
      integer chanmax, ptsmax
      parameter (chanmax=165,ptsmax=331)
c$$$      double precision a(nch,npts)
c$$$      double precision psi(nch,nch,npts)
      double precision a(chanmax,ptsmax)
      double precision psi(chanmax,chanmax,ptsmax)

      double precision r
      integer i,j,k
      integer choice

      if(choice.eq.6) then
         do j = 1,nch
            do i = 1,npts
               read(10,*) r, (psi(j,k,i),k = 1,nch)
            end do
         end do
      else
         do j = 1,nch
            do i = 1,npts
               do k = 1,nch
                  psi(j,k,i) = 0.d0
               end do
               psi(j,j,i) = a(j,i)
            end do
         end do
      endif

      return
      end
C *-/

C *-
      subroutine imats(i1,i2,q1,q2,nch,npts)
      implicit none

      integer nch,npts
CARR
      integer chanmax, ptsmax
      parameter(chanmax=165,ptsmax=331)
c$$$      double precision i1(nch,nch,npts),i2(nch,nch,npts)
c$$$      double precision q1(nch,nch,npts),q2(nch,nch,npts)
      double precision i1(chanmax,chanmax,ptsmax)
      double precision i2(chanmax,chanmax,ptsmax)
      double precision q1(chanmax,chanmax,ptsmax)
      double precision q2(chanmax,chanmax,ptsmax)

      integer i,ip,j,k

      do i = 1,nch
         do j = 1,nch
            i1(j,i,1) = q1(j,i,1)
            i2(j,i,npts) = q2(j,i,npts)
         end do
      end do
      do i = 1,npts - 1
         ip = npts - i
         do j = 1,nch
            do k = 1,nch
               i1(j,k,i+1) = i1(j,k,i) + q1(j,k,i+1)
               i2(j,k,ip) = i2(j,k,ip+1) + q2(j,k,ip)
            end do
         end do
      end do
      
      return
      end

C *-
C * calculate the K matrix from the R matrix
C * use, say, the formula for the K-matrix from Lane's review
C * psi ~ [sin(kr) + cos(kr) * (k**-.5) * K * (k**.5)] * A
C * psi' ~ [k * cos(kr) - sin(kr) * (k**.5) * K * (k**.5)] * A
C * and 
C * Rmat = psi * ([psi']**-1)
C * then solve for K
C * Note that we had better be in the asymptotic region 'cause
C * were are using sins and cosines

      Subroutine rkmat(kmatrix,k,r,l0,rmatn,chin,nch)
      implicit none
      integer i,j,ii

      integer nch
CARR
      integer chanmax
      parameter(chanmax=165)
c$$$      double precision kmatrix(nch,nch),k(nch),r
c$$$      double precision rmatn(nch,nch)
c$$$      integer chin(nch,2)
      double precision kmatrix(chanmax,chanmax),k(chanmax),r
      double precision rmatn(chanmax,chanmax)
      integer chin(chanmax,2)

      double precision rbes,rbesp,rneu,rneup
      double precision dummy(nch,nch)
      double precision numer(nch,nch)
      double precision denom(nch,nch)

C * to carry out LU decomposition and create an inverse
      integer indx(nch)
      double precision temp(nch,nch)

      integer l0,l,lp

C * rmatn now contains the untransformed ratio of the wfn to
C * the deriv of the wfn
      write(6,*)
      write(6,*) 'This is the R-Matrix:'
      write(6,*)
      do i = 1,nch
         write(6,501) (rmatn(i,j), j=1,nch)
      end do
 501  format(6(1pe12.3))

      write(6,*)
      write(6,*) 'k=', (k(i), i=1,nch)
      write(6,*)
      write(6,*) 'r=', r
      write(6,*) 'l0=', l0
      write(6,*) 'end of R-Matrix stuff'

      do i = 1,nch
         l = 2*(chin(i,2)-1) + l0
         do j = 1, nch
            lp = 2*(chin(j,2)-1) + l0
            numer(i,j) = rmatn(i,j)*rbesp(k(j),r,lp)
            dummy(i,j) = rmatn(i,j)*rneup(k(j),r,lp)
         end do
         numer(i,i) = numer(i,i) - rbes(k(i),r,l)
         dummy(i,i) = dummy(i,i) - rneu(k(i),r,l)
      end do

      call nrinverse(dummy,denom,temp,nch,nch,indx)

      do i = 1,nch
         do j = 1,nch
            kmatrix(i,j) = 0.d0
            do ii = 1,nch
               kmatrix(i,j) = kmatrix(i,j) + denom(i,ii)*numer(ii,j)
            end do
            kmatrix(i,j) = kmatrix(i,j)*dsqrt(k(i)/k(j))
         end do
      end do
      write(6,*) 'K matrix at the asymptotic boundary'
      do i = 1,nch
         write(6,500)  (kmatrix(i,j),j = 1,nch)
         write(20,500) (kmatrix(i,j),j = 1,nch)
      end do
 500  format(6(1pe12.4))
      return
      end

C *-
      subroutine kmatrix(asympsi,kmat,k,g1,g2,nch,npts)
      implicit none

      integer nch,npts
CARR
      integer chanmax,ptsmax
      parameter (chanmax=165,ptsmax=331)
c$$$      double precision asympsi(nch,nch)
c$$$      double precision kmat(nch,nch)
c$$$      double precision k(nch),g1(nch,npts),g2(nch,npts)
      double precision asympsi(chanmax,chanmax) 
      double precision kmat(chanmax,chanmax)
      double precision k(chanmax),g1(chanmax,ptsmax),g2(chanmax,ptsmax)

      integer i,j

      write(6,*) 'In Routine kmatrix'
      write(6,*) 'nch=',nch,' npts=',npts
      write(6,*) 'K matrix'

      do i = 1,nch
         do j = 1,nch
            if(i.eq.j) then
               kmat(i,j) = (asympsi(i,j)-g1(i,npts))/(k(i)*g2(i,npts))
            else
               kmat(i,j) = asympsi(i,j)/(dsqrt(k(i)*k(j))*g2(i,npts))
            endif
         end do
         write(6,500) (kmat(i,j),j=1,nch)
      end do

      return
 500  format(6(1pe12.4))
      end

C *-
C * This routine uses the K-Matrix to calculate the eigenphase shifts 
C * and eigenphase sums.  The K-matrix is assumed to be symmetric
C * so that the eigenvalues are real.

      Subroutine eigfas(kmat,kmattemp,eigshft,eigsum,eigvec,nch)
      implicit none
      integer i,j,nch,nrot
CARR
      integer chanmax
      parameter(chanmax=165)
      double precision kmattemp(chanmax,chanmax)
      double precision kmat(chanmax,chanmax),eigshft(chanmax)
      double precision eigsum,eigvec(chanmax,chanmax)

      do i = 1,nch
         do j = 1,nch
            kmattemp(i,j) = kmat(i,j)
         end do
      end do
C     call jacobi(kmattemp,nch,nch,eigshft,eigvec,nrot)
      call jacobi(kmattemp,nch,chanmax,eigshft,eigvec,nrot)
      do i = 1,nch
         eigshft(i) = datan(eigshft(i))
      end do
      eigsum = 0.d0
      do i = 1,nch
         eigsum = eigsum + eigshft(i)
      end do
      return
      end

C *-
C * calculate the Clebsch Gordon elements to be used
C * when calculating the exchange potential

      Subroutine clebexch(clebx1,clebx2,clebx3,clebx4,clebx5,symlam,l0
     $     ,nlproj,mlproj,l0proj,nbound,nexdim)
      implicit none
      integer nexdim,nbound
      integer symlam,l0
      integer nlproj(nbound),l0proj(nbound),mlproj(nbound)

      integer boundmax,chanmax,exlmmax,nlmomax
      parameter(chanmax=165,exlmmax=60,boundmax=6,nlmomax=16)
      double precision clebx1(boundmax,chanmax,exlmmax,nlmomax)
      double precision clebx2(boundmax,chanmax,exlmmax,nlmomax)
      double precision clebx3(boundmax,chanmax,exlmmax,nlmomax)
      double precision clebx4(boundmax,chanmax,exlmmax,nlmomax)
      double precision clebx5(boundmax,chanmax,exlmmax,nlmomax)

      integer m,mi,ibnd,i,j,k,l,lam,lp,lamlb,lamub

C * connect to the Cleb Gordon routines:
      double precision faclog(500)
      double precision rac,u9
      integer ia,ib,ic,id,ie,if,l9(9),kselc
      common/faclg/faclog,rac,u9,ia,ib,ic,id,ie,if,l9,kselc

      m = symlam
      do ibnd = 1,nbound
         mi = mlproj(ibnd)
         do i = 1,nexdim
            l = 2*(i-1)+l0
            do k = 1, nlproj(ibnd)
               lp = 2*(k-1) + l0proj(ibnd)
               lamlb = iabs(l - lp)
               lamub = l + lp
               if (lamub .gt. exlmmax) then
                  write(0,*) "exlmmax = ",exlmmax ," and should be >/="
     $                 ,lamub
                  stop ' INCREASE exlmmax in source!'
               endif
               do lam = lamlb, lamub
                  clebx1(ibnd,i,lam,k) = 0.0d0
                  clebx2(ibnd,i,lam,k) = 0.0d0
                  clebx3(ibnd,i,lam,k) = 0.0d0
                  clebx4(ibnd,i,lam,k) = 0.0d0
                  clebx5(ibnd,i,lam,k) = 0.0d0

                  if (mod(l+lp+lam,2) .eq. 0) then
                     call cleb(l,lam,lp,-m,m-mi,-mi)
                     clebx1(ibnd,i,lam,k) = rac
                     if (mi .eq. 1) then
                        call cleb(l,lam,lp,-m,m+mi,mi)
                        clebx4(ibnd,i,lam,k) = rac
                     else if (mi .eq. 0) then
                        clebx4(ibnd,i,lam,k) = 0.0
                     endif
                     call cleb(l,lam,lp,0,0,0)
                     clebx2(ibnd,i,lam,k) = rac
                     call cleb(l,lam,lp,m,mi-m,mi)
                     clebx3(ibnd,i,lam,k) = rac
                     if (mi .eq. 1) then
                        call cleb(l,lam,lp,m,-mi-m,-mi)
                        clebx5(ibnd,i,lam,k) = rac
                     else if (mi .eq. 0) then
                        clebx5(ibnd,i,lam,k) = 0.0
                     endif
                  endif 
               end do
            end do
         end do
      end do
      return
      end

C *-
C * this routine calculates the exchange kernel an writes it to
C * a file.  It will be used later to calculate the
C * exchange potential

      Subroutine Kerncalc(sphrj,kernel,clebx1,clebx2,clebx3,clebx4
     $     ,clebx5,nlproj,l0proj,exunit,l0,nexdim,nbound,rgs,nptsex)
      implicit none
      integer nptsex
      integer exunit
      integer nexdim,nbound
      integer nlproj(nbound),l0proj(nbound)
      double precision rgs(nptsex)
      double precision kernel(nexdim,nexdim,nptsex)
      double precision sphrj(nbound,nexdim,nptsex)
      integer l0
      integer i,j,k,il,ilp
      integer l,lp,l2p,l3p
      integer ipts,jpts,ibnd
      double precision ralfa,rbeta,rless,rgrtr,rfact

      integer boundmax,chanmax,exlmmax,nlmomax
      parameter(chanmax=165,exlmmax=60,boundmax=6,nlmomax=16)
      double precision clebx1(boundmax,chanmax,exlmmax,nlmomax)
      double precision clebx2(boundmax,chanmax,exlmmax,nlmomax)
      double precision clebx3(boundmax,chanmax,exlmmax,nlmomax)
      double precision clebx4(boundmax,chanmax,exlmmax,nlmomax)
      double precision clebx5(boundmax,chanmax,exlmmax,nlmomax)

      double precision numer
      double precision lamterm,l2term,l3term,ibndterm
      double precision clebtemp
      integer lam,lamlb,lamub

C * Open the file that we will write the exchange kernel to
C *      open(unit=exunit,file=kernfil,form='unformatted')

      do ipts = 1, nptsex
         rbeta = rgs(ipts)
         do jpts = 1, nptsex
            ralfa = rgs(jpts)
            if(ralfa.gt.rbeta)then
               rgrtr = ralfa
               rless = rbeta
            else
               rgrtr = rbeta
               rless = ralfa
            endif
            do il = 1,nexdim
               l = 2*(il-1)+l0
               do ilp = 1, nexdim
                  lp = 2*(ilp-1) + l0
                  numer = dsqrt((2.d0*l+1.d0)*(2.d0*lp+1.d0))
                  
                  ibndterm = 0.0d0
                  do ibnd = 1, nbound
                     
                     l2term = 0.d0
                     do i = 1,nlproj(ibnd)
                        l2p = 2*(i-1) + l0proj(ibnd)

                        l3term = 0.d0
                        do j = 1,nlproj(ibnd)
                           l3p = 2*(j-1)+l0proj(ibnd)

                           lamterm = 0.d0
                           lamlb = max(iabs(l-l2p), iabs(lp-l3p))
                           lamub = min(l+l2p, lp+l3p)
                           do lam = lamlb, lamub
                              if (mod(l+l2p+lam,2) .eq. 0 .and.
     $                             mod(lp+l3p+lam,2) .eq. 0) then 
                                 rfact = rless**(lam)/rgrtr**(lam+1)
                                 clebtemp = (clebx1(ibnd,il,lam,i)
     $                                *clebx3(ibnd,ilp,lam,j)
     $                                +clebx4(ibnd,il,lam,i)
     $                                *clebx5(ibnd,ilp,lam,j))
     $                                *clebx2(ibnd,il,lam,i)
     $                                *clebx2(ibnd,ilp,lam,j)

                                 lamterm = lamterm + clebtemp*rfact
                              endif
                           end do

                           l3term = l3term+lamterm*sphrj(ibnd,j,jpts)
     $                          /dsqrt(2.d0*l3p+1.d0)
                        end do
                        l2term = l2term+l3term*sphrj(ibnd,i,ipts)/
     $                       dsqrt(2.d0*l2p+1.d0)
                     end do
                     ibndterm = ibndterm + l2term
                  end do
                  kernel(il,ilp,jpts) = ibndterm*numer
               end do
            end do
         end do
         write(exunit) kernel
      end do

      return
      end
C *-/

C *-
C * This subrpoutine checks to see if the exchange kernel
C * size agrees with the size of the scattering matrix
C * (i.e. so they have the same number of vib channels
C * and partial waves)
C * If there is a disagreement, this routine resizes the exchange
C * kernel matrix then calls the KERNPSI matrix to carry out the
C * multiplication.

C * I am not sure this is a) efficient or b) an appropriate
C * way for allowing convergence studies to be done from
C * one set of exchange kernels

      subroutine getkpsi(kernel,psi,kpsi,dum1,step,nstep,nvib,nvibx
     $     ,nvibxin,npw,npwavx,npwavxin,nch,nchex,nexdim,npts,nreg
     $     ,nregex,nptsex,exunit,wtt,ity)
      implicit none
      integer nexdim,nch,nchex
      integer nptsex,nregex
      integer nvib,nvibx,nvibxin
      integer npw,npwavx,npwavxin
      integer npts,nreg
      double precision wtt(npts)
      integer ity

C * EXUNIT: the unit to read the Kernel from
      integer exunit

      double precision kernel(nexdim,nexdim,nptsex)
      double precision dum1(nch,nch,nptsex)
CARR
      integer chanmax,ptsmax
      parameter(chanmax=165,ptsmax=331)
c$$$      double precision kpsi(nch,nch)
c$$$      double precision psi(nch,nch,npts)
      double precision kpsi(chanmax,chanmax)
      double precision psi(chanmax,chanmax,ptsmax)

      double precision step(nreg)
      integer nstep(nreg)

C * note the following definitionsfrom the main routine
C *      nexdim = nvibxin*npwavxin
C *       nchex = nvibx*npwavx
C *         nch = nvib*npw
      read(exunit) kernel

      if(nvib.eq.nvibx.and.nvibx.eq.nvibxin.and.
     $     npw.eq.npwavx.and.npwavx.eq.npwavxin) then
         call kernpsi(kernel,psi,kpsi,step,nstep,nch,npts,nreg,nregex
     $        ,nptsex,nchex,nexdim,wtt,ity)
      else

C * THIS LOOKS WRONG----TAKE ANOTHER LOOK AT IT
c * IT IS FOR CHANGING DIMENSIONS OF THE KERNEL
         call movemat(kernel,dum1,nvibxin,npwavxin,nvibx,npwavx
     $        ,nptsex)
         call movemat(dum1,kernel,nvibx,npwavx,nvib,npw,nptsex)
         call kernpsi(kernel,psi,kpsi,step,nstep,nch,npts,nreg,nregex
     $        , nptsex,nch,nch,wtt,ity)
      endif
      return
      end

C *-
C * This is a subroutine to move the relevant parts of the
C * matrix A into the matrix B.  This routine is for
C * doing calculations that use different numbers of
C * exchange channels without having to re-calculate
C * the kernel.

      subroutine movemat(a,b,na1,na2,nb1,nb2,npts)
      implicit none
      integer na1,na2,nb1,nb2,npts
      double precision a(na1*na2,na1*na2,npts)
      double precision b(nb1*nb2,nb1*nb2,npts)
      integer ipts,i1,i2,j1,j2,ia,ib,ja,jb
      do ipts=1,npts
         ib=0
         do i1=1,nb1
            do i2=1,nb2
               ib=ib+1
               jb=0
               do j1=1,nb1
                  do j2=1,nb2
                     jb=jb+1
                     if(i1.le.na1.and.i2.le.na2.and.
     $                    j1.le.na1.and.j2.le.na2) then
                        ia=(i1-1)*na2+i2
                        ja=(j1-1)*na2+j2
                        b(ib,jb,ipts)=a(ia,ja,ipts)
                     else
                        b(ib,jb,ipts)=0.d0
                     endif
                  end do
               end do
            end do
         end do
      end do
      return
      end

C *-
C * This routine calculates the product of the kernel matrix with the
C * wavefunction matrix... the result is integrated over r and 
C * inserted in the Q matrix...

      subroutine Kernpsi(kernel,psi,kpsi,step,nstep,nch,npts,nreg
     $     ,nregex,nptsex,nchex,nexdim,wtt,ity)
      implicit none

C *  NCHEX: the number of exchange channels included in the scattering
C *         calculation
      integer nregex,nptsex
      integer nchex,nexdim
      integer nch,npts,nreg
      integer i,j,m,istep,ipts,ireg,nstep(nreg)
      double precision kernel(nexdim,nexdim,nptsex)
      double precision step(nreg),temp
      double precision wtt(npts)
      integer ity
CARR
      integer chanmax,ptsmax
      parameter(chanmax=165,ptsmax=331)
c$$$      double precision psi(nch,nch,npts),kpsi(nch,nch)
      double precision psi(chanmax,chanmax,ptsmax)
      double precision kpsi(chanmax,chanmax)

      do i=1,nch
         do j=1,nch
            kpsi(i,j)=0.d0
         end do
      end do

C * Gauss-Legendre quadrature use different r points
      if (ity .eq. 0) then
         do i=1,nchex
            do j=1,nchex
               kpsi(i,j)=0.d0
               do m=1,nchex
                  ipts=0
                  temp=kernel(i,m,1)*psi(m,j,1)*step(1)/2.d0
                  do ireg=1,nregex
                     do istep=1,nstep(ireg)-1
                        ipts=ipts+1
C * Hao Feng 
                        if (ipts .ne. 1) then
                           temp=temp+kernel(i,m,ipts)*psi(m,j,ipts)
     $                          *step(ireg)
                        endif
C *-/
                     end do
                     ipts=ipts+1
                     if(ireg.lt.nregex) then
                        temp=temp+kernel(i,m,ipts)*psi(m,j,ipts)*
     $                       (step(ireg)+step(ireg+1))/2.d0
                     else
                        temp=temp+kernel(i,m,ipts)*psi(m,j,ipts)*
     $                       step(ireg)/2.d0
                     endif
                  end do
                  kpsi(i,j)=kpsi(i,j)+temp
               end do
            end do
         end do

      else if (ity .eq. 1 .or. ity .eq. 2) then
         do i = 1, nchex
            do j = 1, nchex
               kpsi(i,j) = 0.d0
               do m = 1, nchex
                  temp = 0.d0
                  do ipts  =  1, nptsex
                     temp = temp + kernel(i,m,ipts)*psi(m,j,ipts)
     $                    *wtt(ipts)
                  end do
                  kpsi(i,j) = kpsi(i,j) + temp
               end do
            end do
         end do
      endif

      do i=nchex+1,nch
         do j=1,nch
            kpsi(i,j)=0.d0
         end do
      end do

      do i=1,nch
         do j=nchex+1,nch
            kpsi(i,j)=0.d0
         end do
      end do
      return
      end

C *-
C * this routine reads in the spherical projections of the Target MO's
C * as from the output of the ALAM code.

      Subroutine sphjread(sphrj,nlproj,nbound,nptspot,nptsex,rgs
     $     ,nexdim)
      implicit none
      integer nexdim,nptspot,nbound,nptsex
      integer i,j,k,nlproj(nbound)
      double precision sphrj(nbound,nexdim,nptsex),r(nptspot)
      double precision sj(nptspot),sy(nptsex)
      double precision rgs(nptsex),csplin(nptspot)

C * read in the r mesh points from the sphj file
C * go to npts+1 because sphj file has point at r=0.0
      write(6,*) 'in sphrj'
      write(6,*) 'nbound    nlproj(1)         nexdim         nptspot'
      write(6,*) nbound,nlproj(1),nexdim,nptspot

      read(3) r(1)
      write (60,*) r(1)
      do i = 1, nptspot 
         read(3) r(i)
         WRITE (60,*) r(i)
      end do

C * read in the spherical projs from the ALAM output file
      do i = 1, nbound
         do j = 1, nlproj(i)

C * as above, we read in the points at r=0.0 but we don't use it since
C * we don't have it on the potential mesh
            read(3) sj(1)
            do k = 1, nptspot
               read(3) sj(k)
               write(61,*) r(k),sj(k)
            end do
C * interpolate sphrj(nptsexch) through original MO coefficients sj(npts)
            call spline(r,sj,nptspot,1.d30,1.d30,csplin)
            do k = 1, nptsex
               call splint(r,sj,csplin,nptspot,rgs(k),sy(k))
               sphrj(i,j,k) = sy(k)
               write(62,*) rgs(k),sphrj(i,j,k)
            enddo
         end do
      end do

 500  format(6(1pe12.3))
      return
      end

C *-
      subroutine mymatmul(m1,m2,m3,nsize,ndim)
      implicit none
      integer nsize,ndim
      double precision m1(ndim,ndim),m2(ndim,ndim),m3(ndim,ndim)
      integer i,j,k
      do i = 1,nsize
         do j = 1,nsize
            m3(i,j) = 0.d0
            do k = 1,nsize
               m3(i,j) = m3(i,j)+m1(i,k)*m2(k,j)
            end do
         end do
      end do
      return
      end

C *-
      subroutine matprdm(m1,m2,m3,nch,npts)
      implicit none
      integer nch,npts
CARR
      integer chanmax,ptsmax
      parameter(chanmax=165,ptsmax=331)
c$$$      double precision m1(nch,nch,npts),m2(nch,nch,npts)
c$$$      double precision m3(nch,nch,npts)
      double precision m1(chanmax,chanmax,ptsmax)
      double precision m2(chanmax,chanmax,ptsmax)
      double precision m3(chanmax,chanmax,ptsmax)


      integer i,j,k,l

      do i = 1,npts
         do j = 1,nch
            do k = 1,nch
               m3(j,k,i) = 0.d0
               do l = 1,nch
                  m3(j,k,i) = m3(j,k,i)+m1(j,l,i)*m2(l,k,i)
               end do
            end do
         end do
      end do
      return
      end

C *-
      subroutine gramschm(xhi,xhifin,psifin,psi,numv,step,nstep,nreg,nch
     $     ,npts,wtt,niter,ity)
      implicit none

C * this subroutine takes nch vectors which are in xhi 
C * and orthonormalizes them with the set in psi 
C * and then adds them to the set psi and increases the
C * number of vectors (numv) to numv + nch
      integer nch,npts,niter,ity,nreg
      integer nstep(nreg)
      double precision step(nreg)
      integer istep,ireg,ipts
      double precision wtt(npts)
CARR
      integer chanmax,ptsmax,itermax
      parameter (chanmax=165,ptsmax=331,itermax=10) 
c$$$      double precision xhi(nch,nch,npts),psifin(niter,nch,nch,npts)
c$$$      double precision xhifin(niter,nch,nch,npts)
c$$$      double precision psi(nch,nch,npts)
      double precision xhi(chanmax,chanmax,ptsmax)
      double precision psifin(itermax,chanmax,chanmax,ptsmax)
      double precision xhifin(itermax,chanmax,chanmax,ptsmax)
      double precision psi(chanmax,chanmax,ptsmax)

      double precision norm,dotprd
      double precision sqtnorm,minnorm
      integer i,j,k,l,m,iters,numv

C * i is the index for the vectors in xhi
C * j is the index for the vectors in psifin
      minnorm = 1.d-30
      iters = numv

C * orthogonalize xhi with the current basis set
C * then we'll normalize it and add it to the set
      do i = 1,nch
         do j = 1,nch
            do k = 1,npts
               psi(j,i,k) = xhi(j,i,k)
            end do
         end do
      end do

C * We now have (in PSI) a new set of NCH vectors to orthonormalize
C * to the set we have accumulated 

C * Recall we get NCH new vectors to orthonormalize with each
C * iteration so at the end we have NITER * NCH orthonormal vectors
      do i = 1,nch
         do j = 1,iters-1
            do k = 1,nch
               dotprd = 0.d0
               do l = 1,nch
                  if (ity .eq. 0) then
                     ipts = 0
                     do ireg = 1,nreg
                        do istep = 1,nstep(ireg)
                           ipts = ipts+1
                           dotprd = dotprd + psifin(j,k,l,ipts)*xhi(l,i
     $                          ,ipts)*step(ireg)
                        end do
                     end do
                  else if (ity .eq. 1 .or. ity .eq. 2) then
                     do ipts = 1, npts
                        dotprd = dotprd + psifin(j,k,l,ipts)*xhi(l,i
     $                       ,ipts)*wtt(ipts)
                     end do
                  endif
               end do
               do l = 1,nch
                  do m = 1,npts
                     psi(l,i,m)  =  psi(l,i,m) - dotprd*psifin(j,k,l,m)
                  end do
               end do
            end do
         end do
      end do

C * Above we orthonormalized the set of NCH vectors in this
C * iteration to the accumulated set.

C * Below we orthonormalize the NCH iterate vectors to each other
      do i = 1,nch
         do j = 1,i-1
            dotprd = 0.d0
            do k = 1,nch
               if (ity .eq. 0) then
                  ipts = 0
                  do ireg = 1,nreg
                     do istep = 1,nstep(ireg)
                        ipts = ipts+1
                        dotprd = dotprd+psi(k,i,ipts)*psi(k,j,ipts)*
     $                       step(ireg)
                     end do
                  end do
               else if (ity .eq. 1 .or. ity .eq. 2) then
                  do ipts = 1, npts
                     dotprd = dotprd+psi(k,i,ipts)*psi(k,j,ipts)
     $                    *wtt(ipts)
                  end do
               endif
            end do
            do k = 1,nch
               do l = 1,npts
                  psi(k,i,l) = psi(k,i,l)-dotprd*psi(k,j,l)
               end do
            end do
         end do
         norm = 0.d0

C * normalize the vector-- if this is the first vector
C * then the subroutine falls to here first
         do j = 1,nch
            if (ity .eq. 0) then
               ipts = 0
               do ireg = 1,nreg
                  do istep = 1,nstep(ireg)
                     ipts = ipts+1
                     norm = norm+psi(j,i,ipts)*psi(j,i,ipts) *step(ireg)
                  end do
               end do
            else if (ity .eq. 1 .or. ity .eq. 2) then
               do ipts = 1, npts
                  norm = norm+psi(j,i,ipts)*psi(j,i,ipts)*wtt(ipts)
               end do
            endif
         end do
         write(7,*) 'the norm is',norm
         if(norm.lt.minnorm) then
            write(7,*) 'another vector less than',minnorm
            sqtnorm = 0.d0
         else
            sqtnorm = 1.d0/dsqrt(norm)
         endif

C * normalize xhi --> psi and psifin
C * we store BOTH the orthonormal and original vectors for later use:
         do k = 1,nch
            do l = 1,npts
               xhifin(iters,i,k,l) = xhi(k,i,l)
               psifin(iters,i,k,l) = psi(k,i,l)*sqtnorm
               psi(k,i,l) = psifin(iters,i,k,l)
            end do
         end do
      end do

      return
      end
C *-/

C *-
      subroutine newxhi(xhi,g1,g2,i1,i2,nch,npts,chin)
      implicit none

C * take the i matrices, the g matrices and make the next iterates
C * the xhi's
C * this is right out of the Linear Algebraic papers
CARR
      integer chanmax,ptsmax
      parameter (chanmax=165,ptsmax=331) 
c$$$      integer nch,npts
c$$$      double precision xhi(nch,nch,npts)
c$$$      double precision g1(nch,npts),g2(nch,npts),b(nch,npts)
c$$$      double precision i1(nch,nch,npts),i2(nch,nch,npts)
c$$$      integer chin(nch,2)
      double precision xhi(chanmax,chanmax,ptsmax)
      double precision g1(chanmax,ptsmax)
      double precision g2(chanmax,ptsmax)
      double precision i1(chanmax,chanmax,ptsmax)
      double precision i2(chanmax,chanmax,ptsmax)
      integer nch,npts,chin(chanmax,2)

      integer i,j,k

      do j = 1,nch
         do k = 1,nch
            do i = 1,npts - 1
               xhi(j,k,i) = g1(j,i)*i2(j,k,i+1)+g2(j,i)*i1(j,k,i)
            end do
            xhi(j,k,npts) = g2(j,npts)*i1(j,k,npts)
         end do
      end do

      return
      end

C *-
      subroutine psifinal(psifin,psi,acoef,asympsi,nch,npts,niter)
      implicit none

C * this routine creates the solution wavefunction psi from
C * the orthonormal basis psifin and the coefficients a
      integer nch,npts,niter
CARR
      integer chanmax,ptsmax,itermax
      parameter (chanmax=165,ptsmax=331,itermax=10) 
c$$$      double precision psifin(niter,nch,nch,npts)
c$$$      double precision psi(nch,nch,npts)
c$$$      double precision acoef(nch*niter,nch)
c$$$      double precision asympsi(nch,nch)
      double precision psifin(itermax,chanmax,chanmax,ptsmax)
      double precision psi(chanmax,chanmax,ptsmax)
      double precision acoef(chanmax*itermax,chanmax)
      double precision asympsi(chanmax,chanmax)

      integer i,j,jj,k,l

C * We already have the coefficients for the basis vectors
C * So here it is just a matter of multiplying through to 
C * get the final wavefunction for the problem.
      do i = 1,nch
         do k = 1,nch
            do l = 1,npts
               psi(k,i,l) = 0.d0
               do j = 1,niter-1
                  do jj = 1,nch
                     psi(k,i,l) = psi(k,i,l)+acoef(nch*(j-1)+jj,i)
     $                    *psifin(j,jj,k,l)
                  end do
               end do
            end do
         end do
      end do

C * now store the final (radial) points of the
C * wavefunction (in the asymptotic regin
C * for k matrix extraction) in a separate matrix
      do i = 1,nch
         do j = 1,nch
            asympsi(i,j) = psi(i,j,npts)
         end do
      end do
      return
      end

C *-
      subroutine printwfn(psi,mprint,nch,step,nstep,nreg,npts)
      implicit none
      integer nch,npts,mprint
      integer nreg,nstep(nreg)
      double precision step(nreg),r
      double precision psi(nch,nch,npts)
      integer i,j,k,l,m
      do j=1,mprint
         r=0.0d0
         m=0
         do i=1,nreg
            do l=1,nstep(i)
               r=r+step(i)
               m=m+1
               write(81,500)r,(psi(j,k,m),k=1,mprint)
            end do
         end do
      end do
 500  format(f7.2,7(1pe12.4))
      return
      end

      subroutine qget3(q1,q2,g1,g2,v,psi,nstep,step,nreg,nch,npts,wtt
     $     ,ity)
      implicit none

      integer nch,npts,nreg
      integer nstep(nreg)
      double precision step(nreg)
      double precision wtt(npts)
      integer ity
CARR
      integer chanmax,ptsmax
      parameter(chanmax=165,ptsmax=331)
c$$$      double precision q1(nch,nch,npts),q2(nch,nch,npts)
c$$$      double precision g1(nch,npts),g2(nch,npts)
c$$$      double precision v(nch,nch,npts)
c$$$      double precision psi(nch,nch,npts)
      double precision q1(chanmax,chanmax,ptsmax)
      double precision q2(chanmax,chanmax,ptsmax)
      double precision g1(chanmax,ptsmax)
      double precision g2(chanmax,ptsmax)
      double precision v(chanmax,chanmax,ptsmax)
      double precision psi(chanmax,chanmax,ptsmax)

      integer i,j,k,ireg,ipts

C * This routine calculates the Q matrices for the "direct" part
C * of the LA method.  
      call matprdm(v,psi,q1,nch,npts)

C * Gauss-Legendre quadratures use different r points
      if (ity .eq. 0) then
         ipts=0
         do ireg=1,nreg
            do i=1,nstep(ireg)-1
               ipts=ipts+1
               do j=1,nch
                  do k=1,nch
                     q2(j,k,ipts)=(g2(j,ipts)*q1(j,k,ipts))*step(ireg)
                     q1(j,k,ipts)=(g1(j,ipts)*q1(j,k,ipts))*step(ireg)
                  end do
               end do
            end do
            ipts=ipts+1
            if(ireg.lt.nreg)then
               do j=1,nch
                  do k=1,nch
                     q2(j,k,ipts)=g2(j,ipts)*q1(j,k,ipts)*(step(ireg)
     $                    +step(ireg+1))/2.d0
                     q1(j,k,ipts)=g1(j,ipts)*q1(j,k,ipts)*(step(ireg)
     $                    +step(ireg+1))/2.d0
                  end do
               end do
            else
               do j=1,nch
                  do k=1,nch
                     q2(j,k,ipts)=g2(j,ipts)*q1(j,k,ipts)*step(ireg)
     $                    /2.d0
                     q1(j,k,ipts)=g1(j,ipts)*q1(j,k,ipts)*step(ireg)
     $                    /2.d0
                  end do
               end do
            endif
         end do
      
      else if (ity .eq. 1 .or. ity .eq. 2) then
         do ipts = 1, npts
            do j = 1, nch
               do k = 1, nch
                  q2(j,k,ipts) = (g2(j,ipts)*q1(j,k,ipts))*wtt(ipts)
                  q1(j,k,ipts) = (g1(j,ipts)*q1(j,k,ipts))*wtt(ipts)
               end do
            end do
         end do
      endif

      do j=1,nch
         do k=1,nch
            q2(j,k,1)=q2(j,k,1)*1.5d0
            q1(j,k,1)=q1(j,k,1)*1.5d0
         end do
      end do
      return
      end

C *-
C * this routine does the same as the qget3 routine, but here
C * we include the exchange potential

      subroutine qget3x(q1,q2,g1,g2,v,kernel,dummy,exunit,psi,kpsi,
     $     nstep,step,nreg,nregex,nptsex,nch,npts,nchex,nexdim,nvib
     $     ,nvibx,nvibxin,npw,npwavx,npwavxin,wtt,ity)
      implicit none

      integer nvib,nvibx,nvibxin,npw,npwavx,npwavxin
      integer exunit,nch,nchex,nexdim
      integer nptsex,nregex,npts,nreg
      integer nstep(nreg)
      double precision step(nreg)
      double precision wtt(npts)
      integer ity
CARR
      integer chanmax,ptsmax
      parameter(chanmax=165,ptsmax=331)
c$$$      double precision q1(nch,nch,npts),q2(nch,nch,npts)
c$$$      double precision g1(nch,npts),g2(nch,npts)
c$$$      double precision v(nch,nch,npts)
c$$$      double precision psi(nch,nch,npts)
c$$$      double precision kpsi(nch,nch)
      double precision q1(chanmax,chanmax,ptsmax)
      double precision q2(chanmax,chanmax,ptsmax)
      double precision g1(chanmax,ptsmax),g2(chanmax,ptsmax)
      double precision v(chanmax,chanmax,ptsmax)
      double precision psi(chanmax,chanmax,ptsmax)
      double precision kpsi(chanmax,chanmax)

      double precision kernel(nexdim,nexdim,nptsex)
      double precision dummy(nchex,nchex,nptsex)
      integer i,j,k,ireg,ipts


C * This routine calculates the Q matrices for the "direct" part
C * of the LA method.  
      call matprdm(v,psi,q1,nchex,npts)

C * Gauss-Legendre quadrature, use different r points
      if (ity .eq. 0) then
         ipts = 0
         write(6,*) 'Exchker * wfn'

         do ireg = 1,nreg
            do i = 1,nstep(ireg)-1
               ipts = ipts + 1
               if(ipts .le. nptsex)then
                  call getkpsi(kernel,psi,kpsi,dummy,step,nstep,nvib
     $                 ,nvibx,nvibxin,npw,npwavx,npwavxin,nch,nchex
     $                 ,nexdim,npts,nreg,nregex,nptsex,exunit,wtt,ity)
               else
                  do j = 1,nchex
                     do k = 1,nchex
                        kpsi(j,k) = 0.d0
                     end do
                  end do
               endif

               do j = 1,nch
                  do k = 1,nch
                     q2(j,k,ipts) = (g2(j,ipts)*(q1(j,k,ipts)-2.d0
     $                    *kpsi(j,k)))*step(ireg)
                     q1(j,k,ipts) = (g1(j,ipts)*(q1(j,k,ipts)-2.d0
     $                    *kpsi(j,k)))*step(ireg)
                  end do
               end do
            end do
            ipts = ipts + 1
            if(ipts.le.nptsex) then
               call getkpsi(kernel,psi,kpsi,dummy,step,nstep,nvib ,nvibx
     $              ,nvibxin,npw,npwavx,npwavxin,nch,nchex ,nexdim,npts
     $              ,nreg,nregex,nptsex,exunit,wtt,ity)
            else
               do j = 1,nch
                  do k = 1,nch
                     kpsi(j,k) = 0.d0
                  end do
               end do
            endif
            if(ireg.lt.nreg)then
               do j = 1,nch
                  do k = 1,nch
                     q2(j,k,ipts) = (g2(j,ipts)*(q1(j,k,ipts)-2.d0
     $                    *kpsi(j,k)))*(step(ireg)+step(ireg+1))/2.d0
                     q1(j,k,ipts) = (g1(j,ipts)*(q1(j,k,ipts)-2.d0
     $                    *kpsi(j,k)))*(step(ireg)+step(ireg+1))/2.d0
                  end do
               end do
            else
               do j = 1,nch
                  do k = 1,nch
                     q2(j,k,ipts) = (g2(j,ipts)*(q1(j,k,ipts)-2.d0
     $                    *kpsi(j,k)))*step(ireg)/2.d0
                     q1(j,k,ipts) = (g1(j,ipts)*(q1(j,k,ipts)-2.d0
     $                    *kpsi(j,k)))*step(ireg)/2.d0
                  end do
               end do
            endif
         end do

      else if (ity .eq. 1 .or. ity .eq. 2) then
         do ipts = 1, nptsex
            call getkpsi(kernel,psi,kpsi,dummy,step,nstep,nvib ,nvibx
     $           ,nvibxin,npw,npwavx,npwavxin,nch,nchex ,nexdim,npts
     $           ,nreg,nregex,nptsex,exunit,wtt,ity)
            do j = 1,nch
               do k = 1,nch
                  q2(j,k,ipts) = (g2(j,ipts)*(q1(j,k,ipts)-2.d0 *kpsi(j
     $                 ,k)))*wtt(ipts)
                  q1(j,k,ipts) = (g1(j,ipts)*(q1(j,k,ipts)-2.d0 *kpsi(j
     $                 ,k)))*wtt(ipts)
               end do
            end do
         end do
      endif

      do j = 1,nch
         do k = 1,nch
            q2(j,k,1) = q2(j,k,1)*1.5d0
            q1(j,k,1) = q1(j,k,1)*1.5d0
         end do
      end do

      return
      end

C *-
      subroutine tmatrix(kmat,tmati,tmatr,temp,nch)
      implicit none

      integer nch

C * numerical recipes inverse routines
CARR
      integer chanmax
      parameter(chanmax=165)
c$$$      double precision temp(nch,nch)
c$$$      double precision kmat(nch,nch)
c$$$      double precision tmatr(nch,nch),tmati(nch,nch)
c$$$      integer indx(nch)
      double precision temp(chanmax,chanmax)
      double precision kmat(chanmax,chanmax)
      double precision tmatr(chanmax,chanmax),tmati(chanmax,chanmax)
      integer indx(chanmax)

      integer i,j

C     call mymatmul(kmat,kmat,tmati,nch,nch)
      call mymatmul(kmat,kmat,tmati,nch,chanmax)
      do i = 1,nch
         tmati(i,i) = tmati(i,i)+1.0d0
      end do

C     call nrinverse(tmati,tmatr,temp,nch,nch,indx)
C     call mymatmul(kmat,tmatr,tmati,nch,nch)
      call nrinverse(tmati,tmatr,temp,nch,chanmax,indx)
      call mymatmul(kmat,tmatr,tmati,nch,chanmax)

      do i = 1,nch
         do j = 1,nch
            tmati(i,j) = -2.d0*tmati(i,j)
         end do
      end do
C     call mymatmul(kmat,tmati,tmatr,nch,nch)
      call mymatmul(kmat,tmati,tmatr,nch,chanmax)
      do i = 1,nch
         do j = 1,nch
            tmatr(i,j) = -tmatr(i,j)
         end do
      end do

      return
      end

C ********************************************************************
      SUBROUTINE CLEB(j1,j2,j,m1,m2,m)
C * -----------------------------------------------------------------
C *     REFERENCE   T. TAMURA   COMPUTER PHYS. COMM. VOL. 1, P. 337
C *     ERRATUM  COMPUTER PHYS. COMM. VOL. 2, P. 174

C *     CALLING SEQUENCE...
C *     CALL FACSET ONCE BEFORE ANY CALLS TO CLEB, RACAH, OR NINEJ
C *     SET IA=2*J1, IB=2*J2, IC=2*J
C *     SET ID=2*M1, IE=2*M2, IF=2*M
C *     CALL CLEB
C *     THE CLEBSCH-GORDAN COEFFICIENT (J1 M1 J2 M2/ J M) WILL BE
C *     RETURNED IN THE VARIABLE RAC
C *     CALLING PARAMETERS ARE PRESERVED
C * -----------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer j1,j2,j,m1,m2,m
      common /clebG/ icleb
      COMMON/FACLG/ FACLOG(500),RAC,U9,IA,IB,IC,ID,IE,IF,L9(9),KSELC
      icleb=icleb+1
      ia=2*j1
      ib=2*j2
      ic=2*j
      id=2*m1
      ie=2*m2
      if=2*m
      RAC=0.0D+00
      IG4=IA+IB+IC
      IF(KSELC.EQ.0) GO TO 10
      IF(ID+IE.NE.IF) GO TO 70
      IF(MOD(IG4,2).NE.0) GO TO 70
      IF(IA+IB-IC.LT.0.OR.IC-IABS(IA-IB).LT.0) GO TO 70
      IF(MIN0(IA-IABS(ID),IB-IABS(IE),IC-IABS(IF)).LT.0) GO TO 70
      IF(MOD(IB+IE,2).NE.0.OR.MOD(IC+IF,2).NE.0) GO TO 70
   10 IF(IA.EQ.0.OR.IB.EQ.0) GO TO 20
      IF(IC.GT.0) GO TO 40
      IF(IC.EQ.0) GO TO 30
      IF(IC.LT.0) GO TO 70
   20 RAC=1.0D+00
      GO TO 70
   30 FB=IB+1
      RAC=((-1.0D+00)**((IA-ID)/2))/DSQRT(FB)
      GO TO 70
   40 IF(ID.NE.0.OR.IE.NE.0) GO TO 50
      IG2=IG4/2
      IF(MOD(IG2,2).NE.0) GO TO 70
      IG=IG2/2
      I1=(IA+IB-IC)/2+1
      I2=(IC+IA-IB)/2+1
      I3=(IC+IB-IA)/2+1
      I4=IG2+2
      I5=IG+1
      I6=(IG2-IA)/2+1
      I7=(IG2-IB)/2+1
      I8=(IG2-IC)/2+1
      F1=DEXP(0.5D+00*(FACLOG(I1)+FACLOG(I2)+FACLOG(I3)-FACLOG(I4))
     $     +(FACLOG(I5)-FACLOG(I6)-FACLOG(I7)-FACLOG(I8)))
      F2=IC+1
      F2=DSQRT(F2)
      S1=1-2*MOD((IG2+IC)/2,2)
      RAC=S1*F1*F2
      GO TO 70
 50   FC2=IC+1
      IABCP=IG4/2+1
      IABC=IABCP-IC
      ICAB=IABCP-IB
      IBCA=IABCP-IA
      IAPD=(IA+ID)/2+1
      IAMD=IAPD-ID
      IBPE=(IB+IE)/2+1
      IBME=IBPE-IE
      ICPF=(IC+IF)/2+1
      ICMF=ICPF-IF
      SQFCLG=0.5D+00*(DLOG(FC2)-FACLOG(IABCP+1)+FACLOG(IABC)
     $     +FACLOG(ICAB)+FACLOG(IBCA)+FACLOG(IAPD)+FACLOG(IAMD)
     $     +FACLOG(IBPE)+FACLOG(IBME)+FACLOG(ICPF)+FACLOG(ICMF))
      NZMIC2=(IB-IC-ID)/2
      NZMIC3=(IA-IC+IE)/2
      NZMI=MAX0(0,NZMIC2,NZMIC3)+1
      NZMX=MIN0(IABC,IAMD,IBPE)
      S1=(-1.0D+00)**(NZMI-1)
      DO NZ=NZMI,NZMX
         NZM1=NZ-1
         NZT1=IABC-NZM1
         NZT2=IAMD-NZM1
         NZT3=IBPE-NZM1
         NZT4=NZ-NZMIC2
         NZT5=NZ-NZMIC3
         TERMLG=SQFCLG-FACLOG(NZ)-FACLOG(NZT1)-FACLOG(NZT2)
     $        -FACLOG(NZT3)-FACLOG(NZT4)-FACLOG(NZT5)
         SSTERM=S1*DEXP(TERMLG)
         RAC=RAC+SSTERM
         S1=-S1
      ENDDO
      IF(DABS(RAC).LT.1.0D-12) RAC=0.0D+00
 70   RETURN
      END

C *-
      SUBROUTINE FACSET
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/FACLG/ FACLOG(500),RAC,U9,IA,IB,IC,ID,IE,IF,L9(9),KSELC
      FACLOG(1)=0.0D+00
      FACLOG(2)=0.0D+00
      F1=1.0D+00
      DO N=3,500
         F1=F1+1.0D+00
         FACLOG(N)=FACLOG(N-1)+DLOG(F1)
      ENDDO
      RETURN
      END

C *-
C * this function calculates the nuclear part of the interaction potential
C * use it for the higher vlambdas where the electron cloud doesn't
C * contribute

C * NOTE: this routine works only with the IRR=1 option set.
C * That is, it works only when doing a Rigid Rotator or VIBAV
C * calculation.

      double precision function vnuc(r,lam,ra,rb,qa,qb)
      implicit none
      double precision r,ra,rb,qa,qb
      integer lam
      double precision factor
      
      factor = -1.0d0

      ra=abs(ra)
      rb=abs(rb)
      vnuc=-qa*(min(r,ra)**lam)/(max(r,ra)**(lam+1))*factor**lam
      vnuc=vnuc-qb*(min(r,rb)**lam)/(max(r,rb)**(lam+1))
      return
      end

C *-
C * this routine reads the FORMATTED input Rigid Rotator or VIBAV
C * LOCAL potential.  It also calculates higher Vlambdas (if needed)
C * from the analytic nuclear terms (see function VNUC)

      subroutine potread(vlam,stepimesh,alpha0,alpha2,q,nlamdain,nlam
     $     ,nptsimesh,infil,nregimesh,nstepimesh)
      implicit none

C  * NLAMDAIN in the number of Legendre terms included in the input
C * potential
C * NLAM: the total number of Legendre terms requested in the 
C *       short range part of the local part of the potential
      double precision alpha0,alpha2,q
      integer i,j,nlamdain,lamda,nlam,npts
      integer nregimesh,nptsimesh

CARR
      integer potptsmax,lammax,vibmax,nrgmax,cntrmax
      parameter(potptsmax=330,lammax=21,vibmax=15,cntrmax=2,nrgmax=100)
c$$$      integer nstepimesh(nregimesh)
c$$$      double precision vlam(nlam,1,1,nptsimesh),vnuc
c$$$      double precision stepimesh(nregimesh)
      integer nstepimesh(nrgmax)
      double precision vlam(lammax,vibmax,vibmax,potptsmax),vnuc
      double precision stepimesh(nrgmax)

      double precision rstrt(nregimesh),rstop(nregimesh)
      double precision rstep(nregimesh)
      integer nptsreg(nregimesh)
      character*8 infil
      integer ind
      character*80 potitle
      integer nreg,ncntr
      double precision charge(cntrmax),nucpos(cntrmax)

      integer ireg,ipts,lam,ptscount
      double precision r
      ind=index(infil,' ')
      if(ind.eq.0) ind=9

      open(unit=9,file=infil(1:ind-1),form='formatted', status =
     $     'unknown')

      write(6,*)'Routine Potread'
      read(9,500) potitle
      write(6,500) potitle
 500  format(a80)
      read(9,*) nreg
      write(6,501) nreg
 501  format(' Number of Regs in potential Mesh is:',i5)
      write(6,502)
      do i=1,nreg
         read(9,*) rstrt(i),rstop(i),rstep(i)
         nptsreg(i) = int((rstop(i)-rstrt(i)+rstep(i)/100.d0)/rstep(i))
         write(6,503) rstrt(i),rstop(i),rstep(i),nptsreg(i)
      end do
 502  format(' Input Potential Mesh',/,'    Rstart    Rstop    Rstep')
 503  format(3f10.5,i5)
      read(9,*) ncntr
      write(6,504) ncntr
 504  format(' Number of Nuclear Centers: ',i4)
      write(6,505)
      do i=1,ncntr
         read(9,*) charge(i),nucpos(i)
         write(6,506) charge(i),nucpos(i)
      end do
 505  format('Charge and Positions of Nuclear Centers:')
 506  format(2f10.5)
      do i=1,nlamdain
         read(9,*) lamda,npts
         write(6,507) lamda,npts
         do j=1,npts
            read(9,508) r,vlam(i,1,1,j)
         end do
      end do
      write(6,*)

C * If we want to do the LA method beyond the range of the Short-Range 
C * potential, calculate the long range potential on the integration 
C * mesh. this assumes that the SR pot mesh and integration mesh are
C * ***THE SAME****
      ipts=npts
      do ireg=nreg+1,nregimesh
         do j=1,nstepimesh(ireg)
            ipts=ipts+1
            r=r+stepimesh(ireg)
            do i=1,nlamdain
               if(i.eq.1)then
                  vlam(i,1,1,ipts)=-alpha0/(2.d0*r**4.d0)
               else if(i.eq.2)then
                  vlam(i,1,1,ipts)=-alpha2/(2.d0*r**4.d0)-q/(r*r*r)
               else
                  vlam(i,1,1,ipts)=0.d0
               endif
            end do 
            write(6,510) r,vlam(1,1,1,ipts),vlam(2,1,1,ipts)
 510        format(f7.3,2(1pe12.3))
         end do 
      end do 
C * end of long-range potential stuff

      if(nlam.gt.nlamdain) then
         write(6,*)' Read ',nlamdain,' Legendre terms from'
         write(6,*)' input potential file.  Extending to ',nlam
         write(6,*)' terms with analytic e-nuclear terms'
         write(6,*)
      endif

      do i=nlamdain+1,nlam
         r=0.d0
         lam=2*(i-1)
         write(52,507)lam,npts
         ptscount=0
         do ireg=1,nreg
            do ipts=1,nptsreg(ireg)
               ptscount=ptscount+1
               r=r+rstep(ireg)
               vlam(i,1,1,ptscount)=vnuc(r,lam,nucpos(1),nucpos(2),
     $              charge(1),charge(2))
               write(52,508)r,vlam(i,1,1,ptscount)
            end do
         end do
      end do
 507  format(' Lamda=',i3,' Npts=',i4)
C508  format(f10.5,d23.16)
 508  format(f10.5,d24.16)

C * Now we have read in the Vlamdas and stored the in the array Vlam
      write(6,509)
 509  format('End of Potential Data')
      close(9)
      return
      end

C *-
C * This routine calculates the angular coupling coefficients for
C * the matrix elements of the local potential
C * essentially it amounts to a sum over lambdas and
C * products with Clebsch Gordons -- See Morrison and Collins PRA 1978

      subroutine fmatrix(fmat,chin,l0,nlamdain,symlam,npw,nch)
      implicit none
      integer l0,nlamdain,symlam,npw,nch

CARR
      integer pwavemax,lammax,chanmax
      parameter(chanmax=165,pwavemax=11,lammax=21)
c$$$      double precision fmat(npw,npw,nlamdain),c1,c2
c$$$      integer chin(nch,2)
      double precision fmat(pwavemax,pwavemax,lammax),c1,c2
      integer chin(chanmax,2)

      integer l,lp,lam,i,j,k
      double precision faclog(500)
      double precision rac,u9
      integer ia,ib,ic,id,ie,if,l9(9),kselc
      common/faclg/faclog,rac,u9,ia,ib,ic,id,ie,if,l9,kselc

C * set up the factorial table
      do i=1,npw
         l=2*(chin(i,2)-1)+l0
         do j=1,npw
            lp=2*(chin(j,2)-1)+l0
            do k=1,nlamdain
               lam=2*(k-1)
               call cleb(lp,lam,l,symlam,0,symlam)
               c1=rac
               if(symlam.eq.0)then
                  c2=c1
               else
                  call cleb(lp,lam,l,0,0,0)
                  c2=rac
               endif
               fmat(i,j,k)=dsqrt((2.d0*lp+1.d0)/(2.d0*l+1.d0))*c1*c2
            end do
         end do
      end do
      return
      end

C *-
C * Now the goal is to extend the code to do e-H2 scattering...
C * the code will be hooked onto the end of the Linear Algebraic code
C * to propagate the R-Matrix out into the asymptotic region

C * 7-30-04 A. Feldt fix - changed nvibin to vibdim to reflect the
C *         change in the value passed to it.  This fixes
C *         the dimensioning error that had obtained if nvibin < vibmax
C *         alpha0, alpha2 and q now properly use this size

      subroutine rprop(rmatp,rstart,rasym,step,k,l0,alpha0,alpha2,q
     $     ,fmat,nlamdain,nlams,chin,vibdim,nch,npw)
      implicit none

      integer i,l0,nch,npts,nlams,nlamdain
      integer vibdim,npw

      integer chanmax, pwavemax, vibmax, lammax
      parameter(chanmax=165,pwavemax=11,vibmax=15,lammax=21)
      double precision fmat(pwavemax,pwavemax,lammax)
      double precision alpha0(vibmax,vibmax),alpha2(vibmax,vibmax)
      double precision q(vibmax,vibmax) 
      double precision rstart,rasym,k(chanmax)
      double precision rmatp(chanmax,chanmax)
      integer chin(chanmax,2)
      double precision rmatn(chanmax,chanmax)
      double precision rptrans(chanmax,chanmax),transp(nch,nch)

      double precision v(nch,nch)
      double precision trans(nch,nch),t1(nch),t2(nch)
      double precision lamsq(nch),abslam(nch),habslam(nch)

      double precision step,r
      integer j,ik

      write(6,*)
      write(6,*) ' In R-Prop'
      write(6,*) 'rstart      rasym      step     nlams     '
      write(6,*) rstart,rasym,step,nlams
      write(6,*) 'nch = ',nch
      call rinit(trans,transp,nch)

C * Set up the initial conditions:
C * now look at the next r point:
C * olen changed the below to nint( ) from an int( )
c * this seems to be necessary for linux
      npts = nint((rasym-rstart)/step)
c * end olen modification

      write(6,*) 'npts = ', npts, ' l0 = ', l0
      do i = 1,npts
C * we want the potential at the midpoints of the grid for the
C * R-matrix propagator
         r = rstart + (i-1)*step + step/2.d0
         call lrpot(v,r,l0,alpha0,alpha2,q,vibdim,fmat,chin,nlamdain
     $        ,nlams,nch,npw)
         call lamdiag(v,k,lamsq,abslam,habslam,trans,nch,step)
         call tget(lamsq,habslam,abslam,t1,t2,nch)
         call rtrans(rmatp,rptrans,trans,transp,nch)
         call rnext(rmatn,rmatp,rptrans,t1,t2,nch)
      end do
      call unitrans(trans,rmatp,nch,0)
      return
      end

C *-
C * some initialization for the R-matrix propagator
      subroutine rinit(trans,transp,nch)
      implicit none
      integer nch,i,j
      double precision trans(nch,nch),transp(nch,nch)
      do i = 1,nch
         do j = 1,nch
            trans(i,j) = 0.d0
            transp(i,j) = 0.d0
         end do
         trans(i,i) = 1.d0
         transp(i,i) = 1.d0
      end do
      return
      end

C *-
C * diagonalize the potential minus the energy
      Subroutine lamdiag(v,k,lamsq,abslam,habslam,trans,nch,step)
      implicit none
      integer nrot
      integer i,j,nch
      double precision v(nch,nch),k(nch),trans(nch,nch)
      double precision lamsq(nch),abslam(nch),habslam(nch)
      double precision step,sum
      do i = 1,nch
         v(i,i) = v(i,i)-k(i)*k(i)
      end do

C * devcsf gets the evals and evects of a real symmetric matrix
C * so that it returns real eigenvectors
C * 8-28-95 WKT  call devcsf(nch,v,nch,lamsq,trans,nch)
      call jacobi(v,nch,nch,lamsq,trans,nrot)
      do j = 1,nch
         sum = 0.d0
         do i = 1,nch
            sum = sum+trans(i,j)*trans(i,j)
         end do
         sum = sqrt(sum)
         do i = 1,nch
            trans(i,j) = trans(i,j)/sum
         end do
      end do
      do i = 1,nch
         abslam(i) = sqrt(abs(lamsq(i)))
         habslam(i) = step*abslam(i)
      end do
      return
      end

C *-
C * Calculate the R matrix at the next box using its value at the prev box

      Subroutine rnext(rmatn,rmatp,rptrans,t1,t2,nch)
      implicit none

      integer nch,i,j
CARR
      integer chanmax
      parameter(chanmax=165)
c$$$      integer indx(nch)
c$$$      double precision rmatp(nch,nch)
c$$$      double precision rmatn(nch,nch),rptrans(nch,nch)
      integer indx(chanmax)
      double precision rmatp(chanmax,chanmax)
      double precision rmatn(chanmax,chanmax),rptrans(chanmax,chanmax)

       double precision t1(nch),t2(nch)

      do i = 1,nch
         rptrans(i,i) = rptrans(i,i)+t1(i)
      end do

C * here use rmatp as a dummy array for the inversion routine
C * the inverse of rptrans ends up in rmatn
C     call nrinverse(rptrans,rmatn,rmatp,nch,nch,indx)
      call nrinverse(rptrans,rmatn,rmatp,nch,chanmax,indx)

      do i = 1,nch
         do j = 1,nch
            rmatn(i,j) = -t2(i)*rmatn(i,j)*t2(j)
            rmatp(i,j) = rmatn(i,j)
         end do
         rmatn(i,i) = rmatn(i,i)+t1(i)
         rmatp(i,i) = rmatn(i,i)
      end do
      return
      end

C *-
C * Get the values for the previous R matrix
      Subroutine tget(lamsq,habslam,abslam,t1,t2,nch)
      implicit none
      integer nch,i
      double precision abslam(nch),habslam(nch),lamsq(nch)
      double precision t1(nch),t2(nch)
      do i = 1,nch
         if(lamsq(i).gt.0.0d0) then
            t1(i) = 1.d0/(tanh(habslam(i))*abslam(i))
            t2(i) = 1.d0/(sinh(habslam(i))*abslam(i))
         else
            t1(i) = -1.d0/(tan(habslam(i))*abslam(i))
            t2(i) = -1.d0/(sin(habslam(i))*abslam(i))
         endif
      end do
      return
      end

C *-
C * Untransform the R-matrix from the previous mesh point and transform 
C * it onto the current point.  The transformation matrices are the 
C * matrices used to diagonalize the interaction potential.

      Subroutine rtrans(rmatp,rptrans,trans,transp,nch)
      implicit none

      integer i,j,nch
CARR
      integer chanmax
      parameter(chanmax=165)
c$$$      double precision rmatp(nch,nch),rptrans(nch,nch)
      double precision rmatp(chanmax,chanmax),rptrans(chanmax,chanmax)

      double precision trans(nch,nch),transp(nch,nch)

      call unitrans(transp,rmatp,nch,0)
      call unitrans(trans,rmatp,nch,1)

C * now put the correct elements in transp for the next loop
      do i = 1,nch
         do j = 1,nch
            transp(i,j) = trans(i,j)
            rptrans(i,j) = rmatp(i,j)
         end do
      end do
      return
      end

C *-
C * unitary transformation 
C * if iswitch = 0 you get T A T(trans)
C * if iswitch = 1 you get T(trans) A T

      Subroutine Unitrans(t,a,nch,iswitch)
      implicit none
      integer nch,i,j,k,l,iswitch
CARR
      integer chanmax
      parameter(chanmax=165)
c$$$      double precision t(nch,nch),a(nch,nch)
      double precision t(nch,nch),a(chanmax,chanmax)

      double precision dummy(nch,nch),tp(nch,nch)

      if(iswitch.eq.0) then
         do i = 1,nch
            do j = 1,nch
               tp(i,j) = t(j,i)
            end do
         end do
      else if (iswitch.eq.1) then
         do i = 1,nch
            do j = 1,nch
               tp(i,j) = t(i,j)
            end do
         end do
      endif
      do i = 1,nch
         do j = 1,nch
            dummy(i,j) = 0.0
            do k = 1,nch
               do l = 1,nch
                  dummy(i,j) = dummy(i,j)+tp(k,i)*a(k,l)*tp(l,j)
               end do
            end do
         end do
      end do
      do i = 1,nch
         do j = 1,nch
            a(i,j) = dummy(i,j)
         end do
      end do
      return
      end

C *-
C * This Subroutine sets up the channel indexes.  The idea is to have an
C * array that tells the codes the vibrational level and partial wave
C * index for a particular channel index.  The advantage of this 
C * method is that it allows the matrices to remain two-dimensional
C * which means that in adding a degree of freedom to the problem
C * one just changes this routine and the potential routines, but 
C * leaves the integration routines unscathed.

      subroutine chindex(chin,npw,nch)
      implicit none
      integer i,npw,nch
CARR
      integer chanmax
      parameter(chanmax=165)
c$$$      integer chin(nch,2)
      integer chin(chanmax,2)

      write(6,*) 'These are the channel labels'
      do i=1,nch
         chin(i,1)=(i-1)/npw+1
         chin(i,2)=mod(i,npw)
         if(chin(i,2).eq.0) chin(i,2)=npw
         write(6,*)i,chin(i,1),chin(i,2)
      end do
      return
      end

C *-
C * Calculate the potential matrix elements for the long-range
C * part of the potential
C * The routine hass been modified to take in vibrational matrix
C * elements of the long range moments to allow vibrational
C * excitation calculations to be carried out.

C * 7-30-04 A. Feldt fix - alpha0, alpha2 and q must be dimensioned
C *         as they were in subroutine indata to allow smaller nvibin

      Subroutine lrpot(v,r,l0,alpha0,alpha2,q,vibdim,fmat,chin,nlamdain
     $     ,nlams,nch,npw)
      implicit none

C * VIBDIM: the number of vib channels for the multipole
C *         moments and polarizabilities as they were *dimensioned* 
C *         when READ IN
      integer vibdim

      integer nch,nlamdain,nlams,npw
CARR 
      integer chanmax,vibmax,pwavemax,lammax
      parameter(chanmax=165,vibmax=15,pwavemax=11,lammax=21)
c$$$      integer chin(nch,2)
c$$$      double precision alpha0(vibdim,vibdim),alpha2(vibdim,vibdim)
c$$$      double precision q(vibdim,vibdim),fmat(npw,npw,nlamdain)
      integer chin(chanmax,2)
      double precision alpha0(vibmax,vibmax),alpha2(vibmax,vibmax)
      double precision q(vibmax,vibmax),fmat(pwavemax,pwavemax,lammax)

      double precision v(nch,nch),r,vlam(vibdim,vibdim)

      integer i,j,k,l,l0,iv,jv,il,jl

      do i = 1,nch
         iv = chin(i,1)
         il = chin(i,2)
         l = 2*(il-1)+l0
         do j = 1,nch
            jv = chin(j,1)
            jl = chin(j,2)
            v(i,j) = 0.d0
            do k = 1,nlams
               if(k.eq.1)then
                  vlam(iv,jv) = -alpha0(iv,jv)/(2.d0*r**4.d0)
               else if (k.eq.2)then
                  vlam(iv,jv) = (-alpha2(iv,jv)/(2.d0*r)-q(iv,jv))/(r
     $                 **3.d0)

C * fix by A. Feldt to allow more lambdas, the angular momentum barrier
C * only now
               else
                  vlam(iv,jv) = 0.d0
C * end fix
               endif
               v(i,j) = v(i,j)+2.d0*fmat(il,jl,k)*vlam(iv,jv)
            end do
         end do
         v(i,i) = v(i,i)+l*(l+1)/(r*r)
      end do
      return
      end

C *-
C * this subroutine calculates the potential matrix elements
C * from the legendre projections.  The routine also works
C * if there is vibrational excitation

      subroutine pot(vll,fmat,vlam,chin,nlamdain,nlamloc,nvibin,nch
     $     ,nptspot,npts,rpot,rgs,npw)
      implicit none
      integer nlamloc
      integer nlamdain,nvibin,nch,nptspot,npts,ipts,npw
CARR
      integer chanmax,potptsmax,ptsmax,pwavemax,lammax,vibmax
      parameter(chanmax=165,potptsmax=330,ptsmax=331,pwavemax=11,lammax
     $     =21,vibmax=15)
c$$$      double precision vll(nch,nch,npts),fmat(npw,npw,nlamdain)
c$$$      double precision vlam(nlamdain,nvibin,nvibin,npts)
c$$$      integer chin(nch,2)
      double precision vll(chanmax,chanmax,ptsmax)
      double precision fmat(pwavemax,pwavemax,lammax)
      double precision vlam(lammax,vibmax,vibmax,potptsmax)
      integer chin(chanmax,2)

      integer i,j,k
      double precision rpot(nptspot),rgs(npts)
      double precision vt1(nptspot),vt2(npts),vt3(lammax,npts)
      double precision csplin(nptspot)

      do i=1,nch
         do j=1,nch
C * interpolate vlam of npts points through nptspot points
            do k = 1, nlamloc
               do ipts = 1, nptspot
                  vt1(ipts) = vlam(k,chin(i,1),chin(j,1),ipts)
               enddo
               call spline(rpot,vt1,nptspot,1.d030,1.d30,csplin)
               do ipts = 1, npts
                  call splint(rpot,vt1,csplin,nptspot,rgs(ipts)
     $                 ,vt2(ipts))
                  vt3(k,ipts) = vt2(ipts)
               enddo
            enddo
            do ipts = 1, npts
               vll(i,j,ipts) = 0.d0
               do k = 1, nlamloc
                  vll(i,j,ipts) = vll(i,j,ipts) + 2.d0*fmat(chin(i,2)
     $                 ,chin(j,2),k)*vt3(k,ipts)
               end do
            end do
         end do
      end do
      return
      end
C *-/

C *-
C * this is the routine that gets the bessel and neuman
C * funtions used on construction of the green's function.
C * The routine has been modified to work for vibrational
C * excitation as well as the Rigid Rotator problem.

      subroutine bessl(irmat,g1,g2,b,k,l0,rgs,chin,nch,npts,nvopen)
      implicit none

      integer npts,nch

      integer irmat,nvopen
      double precision rgs(npts),rmax
CARR
      integer chanmax,ptsmax
      parameter(chanmax=165,ptsmax=331)
c$$$      double precision g1(nch,npts),g2(nch,npts),b(nch,npts)
c$$$      integer chin(nch,2)
      double precision g1(chanmax,ptsmax)
      double precision g2(chanmax,ptsmax)
      double precision  b(chanmax,ptsmax)
      integer chin(chanmax,2)

      double precision k(nch)
      double precision cl(nch)
      double precision rbes,rbesp,rneu,rneup
      double precision irbes,irbesp,rhank,rhankp
      double precision r
      integer ipts
      integer l0,j,l

      do j=1,nch
         l=2*(chin(j,2)-1)+l0
         rmax = rgs(npts)
         if(chin(j,1).gt.nvopen)then
            cl(j) = rhankp(k(j),rmax,l)/irbesp(k(j),rmax,l)
         else
            cl(j) = rneup(k(j),rmax,l)/rbesp(k(j),rmax,l)
         endif
         do ipts = 1, npts
            r = rgs(ipts)
            if(chin(j,1).gt.nvopen) then
               g1(j,ipts) = irbes(k(j),r,l)
               b(j,ipts)  = g1(j,ipts)
               g2(j,ipts) = -rhank(k(j),r,l)/k(j)
               if(irmat.eq.1) then
                  g2(j,ipts) = g2(j,ipts) + cl(j)*g1(j,ipts)/k(j)
               endif
            else
               g1(j,ipts) = rbes(k(j),r,l)
               b(j,ipts)  = g1(j,ipts)
               g2(j,ipts) = -rneu(k(j),r,l)/k(j)

               if(irmat.eq.1) then
                  g2(j,ipts) = g2(j,ipts) + cl(j)*g1(j,ipts)/k(j)
               endif
            endif
C * in accordance with the signs of local potential and exchange kernel 
            g2(j,ipts) = -g2(j,ipts)
C *
         end do
         if(irmat.eq.1) then
            do ipts = 1, npts
               b(j,ipts) = b(j,ipts)*g2(j,npts)
C * in accordance with the signs of local potential and exchange kernel 
               b(j,ipts) = -b(j,ipts)
C *
            end do
         endif
      end do

      return
      end
C *-/

C *-
C * This subroutine calculates the channels energies when
C * vibrational excitation is allowed. 

C * k^2 = 2*[E - (e_v - e_v0)]
C * Now v0 = 0 ONLY! (Hao Feng)

      Subroutine chener(energy,k,chin,nvib,nch,nvopen,iswit,we,wexe)
      implicit none
      integer vibmax

      integer i,nch,nvib,nvopen
      double precision energy,k(nch)

C      double precision viben(vibmax)
      double precision viben(nvib)
CARR
      integer chanmax
      parameter(chanmax=165)
c$$$      integer chin(nch,2)
      integer chin(chanmax,2)

      double precision cmtoau
      integer iswit
      double precision we, wexe

C * ==============================================================
C * The constants were taken from Herzberg (1950)C * D2 constants
C *
C *   we=3115.5d0
C *   wexe=61.8d0

C * H2 constants
C *   we=4395.2d0
C *   wexe=117.9d0

C * N2 constants
C *   we=2358.57d0
C *   wexe=14.324d0

C * DEGENERATE BFVCC!!!
C *   we=0.0d0
C *   wexe=0.d0
C * ==============================================================

C * convert from cm-1 to Hartrees
      cmtoau=1.239841d-4/(2.d0*13.6058d0)
      if (iswit .eq. 1) then
         we = we*cmtoau
         wexe = wexe*cmtoau
         iswit = 0
      endif

      viben(1)=0.d0
      nvopen=0
      do i=1,nvib-1
C        viben(i+1)=viben(i)+cmtoau*we-cmtoau*wexe*2.d0*i
         viben(i+1) = viben(i) + we - wexe*2.d0*i
      end do

      do i=1,nvib
         if((energy-viben(i)).lt.0.d0.and.nvopen.eq.0) nvopen=i-1
      end do

      if(nvopen.eq.0) nvopen = nvib
      write(6,*)'The number of vib channels is:',nvib
      write(6,*)'The number of open channels is:',nvopen
      do i=1,nch
         k(i)=dsqrt(2.d0*dabs(energy-viben(chin(i,1))))
      end do
      return
      end

C *-
C * This subroutine 
C * elements of the Legendre projections of the local part
C * of the static potentisl from unit 9.
C * The vibrational matrix elements are calculated using the 
C * WLAM program.

      subroutine vpotread(vlamvib,infil,nlamdain,nvibin,nofr,nptspot)
      implicit none
      integer nlamdain,nvibin,nofr,nptspot
CARR
      integer lammax,vibmax,potptsmax
      parameter(lammax=21,vibmax=15,potptsmax=330)
c$$$      double precision vlamvib(nlamdain,nvibin,nvibin,nptspot)
      double precision vlamvib(lammax,vibmax,vibmax,potptsmax)

C * set up stuff to do read of old wlam input:
      integer iwlampot

C     double precision temp(10,6,6,8)
      double precision temp(nofr,nvibin,nvibin,nlamdain)

      integer i,j,k,l,m
      integer ncount,ind
      character*8 infil
      ind=index(infil,' ')
      write(6,*)
      write(6,*)' in the routine vpotread:'
      write(6,*)
      write(6,*)'infil   nlamdain   nvibin   nofr   npts:'
      write(6,*)infil,nlamdain,nvibin,nofr,nptspot
      write(6,*)
      if(ind.eq.0)ind=9

      open(unit=9,file=infil(1:ind-1),form='unformatted', status
     $     ='unknown')

      iwlampot = 0

      if (iwlampot.eq.0) then
         ncount = nptspot/nofr
         do i = 1, ncount
            read(9,end=999,err=999) ((((vlamvib(j,k,l,(i-1)*nofr+m),m=1
     $           ,nofr), l=1,nvibin),k=1,nvibin),j=1,nlamdain)
         end do
      else
         do i=1,nptspot,nofr
            read(9)temp
            do j=1,nofr
               do k=1,nvibin
                  do l=1,nvibin
                     do m=1,nlamdain
                        vlamvib(m,k,l,(i-1)+j)=temp(j,k,l,m)
                     end do
                  end do
               end do
            end do
         end do
      endif

      write(6,*)
      write(6,*)'Printing the potential at the end of the'
      write(6,*)' LA region'
      do m=1,nlamdain
         write(6,*)'Lambda number:',m
         write(6,1501)(l,l=1,nvibin)
         do k=1,nvibin
            write(6,1500)k,(vlamvib(m,k,l,nptspot),l=1,nvibin)
         end do
      end do
 1500 format(i3,6(1pe12.4))
 1501 format(3x,6(i6,6x))
      write(6,*)' end of vpotread routine'
      write(6,*)
      close(9)
      return
 999  write(6,*)'i= ',i,' j=',j,' k=',k,' l=',l
      stop
      end

C *-
C * this routine truncates the R-Matrix before propagation into
C * the asymptotic region

C * NCH: the total number pw channels in the problem BEFORE truncation

      subroutine rtrunc(asympsi,nvtrunc,nltrunc,nchtrunc,dummy,chin
     $     ,nch)
      implicit none
      integer nvtrunc,nltrunc,nch,nchtrunc
CARR
      integer chanmax
      parameter(chanmax=165)
c$$$      double precision asympsi(nch,nch),dummy(nchtrunc,nchtrunc)
c$$$      integer chin(nch,2)
      double precision asympsi(chanmax,chanmax),dummy(chanmax,chanmax)
      integer chin(chanmax,2)

      integer i,j,itrun,jtrun

      itrun = 0
      do i = 1,nch
         if(chin(i,1).le.nvtrunc.and.chin(i,2).le.nltrunc) then
            itrun = itrun+1
            jtrun = 0
            do j = 1,nch
               if(chin(j,1).le.nvtrunc.and.chin(j,2).le.nltrunc) then
                  jtrun = jtrun+1
                  dummy(itrun,jtrun) = asympsi(i,j)
               endif
            end do
         end if
      end do
      return
      end

C *-
C * this routine combines the numerical recipes routines
C * LUDCMP, and LUBKSB to get the inverse of a matrix.
C * I have modified the routines to do double precision

      subroutine nrinverse(a,b,temp,n,np,indx)
      implicit none
      integer n,np,indx(np),i,j
      double precision a(np,np),temp(np,np),b(np,np)
      double precision d
      double precision c(n)

      do i = 1,n
         do j = 1,n
            temp(i,j) = a(i,j)
            b(i,j) = 0.d0
         end do
         b(i,i) = 1.d0
      end do
      call ludcmp(temp,n,np,indx,d)
      do j = 1,n
         do i = 1, n
            c(i) = b(i,j)
         enddo

C        call lubksb(temp,n,np,indx,b(1,j))
         call lubksb(temp,n,np,indx,c)

         do i = 1, n
            b(i,j) = c(i)
         enddo

      end do
      return
      end

C *-
C * this routine carries out a LU decomposition of
C * a matrix 'a' and overwrites 'a' in the process.
C * See the numerical recipes book p. 38

      SUBROUTINE ludcmp(a,n,np,indx,d)
      implicit none
      INTEGER n,np,indx(n),NMAX
      DOUBLE PRECISION d,a(np,np),TINY
      PARAMETER (NMAX=5000,TINY=1.0e-20)
      INTEGER i,imax,j,k
      DOUBLE PRECISION aamax,dum,sum,vv(NMAX)
      d=1.d0
      do i=1,n
         aamax=0.d0
         do j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
         end do
         if (aamax.eq.0.d0) stop 'singular matrix in ludcmp'
         vv(i)=1.d0/aamax
      end do
      do j=1,n
         do i=1,j-1
            sum=a(i,j)
            do k=1,i-1
               sum=sum-a(i,k)*a(k,j)
            end do
            a(i,j)=sum
         end do
         aamax=0.d0
         do i=j,n
            sum=a(i,j)
            do k=1,j-1
               sum=sum-a(i,k)*a(k,j)
            end do
            a(i,j)=sum
            dum=vv(i)*abs(sum)
            if (dum.ge.aamax) then
               imax=i
               aamax=dum
            endif
         end do
         if (j.ne.imax)then
            do k=1,n
               dum=a(imax,k)
               a(imax,k)=a(j,k)
               a(j,k)=dum
            end do
            d=-d
            vv(imax)=vv(j)
         endif
         indx(j)=imax
         if(a(j,j).eq.0.d0)a(j,j)=TINY
         if(j.ne.n)then
            dum=1.d0/a(j,j)
            do i=j+1,n
               a(i,j)=a(i,j)*dum
            end do
         endif
      end do
      return
      END

C *-
C * numerical recipes routine to calculate the
C * inverse of a matrix.  This routine
C * requires that ludcmp be run first

      SUBROUTINE lubksb(a,n,np,indx,b)
      implicit none
      INTEGER n,np,indx(n)
      DOUBLE PRECISION a(np,np),b(n)
      INTEGER i,ii,j,ll
      DOUBLE PRECISION sum
      ii=0
      do i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         if (ii.ne.0)then
            do j=ii,i-1
               sum=sum-a(i,j)*b(j)
            end do
         else if (sum.ne.0.d0) then
            ii=i
         endif
         b(i)=sum
      end do
      do i=n,1,-1
         sum=b(i)
         do j=i+1,n
            sum=sum-a(i,j)*b(j)
         end do
         b(i)=sum/a(i,i)
      end do
      return
      END

C *-
C * this subroutine diagonalizes the real, symmetric
C * matrix 'a', the output eigenvectors are in v
C * the eigenvalues are in d, and the upper
C * triangle of 'a' is destroyed

      SUBROUTINE jacobi(a,n,np,d,v,nrot)
      implicit none
      INTEGER n,np,nrot,NMAX
      DOUBLE PRECISION a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=5000)
      INTEGER i,ip,iq,j
      DOUBLE PRECISION c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do ip=1,n
         do iq=1,n
            v(ip,iq)=0.d0
         end do
         v(ip,ip)=1.d0
      end do
      do ip=1,n
         b(ip)=a(ip,ip)
         d(ip)=b(ip)
         z(ip)=0.d0
      end do
      nrot=0
      do i=1,500
         sm=0.d0
         do ip=1,n-1
            do iq=ip+1,n
               sm=sm+abs(a(ip,iq))
            end do
         end do
         if(sm.eq.0.d0)return
         if(i.lt.4)then
            tresh=0.2d0*sm/n**2
         else
            tresh=0.d0
         endif
         do ip=1,n-1
            do iq=ip+1,n
               g=100.d0*abs(a(ip,iq))
               if((i.gt.4).and.(abs(d(ip))+
     $              g.eq.abs(d(ip))).and.(abs(d(iq))
     $              +g.eq.abs(d(iq))))then
                  a(ip,iq)=0.d0
               else if(abs(a(ip,iq)).gt.tresh)then
                  h=d(iq)-d(ip)
                  if(abs(h)+g.eq.abs(h))then
                     t=a(ip,iq)/h
                  else
                     theta=0.5d0*h/a(ip,iq)
                     t=1.d0/(abs(theta)+sqrt(1.d0+theta**2))
                     if(theta.lt.0.d0)t=-t
                  endif
                  c=1.d0/sqrt(1.d0+t**2)
                  s=t*c
                  tau=s/(1.d0+c)
                  h=t*a(ip,iq)
                  z(ip)=z(ip)-h
                  z(iq)=z(iq)+h
                  d(ip)=d(ip)-h
                  d(iq)=d(iq)+h
                  a(ip,iq)=0.d0
                  do j=1,ip-1
                     g=a(j,ip)
                     h=a(j,iq)
                     a(j,ip)=g-s*(h+g*tau)
                     a(j,iq)=h+s*(g-h*tau)
                  end do
                  do j=ip+1,iq-1
                     g=a(ip,j)
                     h=a(j,iq)
                     a(ip,j)=g-s*(h+g*tau)
                     a(j,iq)=h+s*(g-h*tau)
                  end do
                  do j=iq+1,n
                     g=a(ip,j)
                     h=a(iq,j)
                     a(ip,j)=g-s*(h+g*tau)
                     a(iq,j)=h+s*(g-h*tau)
                  end do
                  do j=1,n
                     g=v(j,ip)
                     h=v(j,iq)
                     v(j,ip)=g-s*(h+g*tau)
                     v(j,iq)=h+s*(g-h*tau)
                  end do
                  nrot=nrot+1
               endif
            end do
         end do
         do ip=1,n
            b(ip)=b(ip)+z(ip)
            d(ip)=b(ip)
            z(ip)=0.d0
         end do
      end do
      stop 'too many iterations in jacobi'
      return
      END

C *-
C * the next roughly 200 lines are for numerical recipes bessel routines

      SUBROUTINE beschb(x,gam1,gam2,gampl,gammi)
      implicit none
      INTEGER NUSE1,NUSE2
      DOUBLE PRECISION gam1,gam2,gammi,gampl,x
      PARAMETER (NUSE1=7,NUSE2=8)
      DOUBLE PRECISION xx,c1(7),c2(8),chebev
      SAVE c1,c2
      DATA c1/-1.142022680371172d0,6.516511267076d-3,3.08709017308d-4,
     $     -3.470626964d-6,6.943764d-9,3.6780d-11,-1.36d-13/
      DATA c2/1.843740587300906d0,-.076852840844786d0,1.271927136655d-3,
     $     -4.971736704d-6,-3.3126120d-8,2.42310d-10,-1.70d-13,-1.d-15/
      xx=8.d0*x*x-1.d0
      gam1=chebev(-1.d0,1.d0,c1,NUSE1,xx)
      gam2=chebev(-1.d0,1.d0,c2,NUSE2,xx)
      gampl=gam2-x*gam1
      gammi=gam2+x*gam1
      return
      END

      SUBROUTINE bessik(x,xnu,ri,rk,rip,rkp)
      implicit none
      INTEGER MAXIT
      DOUBLE PRECISION ri,rip,rk,rkp,x,xnu,XMIN
      DOUBLE PRECISION EPS,FPMIN,PI
      PARAMETER (eps=1.d-16,FPMIN=1.d-300,MAXIT=10000,XMIN=2.d0, PI
     $     =3.141592653589793d0)
      INTEGER i,l,nl
      DOUBLE PRECISION a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,
     $     gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu
     $     ,rip1, ripl,ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2
     $     ,xmu,xmu2
      if(x.le.0.d0.or.xnu.lt.0.d0) stop 'bad arguments in bessik'
      nl=int(xnu+.5d0)
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      h=xnu*xi
      if(h.lt.FPMIN)h=FPMIN
      b=xi2*xnu
      d=0.d0
      c=h
      do i=1,MAXIT
         b=b+xi2
         d=1.d0/(b+d)
         c=b+1.d0/c
         del=c*d
         h=del*h
         if(abs(del-1.d0).lt.EPS)goto 1
      end do
      stop 'x too large in bessik; try asymptotic expansion'
 1    continue
      ril=FPMIN
      ripl=h*ril
      ril1=ril
      rip1=ripl
      fact=xnu*xi
      do l=nl,1,-1
         ritemp=fact*ril+ripl
         fact=fact-xi
         ripl=fact*ritemp+ril
         ril=ritemp
      end do
      f=ripl/ril
      if(x.lt.XMIN) then
         x2=.5d0*x
         pimu=PI*xmu
         if(abs(pimu).lt.EPS)then
            fact=1.d0
         else
            fact=pimu/sin(pimu)
         endif
         d=-log(x2)
         e=xmu*d
         if(abs(e).lt.EPS)then
            fact2=1.d0
         else
            fact2=sinh(e)/e
         endif
         call beschb(xmu,gam1,gam2,gampl,gammi)
         ff=fact*(gam1*cosh(e)+gam2*fact2*d)
         sum=ff
         e=exp(e)
         p=0.5d0*e/gampl
         q=0.5d0/(e*gammi)
         c=1.d0
         d=x2*x2
         sum1=p
         do i=1,MAXIT
            ff=(i*ff+p+q)/(1.d0*i*i-xmu2)
            c=c*d/i
            p=p/(1.d0*i-xmu)
            q=q/(1.d0*i+xmu)
            del=c*ff
            sum=sum+del
            del1=c*(p-i*ff)
            sum1=sum1+del1
            if(abs(del).lt.abs(sum)*EPS)goto 2
         end do
         stop 'bessk series failed to converge'
 2       continue
         rkmu=sum
         rk1=sum1*xi2
      else
         b=2.d0*(1.d0+x)
         d=1.d0/b
         delh=d
         h=delh
         q1=0.d0
         q2=1.d0
         a1=.25d0-xmu2
         c=a1
         q=c
         a=-a1
         s=1.d0+q*delh
         do i=2,MAXIT
            a=a-2.d0*(1.d0*i-1.d0)
            c=-a*c/i
            qnew=(q1-b*q2)/a
            q1=q2
            q2=qnew
            q=q+c*qnew
            b=b+2.d0
            d=1.d0/(b+a*d)
            delh=(b*d-1.d0)*delh
            h=h+delh
            dels=q*delh
            s=s+dels
            if(abs(dels/s).lt.EPS)goto 3
         end do
         stop 'bessik: failure to converge in cf2'
 3       continue
         h=a1*h
         rkmu=sqrt(PI/(2.d0*x))*exp(-x)/s
         rk1=rkmu*(xmu+x+.5d0-h)*xi
      endif
      rkmup=xmu*xi*rkmu-rk1
      rimu=xi/(f*rkmu-rkmup)
      ri=(rimu*ril1)/ril
      rip=(rimu*rip1)/ril
      do i=1,nl
         rktemp=(xmu+1.d0*i)*xi2*rk1+rkmu
         rkmu=rk1
         rk1=rktemp
      end do
      rk=rkmu
      rkp=xnu*xi*rkmu-rk1
      return
      END

      SUBROUTINE bessjy(x,xnu,rj,ry,rjp,ryp)
      implicit none
      INTEGER MAXIT
      DOUBLE PRECISION rj,rjp,ry,ryp,x,xnu,XMIN
      DOUBLE PRECISION EPS,FPMIN,PI
      PARAMETER (eps=1.d-16,FPMIN=1.d-300,MAXIT=10000,XMIN=2.d0, PI
     $     =3.141592653589793d0)
      INTEGER i,isign,l,nl
      DOUBLE PRECISION a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e,
     $     f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu
     $     ,pimu2,q, r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup
     $     ,rytemp,sum,sum1, temp,w,x2,xi,xi2,xmu,xmu2
      if(x.le.0.d0.or.xnu.lt.0.d0) stop 'bad arguments in bessjy'
      if(x.lt.XMIN)then
         nl=int(xnu+.5d0)
      else
         nl=max(0,int(xnu-x+1.5d0))
      endif
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      w=xi2/PI
      isign=1
      h=xnu*xi
      if(h.lt.FPMIN)h=FPMIN
      b=xi2*xnu
      d=0.d0
      c=h
      do i=1,MAXIT
         b=b+xi2
         d=b-d
         if(abs(d).lt.FPMIN)d=FPMIN
         c=b-1.d0/c
         if(abs(c).lt.FPMIN)c=FPMIN
         d=1.d0/d
         del=c*d
         h=del*h
         if(d.lt.0.d0)isign=-isign
         if(abs(del-1.d0).lt.EPS)goto 1
      end do
      stop 'x too large in bessjy; try asymptotic expansion'
 1    continue
      rjl=isign*FPMIN
      rjpl=h*rjl
      rjl1=rjl
      rjp1=rjpl
      fact=xnu*xi
      do l=nl,1,-1
         rjtemp=fact*rjl+rjpl
         fact=fact-xi
         rjpl=fact*rjtemp-rjl
         rjl=rjtemp
      end do
      if(rjl.eq.0.d0)rjl=EPS
      f=rjpl/rjl
      if(x.lt.XMIN) then
         x2=.5d0*x
         pimu=PI*xmu
         if(abs(pimu).lt.EPS)then
            fact=1.d0
         else
            fact=pimu/sin(pimu)
         endif
         d=-log(x2)
         e=xmu*d
         if(abs(e).lt.EPS)then
            fact2=1.d0
         else
            fact2=sinh(e)/e
         endif
         call beschb(xmu,gam1,gam2,gampl,gammi)
         ff=2.d0/PI*fact*(gam1*cosh(e)+gam2*fact2*d)
         e=exp(e)
         p=e/(gampl*PI)
         q=1.d0/(e*PI*gammi)
         pimu2=0.5d0*pimu
         if(abs(pimu2).lt.EPS)then
            fact3=1.d0
         else
            fact3=sin(pimu2)/pimu2
         endif
         r=PI*pimu2*fact3*fact3
         c=1.d0
         d=-x2*x2
         sum=ff+r*q
         sum1=p
         do i=1,MAXIT
            ff=(i*ff+p+q)/(1.d0*i*i-xmu2)
            c=c*d/i
            p=p/(1.d0*i-xmu)
            q=q/(1.d0*i+xmu)
            del=c*(ff+r*q)
            sum=sum+del
            del1=c*p-i*del
            sum1=sum1+del1
            if(abs(del).lt.(1.d0+abs(sum))*EPS)goto 2
         end do
         stop 'bessy series failed to converge'
 2       continue
         rymu=-sum
         ry1=-sum1*xi2
         rymup=xmu*xi*rymu-ry1
         rjmu=w/(rymup-f*rymu)
      else
         a=.25d0-xmu2
         p=-.5d0*xi
         q=1.d0
         br=2.d0*x
         bi=2.d0
         fact=a*xi/(p*p+q*q)
         cr=br+q*fact
         ci=bi+p*fact
         den=br*br+bi*bi
         dr=br/den
         di=-bi/den
         dlr=cr*dr-ci*di
         dli=cr*di+ci*dr
         temp=p*dlr-q*dli
         q=p*dli+q*dlr
         p=temp
         do i=2,MAXIT
            a=a+2.d0*(1.d0*i-1.d0)
            bi=bi+2.d0
            dr=a*dr+br
            di=a*di+bi
            if(abs(dr)+abs(di).lt.FPMIN)dr=FPMIN
            fact=a/(cr*cr+ci*ci)
            cr=br+cr*fact
            ci=bi-ci*fact
            if(abs(cr)+abs(ci).lt.FPMIN)cr=FPMIN
            den=dr*dr+di*di
            dr=dr/den
            di=-di/den
            dlr=cr*dr-ci*di
            dli=cr*di+ci*dr
            temp=p*dlr-q*dli
            q=p*dli+q*dlr
            p=temp
            if(abs(dlr-1.d0)+abs(dli).lt.EPS)goto 3
         end do
         stop 'cf2 failed in bessjy'
 3       continue
         gam=(p-f)/q
         rjmu=sqrt(w/((p-f)*gam+q))
         rjmu=sign(rjmu,rjl)
         rymu=rjmu*gam
         rymup=rymu*(p+q/gam)
         ry1=xmu*xi*rymu-rymup
      endif
      fact=rjmu/rjl
      rj=rjl1*fact
      rjp=rjp1*fact
      do i=1,nl
         rytemp=(xmu+1.d0*i)*xi2*ry1-rymu
         rymu=ry1
         ry1=rytemp
      end do
      ry=rymu
      ryp=xnu*xi*rymu-ry1
      return
      END

      double precision FUNCTION chebev(a,b,c,m,x)
      implicit none
      INTEGER m
      DOUBLE PRECISION a,b,x,c(m)
      INTEGER j
      DOUBLE PRECISION d,dd,sv,y,y2
      if ((x-a)*(x-b).gt.0.d0) stop 'x not in range in chebev'
      d=0.d0
      dd=0.d0
      y=(2.d0*x-a-b)/(b-a)
      y2=2.d0*y
      do j=m,2,-1
         sv=d
         d=y2*d-dd+c(j)
         dd=sv
      end do
      chebev=y*d-dd+0.5d0*c(1)
      return
      END

C *-
C * Numerical Recipes
C * function to give the derivative of a Ricatti Neumann function
      double precision function rneup(k,r,l)
      implicit none
      double precision ybes,jbes,ybesp,jbesp
      integer l
      double precision k,r
      double precision pi
      pi=4.d0*datan(1.d0)
      call bessjy(k*r,1.d0*l+0.5d0,jbes,ybes,jbesp,ybesp)
      if(l.eq.0) then
         rneup=k*dsin(k*r)
      else
         rneup=dsqrt(pi*r*k/2.d0)*
     +        ((.5d0/r)*ybes+k*ybesp)
      endif
      return
      end

C * function to give Ricatti Neumann function
      double precision function rneu(k,r,l)
      implicit none
      integer l
      double precision ybes,jbes,ybesp,jbesp
      double precision k,r
      double precision pi
      pi=4.d0*datan(1.d0)
      call bessjy(k*r,1.d0*l+.5d0,jbes,ybes,jbesp,ybesp)
      rneu=dsqrt(pi*r*k/2.d0)*ybes
      return
      end

C * function to give the derivative of a Ricatti bessel function
      double precision function rbesp(k,r,l)
      implicit none
      double precision ybes,jbes,ybesp,jbesp
      integer l
      double precision k,r
      double precision pi
      pi=4.d0*datan(1.d0)
      call bessjy(k*r,1.d0*l+.5d0,jbes,ybes,jbesp,ybesp)
      if(l.eq.0) then
         rbesp=k*dcos(k*r)
      else
         rbesp=dsqrt(pi*k*r/2.d0)*((.5d0/r)*jbes+k*jbesp)
      endif
      return
      end

C * function to give Ricatti bessel function
      double precision function rbes(k,r,l)
      implicit none
      double precision ybes,jbes,ybesp,jbesp
      integer l
      double precision k,r
      double precision pi
      pi=4.d0*datan(1.d0)
      call bessjy(k*r,1.d0*l+.5d0,jbes,ybes,jbesp,ybesp)
      rbes=dsqrt(pi*k*r/2.d0)*jbes
      return
      end

C *-
C * function to give the derivative of a Ricatti bessel function
C *                           _                _
C *                         d|   -l+1  ^        |
C * this function produces --|  i      j (ikr)  |
C *                        dr|_         l      _|

      double precision function irbesp(k,r,l)
      implicit none
      double precision ibes,kbes,ibesp,kbesp
      integer l
      double precision k,r
      double precision pi
      pi=4.d0*datan(1.d0)
      call bessik(k*r,1.d0*l+0.5d0,ibes,kbes,ibesp,kbesp)
      if(l.eq.0) then
         irbesp=-k*(dexp(-k*r)+dexp(k*r))/2.d0
      else
         irbesp=-dsqrt(pi*k*r/2.d0)*((.5d0/r)*ibes+k*ibesp)
      endif
      return
      end

C * function to give Ricatti bessel function with an imaginary argument
C *                           -l+1 ^
C * this function produces   i     j (ikr)
C *                                 l

      double precision function irbes(k,r,l)
      implicit none
      double precision ibes,kbes,ibesp,kbesp
      integer l
      double precision k,r
      double precision pi
      pi=4.d0*datan(1.d0)
      call bessik(k*r,1.d0*l+.5d0,ibes,kbes,ibesp,kbesp)
      irbes=-dsqrt(pi*k*r/2.d0)*ibes
      return
      end

C * function to give the derivative of a Ricatti Hankel function
C *                           _               _
C *                         d|   l+1  ^(+)     |
C * this function produces --|  i     h (ikr)  |
C *                        dr|_        l      _|

      double precision function rhankp(k,r,l)
      implicit none
      double precision ibes,kbes,ibesp,kbesp
      integer l
      double precision k,r
      double precision pi
      pi=4.d0*datan(1.d0)
      call bessik(k*r,1.d0*l+.5d0,ibes,kbes,ibesp,kbesp)
      if(l.eq.0) then
         rhankp=-dexp(-k*r)
      else
         rhankp=dsqrt(2.d0*k*r/pi)*((.5d0/r)*kbes+k*kbesp)
      endif
      return
      end

C * function to give Ricatti Hankel function
C *                           l+1  ^(+)
C * this function produces   i     h (ikr)
C *                                 l

      double precision function rhank(k,r,l)
      implicit none
      double precision ibes,kbes,ibesp,kbesp
      integer l
      double precision k,r
      double precision pi
      pi=4.d0*datan(1.d0)
      call bessik(k*r,1.d0*l+.5d0,ibes,kbes,ibesp,kbesp)
      rhank=dsqrt(2.d0*k*r/pi)*kbes
      return
      end
C *-/

C *-
C * judge if the number of partial waves is enough

      subroutine judgenpwav(nbound,nlproj,nexdim)
      implicit none
      
      integer boundmax
      parameter (boundmax=6)
      integer nexdim,nlproj(boundmax)

      integer nbound,i,nlmax,nexdimc,factor

      nlmax = 0
      do i = 1, nbound
         if(nlproj(i) .gt. nlmax) nlmax = nlproj(i)
      enddo

      if(nexdim .lt. nlmax) then
         write (0,*)
         write (0,*) 'nexdim is TOO SMALL and should be greater than'
     $        ,nlmax
         write (0,*) "   nexdim = npwavxin for cal. exchange kernel!"
         write(0,*) 'or nlproj is TOO LARGE and should be less than'
     $        ,nexdim
         stop
      endif

      return
      end
C *-

