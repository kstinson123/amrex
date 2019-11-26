
#include "myfunc.H"
#include "myfunc_F.H"
#include <AMReX_BCUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>

#include <iostream>
#include <AMReX_PlotFileUtil.H>

void SDC_advance(MultiFab& phi_old,
		 MultiFab& phi_new,
		 std::array<MultiFab, AMREX_SPACEDIM>& flux,
		 Real dt,
		 const Geometry& geom,
		 const Vector<BCRec>& bc,
		 MLMG&  mlmg,
		 MLABecLaplacian& mlabec,
		 SDCstruct &SDC, Real a, Real d, Real r,
         std::array<MultiFab,AMREX_SPACEDIM>& face_bcoef,
         std::array<MultiFab,AMREX_SPACEDIM>& prod_stor, Real time, Real epsilon, Real k_freq, Real kappa, MultiFab& bdry_values, MLABecLaplacian& mlabec_BCfill, int Nprob)
{

  /*  This is a multi-implicit SDC example time step for an 
      advection-diffusion-reaction equation of the form

      phi_t = A(phi)+D(phi)+R(phi)
      
      The advection is treated explicilty and the diffusion and reaction implicitly
      and uncoupled

      the constants a,d, and r control the strength of each term
  */
    
    // Need to advance time in here as we have a forcing term: So we make nodal fraction which you use as you step through the array.
    // std::cout << std::setprecision(30), Print() << SetPrecision(30) << …stuff…, Print().SetPrecision(30) << a << std::endl;

  Real qij;
  Real current_time = time;
  const BoxArray &ba=phi_old.boxArray();
  const DistributionMapping &dm=phi_old.DistributionMap();
  // Copy old phi into first SDC node
  MultiFab::Copy(SDC.sol[0],phi_old, 0, 0, 1, 2);

  // Fill the ghost cells of each grid from the other grids
  // includes periodic domain boundaries
    // Fill Dirichlet values
    bdry_values.setVal(0,0); bdry_values.setBndry(0);
    if(Nprob<3){
        for ( MFIter mfi(bdry_values); mfi.isValid(); ++mfi )
        {          const Box& bx = mfi.validbox();
            fill_bdry_values(BL_TO_FORTRAN_BOX(bx),
                             BL_TO_FORTRAN_ANYD(bdry_values[mfi]),
                             geom.CellSize(), geom.ProbLo(), geom.ProbHi(),&current_time, &epsilon,&k_freq, &kappa, &Nprob);
        }
       
        mlabec_BCfill.fourthOrderBCFill(SDC.sol[0],bdry_values);
    }
    // Fill periodic values and interal ghost cells
    SDC.sol[0].FillBoundary(geom.periodicity());
    

  /*  for ( MFIter mfi(bdry_values); mfi.isValid(); ++mfi )
    {          const Box& bx = mfi.validbox();
        
        print_multifab(BL_TO_FORTRAN_ANYD(SDC.sol[0][mfi]));
    }*/
    // Fill non-periodic physical boundaries (doesn't do Dirichlet; probably only one ghost cell).
  // FillDomainBoundary(SDC.sol[0], geom, bc);
  
  //  Compute the first function value (Need boundary filled here at current time explicit).
  int sdc_m=0;
  SDC_feval(flux,geom,bc,SDC,a,d,r,face_bcoef,prod_stor,sdc_m,-1,time, epsilon, k_freq, kappa, Nprob);
    
    
  // Copy first function value to all nodes
  for (int sdc_n = 1; sdc_n < SDC.Nnodes; sdc_n++)
    {
        
      //        current_time = time+dt*nodeFrac[sdc_n];
        current_time = time+dt*SDC.qnodes[sdc_n];
      //MultiFab::Copy(SDC.f[0][sdc_n],SDC.f[0][0], 0, 0, 1, 0);
        // Subsequent Calculation doesn't need BC conditions
     SDC_feval(flux,geom,bc,SDC,a,d,r,face_bcoef,prod_stor,sdc_n,0,current_time, epsilon, k_freq, kappa, Nprob);
    /*    BoxArray bam(geom.Domain());
        MultiFab Marc(bam,DistributionMapping(bam),1,0);
        Marc.copy(SDC.f[0][0]);
        ParmParse pp;
        std::string myName;
        pp.query("myName",myName);
        //BoxArray bam1 = BoxArray(SDC.f[0][0].boxArray()).grow(2); //./diffmultifab2d.gnu.TEST.ex infile1=JUNK infile2=JUNK1
        //MultiFab Marc1(bam1,DistributionMapping(bam1),1,0);       // ~/Amrvis/amrvis2d.gnu.ex -mf JUNK
        //Marc1.setVal(3);                                     
        VisMF::Write(Marc,myName);
        Abort("MARC");*/
      MultiFab::Copy(SDC.f[1][sdc_n],SDC.f[1][0], 0, 0, 1, 0);
      if (SDC.Npieces==3)
	MultiFab::Copy(SDC.f[2][sdc_n],SDC.f[2][0], 0, 0, 1, 0);      
    }


  //  Now do the actual sweeps
  for (int k=1; k <= SDC.Nsweeps; ++k)
    {

      //  Compute RHS integrals
      SDC.SDC_rhs_integrals(dt);
        
      //  Substep over SDC nodes
      for (int sdc_m = 0; sdc_m < SDC.Nnodes-1; sdc_m++)
	{
	  amrex::Print() << "sweep " << k << ", substep " << sdc_m+1 <<"\n";
        
	  // use phi_new as rhs and fill it with terms at this iteration
	  SDC.SDC_rhs_k_plus_one(phi_new,dt,sdc_m);
	  
	  // get the best initial guess for implicit solve
	  MultiFab::Copy(SDC.sol[sdc_m+1],phi_new, 0, 0, 1, 2);
	  for ( MFIter mfi(SDC.sol[sdc_m+1]); mfi.isValid(); ++mfi )
	    {
	      //	      const Box& bx = mfi.validbox();
	      qij = dt*SDC.Qimp[sdc_m][sdc_m+1];
	      SDC.sol[sdc_m+1][mfi].saxpy(qij,SDC.f[1][sdc_m+1][mfi]);
	    }
     // SDC.sol[sdc_m+1].FillBoundary(geom.periodicity());
        // Get Implicit time value
	  //        current_time = time+dt*nodeFrac[sdc_m+1];
        current_time = time+dt*SDC.qnodes[sdc_m+1];
	
        // Fill Dirichlet Values
        if (Nprob<3){
            for ( MFIter mfi(bdry_values); mfi.isValid(); ++mfi )
            {          const Box& bx = mfi.validbox();
                fill_bdry_values(BL_TO_FORTRAN_BOX(bx),
                                 BL_TO_FORTRAN_ANYD(bdry_values[mfi]),
                                 geom.CellSize(), geom.ProbLo(), geom.ProbHi(),&current_time, &epsilon,&k_freq, &kappa, &Nprob);
            }
            mlabec_BCfill.fourthOrderBCFill(SDC.sol[sdc_m+1],bdry_values);
        }
        // Fill periodic values and interal ghost cells
        SDC.sol[sdc_m+1].FillBoundary(geom.periodicity());
        
        
        
	  // Solve for the first implicit part
        
	  SDC_fcomp(phi_new, flux, geom, bc, SDC, mlmg, mlabec,dt,a,d,r,face_bcoef,prod_stor,sdc_m+1,1, current_time, epsilon, k_freq, kappa, mlabec_BCfill, k, Nprob);
     
        
        

	  if (SDC.Npieces==3)
	    {
	      // Build rhs for 2nd solve
	      MultiFab::Copy(phi_new, SDC.sol[sdc_m+1],0, 0, 1, 2);
	      
	      // Add in the part for the 2nd implicit term to rhs
	      SDC.SDC_rhs_misdc(phi_new,dt,sdc_m);
	      
	      // Solve for the second implicit part
	      SDC_fcomp(phi_new, flux, geom, bc, SDC, mlmg, mlabec,dt, a,d,r,face_bcoef,prod_stor,sdc_m+1,2, current_time, epsilon, k_freq, kappa, mlabec_BCfill, k, Nprob);
	    }
	  // Compute the function values at node sdc_m+1
        
        // Shouldn't need the following BC code as the solve should linearly maintain it.
        if (Nprob<3){
            for ( MFIter mfi(bdry_values); mfi.isValid(); ++mfi )
            {          const Box& bx = mfi.validbox();
                fill_bdry_values(BL_TO_FORTRAN_BOX(bx),
                                 BL_TO_FORTRAN_ANYD(bdry_values[mfi]),
                                 geom.CellSize(), geom.ProbLo(), geom.ProbHi(),&current_time, &epsilon,&k_freq, &kappa, &Nprob);
            }
            mlabec_BCfill.fourthOrderBCFill(SDC.sol[sdc_m+1],bdry_values);
        }
        SDC.sol[sdc_m+1].FillBoundary(geom.periodicity());
        //       mlabec.fillSolutionBC(0, SDC.sol[sdc_m+1], &bdry_values);
      //  amrex::Print() << "current time" << current_time <<"\n";
	//  
	  SDC_feval(flux,geom,bc,SDC,a,d,r,face_bcoef,prod_stor,
                sdc_m+1,-1,current_time,epsilon, k_freq, kappa, Nprob);

	} // end SDC substep loop
        
      
        
    }  // end sweeps loop
    ///////////////////////////////////////////////////////
    // PAUSE
    ///////////////////////////////////////////////////////
   // std::cin.get();
    
    
  // Return the last node in phi_new
  MultiFab::Copy(phi_new, SDC.sol[SDC.Nnodes-1], 0, 0, 1, 2);
    
    

}

void SDC_feval(std::array<MultiFab, AMREX_SPACEDIM>& flux,
	       const Geometry& geom,
	       const Vector<BCRec>& bc,
	       SDCstruct &SDC,
	       Real a,Real d,Real r,
           std::array<MultiFab, AMREX_SPACEDIM>& face_bcoef,
           std::array<MultiFab, AMREX_SPACEDIM>& prod_stor,
	       int sdc_m,int npiece, Real time, Real epsilon, Real k_freq, Real kappa, int Nprob)
{
  /*  Evaluate explicitly the rhs terms of the equation at the SDC node "sdc_m".
      The input parameter "npiece" describes which term to do.  
      If npiece = -1, do all the pieces */
   
    
  const BoxArray &ba=SDC.sol[0].boxArray();
  const DistributionMapping &dm=SDC.sol[0].DistributionMap();

  const Box& domain_bx = geom.Domain();
  const Real* dx = geom.CellSize();
  int nlo,nhi;
  if (npiece < 0)
    {
      nlo=0;
      nhi=SDC.Npieces;
    }
  else
    {
      nlo=npiece;
      nhi=npiece+1;
    }

    
    
    
  for ( MFIter mfi(SDC.sol[sdc_m]); mfi.isValid(); ++mfi )
    {
      const Box& bx = mfi.validbox();
      for (int n = nlo; n < nhi; n++)    
	{
	  SDC_feval_F(BL_TO_FORTRAN_BOX(bx),
		      BL_TO_FORTRAN_BOX(domain_bx),
		      BL_TO_FORTRAN_ANYD(SDC.sol[sdc_m][mfi]),
		      BL_TO_FORTRAN_ANYD(flux[0][mfi]),
		      BL_TO_FORTRAN_ANYD(flux[1][mfi]),
#if (AMREX_SPACEDIM == 3)   
		      BL_TO_FORTRAN_ANYD(flux[2][mfi]),
#endif		       
		      BL_TO_FORTRAN_ANYD(SDC.f[n][sdc_m][mfi]),
		      dx,&a,&d,&r,
              BL_TO_FORTRAN_ANYD(face_bcoef[0][mfi]),
              BL_TO_FORTRAN_ANYD(face_bcoef[1][mfi]),
              BL_TO_FORTRAN_ANYD(prod_stor[0][mfi]),
              BL_TO_FORTRAN_ANYD(prod_stor[1][mfi]),
              &n, &time, &epsilon, &k_freq, &kappa, &Nprob);
	}
      
    }
}
void SDC_fcomp(MultiFab& rhs,
           std::array<MultiFab, AMREX_SPACEDIM>& flux,
	       const Geometry& geom,
	       const Vector<BCRec>& bc,
	       SDCstruct &SDC,
	       MLMG &mlmg,
	       MLABecLaplacian& mlabec,	      
	       Real dt,Real a,Real d,Real r,
           std::array<MultiFab,AMREX_SPACEDIM>& face_bcoef,
           std::array<MultiFab,AMREX_SPACEDIM>& prod_stor,
	       int sdc_m,int npiece, Real time, Real epsilon, Real k_freq, Real kappa, MLABecLaplacian& mlabec_BCfill, int numsweeps, int Nprob)
{
  /*  Solve implicitly for the implicit terms of the equation at the SDC node "sdc_m".
      The input parameter "npiece" describes which term to do.  */
    ParmParse pp;
    
  const BoxArray &ba=SDC.sol[0].boxArray();
  const DistributionMapping &dm=SDC.sol[0].DistributionMap();

  const Box& domain_bx = geom.Domain();
  const Real* dx = geom.CellSize();
  Real qij;
  Real t;
    
    // relative and absolute tolerances for linear solve
  const Real tol_rel = 1.e-12;
  const Real tol_abs = 0.0;
  const Real tol_res = 1.e-12;    // Tolerance on residual
  Real resnorm = 1.e10;    // Tolerance on residual
    Real zeroReal = 0.0;
    Real corrnorm;
  // Make some space for iteration stuff
  MultiFab corr(ba, dm, 1, 2);
  MultiFab resid(ba, dm, 1, 2);
  MultiFab eval_storage(ba, dm, 1, 0);
    
 // Temp storage for checking code
    MultiFab temp_resid(ba, dm, 1, 2);
    MultiFab temp_corr(ba, dm, 1, 2);
    MultiFab temp_fab(ba, dm, 1, 2);
    MultiFab temp_zero(ba, dm, 1, 1);
    MultiFab temp_err(ba, dm, 1, 2);
    MultiFab temp_noghost(ba, dm, 1, 0);
    MultiFab temp_oneghost(ba, dm, 1, 1);
    temp_zero.setVal(0.0); temp_zero.setBndry(0);
    
    
    ///////////////////////////////////////////////////////////
    // Set boundary values
    ///////////////////////////////////////////////////////////
    
    MultiFab bdry_values(ba, dm, 1, 1);
    bdry_values.setVal(0.0); bdry_values.setBndry(0);
    if (Nprob < 3){
        for ( MFIter mfi(bdry_values); mfi.isValid(); ++mfi )
        {          const Box& bx = mfi.validbox();
            fill_bdry_values(BL_TO_FORTRAN_BOX(bx),
                             BL_TO_FORTRAN_ANYD(bdry_values[mfi]),
                             geom.CellSize(), geom.ProbLo(), geom.ProbHi(),&time, &epsilon,&k_freq, &kappa, &Nprob);
        }
        mlabec_BCfill.fourthOrderBCFill(SDC.sol[sdc_m],bdry_values);
    }
    
    SDC.sol[sdc_m].FillBoundary(geom.periodicity());
    
  if (npiece == 1)  
    {
        
            
        
        
        
        //////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////
        
        // Do diffusion solve
    
        // Fill the ghost cells of each grid from the other grids
        // includes periodic domain boundaries
    //    rhs.FillBoundary(geom.periodicity());
    //    SDC.sol[sdc_m].FillBoundary(geom.periodicity());
        
        // Fill non-periodic physical boundaries
       // FillDomainBoundary(rhs, geom, bc);
       // FillDomainBoundary(SDC.sol[sdc_m], geom, bc);
        
        //  Set diffusion scalar in solve
        qij = dt*SDC.Qimp[sdc_m-1][sdc_m];
	Real ascalar = 1.0;
        mlabec.setScalars(ascalar, d*qij);
        
        // set the boundary conditions
      //  mlabec.setLevelBC(0, &rhs);
      //  mlabec.setLevelBC(0, &SDC.sol[sdc_m]);
        // set level BC to 0
        int resk=0;
        
        int maxresk;
        if (Nprob<3){
            maxresk = 1500;
        }
        else {
            maxresk = 10;
        }
        pp.query("maxresk",maxresk);
        int numGS=10;
        if (Nprob<3){
            numGS = 10;
        }
        else {
            numGS = 1000;
        }
        pp.query("numGS",numGS);
        
        while ((resnorm > tol_res) & (resk <=maxresk))
        {
           
            for ( MFIter mfi(SDC.sol[sdc_m]); mfi.isValid(); ++mfi )
            {
                const Box& bx = mfi.validbox();
                SDC_feval_F(BL_TO_FORTRAN_BOX(bx),
                            BL_TO_FORTRAN_BOX(domain_bx),
                            BL_TO_FORTRAN_ANYD(SDC.sol[sdc_m][mfi]),
                            BL_TO_FORTRAN_ANYD(flux[0][mfi]),
                            BL_TO_FORTRAN_ANYD(flux[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                            BL_TO_FORTRAN_ANYD(flux[2][mfi]),
#endif
                            BL_TO_FORTRAN_ANYD(eval_storage[mfi]),
                            dx,&a,&d,&r,
                            BL_TO_FORTRAN_ANYD(face_bcoef[0][mfi]),
                            BL_TO_FORTRAN_ANYD(face_bcoef[1][mfi]),
                            BL_TO_FORTRAN_ANYD(prod_stor[0][mfi]),
                            BL_TO_FORTRAN_ANYD(prod_stor[1][mfi]),
                            &npiece, &zeroReal, &zeroReal, &zeroReal, &zeroReal, &Nprob);
                
            }
            //Rescale eval_storage to make resid.
            
            
            resid.setVal(0.0);
            MultiFab::Saxpy(resid,qij,eval_storage,0,0,1,0);
            MultiFab::Saxpy(resid,1.0,rhs,0,0,1,0);
            MultiFab::Saxpy(resid,-1.0,SDC.sol[sdc_m],0,0,1,0);
            
            //corrnorm=eval_storage.norm0();
           // amrex::Print() << "iter " << resk << ",  Diffusion operator norm " << corrnorm << "\n";
            
            resnorm=resid.norm0();
            ++resk;
            
           if(resnorm <= tol_res){
               amrex::Print() << "Reached tolerance: iter " << resk << ",  residual norm " << resnorm << "\n";
                
                break;
            }
            
        
            amrex::Print() << "iter " << resk << ",  residual norm " << resnorm << "\n";
            // includes periodic domain boundaries
           // resid.FillBoundary(geom.periodicity());
            
            // Fill non-periodic physical boundaries
           // FillDomainBoundary(resid, geom, bc);
           
            //  Do the multigrid solve
            //mlmg.solve({&SDC.sol[sdc_m]}, {&rhs}, tol_rel, tol_abs);
            //MultiFab::Copy(corr,SDC.sol[sdc_m], 0, 0, 1, 0);
            // set the boundary conditions, which are homogeneous
 /*
           // mlabec.setLevelBC(0, &resid);
            temp_noghost.setVal(0.5);
            temp_oneghost.setVal(.33333);
            temp_fab.setVal(2.0);
            
       //     MultiFab::Saxpy(temp_fab,1.0,temp_oneghost,0,0,1,1);
            
            
            
            BoxArray bam(geom.Domain());
            MultiFab Marc(bam,DistributionMapping(bam),1,2);
            Marc.copy(temp_fab);
            ParmParse pp;
            std::string myName;
            pp.query("myName",myName);
            //BoxArray bam1 = BoxArray(SDC.f[0][0].boxArray()).grow(2); //./diffmultifab2d.gnu.TEST.ex infile1=JUNK infile2=JUNK1
            //MultiFab Marc1(bam1,DistributionMapping(bam1),1,0);       // ~/Amrvis/amrvis2d.gnu.ex -mf JUNK
            //Marc1.setVal(3);
            VisMF::Write(Marc,myName);
            Abort("MARC");
    */
            
            
            
            
            
            
            
            //mlmg.setFixedIter(3);
            mlabec.setLevelBC(0, &temp_zero);
    //        mlmg.solve({&corr}, {&resid}, tol_rel, tol_abs);
            
            corr.setVal(0.0);
       
            /////////////////////////////////////////////////////////////////
            // SMOOTHER
            /////////////////////////////////////////////////////////////////
            for(int g = 1; g<=numGS;g++){
	      //                if (Nprob<3){mlabec.fourthOrderBCFill(corr,temp_zero); }
	      if (Nprob<3){mlabec_BCfill.fourthOrderBCFill(corr,temp_zero); }
             corr.FillBoundary(geom.periodicity());
            
	     mlmg.setFixedIter(3);                
	     mlabec.prepareForSolve();
	     mlmg.solve({&corr}, {&resid}, tol_rel, tol_abs);	     
	     //mlabec.Fsmooth(0, 0, corr , resid, 0);
            }
            //////////////////////////////////////////////////////////////////
            ////////////////////////////
            // Box check
            //////////////////////////////
            
            ////////////////////////////
            // Box check
            //////////////////////////////
            
            
            
         //   mlabec_BCfill.fourthOrderBCFill(SDC.sol[sdc_m],bdry_values);
         //   SDC.sol[sdc_m].FillBoundary(geom.periodicity());
            
          //  mlabec_BCfill.fourthOrderBCFill(corr,temp_zero);
            
            MultiFab::Saxpy(SDC.sol[sdc_m],1.0,corr,0,0,1,0);
            if (Nprob<3){mlabec_BCfill.fourthOrderBCFill(SDC.sol[sdc_m],bdry_values);}
            SDC.sol[sdc_m].FillBoundary(geom.periodicity());
           
            /////////////////////////////////////////////////////////////////////////////
            // We plot the error after each smooth to see when perturbation arises.
            /////////////////////////////////////////////////////////////////////////////
        /*    int n = 10*numsweeps + sdc_m;
            if (1 > 0 )
             {
                 temp_err.setVal(1000.0);
             MultiFab::Copy(temp_err, SDC.sol[sdc_m], 0, 0, 1, 2);
                 
                 for ( MFIter mfi(temp_err); mfi.isValid(); ++mfi )
                 {
                     const Box& bx = mfi.validbox();
                     err_phi(BL_TO_FORTRAN_BOX(bx),
                             BL_TO_FORTRAN_ANYD(temp_err[mfi]),
                             geom.CellSize(), geom.ProbLo(), geom.ProbHi(),&a,&d,&r,&time, &epsilon,&k_freq, &kappa);
                 }
                 
             const std::string& pltfile = amrex::Concatenate("plt",n,5);
             WriteSingleLevelPlotfile(pltfile, temp_err, {"phi"}, geom, time, n);
             }*/
             
            
             
           
            
            
            
        }
        
        
  //////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
        
    }
  else
    {  // Do reaction solve  y - qij*y*(1-y)*(y-1/2) = rhs

      //  make a flag to change how the reaction is done
      int nflag=1;  // Lazy approximation

      qij = r*dt*SDC.Qimp[sdc_m-1][sdc_m];	            
      for ( MFIter mfi(SDC.sol[sdc_m]); mfi.isValid(); ++mfi )
	{
	  const Box& bx = mfi.validbox();
	  SDC_fcomp_reaction_F(BL_TO_FORTRAN_BOX(bx),
			       BL_TO_FORTRAN_BOX(domain_bx),
			       BL_TO_FORTRAN_ANYD(SDC.sol[sdc_m][mfi]),
			       BL_TO_FORTRAN_ANYD(rhs[mfi]),		      
			       BL_TO_FORTRAN_ANYD(SDC.f[2][sdc_m][mfi]),
			       &qij,&nflag); 
      }
      
    }

}




