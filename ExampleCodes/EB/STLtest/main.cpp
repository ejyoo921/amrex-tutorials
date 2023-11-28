#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EB_utils.H>
#include <AMReX_EB_STL_utils.H>

// EY:
#include <AMReX_Print.H>
#include <chrono> // Timing tool 


using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {
        int nghost = 1;
        int max_grid_size=64;

        std::string stl_fname;

        Vector<Real> plo;
        Vector<Real> phi;
        Vector<int> ncells;
        Real dx[3];

        ParmParse pp;
        pp.getarr("prob_lo",plo);
        pp.getarr("prob_hi",phi);
        pp.getarr("ncells",ncells);
        pp.get("stl_file",stl_fname);
        pp.query("max_grid_size",max_grid_size);

        RealBox real_box({AMREX_D_DECL(plo[0], plo[1], plo[2])},
                {AMREX_D_DECL(phi[0], phi[1], phi[2])});

        Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};

        IntVect domain_lo(AMREX_D_DECL(0,0,0));
        IntVect domain_hi(AMREX_D_DECL(ncells[0]-1,ncells[1]-1,ncells[2]-1));

        dx[0]=(phi[0]-plo[0])/ncells[0];
        dx[1]=(phi[1]-plo[1])/ncells[1];
        dx[2]=(phi[2]-plo[2])/ncells[2];

        Box domain(domain_lo, domain_hi);
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);

        Geometry geom(domain,real_box,CoordSys::cartesian,is_periodic);
        DistributionMapping dm(ba);

        std::string pltfile;

        /* //---- Marker fill
        MultiFab marker;
        BoxArray nodal_ba = amrex::convert(ba, IntVect::TheNodeVector());
        marker.define(nodal_ba, dm, 1, nghost);

        STLtools stlobj; 

        Real scale = 1.0;
        Array<Real,3> const& center = {0,0,0};
        int reverse_normal = 0;

        stlobj.read_stl_file(stl_fname, scale, center, reverse_normal);
        stlobj.fill(marker,{0,0,0}, geom, -1.0, 1.0); // Default outside = 1, inside = -1, bounday = 0.

        // // write plot file
        WriteSingleLevelPlotfile("plt", marker, {"marker"}, geom, 0.0, 0);
        */

        

        /*---------// From Weiqun--------------------------------------------------------- 
        int required_coarsening_level = 0; // typically the same as the max AMR level index
        int max_coarsening_level = 100;    // typically a huge number so MG coarsens as much as possible
        // build a simple geometry using the "eb2." parameters in the inputs file
        EB2::Build(geom, required_coarsening_level, max_coarsening_level);

        std::string plot_file{"plt"};
        auto const& factory = makeEBFabFactory(geom, ba, dm, {1,1,1}, EBSupport::full);
        MultiFab const& vfrc = factory->getVolFrac();
        amrex::WriteMLMF(plot_file, {&vfrc}, {geom});
        ---------------------------------------------------------------------------------- */

        int required_coarsening_level = 0; // typically the same as the max AMR level index
        int max_coarsening_level = 100;    // typically a huge number so MG coarsens as much as possible
        // build a simple geometry using the "eb2." parameters in the inputs file

        // EY: Timing for EB2::Build
        auto t0 = std::chrono::high_resolution_clock::now();    
        EB2::Build(geom, required_coarsening_level, max_coarsening_level);

        auto t1 = std::chrono::high_resolution_clock::now();
        auto dt = 1.e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
        amrex::Print() << "EB2::Build time = " << dt << "(s)" << "\n";

        
        auto const& factory = makeEBFabFactory(geom, ba, dm, {1,1,1}, EBSupport::full);
        MultiFab const& vfrc = factory->getVolFrac();
        amrex::WriteMLMF("plt", {&vfrc}, {geom});
    }

    amrex::Print() << "Exit now" << "\n";
    amrex::Finalize();
}
