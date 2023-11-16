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

        MultiFab marker;
        BoxArray nodal_ba = amrex::convert(ba, IntVect::TheNodeVector());
        marker.define(nodal_ba, dm, 1, nghost);

        STLtools stlobj; 

        Real scale = 1.0;
        Array<Real,3> const& center = {0,0,0};
        int reverse_normal = 0;

        stlobj.read_stl_file(stl_fname, scale, center, reverse_normal);
        stlobj.fill(marker,{0,0,0},geom, -1.0, 1.0); // Default outside = 1, inside = -1, bounday = 0.

        // // write plot file
        std::string pltfile;
        WriteSingleLevelPlotfile("plt", marker, {"marker"}, geom, 0.0, 0);

        /* old one --------------------------------------------------------
        // EY: what are these apx, apy, apz doing here?
        MultiFab marker,apx,apy,apz;
        // pp.getarr("outside_point",pointoutside);
        BoxArray nodal_ba = amrex::convert(ba, IntVect::TheNodeVector());
        marker.define(nodal_ba, dm, 1, nghost);
        // EY: read_ascii_stl_file is not public
        stlobj.read_ascii_stl_file(stl_fname); 
        // EY: error says we need 4 arguments

        int reverse_normal = 0;
        stlobj.read_stl_file(stl_fname, scale, center, reverse_normal);

        Real plo_arr[]={plo[0],plo[1],plo[2]}; // EY what do we do with this?
        Real po_arr[]={pointoutside[0],pointoutside[1],pointoutside[2]};

        EY: wrong name
        stlobj.stl_to_markerfab(marker,geom,po_arr);
        EY: nghost needs to be IntVect form
        stlobj.fill(marker,{0,0,0},geom);
        EY: what is the difference?
        int box_type_num;
        box_type_num = stlobj.getBoxType (domain, geom, RunOn::Gpu);
        amrex::PrintToFile("box_type") << "domain  \n" << box_type_num;

        marker.FillBoundary(geom.periodicity());
        const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
        const EB2::Level& eb_level = eb_is.getLevel(geom);

        // options are basic, volume, or full
        EBSupport ebs = EBSupport::full;
        EBSupport ebs = EBSupport::volume; // EY

        // number of ghost cells for each of the 3 EBSupport types
        Vector<int> ng_ebs = {2,2,2};

        // This object provides access to the EB database in the format of basic AMReX objects
        // such as BaseFab, FArrayBox, FabArray, and MultiFab
        EBFArrayBoxFactory factory(eb_level, geom, ba, dm, ng_ebs, ebs);

        // write plot file
        std::string pltfile;
        pltfile = "plt";
        WriteSingleLevelPlotfile(pltfile, marker, {"marker"}, geom, 0.0, 0);

        // From Weiqun 
        int required_coarsening_level = 0; // typically the same as the max AMR level index
        int max_coarsening_level = 100;    // typically a huge number so MG coarsens as much as possible
        // build a simple geometry using the "eb2." parameters in the inputs file
        EB2::Build(geom, required_coarsening_level, max_coarsening_level);

        std::string plot_file{"plt"};
        auto const& factory = makeEBFabFactory(geom, ba, dm, {1,1,1}, EBSupport::full);
        MultiFab const& vfrc = factory->getVolFrac();
        amrex::WriteMLMF(plot_file, {&vfrc}, {geom});
        ---------------------------------------------------------------------------------- */
    }
    amrex::Print() << "Exit now" << "\n";
    amrex::Finalize();
}
