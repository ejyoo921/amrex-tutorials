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

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {
        int nghost = 1;
        int max_grid_size=64;
        // EY: what are these apx, apy, apz doing here?
        // MultiFab marker,apx,apy,apz;
        MultiFab marker;
        std::string stl_fname;

        Vector<Real> plo;
        Vector<Real> phi;
        Vector<int> ncells;
        Vector<Real> pointoutside;
        Real dx[3];

        ParmParse pp;
        pp.getarr("prob_lo",plo);
        pp.getarr("prob_hi",phi);
        pp.getarr("ncells",ncells);
        pp.get("stl_file",stl_fname);
        // pp.getarr("outside_point",pointoutside);
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
        BoxArray nodal_ba = amrex::convert(ba, IntVect::TheNodeVector());

        marker.define(nodal_ba, dm, 1, nghost);


        STLtools stlobj; 

        // EY: read_ascii_stl_file is not public
        // stlobj.read_ascii_stl_file(stl_fname); 
        // EY: error says we need 4 arguments
        Real scale = 1.0;
        Array<Real,3> const& center = {0,0,0};
        int reverse_normal = 0;
        stlobj.read_stl_file(stl_fname, scale, center, reverse_normal);

        // Real plo_arr[]={plo[0],plo[1],plo[2]}; // EY what do we do with this?
        // Real po_arr[]={pointoutside[0],pointoutside[1],pointoutside[2]};

        // EY: wrong name
        // stlobj.stl_to_markerfab(marker,geom,po_arr);
        // EY: nghost needs to be IntVect form
        stlobj.fill(marker,{0,0,0},geom);

        marker.FillBoundary(geom.periodicity());

        //write plot file
        std::string pltfile;
        pltfile = "plt";
        WriteSingleLevelPlotfile(pltfile, marker, {"marker"}, geom, 0.0, 0);
    }

    amrex::Finalize();
}
