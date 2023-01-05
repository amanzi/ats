
#ifndef ATS_CLM_INTERFACE_HH_
#define ATS_CLM_INTERFACE_HH_

#include <cstdint>
#include <vector>

#include "Epetra_MultiVector.h"

#define NUM_LC_CLASSES 18

namespace ATS {
namespace ELMKernels {

  class ELMKernelDriver {

    std::shared_ptr<ELMStateType> ELMState_;

    // time-invariant, no spatial aspect
    std::shared_ptr<ELM::SnicarData<ViewD1, ViewD2, ViewD3>> snicar_data_;
    std::shared_ptr<ELM::SnwRdsTable<ViewD3>> snw_rds_table_;
    std::shared_ptr<ELM::PFTData<ViewD1>> pft_data_;

    // time variable, but no spatial aspect, unless domain very large
    // holds 12 monthly values for aerosol dep
    // interpolates input file using lat and lon
    std::shared_ptr<ELM::AerosolDataManager<ViewD1>> aerosol_data_;

    // time and space variable - must be time-interpolated at every step
    std::shared_ptr<ELM::AtmForcObjects> atm_forcing_;
    std::shared_ptr<ELM::PhenologyDataManager<ViewD2>> phen_data_;

    // derived from aerosol_data and snow physics
    std::shared_ptr<ELM::AerosolMasses<ViewD2>> aerosol_masses_;
    std::shared_ptr<ELM::AerosolConcentrations<ViewD2>> aerosol_concentrations_;

    // the path/to/file for phenology and atm forcing inputs
    std::string fname_surfdata_;
    std::string fname_forc_;

    // TODO better 
    // number of columns, cell per column
    int ncols_{0};
    int cell_per_col{0};

    // TODO should simplify this
    // don't need 2D domain decomp with ATS
    // will also likely have ATS do the file reads
    DomainDecomposition<2> dd_;

    // for phenology input data - should hide this in phen_data
    std::unordered_map<std::string, h_ViewD2> host_phen_views_;

    ELM::Utils::Date start_;


// mapping between ATS and ELM
// if ATS list of local nodes for each rank isn't contiguously ordered
// use a std::map<size_t, size_t>
// where map[0...ncols]  = ATS node id
// so ELM_var(i) = ATS_var(map[i])
// this should be fast to iterate?

// SerialDenseVec?
// soil temp solve?
// mesh definition - pass in mesh_, or just vector of interfaces (zisoi) to get dz, zsoi, zisoi
// slope aspect

// need:
// Land
// lat
// lon
// lat_r
// lon_r
// dewmx
// irrig_rate
// n_irrig_steps_left
// oldfflag
// veg_active
// do_capsnow
// topo_slope
// topo_std
// t_h2osfc
// tsoi
// altmax_indx
// altmax_lastyear_indx
// t10
// t_veg
// dz      |
// zsoi    | - assigned together
// zisoi   |
// h2osoi_ice |
// h2osoi_liq | - together?
// h2osoi_vol |

// need to make phenology_data and helpers inot it's own call
// make AtmForc read start date from the files, instead of providing
// should make a decision on clock sync
// put coszen into functions, also max_dayl and dayl




  public:

    int init_lsm(int ncols, const MPI_Comm& comm, std::string_view fsurfdata, 
                 std::string_view fsnicar, std::string_view fforc,
                 std::string_view fparam, std::string_view faerosol,
                 std::string_view fsnowage) {
      // assign filenames
      fname_surfdata_{fsurfdata};
      fname_forc_{fforc};

      ncols_ = ncols;
      cell_per_col = nlevgrnd; // hardwired in various places within ELMKernels
                               // may break if not set at 15

      // assert a bunch to make sure everything lines up

      // need to figure out mapping
      // implement better way of defining local/global problem
      // knowledge of global problem is currently required for 
      // pnetcdf usage - if ATS reads, we don't need it
      dd_ = ELM::Utils::create_domain_decomposition_2D(
          { ncols, ncols },
          { 1, 1 },
          { 0, 0 });
      dd_.comm = comm;
      
      // bloated state struct - will move some of this data
      ELMState_ = std::make_shared<ELMStateType>(ncols);

      // time invariant constants
      // without spatial awareness 
      snicar_data_ = std::make_shared<ELM::SnicarData<ViewD1, ViewD2, ViewD3>>();
      snw_rds_table_ = std::make_shared<ELM::SnwRdsTable<ViewD3>>();
      pft_data_ = std::make_shared<ELM::PFTData<ViewD1>>();

      // one year of aerosol deposition from nearest gridcell
      // currently invariant after this read
      aerosol_data_ = std::make_shared<ELM::AerosolDataManager<ViewD1>>();

      ELM::initialize_kokkos_elm(*ELMState_, *snicar_data_, *snw_rds_table_, *pft_data_,
                                 *aerosol_data_, dd_, fname_surfdata_, fname_param_,
                                 fname_snicar_, fname_snowage_, fname_aerosol_);
      return 0;
    }


    // initialization of state fields
    // must deep copy all of these into kokkos

    // can also be provided on cell-by-cell basis,
    // which will be necessary for restarts
    int init_snow(double snow_depth, double snow_temp) {
      return 0;
    }
    int init_snow(const Epetra_MultiVector& snow_depth) {
      return 0;
    } // snow temp needs to be taken care of here as well

    // should probably be run after snl calculated
    int init_temperature(double const_temp) {
      // deep copy to State view
      assign(ELMState_->t_h2osfc, const_temp);
      // copy constant value into subsurface only
      // snow layers should be initialized as 0.0
      auto h_tsoi = Kokkos::create_mirror_view(ELMState_->t_soisno);
      for (int i = 0; i < ncols_; ++i) {
        for (int j = nlevgrnd; j < nlevsno + nlevgrnd; ++j) {
          h_tsoi(i, j) = const_temp;
        }
      }
      assign(ELMState_->t_soisno, h_tsoi);
      return 0;
    }

    // should make a decision on clock sync
    int init_start_time(int year, int month, int day) {
      start_ = ELM::Utils::Date(year, month, day);
      return 0;
    }

    // just pass in mesh_? or zisoi?
    int init_mesh(const Epetra_MultiVector& dz) {
      return 0;
    } // should be constant with 1 column worth of dz, or zisoi






    int set_energy_state(const Epetra_MultiVector& soil_temp) {

      // place soil_temp from ATS into t_soisno once mapping decided upon
      return 0;
    }

    int set_water_state(const Epetra_MultiVector& saturation, const Epetra_MultiVector& porosity) {
      // place saturation from ATS into watsat once mapping decided upon
      // h2osoi_liq, h2osoi_ice?, watsat
      return 0;
    }


    int transfer_function_VG_to_CH() {
      // we're ignoring pressure, so this probably isn't needed
      return 0;
    }


    int setup_lsm() {

      int atm_nsteps = 101; // how many timesteps to read?
      const auto fstart = ELM::Utils::Date(1985, 1, 1); // make atm_forcing read this - may already be done

      // need to interpolate forcing at every step and 
      // reread from file periodically - happens automatically when forcing step (atm_nsteps-1)/atm_nsteps has been consumed 
      atm_forcing_ = std::make_shared<ELM::AtmForcObjects>(fname_forc_, fstart, atm_nsteps, ncols);

      // phenology data manager
      // make host mirrors - need to be persistent
      phen_data_ = std::make_shared<ELM::PhenologyDataManager<ViewD2>>(dd, ncells, 17);
      host_phen_views_ = get_phen_host_views(*phen_data);

      // containers for aerosol deposition and concentration within snowpack layers
      aerosol_masses_ = std::make_shared<ELM::AerosolMasses<ViewD2>>(ncells);
      aerosol_concentrations_ = std::make_shared<ELM::AerosolConcentrations<ViewD2>>(ncells);
    }



    int get_water_flux_diagnostic() {
      // get some/all of the optional qflx outputs
      // need to calc qflx_infiltration
      // with qflx_top_soil, qflx_snomelt, etc
      // from SoilHydrologyMod.F90
      // qflx_snow_grnd
      // qflx_rain_grnd
      // qflx_snow_melt
      // qflx_evap_tot
      // qflx_evap_soil
      // qflx_evap_veg
      // qflx_tran_veg
      // qflx_irr
      // irr_flag
      return 0;
    }

    int get_water_flux_total() {
      // get the accumulated qflx outputs
      // convert subsurf flux mm/s to 1/s -- mm to m 1e3
      // keep sfc flux in L/T
      return 0;
    }

    int get_energy_flux_diagnostic() {
      // get some/all of the optional eflx outputs
      //eflx_sh_grnd
      //eflx_sh_snow
      //eflx_sh_soil
      //eflx_sh_h2osfc
      return 0;
    }

    int get_energy_flux_total() {
      // get the accumulated eflx outputs
      // eflx_lh_tot
      // eflx_sh_tot
      // eflx_lwrad_out - or eflx_lwrad_net??
      // eflx_soil
      return 0;
    }


    int advance_time(double time, double dt) {
      // call 2 functions - prepare init_timestep() and run_lsm_physics
      return 0;
    }


  };


} // namespace ATS
} // namespace ELMKernels



#endif
