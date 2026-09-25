# =====================================================================
# Born cross section of e+e- -> K_S0 K±pi∓ (K_S0 -> pi+ pi-)
# at the 15 BESIII 703 scan points (sqrt(s) ~ 3.808 - 4.600 GeV),
# using real data only + exclusive ConExc signal MC.
# =====================================================================

### Dataset description ###
# 15 BESIII 703 scan points (XYZ region, continuum)
data_3810 = DatasetManager.real_data.find("703_3810")
data_3900 = DatasetManager.real_data.find("703_3900")
data_4090 = DatasetManager.real_data.find("703_4090")
data_4190 = DatasetManager.real_data.find("703_4190")
data_4210 = DatasetManager.real_data.find("703_4210")
data_4220 = DatasetManager.real_data.find("703_4220")
data_4230 = DatasetManager.real_data.find("703_4230")
data_4245 = DatasetManager.real_data.find("703_4245")
data_4310 = DatasetManager.real_data.find("703_4310")
data_4390 = DatasetManager.real_data.find("703_4390")
data_4420 = DatasetManager.real_data.find("703_4420")
data_4470 = DatasetManager.real_data.find("703_4470")
data_4530 = DatasetManager.real_data.find("703_4530")
data_4575 = DatasetManager.real_data.find("703_4575")
data_4600 = DatasetManager.real_data.find("703_4600")

data_points = [data_3810, data_3900, data_4090, data_4190, data_4210,
               data_4220, data_4230, data_4245, data_4310, data_4390,
               data_4420, data_4470, data_4530, data_4575, data_4600]

# ConExc decay card: models the ISR / vacuum-polarisation corrections of the
# measured Born cross section (continuum -> vpho). ConExc equally populates
# K_S0 K+ pi- and its charge conjugate K_S0 K- pi+. No "Particle vpho" line is
# written: for a multi-energy scan the DSL injects Particle vpho <ECMS> 0.0 at
# each energy point.
decay_card_kpi = <<~DECAYCARD
    Decay vpho
    0.500 K_S0 K+ pi- ConExc;
    0.500 K_S0 K- pi+ ConExc;
    Enddecay

    Decay K_S0
    1.000 pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC: 100k events at each of the 15 energy points.
exMC_kpi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "exmc_ks0kpi_scan"   # auto-suffixed per energy point
  config.events        = 100000
  config.decay_card    = decay_card_kpi
  config.cross_section = :default
end

### Event selection (BOSS) ###
# Mode I: e+e- -> K_S0 K+ pi-
alg_name_kp = "Ks0KpiKplus"
alg_kp = Algorithm.new(alg_name_kp)
alg_kp.set_header(["#{alg_name_kp}Alg/#{alg_name_kp}.h"])
      .set_constant({"ECMS" => [:double, 4.260]})  # nominal; ConExc injects the per-point ECMS

# Mode II: e+e- -> K_S0 K- pi+ (charge conjugate)
alg_name_km = "Ks0KpiKminus"
alg_km = Algorithm.new(alg_name_km)
alg_km.set_header(["#{alg_name_km}Alg/#{alg_name_km}.h"])
      .set_constant({"ECMS" => [:double, 4.260]})

# Common selection: charged tracks, photons, PID, K_S0 pion E/p suppression
event_selection_common = Selection.new
    .select_track {
        cos_theta 0.93       # |cos(theta)| < 0.93
        Vz        10.0       # |Vz| < 10 cm
        Vr        1.0        # Vr < 1 cm
        nChrp    "==2"       # exactly 2 positive tracks
        nChrn    "==2"       # exactly 2 negative tracks
        nNet     "==0"       # net charge zero
    }
    .select_photon {
        tdc_emc_start     0      # TDC window start
        tdc_emc_end       14     # TDC window end
        angle_to_track    10.0   # min angle to nearest charged track (deg)
        energyThreshold_b 0.025  # 25 MeV (barrel)
        energyThreshold_e 0.050  # 50 MeV (endcap)
        # no explicit photon multiplicity requirement
    }
    .pid(method: :probability) {          # dE/dx - TOF probability method
        prob_cut 0.001
        identify :kaon, against: [:pion, :proton]
        identify :pion, against: [:kaon, :proton]
    }
    .remove(:pip) { condition "ep_ratio_of(:pip) > 0.8" }  # suppress gamma conversion
    .remove(:pim) { condition "ep_ratio_of(:pim) > 0.8" }

# Mode I selection: reconstruct K_S0 -> pi+ pi- and 4C fit to K_S0 K+ pi-
kp_selection = event_selection_common.dup
    .secondary_vertex_fit([:pip, :pim]) {
        build_virtual_particle(:K_S0).by_minimizing_mass_difference  # closest to K_S0 mass
        remove_used_particle_from_candidate_list
    }
    .kinematic_fit([:K_S0, :kp, :pim]) {
        nominal                    # result of this fit is saved
        constrain_four_momentum    # 4C energy-momentum constraint
        chi2_cut 200               # loose cut; tighter cut optimised in ROOT
    }

# Mode II selection: reconstruct K_S0 -> pi+ pi- and 4C fit to K_S0 K- pi+
km_selection = event_selection_common.dup
    .secondary_vertex_fit([:pip, :pim]) {
        build_virtual_particle(:K_S0).by_minimizing_mass_difference
        remove_used_particle_from_candidate_list
    }
    .kinematic_fit([:K_S0, :km, :pip]) {
        nominal
        constrain_four_momentum
        chi2_cut 200
    }

# Extra K_S0 selection criteria not expressible in the DSL -> notes
alg_kp
  .note(:ks0_mass_window, "K_S0 candidate required |M(pi+pi-) - M(K_S0)| < 0.020 GeV/c^2")
  .note(:decay_length,    "K_S0 decay length significance > 2 sigma")
  .note(:dca_cut,         "K_S0 daughter pion DCA < 25 cm in z and < 20 cm in r-phi")
  .note(:background_veto, "K_S0 sidebands 0.435-0.455 and 0.545-0.565 GeV/c^2 used for background estimation")
  .with_decay_card(decay_card_kpi).apply(kp_selection)

alg_km
  .note(:ks0_mass_window, "K_S0 candidate required |M(pi+pi-) - M(K_S0)| < 0.020 GeV/c^2")
  .note(:decay_length,    "K_S0 decay length significance > 2 sigma")
  .note(:dca_cut,         "K_S0 daughter pion DCA < 25 cm in z and < 20 cm in r-phi")
  .note(:background_veto, "K_S0 sidebands 0.435-0.455 and 0.545-0.565 GeV/c^2 used for background estimation")
  .with_decay_card(decay_card_kpi).apply(km_selection)

# Execute on the 15 real-data points and the 15 corresponding signal-MC samples
root_files_kp = alg_kp.execute_on(data_points + exMC_kpi)
root_files_km = alg_km.execute_on(data_points + exMC_kpi)