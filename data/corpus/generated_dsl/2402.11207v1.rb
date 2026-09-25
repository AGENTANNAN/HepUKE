# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# XYZ scan data, sqrt(s) = 4.130 - 4.700 GeV
# Sample name convention: [BOSS version]_[ECMS in MeV]
scan_data_names = %w[
  705_4130 705_4160 703_4180 703_4200 703_4220 703_4230 703_4260 705_4290
  705_4315 705_4340 703_4360 705_4380 705_4400 703_4420 705_4440 706_4620
  706_4640 706_4660 706_4680 706_4700
]
scan_data = scan_data_names.map { |name| DatasetManager.real_data.find(name) }

# Corresponding inclusive MC samples
# (no inclusive MC sample is published for the 4290 / 4315 / 4340 / 4380 / 4400 / 4440 points)
scan_incMC = %w[
  705_4130 705_4160 703_4180 703_4200 703_4220 703_4230 703_4260 703_4360
  703_4420 706_4620 706_4640 706_4660 706_4680 706_4700
].map { |name| DatasetManager.inclusive_mc.find(name) }

### Decay cards (EvtGen format); top mother psi(4260) for the KKMC convention ###
# Channel 1: e+e- -> p+ p+ pi- anti-d  (anti-deuteron)
decay_card_pppim_dbar = <<~DECAYCARD
    Decay psi(4260)
    1.0000 p+ p+ pi- anti-d    PHSP;
    Enddecay
    End
DECAYCARD

# Channel 2 (charge conjugate): e+e- -> anti-p- anti-p- pi+ d  (deuteron)
decay_card_ppbarpip_d = <<~DECAYCARD
    Decay psi(4260)
    1.0000 anti-p- anti-p- pi+ d    PHSP;
    Enddecay
    End
DECAYCARD

### Exclusive MC: 100k events for each channel at every energy point ###
exMC_pppim_dbar = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_pppim_dbar"   # suffixed per energy point
  config.events        = 100_000
  config.decay_card    = decay_card_pppim_dbar
  config.cross_section = :default
end

exMC_ppbarpip_d = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name   = "exmc_ppbarpip_d"   # suffixed per energy point
  config.events        = 100_000
  config.decay_card    = decay_card_ppbarpip_d
  config.cross_section = :default
end

### Event selection (BOSS) ###
# The two charge-conjugate channels require different final-state charges, so they
# are built as two algorithms that share the same selection chain structure.
alg_pppim_dbar = Algorithm.new("PpPimDbar")        # e+e- -> p+ p+ pi- anti-d
alg_pppim_dbar.set_header(["PpPimDbarAlg/PpPimDbar.h"])
              .set_constant({"ECMS" => [:double, 4.26]})   # reference energy of the scan
              .set_alias({"std::vector<double>" => "Vdouble"})

alg_ppbarpip_d = Algorithm.new("PpbarPipD")        # e+e- -> anti-p- anti-p- pi+ d
alg_ppbarpip_d.set_header(["PpbarPipDAlg/PpbarPipD.h"])
              .set_constant({"ECMS" => [:double, 4.26]})
              .set_alias({"std::vector<double>" => "Vdouble"})

# ---- Channel 1: e+e- -> p+ p+ pi- anti-d (anti-deuteron potentially missed) ----
selection_pppim_dbar = Selection.new
selection_pppim_dbar
  .select_track {                 # charged track quality cuts
      cos_theta 0.93              # |cos(theta)| < 0.93
      Vz        10.0              # |Vz| < 10 cm
      Vr        1.0               # Vr < 1 cm
      nChrp     "==2"             # two positive tracks
      nChrn     "==1"             # one negative track
      nNet      "==1"             # net charge +1
  }
  .pid(method: :probability) {    # PID with the probability method
      prob_cut 0.001              # PID probability > 0.001
      identify :proton, against: [:kaon, :pion]    # p+ / pbar separated from K and pi
      identify :pion,   against: [:kaon, :proton]  # pi+ / pi- separated from K and p
      nprp "==2"                  # two protons
      npim "==1"                  # one pi-
  }
  .secondary_vertex_fit([:prp, :prp, :pim]) {     # common vertex for the two protons and the pion
      build_virtual_particle(:pppim).by_minimizing_verfit_chi2   # pick the best vertex (chi2_VF)
  }
  .partial_miss([4]) {            # recID 4 = anti-d is left undetected (partial reconstruction)
      require_recoil_mass 1.80, 2.30   # recoil-mass window 1.80 - 2.30 GeV
  }

# ---- Channel 2 (charge conjugate): e+e- -> anti-p- anti-p- pi+ d ----
selection_ppbarpip_d = Selection.new
selection_ppbarpip_d
  .select_track {
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     "==1"             # one positive track
      nChrn     "==2"             # two negative tracks
      nNet      "==-1"            # net charge -1
  }
  .pid(method: :probability) {
      prob_cut 0.001
      identify :proton, against: [:kaon, :pion]
      identify :pion,   against: [:kaon, :proton]
      nprm "==2"                  # two anti-protons
      npip "==1"                  # one pi+
  }
  .secondary_vertex_fit([:prm, :prm, :pip]) {
      build_virtual_particle(:ppbarpip).by_minimizing_verfit_chi2
  }
  .partial_miss([4]) {            # recID 4 = d is left undetected
      require_recoil_mass 1.80, 2.30
  }

# BOSS-side procedures that have no DSL construct
alg_pppim_dbar
  .note(:vertex_fit_chi2, "the secondary vertex fit of the two protons and the pion is
    required to satisfy chi2_VF < 70")
  .note(:vertex_position, "the fitted secondary-vertex position must lie within 0.5 cm of
    (0.12, -0.14) cm in the transverse plane, measured with respect to the interaction point")
  .note(:recoil_pt_cut, "for the three-track (deuteron-missed) topology, pt(recoil) < 0.35 GeV/c
    is applied for events with |cos(theta_recoil)| <= 0.93")
  .note(:deuteron_tof_mass, "four-track case (deuteron detected): the deuteron candidate must have
    M^2 from TOF between 3.2 and 4.1 (GeV/c^2)^2 together with cos(theta') > 0.9")
  .note(:background_veto, "four-track case: proton/anti-proton vetoes are applied to suppress the
    dominant background before the final selection")

alg_ppbarpip_d
  .note(:vertex_fit_chi2, "the secondary vertex fit of the two anti-protons and the pion is
    required to satisfy chi2_VF < 70")
  .note(:vertex_position, "the fitted secondary-vertex position must lie within 0.5 cm of
    (0.12, -0.14) cm in the transverse plane, measured with respect to the interaction point")
  .note(:recoil_pt_cut, "for the three-track (deuteron-missed) topology, pt(recoil) < 0.35 GeV/c
    is applied for events with |cos(theta_recoil)| <= 0.93")
  .note(:deuteron_tof_mass, "four-track case (deuteron detected): the deuteron candidate must have
    M^2 from TOF between 3.2 and 4.1 (GeV/c^2)^2 together with cos(theta') > 0.9")
  .note(:background_veto, "four-track case: proton/anti-proton vetoes are applied to suppress the
    dominant background before the final selection")

# Generate the algorithms for the two charge-conjugate channels
alg_pppim_dbar.with_decay_card(decay_card_pppim_dbar).apply(selection_pppim_dbar)
alg_ppbarpip_d.with_decay_card(decay_card_ppbarpip_d).apply(selection_ppbarpip_d)

# Execute on the scan data, the corresponding inclusive MC, and the exclusive MC of each channel
root_files_pppim_dbar = alg_pppim_dbar.execute_on(scan_data + scan_incMC + exMC_pppim_dbar)
root_files_ppbarpip_d = alg_ppbarpip_d.execute_on(scan_data + scan_incMC + exMC_ppbarpip_d)