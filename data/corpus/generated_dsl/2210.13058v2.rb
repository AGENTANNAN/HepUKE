# ============================================================================
#  BOSS event-selection DSL
#  e+e- -> gamma phi J/psi  (via phi chi_c1,c2 ; exclusive MC through
#  psi(4260) -> phi gamma J/psi)  at 13 energy points 4.600 - 4.951 GeV
#  phi -> K+K-  or  K_S0 K_L0 ;  J/psi -> l+l-  (l = e, mu)
# ============================================================================

### ------------------------------- datasets ------------------------------- ###

# sample-name convention: [BOSS version]_[CMS energy in MeV]
energy_names = %w[
  703_4600  706_4610  706_4620  706_4640  706_4660  706_4680  706_4700
  707_4740  707_4750  707_4780  707_4840  707_4914  707_4946
]

real_data    = energy_names.map { |n| DatasetManager.real_data.find(n) }      # all real data at the 13 points
inclusive_mc = energy_names.map { |n| DatasetManager.inclusive_mc.find(n) }  # inclusive MC at the 13 points

### ------------------------------ decay cards ------------------------------ ###

# psi(4260) -> phi gamma J/psi ; phi -> K+ K- ; J/psi -> e+ e-
decay_card_kk = <<~DECAYCARD
    Decay psi(4260)
    1.000 phi gamma J/psi  PHSP;
    Enddecay

    Decay phi
    1.000 K+ K-  VSS;
    Enddecay

    Decay J/psi
    1.000 e+ e-  PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# psi(4260) -> phi gamma J/psi ; phi -> K_S0 K_L0 ; J/psi -> e+ e-
decay_card_kskl = <<~DECAYCARD
    Decay psi(4260)
    1.000 phi gamma J/psi  PHSP;
    Enddecay

    Decay phi
    1.000 K_S0 K_L0  PHSP;
    Enddecay

    Decay K_S0
    1.000 pi+ pi-  PHSP;
    Enddecay

    Decay J/psi
    1.000 e+ e-  PHOTOS VLL;
    Enddecay

    End
DECAYCARD

### ------------------------------ exclusive MC ------------------------------ ###

# 100k-event signal MC for each phi decay mode, one sample per energy point
# (signal generated over the whole scan via psi(4260) -> phi gamma J/psi)
exmc_kk = DatasetManager.create_exclusive_mc_for(real_data) do |config|
  config.sample_name   = "exmc_phiKK_jpsi_ee"
  config.events        = 100_000
  config.decay_card    = decay_card_kk
  config.cross_section = :default
end

exmc_kskl = DatasetManager.create_exclusive_mc_for(real_data) do |config|
  config.sample_name   = "exmc_phiKsKl_jpsi_ee"
  config.events        = 100_000
  config.decay_card    = decay_card_kskl
  config.cross_section = :default
end

### -------- Mode I : phi -> K+K-, exactly 3 tracks (K- constrained missing) -------- ###

alg_name1 = "GammaPhiJpsiKK3Trk"
alg1 = Algorithm.new(alg_name1)
alg1.set_header(["#{alg_name1}Alg/#{alg_name1}.h"])
    .set_constant({"ECMS" => [:double, 4.6]})
    .set_alias({"std::vector<double>" => "Vdouble"})

sel1 = Selection.new
  .select_track {
    cos_theta 0.93       # |cos(theta)| < 0.93
    Vz        10.0       # |Vz| < 10 cm
    Vr        1.0        # Vr < 1 cm
    nTot      "==3"      # exactly 3 charged tracks : K+ l+ l-
  }
  .select_photon {
    tdc_emc_start     0      # EMC TDC window 0 - 700 ns
    tdc_emc_end       14
    angle_to_track    10.0   # > 10 deg to nearest charged track
    energyThreshold_b 0.025  # barrel : 25 MeV
    energyThreshold_e 0.050  # endcap : 50 MeV
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.95,
                                   treat_as_electron_if_energy_above: 1.0
    identify :kaon, against: [:pion]   # K / pi separation
    nkp ">=1"                          # at least one K+
    nlp "==1"                          # l+
    nlm "==1"                          # l-
  }
  .kinematic_fit([:kp, :lp, :lm]) {
    nominal
    miss_track_of :km                                # missing K-, mass constrained to nominal (1C)
    invariant_mass_of(:lp, :lm).within(3.06, 3.14)   # |M(l+l-) - m_J/psi| < 3 sigma
    chi2_cut 200
  }

alg1.note(:muon_hit_depth, "3-track (K+ l+ l-) mode: in addition to the probability PID, a muon MUC hit depth > 40 cm is required for muon identification; the MUC hit-depth variable is not expressible in the DSL and the cut is applied in the ROOT analysis")

alg1.with_decay_card(decay_card_kk).apply(sel1)

### -------- Mode II : phi -> K+K-, exactly 4 tracks (K+ K- l+ l-), 4C -------- ###

alg_name2 = "GammaPhiJpsiKK4Trk"
alg2 = Algorithm.new(alg_name2)
alg2.set_header(["#{alg_name2}Alg/#{alg_name2}.h"])
    .set_constant({"ECMS" => [:double, 4.6]})
    .set_alias({"std::vector<double>" => "Vdouble"})

sel2 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nTot      "==4"      # exactly 4 charged tracks : K+ K- l+ l-
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.95,
                                   treat_as_electron_if_energy_above: 1.0
    identify :kaon, against: [:pion]
    nkp "==1"
    nkm "==1"
    nlp "==1"
    nlm "==1"
  }
  .kinematic_fit([:kp, :km, :lp, :lm]) {
    nominal
    constrain_four_momentum                          # 4C fit
    invariant_mass_of(:lp, :lm).within(3.06, 3.14)   # |M(l+l-) - m_J/psi| < 3 sigma
    chi2_cut 200
  }

alg2.with_decay_card(decay_card_kk).apply(sel2)

### ----- Mode III : phi -> K_S0 K_L0, K_S0 -> pi+pi- (>= 4 tracks) ----- ###

alg_name3 = "GammaPhiJpsiKsKl"
alg3 = Algorithm.new(alg_name3)
alg3.set_header(["#{alg_name3}Alg/#{alg_name3}.h"])
    .set_constant({"ECMS" => [:double, 4.6]})
    .set_alias({"std::vector<double>" => "Vdouble"})

sel3 = Selection.new
  .select_track {
    cos_theta 0.93
    Vz        10.0
    Vr        1.0
    nTot      ">=4"      # at least 4 tracks : pi+ pi- l+ l-
  }
  .select_photon {
    tdc_emc_start     0
    tdc_emc_end       14
    angle_to_track    10.0
    energyThreshold_b 0.025
    energyThreshold_e 0.050
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 0.95,
                                   treat_as_electron_if_energy_above: 1.0
    identify :pion, against: [:kaon]   # K_S0 daughters (K / pi separation)
    nlp "==1"
    nlm "==1"
  }
  .secondary_vertex_fit([:pip, :pim]) {                  # K_S0 -> pi+ pi-
    build_virtual_particle(:K_S0).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  .kinematic_fit([:K_S0, :lp, :lm]) {
    nominal
    miss_track_of :K_L0                              # missing K_L0, mass constrained to nominal (1C)
    invariant_mass_of(:lp, :lm).within(3.06, 3.14)   # |M(l+l-) - m_J/psi| < 3 sigma
    chi2_cut 200
  }

alg3.note(:ks_mass_window, "K_S0 candidates are required in a |M(pi+pi-) - m_K_S0| < 3 sigma window; the half-width is set by the K_S0 signal-MC mass resolution and the window is applied in the ROOT analysis")
alg3.note(:ks_decay_length, "K_S0 candidates are required to have a decay length significance > 2 sigma from the secondary-vertex fit; applied in the ROOT analysis")

alg3.with_decay_card(decay_card_kskl).apply(sel3)

### ----------------------------------- run ----------------------------------- ###

alg1.execute_on(real_data + inclusive_mc + exmc_kk)
alg2.execute_on(real_data + inclusive_mc + exmc_kk)
alg3.execute_on(real_data + inclusive_mc + exmc_kskl)