# =====================================================================
# J/psi(3097): search for the Delta-S = Delta-Q violating decay
#   Xi0 -> Sigma- e+ nu_e
# using a double-tag technique:
#   single tag (anti-Xi0) : anti-Xi0 -> anti-Lambda pi0,
#                           anti-Lambda -> anti-p pi+, pi0 -> gamma gamma
#   double tag (Xi0)      : Xi0 -> Sigma- e+ nu_e, Sigma- -> n pi-
#                           (neutron and neutrino undetected)
# =====================================================================

### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # J/psi(3097) real data
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # matching inclusive MC

# Decay card for the full J/psi -> Xi0 anti-Xi0 chain (EvtGen syntax)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 Xi0 anti-Xi0       PHSP;
    Enddecay

    Decay anti-Xi0
    1.0000 anti-Lambda0 pi0   PHSP;
    Enddecay

    Decay anti-Lambda0
    1.0000 anti-p- pi+        HypWK;
    Enddecay

    Decay pi0
    1.0000 gamma gamma        PHSP;
    Enddecay

    Decay Xi0
    1.0000 Sigma- e+ nu_e     PHSP;
    Enddecay

    Decay Sigma-
    1.0000 n0 pi-             HypWK;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC generated with the complete J/psi -> Xi0 anti-Xi0 chain
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_Xi0Xibar0_dtag"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "Xi0Dtag"
xi0_alg = Algorithm.new(alg_name)
xi0_alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
       .set_constant({"ECMS" => [:double, 3.097]})          # sqrt(s) = 3.097 GeV
       .set_alias({"std::vector<double>" => "Vdouble"})

xi0_selection = Selection.new
xi0_selection
  # ---------------- charged tracks ----------------
  .select_track {
    cos_theta 0.93            # |cos(theta)| < 0.93
    Vz        10.0
    Vr        1.0
  }
  # ---------------- photons ----------------
  .select_photon {
    tdc_emc_start     0       # EMC shower time window [0, 700] ns
    tdc_emc_end       14
    angle_to_track    10.0    # minimum 10 deg separation from any charged track
    energyThreshold_b 0.025   # 25 MeV in the barrel
    energyThreshold_e 0.050   # 50 MeV in the endcap
    nGam              ">=2"   # two photons from the tag-side pi0
  }
  # ---------------- tag side: anti-Lambda -> anti-p pi+ ----------------
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]     # anti-proton
    identify :pion,   against: [:kaon, :proton]   # pi+ (tag) / pi- (signal)
  }
  .secondary_vertex_fit([:prm, :pip]) {           # anti-Lambda -> anti-p pi+
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # ---------------- pi0 -> gamma gamma (Kalman mass-constrained fit) ----------------
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # ---------------- double-tag side: two remaining tracks of opposite charge ----------------
  .remove([:prm <= :chrgn, :pip <= :chrgp])       # drop the tag-side tracks
  .assign({:chrgp => :ep, :chrgn => :pim})        # remaining tracks: e+ (positive), pi- (negative)
  .remove(:pim) { condition "three_momentum_of(:pim) < 0.20 || three_momentum_of(:pim) > 0.38" }  # p(pi-) in (0.20, 0.38)
  .remove(:ep)  { condition "three_momentum_of(:ep) > 0.20" }                                     # p(e+) < 0.20
  # ---------------- final kinematic fit (neutron treated as missing) ----------------
  .kinematic_fit([:Lambda_bar, :pi0, :pim, :ep]) {
    nominal
    miss_track_of(:n_bar)                                       # neutron undetected
    invariant_mass_of(:Lambda_bar, :pi0).within(1.295, 1.335)   # |M(Lambda pi0) - m_Xi0| < 20 MeV/c^2
    constrain_four_momentum
    chi2_cut 200
  }

xi0_alg
  .note(:lambda_selection, "Lambda -> p pi- and anti-Lambda -> anti-p pi+ are both built by secondary vertex fits; only candidates with decay-length significance L/sigma_L > 2 and |M(p pi) - m_Lambda| < 5 MeV/c^2 are retained")
  .note(:pi0_selection, "pi0 -> gamma gamma: candidates with both photons reconstructed in the EMC endcap are rejected; M(gamma gamma) required in (115, 150) MeV/c^2")
  .note(:pid_correction_method, "double-tag-side e+ uses a custom low-momentum electron selector (probability method): L(e) > 0.001 and L(e)/(L(e)+L(pi)+L(K)) > 0.8, chi_dE/dx(e as pi) < -4.5, chi_dE/dx(pi as e) < -2.5, together with p(e+) < 0.20 GeV/c")
  .note(:tag_selection, "the Xi0/anti-Xi0 -> Lambda pi0 tag candidate is the combination with the smallest |M(Lambda pi0) - m_Xi0| inside the 20 MeV/c^2 window")
  .note(:signal_selection, "the double-tag side requires exactly two remaining charged tracks of opposite charge after the tag-side tracks are consumed")
  .note(:missing_particles, "both the neutron (Sigma- -> n pi-) and the neutrino (Xi0 -> Sigma- e+ nu_e) are undetected; the neutron is handled as a missing particle and the neutrino contributes to the missing four-momentum")
  .note(:mbc_fit, "the Delta-S = Delta-Q violating signal yield is extracted from a fit to the M_BC spectrum of the single-tag Xi0/anti-Xi0 -> Lambda pi0 candidates, not from a BOSS-level cut")
  .with_decay_card(decay_card_signal)
  .apply(xi0_selection)

# Execute on real data, inclusive MC and the exclusive signal MC
xi0_alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])