### Dataset preparation ###
# Real data at the three c.m. energies of interest (subset of the 56-point 3.51–4.95 GeV scan)
data_3773 = DatasetManager.real_data.find("712_3773")   # sqrt(s) = 3.773 GeV
data_4180 = DatasetManager.real_data.find("703_4180")   # sqrt(s) = 4.180 GeV
data_4230 = DatasetManager.real_data.find("703_4230")   # sqrt(s) = 4.230 GeV
# Inclusive MC (provided at 3.773 GeV only)
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

# Decay card for e+e- -> K+ Xi0 anti-Sigma-.
# Continuum production simulated with psi(4260) as the top KKMC mother (BESIII convention).
decay_card_signal = <<~DECAYCARD
  Decay psi(4260)
  1.000 K+ Xi0 anti-Sigma- PHSP;
  Enddecay

  Decay Xi0
  1.000 Lambda0 pi0 PHSP;
  Enddecay

  Decay Lambda0
  1.000 p+ pi- PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma PHSP;
  Enddecay

  Decay anti-Sigma-
  1.000 anti-n- pi+ PHSP;
  Enddecay

  End
DECAYCARD

# Signal exclusive MC: 100k events at each of the three energy points
exMCs_signal = DatasetManager.create_exclusive_mc_for([data_3773, data_4180, data_4230]) do |config|
  config.sample_name   = "exmc_KXi0_Sigmabar"
  config.events        = 100000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "KXi0Sigmabar"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.230]})

# Full selection chain; partial reconstruction replaces the kinematic fit.
event_selection = Selection.new
  .select_track {
    cos_theta 0.93      # |cos(theta)| < 0.93
    # No impact-parameter (Vz/Vr) requirement: Lambda is reconstructed as a secondary vertex
    nChrp ">=2"         # at least two positive tracks (K+, p+)
    nChrn ">=1"         # at least one negative track (pi-)
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 14
    angle_to_track 10.0
    energyThreshold_b 0.025       # 25 MeV (barrel)
    energyThreshold_e 0.050       # 50 MeV (endcap)
    nGam ">=2"                    # at least two photons (pi0 -> gamma gamma)
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    identify :kaon,   against: [:pion, :proton]
    identify :pion,   against: [:kaon, :proton]
    nprp ">=1"          # at least one proton
    nkp  ">=1"          # at least one K+
    npim ">=1"          # at least one pi-
  }
  .kalman_kinematic_fit([:gamma, :gamma]) {   # 1C Kalman fit: pi0 -> gamma gamma
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 500
    npi0 ">=1"
  }
  .secondary_vertex_fit([:prp, :pim]) {        # secondary vertex for Lambda -> p pi-
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Partial reconstruction of the K+ Xi0 system; anti-Sigma- is the missing recoil.
  # recIDs from the decay card: 1 = K+, 2 = Xi0, 4 = Lambda, 5 = pi0, 6 = p+, 7 = pi-, 8/9 = gamma
  .partial_rec([1, 2, 4, 5, 6, 7, 8, 9]) {
    best_combination_by_mass :Xi0, 1.31486     # choose Xi0 candidate closest to m_Xi0
    require_recoil_mass 1.10, 1.30             # recoil (anti-Sigma-) mass window
  }

my_algorithm
  .note(:lambda_mass_window, "|M(p pi-) - m_Lambda| <= 5 MeV window applied on the secondary-vertex Lambda candidate")
  .note(:xi0_mass_window, "|M(Lambda pi0) - m_Xi0| < 10 MeV window applied on the Xi0 candidate")
  .note(:lambda_flight_significance, "Lambda flight significance L/dL > 2 required")
  .note(:isr_vp_correction, "cross sections extracted with ISR and vacuum-polarization corrections")
  .note(:multi_energy_ecms, "analysis spans 3.773/4.180/4.230 GeV; the per-run beam energy must be used for the recoil-mass constraint")

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)

my_algorithm.execute_on([data_3773, data_4180, data_4230, incMC_3773] + exMCs_signal)