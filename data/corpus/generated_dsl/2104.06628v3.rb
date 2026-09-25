# J/psi -> D- e+ nu_e , D- -> K+ pi- pi-   @ sqrt(s) = 3.097 GeV
### Dataset preparation ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # 3.097 GeV real data (J/psi)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")  # matching inclusive MC

# Signal decay card: semileptonic J/psi -> D- e+ nu_e treated with VLL,
# D- -> K+ pi- pi- generated with phase space.
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 D- e+ nu_e    VLL;
    Enddecay

    Decay D-
    1.0000 K+ pi- pi-    PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive signal MC (500k events)
exMC_signal = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_jpsi_Denu_Kpipi"
    config.related_dataset = jpsi_data
    config.events          = 500000
    config.decay_card      = decay_card_signal
    config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiDenue"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # COMPASS: c.m. energy 3.097 GeV

event_selection = Selection.new
event_selection.select_track {                 # exactly 4 charged tracks, net charge 0
                  cos_theta 0.93               # |cos(theta)| < 0.93
                  Vz        10.0               # |Vz| < 10 cm
                  Vr        1.0                # Vr < 1 cm
                  nChrp     "==2"              # 2 positive tracks  (K+, e+)
                  nChrn     "==2"              # 2 negative tracks  (pi-, pi-)
                  nNet      "==0"              # net charge zero
                }
               .select_photon {                # good photon candidates
                  tdc_emc_start      0         # EMC time window 0-700 ns
                  tdc_emc_end        14
                  energyThreshold_b  0.025     # > 25 MeV in the barrel
                  energyThreshold_e  0.050     # > 50 MeV in the endcap
                  angle_to_track     10.0      # >= 10 deg from any charged track (K+, pi-, e+)
                }
               .pid(method: :probability) {    # probability-method PID
                  prob_cut 0.001               # probability > 0.001
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6   # e+ (high-momentum lepton)
                  identify :kaon, against: [:pion, :proton]   # K+ / K-  (likelihood separation vs pi, p)
                  identify :pion, against: [:kaon, :proton]   # pi+ / pi-
                  nkp   "==1"                  # exactly one K+
                  npim  "==2"                  # exactly two pi-
                  nlp   "==1"                  # exactly one e+
                }
                # 1C kinematic fit: constrain the K+ pi- pi- invariant mass to the nominal D- mass.
                # No 4C fit is used because the neutrino escapes detection.
               .kalman_kinematic_fit([:kp, :pim, :pim]) {
                  invariant_mass_of(:kp, :pim, :pim).constrain_to_nominal_mass_of(:Dm)  # 1C mass constraint
                  invariant_mass_of(:kp, :pim, :pim).within(1.85, 1.89)                 # D- mass window (GeV/c^2)
                  chi2_cut 10                                                           # chi2 < 10
               }
                # Partial reconstruction of D- e+ : the untagged nu_e (recID 3) is inferred from recoil.
                # This replaces the (absent) 4C kinematic fit for the missing-neutrino topology.
               .partial_miss([3])

my_algorithm
  .note(:pid_correction_method, "positron identified with P_e > 0.001 and P_e/(P_pi + P_K) > 1/4 by the probability method; the E/p window 0.85-1.05 is applied at ROOT level")
  .note(:background_veto, "total energy of good photons required < 0.2 GeV to suppress extra-photon backgrounds")
  .note(:missing_momentum_cut, "partial reconstruction of D- e+ plus missing nu_e; requirement |p_miss| > 50 MeV/c and the signal-yield extraction from U_miss = E_miss - c|p_miss| are applied at the ROOT level")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Run over real data, inclusive MC and signal exclusive MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])