# Core DSL classes and dependencies will be loaded automatically at execution
### Dataset description ###
# Real data and inclusive MC at the energy points that provide both samples
data_4127 = DatasetManager.real_data.find("705_4130")   # ~4.127 GeV
data_4226 = DatasetManager.real_data.find("703_4230")   # 4.226 GeV
data_4436 = DatasetManager.real_data.find("705_4440")   # ~4.436 GeV
data_4600 = DatasetManager.real_data.find("703_4600")   # 4.600 GeV

incMC_4127 = DatasetManager.inclusive_mc.find("705_4130")   # inclusive MC at ~4.127 GeV
incMC_4226 = DatasetManager.inclusive_mc.find("703_4230")   # inclusive MC at 4.226 GeV
incMC_4436 = DatasetManager.inclusive_mc.find("705_4440")   # inclusive MC at ~4.436 GeV
incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")   # inclusive MC at 4.600 GeV

# Decay card for e+e- -> K+ K- J/psi, J/psi -> e+e-
# The K+ K- system carries the weighted PHSP, f0(980) and f2(1270) intermediate components
decay_card_ee = <<~DECAYCARD
    Decay psi(4260)
    0.340  K+  K-  J/psi               PHSP;
    0.330  f_0(980)  J/psi             PHSP;
    0.330  f_2(1270)  J/psi            PHSP;
    Enddecay

    Decay f_0(980)
    1.000  K+  K-                      PHSP;
    Enddecay

    Decay f_2(1270)
    1.000  K+  K-                      PHSP;
    Enddecay

    Decay J/psi
    1.000  e+  e-                      PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Decay card for e+e- -> K+ K- J/psi, J/psi -> mu+mu-
decay_card_mumu = <<~DECAYCARD
    Decay psi(4260)
    0.340  K+  K-  J/psi               PHSP;
    0.330  f_0(980)  J/psi             PHSP;
    0.330  f_2(1270)  J/psi            PHSP;
    Enddecay

    Decay f_0(980)
    1.000  K+  K-                      PHSP;
    Enddecay

    Decay f_2(1270)
    1.000  K+  K-                      PHSP;
    Enddecay

    Decay J/psi
    1.000  mu+  mu-                    PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# 300k-event exclusive MC for each lepton mode at 4.226, 4.436 and 4.600 GeV
exMC_energy_points = [data_4226, data_4436, data_4600]

exMC_ee = DatasetManager.create_exclusive_mc_for(exMC_energy_points) do |config|
  config.sample_name   = "exmc_kkjpsi_ee"
  config.events        = 300000
  config.decay_card    = decay_card_ee
  config.cross_section = :default
end

exMC_mumu = DatasetManager.create_exclusive_mc_for(exMC_energy_points) do |config|
  config.sample_name   = "exmc_kkjpsi_mumu"
  config.events        = 300000
  config.decay_card    = decay_card_mumu
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name  = "KKJpsi"
algorithm = Algorithm.new(alg_name)
algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
         .set_constant({"ECMS" => [:double, 4.226]})        # reference CMS energy (per-run beam energy used at run time)
         .set_alias({"std::vector<double>" => "Vdouble"})

# Common selection for both lepton modes (J/psi -> e+e- and J/psi -> mu+mu-)
event_selection = Selection.new
event_selection
    .select_track {                     # Charged track selection
      cos_theta  0.93                   # |cos(theta)| < 0.93
      Vz         10.0                   # |Vz| < 10 cm
      Vr         1.0                    # Vr < 1 cm
      nTot       ">=3"                  # at least three charged tracks
    }
    .pid(method: :probability) {        # Particle identification (probability method)
      prob_cut 0.001                    # C.L. > 0.001
      # High-momentum tracks (p > 1.0 GeV) are treated as leptons (e/mu separated
      # by EMC energy); applied together with the kaon identification below
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6
      identify :kaon, against: [:pion]  # K+ selected against pions with C.L.(K) > 0.001
      nkp ">=1"                         # at least one K+
      nlp "==1"                         # one charged lepton l+
      nlm "==1"                         # one charged lepton l-
    }
    # Nominal kinematic fit on the l+ l- K+ K- system. The K- is not detected and is
    # inferred through the missing-track 1C constraint to the kaon mass, then the
    # four-momentum of the system is constrained to the CMS. Candidate with the
    # smallest chi2 is chosen automatically.
    .kinematic_fit([:lp, :lm, :kp, :km]) {
      nominal
      miss_track_of(:km)
      constrain_four_momentum
      chi2_cut 200
    }
    # Partial-reconstruction hypothesis: combined vertex + 1C quality variable
    .kinematic_fit([:lp, :lm, :kp, :km]) {
      miss_track_of(:km)
      vertex_fit([0, 1])
      chi2_cut 20
    }

algorithm
    .note(:pid_correction_method, "high-momentum lepton PID: electrons required E/p > 0.8 in the EMC, and muons required MUC penetration depth > 40 cm for at least one muon; these detector-level conditions are not expressible in the DSL and are applied on top of identify_high_momentum_leptons")
    .note(:background_veto, "radiative Bhabha background removed by requiring the opening-angle cos(theta) < 0.98 for every opposite-charge track pair")
    .note(:efficiency_curve, "the selection is run at 28 centre-of-mass energies (4.127-4.600 GeV); the beam energy is read per run and the ECMS constant above is only a reference value")
    .with_decay_card(decay_card_ee).apply(event_selection)

# Execute on real data, inclusive MC and the exclusive MC samples
root_files = algorithm.execute_on(
  [data_4127, data_4226, data_4436, data_4600,
   incMC_4127, incMC_4226, incMC_4436, incMC_4600] + exMC_ee + exMC_mumu
)