# Paper: 2504.05584v2 — Transverse polarization and psionic form factors of Lambda hyperon at 3.773 GeV
# e+e- -> Lambda Lambda_bar -> p p_bar pi+ pi- at psi(3770) (3.773 GeV, 20.3 fb^-1)
# Standard selection with secondary vertex fits + 4C kinematic fit

### Dataset preparation ###
data_3773 = DatasetManager.real_data.find("712_3773")
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")

decay_card_signal = <<~DECAYCARD
    Decay psi(3770)
    1.000 Lambda Lambda_bar PHSP;
    Enddecay

    Decay Lambda
    1.000 p+ pi- HypWK;
    Enddecay

    Decay Lambda_bar
    1.000 anti-p- pi+ HypWK;
    Enddecay

    End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "sig_LambdaLambdabar"
  config.related_dataset = data_3773
  config.events = 500_000
  config.decay_card = decay_card_signal
  config.cross_section = :default
end

### Event selection ###
alg = Algorithm.new("LambdaLambdabar")
alg.set_header(["LambdaLambdabarAlg/LambdaLambdabar.h"])
    .set_constant({"ECMS" => [:double, 3.773]})

event_selection = Selection.new
event_selection
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"    # p+ + pi+
    nChrn ">=2"    # anti-p- + pi-
    nNet "==0"
  }
  # PID: identify protons and anti-protons
  .pid(method: :probability) {
    prob_cut 0.001
    identify :proton, against: [:kaon, :pion]
    nprp ">=1"
    nprm ">=1"
  }
  .remove([:prp <= :chrgp, :prm <= :chrgn])
  # Remaining tracks assigned as pions
  .assign({:chrgp => :pip, :chrgn => :pim})
  # Soft pion removal: p > 0.6 GeV/c (inexpressible in DSL)
  # Reconstruct Lambda -> p pi-
  .secondary_vertex_fit([:prp, :pim]) {
    build_virtual_particle(:Lambda).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # Reconstruct Lambda_bar -> anti-p pi+
  .secondary_vertex_fit([:prm, :pip]) {
    build_virtual_particle(:Lambda_bar).by_minimizing_mass_difference
    remove_used_particle_from_candidate_list
  }
  # 4C kinematic fit: Lambda + Lambda_bar
  .kinematic_fit([:Lambda, :Lambda_bar]) {
    nominal
    constrain_four_momentum
    chi2_cut 200       # loose cut; tight cut (chi2 < 100) applied in ROOT
  }

# Track p > 0.6 GeV/c cut applied in ROOT (or via for_each)
alg.note(:track_momentum_cut, "Charged track momentum p > 0.6 GeV/c applied in ROOT")

# Lambda mass window |M - m_Lambda| < 5 MeV/c^2 applied in ROOT
alg.note(:mass_windows, "Lambda mass window |M - m_Lambda| < 5 MeV/c^2 applied in ROOT")

alg.with_decay_card(decay_card_signal).apply(event_selection)
alg.execute_on([data_3773, incMC_3773, exMC_signal])