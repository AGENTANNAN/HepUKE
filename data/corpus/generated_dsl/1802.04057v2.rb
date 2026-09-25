# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")     # psi(3686) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # corresponding psi(3686) inclusive MC

# Decay card for the signal process (EvtGen format):
#   psi(3686) -> Lambda_c+ pbar e+ e-,  Lambda_c+ -> p K- pi+
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.0000 Lambda_c+ anti-p- e+ e-      PHSP;
    Enddecay

    Decay Lambda_c+
    1.0000 p+ K- pi+                    PHSP;
    Enddecay

    End
DECAYCARD

# 100k-event exclusive MC for the signal mode
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_psip_lambdac_pbar_ee"
  config.related_dataset = psip_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "PsipToLambdaCpbarEE"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})

event_selection = Selection.new
event_selection
  .select_track {
      cos_theta 0.93     # |cos(theta)| < 0.93
      Vz        10.0     # |Vz| < 10 cm
      Vr        1.0      # Vr < 1 cm
      nChrp     "==3"    # exactly six charged tracks: 3 positive ...
      nChrn     "==3"    # ... and 3 negative
      nNet      "==0"    # net charge zero
  }
  .pid(method: :probability) {
      prob_cut 0.001
      # high-momentum tracks (p > 1.0 GeV/c) treated as leptons:
      # EMC energy > 0.6 GeV -> electron, otherwise muon
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                      treat_as_electron_if_energy_above: 0.6
      identify :proton, against: [:kaon, :pion]    # p and pbar
      identify :kaon,   against: [:pion, :proton]  # pi/K separation
      identify :pion,   against: [:kaon, :proton]  # pi/K separation
      nprp ">=1"   # at least two protons in total (p and pbar)
      nprm ">=1"
      nkm  "==1"   # one K-
      npip "==1"   # one pi+
      nlp  "==1"   # one l+
      nlm  "==1"   # one l-
  }
  .kinematic_fit([:prp, :prm, :km, :pip, :lp, :lm]) {
      nominal
      vertex_fit([0, 2, 3])    # common vertex for the Lambda_c+ daughters p(0) K-(2) pi+(3)
      constrain_four_momentum  # 4C fit on the six charged particles
      chi2_cut 200
      # anti-Lambda veto: reject M(pbar pi+) > 1.13 GeV/c^2
      invariant_mass_of(:prm, :pip).within(0.0, 1.13)
      # Lambda_c+ reconstruction in the p K- pi+ mass window
      invariant_mass_of(:prp, :km, :pip).within(2.25, 2.32)
  }

my_algorithm
  .note(:vertex_fit_convergence, "the vertex fit of the Lambda_c+ daughters (p K- pi+) is required to converge before the 4C kinematic fit; enforced through the vertex_fit constraint inside the kinematic_fit block")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Execute the algorithm on real data, inclusive MC and the signal exclusive MC
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal])