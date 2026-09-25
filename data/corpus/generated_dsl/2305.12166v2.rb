# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# Six low-energy R-scan points (BOSS 7.1.3). Sample name = [BOSS]_[CMS energy in MeV].
scan_points = ["713_2309", "713_2386", "713_2396", "713_2500", "713_2644", "713_2646"]
data_points = scan_points.map { |s| DatasetManager.real_data.find(s) }        # real data at each energy
incMC_points = scan_points.map { |s| DatasetManager.inclusive_mc.find(s) }    # matching inclusive MC

# ConExc decay card for e+e- -> Delta++ Delta--  (continuum / R-scan, ISR Born cross section).
# The literal token `ConExc` makes the DSL use the no-KKMC simulation template and inject
# `Particle vpho <ECMS> 0.0` per scan point (so no Particle vpho line is written here).
decay_card_signal = <<~DECAYCARD
    Decay vpho
    1.0000 Delta++ anti-Delta-- ConExc;
    Enddecay

    Decay Delta++
    1.0000 p+ pi+ PHSP;
    Enddecay

    Decay anti-Delta--
    1.0000 anti-p- pi- PHSP;
    Enddecay

    End
DECAYCARD

# 500k ConExc signal events at each of the six energy points (same signal MC / decay card)
exMCs_signal = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name = "exmc_deltadeltabar"
  config.events = 500_000
  config.decay_card = decay_card_signal
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "DeltaDeltaBar"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
# ECMS is injected per energy point by the ConExc template, so no fixed ECMS constant is set.

event_selection = Selection.new
  .select_track {          # exactly four charged tracks, net charge zero
    cos_theta 0.93         # |cos(theta)| < 0.93
    Vz 10.0                # |Vz| < 10 cm
    Vr 1.0                 # Vr < 1 cm
    nChrp "==2"            # two positive tracks
    nChrn "==2"            # two negative tracks
    nNet "==0"             # net charge zero
  }
  .pid(method: :probability) {   # probability PID
    prob_cut 0.001               # prob > 0.001
    identify :proton, against: [:kaon, :pion]   # p+ and pbar, vs K and pi
    identify :pion, against: [:kaon]            # pi+ and pi-, vs K
    nprp ">=1"                   # one proton candidate
    nprm ">=1"                   # one anti-proton candidate
    npip ">=1"                   # one pi+ candidate
    npim ">=1"                   # one pi- candidate
  }
  # Tracks are taken as p, pbar, pi+, pi- and fitted under the p pbar pi+ pi- hypothesis
  .kinematic_fit([:prp, :prm, :pip, :pim]) {
    nominal                  # nominal 4C fit
    constrain_four_momentum  # four-momentum conservation
    chi2_cut 50              # chi^2 < 50
  }

alg.with_decay_card(decay_card_signal).apply(event_selection)
root_files = alg.execute_on(data_points + incMC_points + exMCs_signal)