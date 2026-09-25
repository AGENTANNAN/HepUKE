# Core DSL classes and dependencies will be loaded automatically at execution

### Dataset description ###
# BESIII R-scan data (BOSS release 713) covering 2.100-3.080 GeV at 20 c.m. energies.
# The representative point used here is 2.125 GeV.
rscan_data  = DatasetManager.real_data.find("713_2125")
rscan_incMC = DatasetManager.inclusive_mc.find("713_2125")

# ConExc decay card for e+e- -> K+K-K+K-.
# ConExc mode 16 (K+K-K+K-, phase space) models ISR up to second order and the
# measured sigma_0(sqrt(s)). The DSL auto-detects the literal token "ConExc",
# switches to the no-KKMC simulation template, and injects `Particle vpho <ECMS> 0.0`
# per energy point - so `Particle vpho` is intentionally omitted here.
decay_card_4K = <<~DECAYCARD
    Decay vpho
    1.000 ConExc 16;
    Enddecay
    End
DECAYCARD

# 100k-event exclusive MC for the 4K final state
exMC_4K = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_rscan_4K_conexc16"
  config.related_dataset = rscan_data
  config.events          = 100000
  config.decay_card      = decay_card_4K
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# ---------------- Mode I: e+e- -> K+K-K+K- (4C kinematic fit) ----------------
alg_name_4K = "Rscan4K"
alg_4K = Algorithm.new(alg_name_4K)
alg_4K.set_header(["#{alg_name_4K}Alg/#{alg_name_4K}.h"])
      .set_constant({"ECMS" => [:double, 2.125]})
      .set_alias({"std::vector<double>" => "Vdouble"})

sel_4K = Selection.new
sel_4K.select_track {
         cos_theta 0.93
         Vz        100.0
         Vr        10.0
         nChrp     ">=2"
         nChrn     ">=2"
         nNet      "==0"
       }
       .pid(method: :probability) {
         prob_cut 0.001
         identify :kaon, against: [:pion, :proton]   # K+ and K- (charge-conjugation shorthand)
         nkp "==2"
         nkm "==2"
       }
       # QED suppression: identified kaons must have p < 0.8 p_beam (= 0.8 * ECMS/2)
       .remove(:kp) { condition "three_momentum_of(:kp) > 0.8 * ECMS / 2.0" }
       .remove(:km) { condition "three_momentum_of(:km) > 0.8 * ECMS / 2.0" }
       # 4C fit of the four kaons to the CMS four-momentum
       .kinematic_fit([:kp, :kp, :km, :km]) {
         nominal
         constrain_four_momentum
         chi2_cut 200
       }

alg_4K.with_decay_card(decay_card_4K).apply(sel_4K)

# ------- Mode II: e+e- -> phi K+K- (phi -> K+K-), 1C fit with one missing kaon -------
alg_name_phi = "RscanPhiKK"
alg_phi = Algorithm.new(alg_name_phi)
alg_phi.set_header(["#{alg_name_phi}Alg/#{alg_name_phi}.h"])
       .set_constant({"ECMS" => [:double, 2.125]})
       .set_alias({"std::vector<double>" => "Vdouble"})

sel_phi = Selection.new
sel_phi.select_track {
          cos_theta 0.93
          Vz        100.0
          Vr        10.0
          nChrp     ">=1"
          nChrn     ">=1"
        }
        .pid(method: :probability) {
          prob_cut 0.001
          identify :kaon, against: [:pion, :proton]   # K+ and K-
          nkp ">=1"
          nkm ">=1"
        }
        # QED suppression: identified kaons must have p < 0.8 p_beam
        .remove(:kp) { condition "three_momentum_of(:kp) > 0.8 * ECMS / 2.0" }
        .remove(:km) { condition "three_momentum_of(:km) > 0.8 * ECMS / 2.0" }
        # 1C fit: the full K+K-K+K- system with one kaon missing (kaon mass assigned);
        # 4-momentum conservation with one unmeasured massive kaon leaves 1 constraint
        .kinematic_fit([:kp, :kp, :km, :km]) {
          nominal
          miss_track_of :km
          constrain_four_momentum
          chi2_cut 20
        }

alg_phi.with_decay_card(decay_card_4K).apply(sel_phi)

### Execute: both modes run on the same real data, inclusive MC, and 4K exclusive MC ###
root_files_4K  = alg_4K.execute_on([rscan_data, rscan_incMC, exMC_4K])
root_files_phi = alg_phi.execute_on([rscan_data, rscan_incMC, exMC_4K])