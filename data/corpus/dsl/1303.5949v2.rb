# arXiv:1303.5949v2 — Observation of a charged charmoniumlike structure in
# e+e- -> pi+ pi- J/psi at sqrt(s) = 4.260 GeV (BESIII, 525 pb^-1)

### Dataset description ###
data_4260  = DatasetManager.real_data.find("703_4260")      # Real data at 4.260 GeV
incMC_4260 = DatasetManager.inclusive_mc.find("703_4260")   # Corresponding inclusive MC

# Decay card for e+e- -> pi+ pi- J/psi, J/psi -> e+ e-
decay_card_signal_ee = <<~DECAYCARD
    Decay psi(4260)
    1  pi+  pi-  J/psi    PHSP;
    Enddecay

    Decay J/psi
    1  e+  e-    PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Decay card for e+e- -> pi+ pi- J/psi, J/psi -> mu+ mu-
decay_card_signal_mumu = <<~DECAYCARD
    Decay psi(4260)
    1  pi+  pi-  J/psi    PHSP;
    Enddecay

    Decay J/psi
    1  mu+  mu-    PHOTOS VLL;
    Enddecay

    End
DECAYCARD

# Exclusive MC for both J/psi -> l+ l- modes
exMC_signal_ee = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_4260_pipijpsi_ee"
    config.related_dataset = data_4260
    config.events          = 100000
    config.decay_card      = decay_card_signal_ee
    config.cross_section   = :default
end

exMC_signal_mumu = DatasetManager.create_exclusive_mc do |config|
    config.sample_name     = "exmc_4260_pipijpsi_mumu"
    config.related_dataset = data_4260
    config.events          = 100000
    config.decay_card      = decay_card_signal_mumu
    config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "pipiJpsiZc"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 4.26]})
            .set_alias({"std::vector<double>" => "Vdouble"})

# Common selection for the e+e- and mu+mu- modes: 4 charged tracks, net charge
# zero, two of them identified as pions and two as a high-momentum lepton pair.
event_selection = Selection.new
event_selection.select_track {
                    cos_theta       0.93      # |cos(theta)| < 0.93
                    Vz              10.0      # |Vz| < 10 cm along the beam direction
                    Vr              1.0       # Vr < 1 cm in the plane perpendicular to the beam
                    nChrp           "==2"     # exactly 2 positive tracks
                    nChrn           "==2"     # exactly 2 negative tracks
                    nNet            "==0"     # net charge zero
                }
               .pid(method: :probability) {
                    # Pions and leptons are kinematically well separated: tracks with
                    # p > 1.0 GeV/c are treated as leptons, the rest as pions. The EMC
                    # deposit separates electrons (> 1.1 GeV) from muons (< 0.35 GeV).
                    identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                   treat_as_electron_if_energy_above: 0.6
                    identify :pion, against: [:kaon, :proton]
                    npip   "==1"      # one pion+ (the lower-momentum pair)
                    npim   "==1"      # one pion-
                    nlp    "==1"      # one lepton+
                    nlm    "==1"      # one lepton-
               }
               # 4C kinematic fit of the lepton pair and the two pions to the total
               # initial four-momentum of the colliding beams
               .kinematic_fit([:pip, :pim, :lp, :lm]) {
                    nominal
                    constrain_four_momentum   # 4C energy-momentum constraint
                    chi2_cut 60               # chi^2 < 60
               }

# Both J/psi -> e+e- and J/psi -> mu+mu- share the same final-state topology and
# selection, so a single algorithm instance is used.
my_algorithm.with_decay_card(decay_card_signal_ee).apply(event_selection)
root_files = my_algorithm.execute_on([data_4260, incMC_4260,
                                      exMC_signal_ee, exMC_signal_mumu])
