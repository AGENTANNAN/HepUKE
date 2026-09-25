### Dataset description ###
# J/psi dataset at 3.097 GeV (BOSS 708) and its inclusive MC
jpsi_data  = DatasetManager.real_data.find("708_3097")
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Decay card for J/psi -> gamma X(1840), X(1840) -> 3(pi+ pi-) (phase space)
decay_card_signal = <<~DECAYCARD
    Decay J/psi
    1.0000 gamma X(1840) PHSP;
    Enddecay

    Decay X(1840)
    1.0000 pi+ pi- pi+ pi- pi+ pi- PHSP;
    Enddecay

    End
DECAYCARD

# Phase-space exclusive MC for the signal process
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_jpsi_gamma_X1840_3pip3pim"
  config.related_dataset = jpsi_data
  config.events          = 100000
  config.decay_card      = decay_card_signal
  config.cross_section   = :default
end

### Event selection (BOSS) ###
alg_name = "JpsiGammaX1840"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.097]})   # CMS energy 3.097 GeV
            .set_alias({"std::vector<double>" => "Vdouble"})

event_selection = Selection.new
    # Exactly six charged tracks, net charge zero; no |cos(theta)|, Vz, or Vr cuts
    .select_track {
        nChrp "==3"   # three pi+
        nChrn "==3"   # three pi-
        nNet  "==0"   # net charge zero
    }
    # At least one photon; no explicit energy threshold
    .select_photon {
        nGam ">=1"
    }
    # No PID performed: treat all remaining positive tracks as pi+, negatives as pi-
    .assign({:chrgp => :pip, :chrgn => :pim})
    # 4C kinematic fit to the gamma 3(pi+pi-) hypothesis, constrained to the CMS 4-momentum
    .kinematic_fit([:gamma, :pip, :pip, :pip, :pim, :pim, :pim]) {
        nominal                  # nominal fit - this is the fit whose four-momenta are kept
        constrain_four_momentum  # constrain total four-momentum to 3.097 GeV
        chi2_cut 30              # chi2 < 30
    }
    # Competing-hypothesis fit (gamma gamma 3(pi+pi-)); no chi2_cut and no nominal so that
    # only its chi2 is stored for the ROOT-level veto chi2_4C(gamma 3pi) < chi2_4C(gamma gamma 3pi)
    .kinematic_fit([:gamma, :gamma, :pip, :pip, :pip, :pim, :pim, :pim]) {
        constrain_four_momentum
    }

# ROOT-stage vetoes that cannot be expressed as BOSS selection steps
my_algorithm
  .note(:ks0_veto, "ROOT-stage veto: reject events with two or more K_S0 candidates built from secondary-vertex fits on pi+pi- pairs satisfying |M(pi+pi-) - m_K_S0| < 5 MeV")
  .note(:pi0_veto, "ROOT-stage veto: for events with at least two photons require |M(gamma gamma) - m_pi0| > 10 MeV")
  .note(:gg_veto, "ROOT-stage veto: require chi2_4C(gamma 3(pi+pi-)) < chi2_4C(gamma gamma 3(pi+pi-)) to suppress the gamma gamma 3(pi+pi-) background (competing gamma-gamma chi2 stored by the second kinematic fit)")
  .with_decay_card(decay_card_signal)
  .apply(event_selection)

# Run the selection on the real data, inclusive MC, and the signal exclusive MC
root_files = my_algorithm.execute_on([jpsi_data, jpsi_incMC, exMC_signal])