# ============================================================================
# Datasets : psi(3686) real data + inclusive MC (BOSS 709 / 3.686 GeV)
# ============================================================================
psip_data  = DatasetManager.real_data.find("709_3686")     # psi(3686) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")  # inclusive MC at 3.686 GeV

# ============================================================================
# Decay cards : psi(3686) -> pi+ pi- J/psi , J/psi -> 4 leptons
# ============================================================================
decay_card_4e = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi+ pi- J/psi PHSP;
    Enddecay
    Decay J/psi
    1.0000 e+ e- e+ e- PHSP;
    Enddecay
    End
DECAYCARD

decay_card_2e2mu = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi+ pi- J/psi PHSP;
    Enddecay
    Decay J/psi
    1.0000 e+ e- mu+ mu- PHSP;
    Enddecay
    End
DECAYCARD

decay_card_4mu = <<~DECAYCARD
    Decay psi(2S)
    1.0000 pi+ pi- J/psi PHSP;
    Enddecay
    Decay J/psi
    1.0000 mu+ mu- mu+ mu- PHSP;
    Enddecay
    End
DECAYCARD

# ============================================================================
# Exclusive MC : 500k events for each of the three J/psi decay modes
# ============================================================================
exMC_4e = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_pipiJpsi_4e"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_4e
  config.cross_section   = :default
end

exMC_2e2mu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_pipiJpsi_2e2mu"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_2e2mu
  config.cross_section   = :default
end

exMC_4mu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3686_pipiJpsi_4mu"
  config.related_dataset = psip_data
  config.events          = 500000
  config.decay_card      = decay_card_4mu
  config.cross_section   = :default
end

# ============================================================================
# Common selection : charged tracks + PID (shared by the three modes)
# ============================================================================
common_selection = Selection.new
    .select_track {                 # charged-track quality + multiplicity
        cos_theta  0.93             # |cos(theta)| < 0.93
        Vz         10.0             # |Vz| < 10 cm
        Vr         1.0              # Vr < 1 cm
        nTot       ">=6"            # at least six charged tracks
        nNet       "==0"            # net charge zero
    }
    .pid(method: :probability) {    # probability PID
        prob_cut 0.001              # probability > 0.001
        # high-momentum tracks (p>1.0 GeV/c) treated as leptons; lepton with EMC E>0.6 GeV -> electron, else muon
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                       treat_as_electron_if_energy_above: 0.6
        identify :pion, against: [:kaon]   # pi/K separation (pi+ and pi-)
    }
    # soft-pion pi+pi- J/psi tag : keep pions with p < 0.45 GeV/c (drop the hard ones)
    .remove(:pip) { condition "three_momentum_of(:pip) >= 0.45" }
    .remove(:pim) { condition "three_momentum_of(:pim) >= 0.45" }

# ============================================================================
# Mode I : psi(3686) -> pi+ pi- J/psi , J/psi -> e+ e- e+ e-
# ============================================================================
alg_name_I = "PsipToPiPiJpsi4e"
alg_I = Algorithm.new(alg_name_I)
alg_I.set_header(["#{alg_name_I}Alg/#{alg_name_I}.h"])
     .set_constant({"ECMS" => [:double, 3.686]})
     .set_alias({"std::vector<double>" => "Vdouble"})
     .note(:helix_correction,
           "simulated (MC) tracks are helix-parameter corrected before the kinematic fit")
     .note(:background_veto,
           "photon-conversion veto: electron pairs are rejected using the vertex radius R_xy and the
            momentum asymmetry dp = (p_e_high - p_e_low) > 1.0 GeV/c")
     .note(:lepton_pid_selection,
           "additional lepton PID beyond the high-momentum-lepton identification: muon candidates
            require prob(mu) > prob(e), prob(K), EMC energy in [0.1, 0.3] GeV and MUC depth
            requirements; electron candidates require prob(e) > prob(mu), prob(K)")

sel_I = common_selection.dup
    .kinematic_fit([:pip, :pim, :ep, :em, :ep, :em]) do
        nominal
        constrain_four_momentum
        # 5-sigma J/psi recoil-mass window (sigma ~ 4.0 MeV for the 4e mode)
        invariant_mass_of(:ep, :em, :ep, :em).within(3.077, 3.117)
        chi2_cut 200
    end

alg_I.with_decay_card(decay_card_4e).apply(sel_I)
alg_I.execute_on([psip_data, psip_incMC, exMC_4e])

# ============================================================================
# Mode II : psi(3686) -> pi+ pi- J/psi , J/psi -> e+ e- mu+ mu-
# ============================================================================
alg_name_II = "PsipToPiPiJpsi2e2mu"
alg_II = Algorithm.new(alg_name_II)
alg_II.set_header(["#{alg_name_II}Alg/#{alg_name_II}.h"])
      .set_constant({"ECMS" => [:double, 3.686]})
      .set_alias({"std::vector<double>" => "Vdouble"})
      .note(:helix_correction,
            "simulated (MC) tracks are helix-parameter corrected before the kinematic fit")
      .note(:background_veto,
            "photon-conversion veto: electron pairs are rejected using the vertex radius R_xy and the
             momentum asymmetry dp = (p_e_high - p_e_low) > 1.0 GeV/c")
      .note(:lepton_pid_selection,
            "additional lepton PID beyond the high-momentum-lepton identification: muon candidates
             require prob(mu) > prob(e), prob(K), EMC energy in [0.1, 0.3] GeV and MUC depth
             requirements; electron candidates require prob(e) > prob(mu), prob(K)")

sel_II = common_selection.dup
    .kinematic_fit([:pip, :pim, :ep, :em, :mup, :mum]) do
        nominal
        constrain_four_momentum
        # 5-sigma J/psi recoil-mass window (sigma ~ 2.8 MeV for the 2e2mu mode)
        invariant_mass_of(:ep, :em, :mup, :mum).within(3.083, 3.111)
        chi2_cut 200
    end
    # competing hypothesis : the two muons interpreted as pions (e+e-pi+pi-).
    # No chi2_cut and no nominal : the competing chi2 is stored for the ROOT-level
    # veto "signal chi2 < competing chi2".
    .assign({:mup => :pip, :mum => :pim})
    .kinematic_fit([:pip, :pim, :ep, :em, :pip, :pim]) do
        constrain_four_momentum
    end

alg_II.with_decay_card(decay_card_2e2mu).apply(sel_II)
alg_II.execute_on([psip_data, psip_incMC, exMC_2e2mu])

# ============================================================================
# Mode III : psi(3686) -> pi+ pi- J/psi , J/psi -> mu+ mu- mu+ mu-
# ============================================================================
alg_name_III = "PsipToPiPiJpsi4mu"
alg_III = Algorithm.new(alg_name_III)
alg_III.set_header(["#{alg_name_III}Alg/#{alg_name_III}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})
       .set_alias({"std::vector<double>" => "Vdouble"})
       .note(:helix_correction,
             "simulated (MC) tracks are helix-parameter corrected before the kinematic fit")
       .note(:background_veto,
             "photon-conversion veto: electron pairs are rejected using the vertex radius R_xy and the
              momentum asymmetry dp = (p_e_high - p_e_low) > 1.0 GeV/c")
       .note(:lepton_pid_selection,
             "additional lepton PID beyond the high-momentum-lepton identification: muon candidates
              require prob(mu) > prob(e), prob(K), EMC energy in [0.1, 0.3] GeV and MUC depth
              requirements; electron candidates require prob(e) > prob(mu), prob(K)")

sel_III = common_selection.dup
    .kinematic_fit([:pip, :pim, :mup, :mum, :mup, :mum]) do
        nominal
        constrain_four_momentum
        # 5-sigma J/psi recoil-mass window (sigma ~ 2.8 MeV, 2e2mu reference)
        invariant_mass_of(:mup, :mum, :mup, :mum).within(3.083, 3.111)
        chi2_cut 200
    end
    # competing hypothesis : all four muons interpreted as pions (pi+pi-pi+pi-).
    # Competing chi2 stored for the ROOT-level veto (signal chi2 < competing chi2).
    # The two peaking pi+pi-pi+pi- background events are handled at ROOT level.
    .assign({:mup => :pip, :mum => :pim})
    .kinematic_fit([:pip, :pim, :pip, :pim, :pip, :pim]) do
        constrain_four_momentum
    end

alg_III.with_decay_card(decay_card_4mu).apply(sel_III)
alg_III.execute_on([psip_data, psip_incMC, exMC_4mu])