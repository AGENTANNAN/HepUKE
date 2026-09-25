### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")       # psi(3686) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # inclusive MC for psi(3686)

# Decay cards for the three chi_cJ signal channels
# chi_cJ -> pi+ pi- pi0 pi0, dominant intermediate rho+ rho-.
decay_card_chic0 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c0                          PHSP;
    Enddecay
    Decay chi_c0
    1.000 rho+ rho-                             VSS;
    Enddecay
    Decay rho+
    1.000 pi+ pi0                               VSS;
    Enddecay
    Decay rho-
    1.000 pi- pi0                               VSS;
    Enddecay
    Decay pi0
    1.000 gamma gamma                           PHSP;
    Enddecay
    End
DECAYCARD

decay_card_chic1 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1                          PHSP;
    Enddecay
    Decay chi_c1
    1.000 rho+ rho-                             VSS;
    Enddecay
    Decay rho+
    1.000 pi+ pi0                               VSS;
    Enddecay
    Decay rho-
    1.000 pi- pi0                               VSS;
    Enddecay
    Decay pi0
    1.000 gamma gamma                           PHSP;
    Enddecay
    End
DECAYCARD

decay_card_chic2 = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c2                          PHSP;
    Enddecay
    Decay chi_c2
    1.000 rho+ rho-                             VSS;
    Enddecay
    Decay rho+
    1.000 pi+ pi0                               VSS;
    Enddecay
    Decay rho-
    1.000 pi- pi0                               VSS;
    Enddecay
    Decay pi0
    1.000 gamma gamma                           PHSP;
    Enddecay
    End
DECAYCARD

# Exclusive MC samples for chi_c0, chi_c1, chi_c2
exMC_chic0 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name    = "exclusive_chic0_pipi_pi0pi0"
    config.related_dataset = psip_data
    config.events         = 500000
    config.decay_card     = decay_card_chic0
    config.cross_section  = :default
end

exMC_chic1 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name    = "exclusive_chic1_pipi_pi0pi0"
    config.related_dataset = psip_data
    config.events         = 500000
    config.decay_card     = decay_card_chic1
    config.cross_section  = :default
end

exMC_chic2 = DatasetManager.create_exclusive_mc do |config|
    config.sample_name    = "exclusive_chic2_pipi_pi0pi0"
    config.related_dataset = psip_data
    config.events         = 500000
    config.decay_card     = decay_card_chic2
    config.cross_section  = :default
end

### Event selection (shared across chi_c0, chi_c1, chi_c2) ###
alg_name_c0 = "ChiC0PiPiPi0Pi0"
alg_c0 = Algorithm.new(alg_name_c0)
alg_c0.set_header(["#{alg_name_c0}Alg/#{alg_name_c0}.h"])
      .set_constant({"ECMS" => [:double, 3.686]})

alg_name_c1 = "ChiC1PiPiPi0Pi0"
alg_c1 = Algorithm.new(alg_name_c1)
alg_c1.set_header(["#{alg_name_c1}Alg/#{alg_name_c1}.h"])
      .set_constant({"ECMS" => [:double, 3.686]})

alg_name_c2 = "ChiC2PiPiPi0Pi0"
alg_c2 = Algorithm.new(alg_name_c2)
alg_c2.set_header(["#{alg_name_c2}Alg/#{alg_name_c2}.h"])
      .set_constant({"ECMS" => [:double, 3.686]})

# Selection chain: two charged pions + at least five photons (radiative gamma + 2 pi0 -> 4 gamma)
event_selection = Selection.new
event_selection.select_track {
                  cos_theta 0.93     # |cos(theta)| < 0.93
                  Vz        10.0     # |Vz|  < 10 cm
                  Vr        10.0     # |Vxy| < 1 cm (= 10 mm)
                  nChrp     "==1"    # exactly one positive track
                  nChrn     "==1"    # exactly one negative track
                  nNet      "==0"    # zero net charge
                }
               .select_photon {
                  tdc_emc_start     0        # EMC timing start (0 ns)
                  tdc_emc_end       14       # EMC timing end  (~700 ns; 14 * 50 ns)
                  angle_to_track    10.0     # min opening angle to nearest charged track (deg)
                  energyThreshold_b 0.025    # barrel energy threshold
                  energyThreshold_e 0.050    # endcap energy threshold
                  nGam              ">=5"    # at least 5 photons (1 radiative + 4 from 2 pi0)
                }
               .pid(method: :probability) {
                  prob_cut 0.001
                  identify :pion, against: [:kaon, :proton]  # PID: L(pi) > L(K), L(pi) > L(p)
                  npip     "==1"
                  npim     "==1"
                }
               # Reconstruct pi0 candidates from photon pairs (mass constrained to nominal pi0)
               .kalman_kinematic_fit([:gamma, :gamma]) {
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 25
                  npi0 ">=2"                 # need at least two pi0 candidates
               }
               # Nominal 6C fit: 4C energy-momentum conservation + two pi0 mass constraints
               # (pi0 mass constraint already applied via kalman_kinematic_fit above),
               # participants: radiative photon + two charged pions + two pi0
               .kinematic_fit([:gamma, :pip, :pim, :pi0, :pi0]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 200               # loose BOSS cut; tight chi2_6C < 60 applied in ROOT
               }
               # Competing-hypothesis fit (6 photons -> chi2_more, no nominal, no chi2_cut)
               .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :gamma, :gamma, :pip, :pim]) {
                  constrain_four_momentum
               }
               # Competing-hypothesis fit (4 photons -> chi2_less, no nominal, no chi2_cut)
               .kinematic_fit([:gamma, :gamma, :gamma, :gamma, :pip, :pim]) {
                  constrain_four_momentum
               }

# Capture inexpressible BOSS-side procedures (background vetoes / kinematic-fit helix correction)
alg_c0.note(:background_veto,
            "Veto pi0pi0 J/psi: M_recoil(pi0 pi0) > 3.07 GeV/c^2. " \
            "Veto pi+pi- J/psi: 3.09 < M_recoil(pi+ pi-) < 3.10 GeV/c^2. " \
            "Veto eta/J/psi: 3.09 < M(pi+ pi- pi0_high) < 3.18 GeV/c^2 (higher-momentum pi0). " \
            "Veto chi_cJ -> KsKs peaking bkg: 2D window M(pi+pi-) in (0.445, 0.515) GeV/c^2 " \
            "AND M(pi0 pi0) in (0.475, 0.520) GeV/c^2.")
      .note(:helix_correction,
            "Track helix parameter correction applied to charged tracks before kinematic fit; " \
            "systematic estimated by comparing efficiencies with/without correction.")

alg_c1.note(:background_veto,
            "Same set of vetoes as chi_c0 selection: pi0pi0 J/psi recoil, pi+pi- J/psi recoil, " \
            "eta/J/psi in pi+pi-pi0_high, and 2D KsKs veto.")
      .note(:helix_correction,
            "Track helix parameter correction applied before kinematic fit.")

alg_c2.note(:background_veto,
            "Same set of vetoes as chi_c0 selection: pi0pi0 J/psi recoil, pi+pi- J/psi recoil, " \
            "eta/J/psi in pi+pi-pi0_high, and 2D KsKs veto.")
      .note(:helix_correction,
            "Track helix parameter correction applied before kinematic fit.")

# Render selection into C++ for each chi_cJ decay card
alg_c0.with_decay_card(decay_card_chic0).apply(event_selection)
alg_c1.with_decay_card(decay_card_chic1).apply(event_selection)
alg_c2.with_decay_card(decay_card_chic2).apply(event_selection)

# Execute the algorithms on data, inclusive MC and exclusive MC
root_files_c0 = alg_c0.execute_on([psip_data, psip_incMC, exMC_chic0])
root_files_c1 = alg_c1.execute_on([psip_data, psip_incMC, exMC_chic1])
root_files_c2 = alg_c2.execute_on([psip_data, psip_incMC, exMC_chic2])
