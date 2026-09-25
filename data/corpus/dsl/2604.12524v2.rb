### Dataset description ###
psip_data  = DatasetManager.real_data.find("709_3686")       # psi(2S) real data
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")    # inclusive MC for psi(2S)

# Decay card: psi(2S) -> gamma chi_c1, chi_c1 -> pi+ pi- eta', eta' -> gamma pi+ pi-
decay_card_etap_gpipi = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1                          P2GC1;
    Enddecay
    Decay chi_c1
    1.000 pi+ pi- eta'                          PHSP;
    Enddecay
    Decay eta'
    1.000 pi+ pi- gamma                         ETA_DALITZ;
    Enddecay
    End
DECAYCARD

# Decay card: psi(2S) -> gamma chi_c1, chi_c1 -> pi+ pi- eta', eta' -> pi+ pi- eta, eta -> gamma gamma
decay_card_etap_pipieta = <<~DECAYCARD
    Decay psi(2S)
    1.000 gamma chi_c1                          P2GC1;
    Enddecay
    Decay chi_c1
    1.000 pi+ pi- eta'                          PHSP;
    Enddecay
    Decay eta'
    1.000 pi+ pi- eta                           PHSP;
    Enddecay
    Decay eta
    1.000 gamma gamma                           PHSP;
    Enddecay
    End
DECAYCARD

# Exclusive MC samples
exMC_etap_gpipi = DatasetManager.create_exclusive_mc do |config|
    config.sample_name    = "exclusive_chic1_pipi_etap_gpipi"
    config.related_dataset = psip_data
    config.events         = 500000
    config.decay_card     = decay_card_etap_gpipi
    config.cross_section  = :default
end

exMC_etap_pipieta = DatasetManager.create_exclusive_mc do |config|
    config.sample_name    = "exclusive_chic1_pipi_etap_pipieta"
    config.related_dataset = psip_data
    config.events         = 500000
    config.decay_card     = decay_card_etap_pipieta
    config.cross_section  = :default
end

### Event selection: channel 1 - eta' -> gamma pi+ pi- ###
alg_name_ch1 = "ChiC1EtapGPiPi"
alg_ch1 = Algorithm.new(alg_name_ch1)
alg_ch1.set_header(["#{alg_name_ch1}Alg/#{alg_name_ch1}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})

sel_ch1 = Selection.new
sel_ch1.select_track {
          cos_theta 0.93            # |cos(theta)| < 0.93
          Vz        10.0            # |Vz|  < 10 cm
          Vr        10.0            # |Vxy| < 1 cm
          nChrp     "==2"           # exactly 2 positive tracks
          nChrn     "==2"           # exactly 2 negative tracks
          nNet      "==0"           # zero net charge
        }
       .select_photon {
          tdc_emc_start     0
          tdc_emc_end       14      # EMC timing in [0, 700] ns
          angle_to_track    10.0    # opening angle to charged track > 10 deg
          energyThreshold_b 0.100   # >= 100 MeV in barrel
          energyThreshold_e 0.100   # >= 100 MeV in endcap
          nGam              ">=2"   # at least 2 photons
        }
       # All charged tracks assumed to be pions (no explicit PID applied).
       .assign({:chrgp => :pip, :chrgn => :pim})
       # Nominal 4C fit for psi(2S) -> gamma gamma pi+ pi- pi+ pi- (radiative + eta'->gamma pi+pi-)
       .kinematic_fit([:gamma, :gamma, :pip, :pip, :pim, :pim]) {
          nominal
          constrain_four_momentum
          chi2_cut 200              # loose BOSS cut; tight chi2_4C < 40 applied in ROOT
        }

alg_ch1.note(:non_etap_veto,
             "Require exactly one combination with |M(gamma pi+ pi-) - m_eta'| < 45 MeV/c^2 " \
             "and M(gamma pi+ pi- pi+ pi-) > 3.4 GeV/c^2; final eta' window " \
             "|M(gamma pi+ pi-) - m_eta'| < 15 MeV/c^2; chi_c1 window " \
             "|M(gamma pi+ pi- pi+ pi-) - m_chi_c1| < 15 MeV/c^2. " \
             "Require M(pi+ pi-) from eta' > 0.60 GeV/c^2 to suppress non-eta' background. " \
             "Veto pi0/eta contamination: |M(gamma gamma) - m_pi0| < 17 MeV/c^2 " \
             "and |M(gamma pi+ pi-) - m_eta| < 22 MeV/c^2 rejected. " \
             "J/psi vetoes: |M(2(pi+pi-)) - m_Jpsi| > 28 MeV/c^2 and |M_recoil(pi+ pi-) - m_Jpsi| > 6 MeV/c^2.")

alg_ch1.with_decay_card(decay_card_etap_gpipi).apply(sel_ch1)
root_files_ch1 = alg_ch1.execute_on([psip_data, psip_incMC, exMC_etap_gpipi])

### Event selection: channel 2 - eta' -> pi+ pi- eta (eta -> gamma gamma) ###
alg_name_ch2 = "ChiC1EtapPiPiEta"
alg_ch2 = Algorithm.new(alg_name_ch2)
alg_ch2.set_header(["#{alg_name_ch2}Alg/#{alg_name_ch2}.h"])
       .set_constant({"ECMS" => [:double, 3.686]})

sel_ch2 = Selection.new
sel_ch2.select_track {
          cos_theta 0.93
          Vz        10.0
          Vr        10.0
          nChrp     "==2"
          nChrn     "==2"
          nNet      "==0"
        }
       .select_photon {
          tdc_emc_start     0
          tdc_emc_end       14
          angle_to_track    10.0
          energyThreshold_b 0.100
          energyThreshold_e 0.100
          nGam              ">=3"   # radiative photon + 2 from eta -> gamma gamma
        }
       .assign({:chrgp => :pip, :chrgn => :pim})
       # 4C fit under psi(2S) -> gamma gamma gamma pi+ pi- pi+ pi-
       .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pip, :pim, :pim]) {
          constrain_four_momentum
          chi2_cut 200
        }
       # Nominal 5C fit: additional mass constraint on gamma-gamma pair to eta
       .kinematic_fit([:gamma, :gamma, :gamma, :pip, :pip, :pim, :pim]) {
          nominal
          constrain_four_momentum
          invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
          chi2_cut 200              # tight chi2_5C < 40 applied in ROOT
        }

alg_ch2.note(:non_etap_veto,
             "Require exactly one combination with |M(pi+ pi- gamma gamma) - m_eta'| < 21 MeV/c^2 " \
             "and M(gamma gamma pi+ pi- pi+ pi-) > 3.4 GeV/c^2. " \
             "Final eta' window |M(pi+ pi- eta) - m_eta'| < 7 MeV/c^2; " \
             "chi_c1 window |M(eta pi+ pi- pi+ pi-) - m_chi_c1| < 15 MeV/c^2.")

alg_ch2.with_decay_card(decay_card_etap_pipieta).apply(sel_ch2)
root_files_ch2 = alg_ch2.execute_on([psip_data, psip_incMC, exMC_etap_pipieta])
