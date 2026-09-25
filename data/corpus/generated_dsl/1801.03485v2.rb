### Dataset preparation ###
# Real data and inclusive MC at the psi(2S) energy point (3.686 GeV)
psip_data = DatasetManager.real_data.find("709_3686")
psip_incMC = DatasetManager.inclusive_mc.find("709_3686")

# Signal decay card: psi(2S) -> pi0 h_c, h_c -> pi+ pi- J/psi, J/psi -> l+ l-
decay_card_signal = <<~DECAYCARD
    Decay psi(2S)
    1.000 pi0 h_c  PHSP;
    Enddecay

    Decay h_c
    1.000 pi+ pi- J/psi  PHSP;
    Enddecay

    Decay J/psi
    0.5 e+ e-    PHOTOS VLL;
    0.5 mu+ mu-  PHOTOS VLL;
    Enddecay

    Decay pi0
    1.000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Reference decay card: psi(2S) -> eta J/psi, eta -> pi0 pi+ pi-, J/psi -> l+ l-
decay_card_reference = <<~DECAYCARD
    Decay psi(2S)
    1.000 eta J/psi  PHSP;
    Enddecay

    Decay eta
    1.000 pi0 pi+ pi-  PHSP;
    Enddecay

    Decay J/psi
    0.5 e+ e-    PHOTOS VLL;
    0.5 mu+ mu-  PHOTOS VLL;
    Enddecay

    Decay pi0
    1.000 gamma gamma  PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC: 200k events for the signal mode
exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_psip_pi0hc_ll"
  config.related_dataset = psip_data
  config.events = 200000
  config.decay_card = decay_card_signal
  config.cross_section = :default
end

# Exclusive MC: 500k events for the reference mode
exMC_reference = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "exmc_psip_etajpsi_ll"
  config.related_dataset = psip_data
  config.events = 500000
  config.decay_card = decay_card_reference
  config.cross_section = :default
end

### Event selection (BOSS) ###
alg_name = "Pi0hcToLL"
my_algorithm = Algorithm.new(alg_name)
my_algorithm.set_header(["#{alg_name}Alg/#{alg_name}.h"])
            .set_constant({"ECMS" => [:double, 3.686]})
            .set_alias({"std::vector<double>" => "Vdouble"})

# Signal and reference share the identical final state pi0 pi+ pi- l+ l- and identical
# selection, so a single Algorithm instance with one decay card suffices.
event_selection = Selection.new
event_selection.select_track {          # Charged track selection
                  cos_theta 0.93        # |cos(theta)| < 0.93
                  Vz        10.0        # |Vz| < 10 cm
                  Vr        1.0         # Vr < 1 cm in the transverse plane
                  nChrp     ">=2"       # At least 2 positive tracks
                  nChrn     ">=2"       # At least 2 negative tracks
                  nNet      "==0"       # Net charge zero
                }
               .select_photon {         # Photon selection
                  tdc_emc_start     0
                  tdc_emc_end       14
                  energyThreshold_b 0.025   # 25 MeV in the barrel
                  energyThreshold_e 0.050   # 50 MeV in the endcap
                  nGam              ">=2"   # At least 2 photons
                }
               .pid(method: :probability) {   # Particle identification (probability method)
                  prob_cut 0.001
                  identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                                 treat_as_electron_if_energy_above: 0.6  # p>1.0 GeV/c -> lepton; EMC E>0.6 GeV -> e, else mu
                  nlp ">=1"        # At least one l+
                  nlm ">=1"        # At least one l-
                }
               .remove([:lp <= :chrgp, :lm <= :chrgn])  # Drop identified leptons from the charged lists
               .assign({:chrgp => :pip, :chrgn => :pim})  # Remaining good tracks -> pi+ / pi-
               .kalman_kinematic_fit([:gamma, :gamma]) {  # Reconstruct pi0 from two photons (1-C)
                  invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
                  chi2_cut 25
                  npi0 ">=1"     # At least one pi0 candidate
                }
               # 5C kinematic fit: 4C (CMS four-momentum) plus the pi0 mass constraint
               # (carried by reconstructing pi0 as a single fixed-mass participant)
               .kinematic_fit([:pi0, :pip, :pim, :lp, :lm]) {
                  nominal
                  constrain_four_momentum
                  chi2_cut 60
                }

my_algorithm.with_decay_card(decay_card_signal).apply(event_selection)
# Identical selection applied to real data, inclusive MC, and both exclusive-MC samples
root_files = my_algorithm.execute_on([psip_data, psip_incMC, exMC_signal, exMC_reference])