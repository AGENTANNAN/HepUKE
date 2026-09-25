# Observation of Zc(3900)^0 in e+e- -> pi0 pi0 J/psi, J/psi -> l+ l-
# Multiple energy points from 4.190 to 4.420 GeV

### Dataset description ###
ecms_list = [4.190, 4.210, 4.220, 4.230, 4.245, 4.260, 4.310, 4.360, 4.390, 4.420]
sample_names = ["703_4190","703_4210","703_4220","703_4230","703_4245",
                "703_4260","703_4310","703_4360","703_4390","703_4420"]

xyz_data = sample_names.map { |n| DatasetManager.real_data.find(n) }.compact
xyz_incMC = sample_names.map { |n| DatasetManager.inclusive_mc.find(n) }.compact

# Decay card: e+e- -> pi0 Zc(3900)0 -> pi0 pi0 J/psi
decay_card_pi0pi0Jpsi = <<~DECAYCARD
  Decay psi(4260)
  1.000 pi0 pi0 J/psi                    PHSP;
  Enddecay

  Decay J/psi
  0.500 e+ e-                            PHOTOS PHSP;
  0.500 mu+ mu-                          PHOTOS PHSP;
  Enddecay

  Decay pi0
  1.000 gamma gamma                      PHSP;
  Enddecay

  End
DECAYCARD

# Exclusive MC generated at each energy point
exMC_signal = DatasetManager.create_exclusive_mc_for(xyz_data) do |config|
  config.sample_name    = "pi0pi0Jpsi_signal"
  config.events         = 200000
  config.decay_card     = decay_card_pi0pi0Jpsi
  config.cross_section  = :default
end

### Event selection (BOSS) ###
alg_name = "Pi0Pi0Jpsi"
alg = Algorithm.new(alg_name)
alg.set_header(["#{alg_name}Alg/#{alg_name}.h"])
   .set_constant({"ECMS" => [:double, 4.260]})

sel = Selection.new
sel.select_track {
      cos_theta 0.93
      Vz        10.0
      Vr        1.0
      nChrp     "==1"
      nChrn     "==1"
      nNet      "==0"
    }
    .select_photon {
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    5.0
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      nGam              ">=4"
    }
    .pid(method: :probability) {
      prob_cut 0.001
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6
    }
    # Reconstruct pi0 from photon pairs (1C Kalman fit)
    .kalman_kinematic_fit([:gamma, :gamma]) {
      invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
      chi2_cut 200
      npi0     ">=2"
    }
    # Nominal 4C kinematic fit under e+e- -> pi0 pi0 J/psi hypothesis (J/psi -> l+l-)
    .kinematic_fit([:pi0, :pi0, :lp, :lm]) {
      nominal
      constrain_four_momentum
      chi2_cut 200
    }
    # 7C fit: 4C + two pi0 masses + J/psi mass
    .kinematic_fit([:pi0, :pi0, :lp, :lm]) {
      constrain_four_momentum
      invariant_mass_of(:lp, :lm).constrain_to_nominal_mass_of(:jpsi)
    }

alg.note(:pi0_mass_window,
         "Loose pi0 requirement 100 < M(gamma gamma) < 160 MeV/c^2 for candidate selection; " \
         "tighter 120 < M(gg) < 150 MeV/c^2 used later to reject events with more than one " \
         "extra pi0 pi0 combination.")
   .note(:pi0_pair_selection,
         "Select the pi0 pi0 combination with smallest chi^2 = chi2_1C(pi0a)+chi2_1C(pi0b)+chi2_4C, " \
         "requiring that the two pi0s do not share photons.")
   .note(:electron_id,
         "Electron candidates from J/psi -> e+ e- required to satisfy E/p > 0.7.")
   .note(:muon_id,
         "Muon candidates from J/psi -> mu+ mu- required to satisfy E/p < 0.3 and at least one " \
         "muon with hits in more than six MUC layers.")
   .note(:bhabha_veto,
         "Two-track opening angle < 175 deg for any e+/e- with |cos theta| > 0.5, to suppress " \
         "two-photon and Bhabha backgrounds.")
   .note(:jpsi_mass_window,
         "Dilepton invariant mass 2.95 < M(ll) < 3.2 GeV/c^2 to select J/psi candidates.")

alg.with_decay_card(decay_card_pi0pi0Jpsi).apply(sel)
alg.execute_on(xyz_data + xyz_incMC + exMC_signal)
