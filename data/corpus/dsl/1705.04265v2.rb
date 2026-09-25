# BESIII dark-photon search in untagged ISR, sqrt(s) = 3.773 GeV, 2.93 fb^-1
# (arXiv:1705.04265v2)

### Dataset description ###
data_3773  = DatasetManager.real_data.find("712_3773")     # psi(3770) data at 3.773 GeV
incMC_3773 = DatasetManager.inclusive_mc.find("712_3773")  # corresponding inclusive MC

# Signal decay cards (ISR production of a lepton pair plus an undetected photon)
decay_card_mumu = <<~DECAYCARD
    Decay psi(4260)
    1.000  mu+  mu-  gamma   PHSP;
    Enddecay

    End
DECAYCARD

decay_card_ee = <<~DECAYCARD
    Decay psi(4260)
    1.000  e+  e-  gamma   PHSP;
    Enddecay

    End
DECAYCARD

exMC_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_mumu_gamma_isr"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_mumu
  config.cross_section   = :default
end

exMC_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_3773_ee_gamma_isr"
  config.related_dataset = data_3773
  config.events          = 100000
  config.decay_card      = decay_card_ee
  config.cross_section   = :default
end

### Event selection (BOSS) ###
# --- Channel I: e+e- -> mu+mu- gamma_ISR ---
alg_mumu = Algorithm.new("DarkPhotonMuMu")
alg_mumu.set_header(["DarkPhotonMuMuAlg/DarkPhotonMuMu.h"])
        .set_constant({"ECMS" => [:double, 3.773]})

sel_mumu = Selection.new
sel_mumu.select_track do
          cos_theta 0.921    # 0.4 < theta < pi - 0.4 rad is equivalent to |cos(theta)| < 0.921
          Vz        10.0     # |Vz| < 10 cm along the beam axis
          Vr        1.0      # Vr < 1 cm in the transverse plane
          nChrp     "==1"    # exactly one positively charged track
          nChrn     "==1"    # exactly one negatively charged track
          nNet      "==0"    # net charge zero
        end
        # suppress spiraling tracks: transverse momentum above 300 MeV/c for both tracks
        .for_each(:charged) do
          define(:pt) { (px * px + py * py).__pow__(0.5) }
          where { pt < 0.3 }
          remove
        end
        .pid(method: :probability) do
          # muon identification: P(mu) > P(e); the lepton path uses EMC energy vs MUC depth
          identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                          treat_as_electron_if_energy_above: 0.6
          nlp "==1"
          nlm "==1"
        end
        # 1C kinematic fit to e+e- -> mu+mu- gamma_ISR (missing-photon mass constraint)
        .kinematic_fit([:mup, :mum]) do
          nominal
          miss_track_of :gamma
          constrain_four_momentum
          chi2_cut 20
        end

alg_mumu.with_decay_card(decay_card_mumu).apply(sel_mumu)

# --- Channel II: e+e- -> e+e- gamma_ISR ---
alg_ee = Algorithm.new("DarkPhotonEE")
alg_ee.set_header(["DarkPhotonEEAlg/DarkPhotonEE.h"])
      .set_constant({"ECMS" => [:double, 3.773]})

sel_ee = Selection.new
sel_ee.select_track do
        cos_theta 0.921    # 0.4 < theta < pi - 0.4 rad
        Vz        10.0
        Vr        1.0
        nChrp     "==1"
        nChrn     "==1"
        nNet      "==0"
      end
      # suppress spiraling tracks: pt > 300 MeV/c
      .for_each(:charged) do
        define(:pt) { (px * px + py * py).__pow__(0.5) }
        where { pt < 0.3 }
        remove
      end
      .pid(method: :probability) do
        # electron identification: E/p > 0.8, i.e. lepton with a large EMC energy deposit
        identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                        treat_as_electron_if_energy_above: 0.6
        nlp "==1"
        nlm "==1"
      end
      # stronger requirements than the mu+mu- channel (higher non-ISR background)
      .kinematic_fit([:ep, :em]) do
        nominal
        miss_track_of :gamma
        constrain_four_momentum
        chi2_cut 5
      end

alg_ee.with_decay_card(decay_card_ee).apply(sel_ee)

# BOSS-side procedures that have no formal DSL construct
alg_mumu
  .note(:untagged_isr_selection,
        "untagged ISR: the ISR photon is not required to be reconstructed in the EMC; " \
        "the 1C fit infers the missing photon from energy-momentum conservation")
  .note(:background_veto,
        "the polar angle of the missing photon predicted by the 1C fit is required to be " \
        "theta_gamma < 0.1 or > pi - 0.1 rad (mu+mu-), and < 0.05 or > pi - 0.05 rad (e+e-); " \
        "this is a quantity derived from the kinematic fit and is applied downstream")
  .note(:bhabha_correction,
        "the e+e- gamma event yield is corrected for SM Bhabha scattering in bins of " \
        "m(e+e-): phokhara e+e- annihilation events divided by the sum of annihilation " \
        "(phokhara, mu mass replaced by e mass) and Bhabha (babayaga@nlo) events; " \
        "correction factor 2%-8%")

alg_ee
  .note(:untagged_isr_selection,
        "untagged ISR: the ISR photon is not required to be reconstructed in the EMC; " \
        "the 1C fit infers the missing photon from energy-momentum conservation")
  .note(:background_veto,
        "the polar angle of the missing photon predicted by the 1C fit is required to be " \
        "theta_gamma < 0.1 or > pi - 0.1 rad (mu+mu-), and < 0.05 or > pi - 0.05 rad (e+e-); " \
        "this is a quantity derived from the kinematic fit and is applied downstream")
  .note(:bhabha_correction,
        "the e+e- gamma event yield is corrected for SM Bhabha scattering in bins of " \
        "m(e+e-): phokhara e+e- annihilation events divided by the sum of annihilation " \
        "(phokhara, mu mass replaced by e mass) and Bhabha (babayaga@nlo) events; " \
        "correction factor 2%-8%")

### Execution ###
root_files_mumu = alg_mumu.execute_on([data_3773, incMC_3773, exMC_mumu])
root_files_ee   = alg_ee.execute_on([data_3773, incMC_3773, exMC_ee])
