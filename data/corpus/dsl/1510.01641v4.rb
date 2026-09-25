### Dataset description ###
jpsi_data  = DatasetManager.real_data.find("708_3097")     # 225.0M J/psi (2009)
jpsi_incMC = DatasetManager.inclusive_mc.find("708_3097")

# Signal MC: J/psi -> gamma A0, A0 -> mu+ mu-
decay_card_signal = <<~DECAYCARD
  Alias A0 gamma

  Decay J/psi
  1.0000 gamma A0                       PHSP;
  Enddecay

  Decay A0
  1.0000 mu+ mu-                        PHSP;
  Enddecay

  End
DECAYCARD

exMC_signal = DatasetManager.create_exclusive_mc do |config|
  config.sample_name = "Jpsi_gamma_A0_to_mumu"
  config.related_dataset = jpsi_data
  config.events = 200000
  config.decay_card = decay_card_signal
  config.cross_section = :default
end
exMC_signal.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event Selection ###
alg = Algorithm.new("JpsiGammaA0Mumu")
alg.set_header(["JpsiGammaA0MumuAlg/JpsiGammaA0Mumu.h"])
   .set_constant({"ECMS" => [:double, 3.097]})

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
      energyThreshold_b 0.025
      energyThreshold_e 0.050
      tdc_emc_start     0
      tdc_emc_end       14
      angle_to_track    10.0
      nGam              ">=1"
    }
   .pid(method: :probability) {
      prob_cut 0.001
      identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                     treat_as_electron_if_energy_above: 0.6
      nlp ">=1"
      nlm ">=1"
    }
   .kinematic_fit([:mup, :mum, :gamma]) {
      nominal
      constrain_four_momentum
      chi2_cut 200
    }

alg.note(:electron_veto,
         "Suppress electron contamination by requiring E_cal(mu)/p < 0.9c on both charged tracks.")
   .note(:muon_pid,
         "Assign muon mass hypothesis to both tracks; require at least one track " \
         "identified as muon via MUC PID: 0.1 < E_cal(mu) < 0.3 GeV, |Delta t^{TOF}| < 0.26 ns, " \
         "and MUC penetration depth > (-40 + 70*p/GeV) cm for 0.5 <= p <= 1.1 GeV/c or 40 cm for p > 1.1 GeV/c.")
   .note(:photon_angle_mass_dependent,
         "Photon angle to nearest extrapolated track > 20 deg for m(A0) <= 0.3 GeV, > 10 deg for m(A0) > 0.3 GeV.")
   .note(:common_vertex,
         "The two muon candidates are fitted to a common vertex to form the A0 candidate.")
   .note(:chi2_4C_cut,
         "For multiple gamma mu+ mu- candidates, keep the one with minimum 4C chi2; require chi2_4C < 40.")
   .with_decay_card(decay_card_signal)
   .apply(sel)

alg.execute_on([jpsi_data, jpsi_incMC, exMC_signal])
