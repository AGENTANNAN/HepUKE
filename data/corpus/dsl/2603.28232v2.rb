# Lambda_c+ -> n pi+ eta, BESIII single-tag with missing neutron
# Two Lambda_c tag modes are declared; per Rule T1, the two eta reconstruction channels
# (eta -> gamma gamma and eta -> pi+ pi- pi0) use separate TagAnalysis algorithms.
# Energy points 4.600 - 4.843 GeV; DTagTool has 11 anti-Lc- hadronic tag modes.

### Datasets ###
data_4600 = DatasetManager.real_data.find("703_4600")
data_4610 = DatasetManager.real_data.find("706_4610")
data_4620 = DatasetManager.real_data.find("706_4620")
data_4640 = DatasetManager.real_data.find("706_4640")
data_4660 = DatasetManager.real_data.find("706_4660")
data_4680 = DatasetManager.real_data.find("706_4680")
data_4700 = DatasetManager.real_data.find("706_4700")
data_4740 = DatasetManager.real_data.find("707_4740")
data_4750 = DatasetManager.real_data.find("707_4750")
data_4780 = DatasetManager.real_data.find("707_4780")
data_4840 = DatasetManager.real_data.find("707_4840")

incMC_4600 = DatasetManager.inclusive_mc.find("703_4600")
incMC_4620 = DatasetManager.inclusive_mc.find("706_4620")
incMC_4640 = DatasetManager.inclusive_mc.find("706_4640")
incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")
incMC_4740 = DatasetManager.inclusive_mc.find("707_4740")
incMC_4750 = DatasetManager.inclusive_mc.find("707_4750")
incMC_4780 = DatasetManager.inclusive_mc.find("707_4780")
incMC_4840 = DatasetManager.inclusive_mc.find("707_4840")

data_points     = [data_4600, data_4610, data_4620, data_4640, data_4660, data_4680,
                   data_4700, data_4740, data_4750, data_4780, data_4840]
inc_mc_points   = [incMC_4600, incMC_4620, incMC_4640, incMC_4660, incMC_4680,
                   incMC_4700, incMC_4740, incMC_4750, incMC_4780, incMC_4840]

# Decay card for Lc+ -> X pi+ eta (X = n, Lambda, Sigma0) signal; anti-Lc- decays inclusively
decay_card_signal = <<~DECAYCARD
    Decay psi(4260)
    1.0000  Lambda_c+   anti-Lambda_c-              PHSP;
    Enddecay

    Decay Lambda_c+
    0.3333  n0    pi+    eta                        PHSP;
    0.3333  Lambda0    pi+    eta                   PHSP;
    0.3334  Sigma0    pi+    eta                    PHSP;
    Enddecay

    Decay anti-Lambda_c-
    1.0000  anti-p-   K+   pi-                      PHSP;
    Enddecay

    Decay Sigma0
    1.0000  gamma  Lambda0                          PHSP;
    Enddecay

    Decay Lambda0
    1.0000  n0   pi0                                PHSP;
    Enddecay

    Decay eta
    0.3941  gamma gamma                             PHSP;
    0.3251  pi+   pi-   pi0                         PHSP;
    0.2808  pi0   pi0   pi0                         PHSP;
    Enddecay

    Decay pi0
    1.0000  gamma gamma                             PHSP;
    Enddecay

    End
DECAYCARD

exMC_signal_multi = DatasetManager.create_exclusive_mc_for(data_points) do |config|
  config.sample_name   = "Lc_npieta_signal_multiE"
  config.events        = 500000
  config.decay_card    = decay_card_signal
  config.cross_section = :default
end
exMC_signal_multi.each { |m| m.save_to_config(format: :yaml, file_path: 'temp_for_test') }

### Event selection (BOSS) ###
# --- Mode I: eta -> gamma gamma ---
alg_name_gg = "LcNPiEtaGG"
alg_gg = TagAnalysis.new(alg_name_gg)
alg_gg.set_header(["#{alg_name_gg}Alg/#{alg_name_gg}.h"])
      .set_constant({ "ECMS" => [:double, 4.680] })
      .set_alias({ "std::vector<double>" => "Vdouble" })

# Tag side: anti-Lc- reconstructed in 11 hadronic modes
alg_gg.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP, :LambdacPtoKPiP, :LambdacPtoKsPi0P, :LambdacPtoKsPiPiP,
          :LambdacPtoKPiPi0P, :LambdacPtoLambdaPi, :LambdacPtoLambdaPiPi0,
          :LambdacPtoLambdaPiPiPi, :LambdacPtoSigma0Pi, :LambdacPtoSigmaPPi0,
          :LambdacPtoSigmaPPiPi
  t.charm -1
end

# Signal side: one pi+, two photons for eta -> gamma gamma, missing neutron
alg_gg.signal_side do |s|
  s.charged(pip: 1)
  s.require_charge 1
  s.photons 2
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.missing :n0
end

alg_gg.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_gg.note(:signal_track_loose_selection,
            "Signal-side charged track selection uses loose criteria |Vz|<20 cm, |cos(theta)|<0.93 (looser than the tag-side track selection).")
      .note(:eta_gg_mass_window,
            "eta -> gamma gamma candidate requires M(gg) in [0.505, 0.575] GeV/c^2 and 1-C Kalman fit to eta mass with chi^2 < 20.")
      .note(:sigmap_veto_gg,
            "To suppress Lc+ -> Sigma+[-> n pi+] eta peaking background in the eta->gamma gamma channel, veto M_recoil(eta) in (1.13, 1.25) GeV/c^2. Applied in ROOT.")
      .note(:missing_mass_observable,
            "Missing mass M_miss = sqrt(E_miss^2 - |p_miss|^2) is the signal extraction observable; simultaneously extracts X = n, Lambda, Sigma0 components.")
      .note(:deep_learning_selection,
            "After primary event selection, a Transformer-based DNN (Particle Transformer, model-ensemble of 50 networks) classifies events into signal, Lc+ background, non-Lc+ background. Cuts: score(Lc+ bkg) < 0.05 and score(non-Lc+ bkg) < 0.1. This is a post-BOSS ROOT-side operation.")

alg_gg.apply
alg_gg.execute_on(data_points + inc_mc_points + exMC_signal_multi)

# --- Mode II: eta -> pi+ pi- pi0 ---
alg_name_3pi = "LcNPiEta3Pi"
alg_3pi = TagAnalysis.new(alg_name_3pi)
alg_3pi.set_header(["#{alg_name_3pi}Alg/#{alg_name_3pi}.h"])
       .set_constant({ "ECMS" => [:double, 4.680] })
       .set_alias({ "std::vector<double>" => "Vdouble" })

alg_3pi.tag_side(:Lambdac) do |t|
  t.modes :LambdacPtoKsP, :LambdacPtoKPiP, :LambdacPtoKsPi0P, :LambdacPtoKsPiPiP,
          :LambdacPtoKPiPi0P, :LambdacPtoLambdaPi, :LambdacPtoLambdaPiPi0,
          :LambdacPtoLambdaPiPiPi, :LambdacPtoSigma0Pi, :LambdacPtoSigmaPPi0,
          :LambdacPtoSigmaPPiPi
  t.charm -1
end

# Signal side: three charged tracks (pi+ pi+ pi-) plus two photons for pi0, missing neutron
alg_3pi.signal_side do |s|
  s.charged(pip: 2, pim: 1)
  s.require_charge 1
  s.photons 2
  s.min_photon_angle 10.0
  s.min_photon_energy 0.025
  s.missing :n0
end

alg_3pi.fit do |f|
  f.constrain_four_momentum
  f.invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
  f.invariant_mass_of(:pip, :pim, :pi0).constrain_to_nominal_mass_of(:eta)
  f.chi2_cut 200
end

alg_3pi.note(:signal_track_loose_selection,
             "Signal-side charged track selection uses loose criteria |Vz|<20 cm, |cos(theta)|<0.93.")
       .note(:eta_3pi_mass_window,
             "eta -> pi+ pi- pi0 candidate requires M(pi+pi-pi0) in [0.505, 0.575] GeV/c^2 and 1-C fit to eta mass with chi^2 < 20; pi0 -> gamma gamma from M(gg) in [0.115, 0.150] GeV/c^2 with 1-C fit chi^2 < 200.")
       .note(:sigmap_veto_3pi,
             "For eta -> pi+pi-pi0 channel, veto M_recoil(eta) in (1.16, 1.22) GeV/c^2 to suppress Lc+ -> Sigma+[->n pi+] eta. Applied in ROOT.")
       .note(:lambda_3pi_veto,
             "Veto M_recoil(pi+ pi- pi+) in (1.10, 1.13) GeV/c^2 to suppress Lc+ -> Lambda[->n pi0] pi+ pi- pi+. Applied in ROOT.")
       .note(:missing_mass_observable,
             "Missing mass M_miss = sqrt(E_miss^2 - |p_miss|^2) is the signal extraction observable; simultaneously extracts X = n, Lambda, Sigma0.")
       .note(:deep_learning_selection,
             "After primary selection, Transformer-based DNN (ParT ensemble of 50) applied with cuts score(Lc+ bkg) < 0.05, score(non-Lc+ bkg) < 0.1. Handled downstream of BOSS.")

alg_3pi.apply
alg_3pi.execute_on(data_points + inc_mc_points + exMC_signal_multi)
