### Dataset description ###
# Representative c.m. energy point (4681.92 MeV, BOSS 7.0.6) of the 19-point scan
# spanning 4395.38 - 4950.93 MeV
data_4680  = DatasetManager.real_data.find("706_4680")     # 706-1, 4681.92 MeV
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")  # matching inclusive MC

# Decay card: e+e- -> K_S0 K- pi+ J/psi (+ c.c.), J/psi -> e+ e-
decay_card_ee = <<~DECAYCARD
    Decay psi(4260)
    0.5000 K_S0 K- pi+ J/psi    PHSP;
    0.5000 K_S0 K+ pi- J/psi    PHSP;
    Enddecay

    Decay J/psi
    1.0000 e+ e-    PHOTOS VLL;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-    PHSP;
    Enddecay

    End
DECAYCARD

# Decay card: e+e- -> K_S0 K- pi+ J/psi (+ c.c.), J/psi -> mu+ mu-
decay_card_mumu = <<~DECAYCARD
    Decay psi(4260)
    0.5000 K_S0 K- pi+ J/psi    PHSP;
    0.5000 K_S0 K+ pi- J/psi    PHSP;
    Enddecay

    Decay J/psi
    1.0000 mu+ mu-    PHOTOS VLL;
    Enddecay

    Decay K_S0
    1.0000 pi+ pi-    PHSP;
    Enddecay

    End
DECAYCARD

# Signal exclusive MC: 100k events for each of the two J/psi decay modes
exMC_ee = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4680_kskpipijpsi_ee"
  config.related_dataset = data_4680
  config.events          = 100_000
  config.decay_card      = decay_card_ee
  config.cross_section   = :default
end
exMC_mumu = DatasetManager.create_exclusive_mc do |config|
  config.sample_name     = "exmc_4680_kskpipijpsi_mumu"
  config.related_dataset = data_4680
  config.events          = 100_000
  config.decay_card      = decay_card_mumu
  config.cross_section   = :default
end
exMC_ee.save_to_config(format: :yaml, file_path: 'temp_for_test')
exMC_mumu.save_to_config(format: :yaml, file_path: 'temp_for_test')

### Event selection (BOSS) ###

# ================= Mode I: J/psi -> e+ e- =================
alg_name_ee = "KsKPiJpsiEE"
alg_ee = Algorithm.new(alg_name_ee)
alg_ee.set_header(["#{alg_name_ee}Alg/#{alg_name_ee}.h"])
      .set_constant({"ECMS" => [:double, 4.68192]})   # ECMS = 4.68192 GeV
      .set_alias({"std::vector<double>" => "Vdouble"})
      .note(:efficiency_curve, "the full analysis covers 19 c.m. energies from
        4395.38 to 4950.93 MeV; the BOSS selection is run on the representative
        4681.92 MeV data/inclusive-MC point and the signal MC is generated there,
        so the energy dependence of the selection efficiency is not described by
        the single-point BOSS job (handled with the scan samples in the later stage)")

selection_ee = Selection.new
selection_ee.select_track {          # Charged track selection
              cos_theta 0.93         # |cos(theta)| < 0.93
              Vz        10.0         # |Vz| < 10 cm
              Vr        1.0          # Vr < 1 cm
              nChrp     "==3"        # exactly 3 positive tracks
              nChrn     "==3"        # exactly 3 negative tracks
              nNet      "==0"        # net charge zero
            }
            .pid(method: :probability) {   # PID with the probability method
              prob_cut 0.001
              # high-momentum tracks (p > 1.0 GeV) are treated as leptons;
              # in the e+e- mode a lepton with EMC energy > 1.0 GeV is an electron
              identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                             treat_as_electron_if_energy_above: 1.0
              identify :kaon, against: [:pion, :proton]   # pi/K separation
              nkp "==1"        # one K+
              nkm "==0"        # zero K-
              nlp "==1"        # one l+
              nlm "==1"        # one l-
            }
            # remove the identified K+ and leptons, the rest are treated as pions
            .remove([:kp <= :chrgp, :lp <= :chrgp, :lm <= :chrgn])
            .assign({:chrgp => :pip, :chrgn => :pim})
            .secondary_vertex_fit([:pip, :pim]) {   # reconstruct K_S0 -> pi+ pi-
              build_virtual_particle(:K_S0).by_minimizing_mass_difference
              remove_used_particle_from_candidate_list
            }
            # 4C kinematic fit on K_S0 K+ pi- l+ l-
            .kinematic_fit([:K_S0, :kp, :pim, :lp, :lm]) {
              nominal
              constrain_four_momentum
              chi2_cut 200      # loose chi2 cut; tight cut applied in ROOT
            }

alg_ee.with_decay_card(decay_card_ee).apply(selection_ee)
alg_ee.execute_on([data_4680, incMC_4680, exMC_ee])

# ================= Mode II: J/psi -> mu+ mu- =================
alg_name_mumu = "KsKPiJpsiMuMu"
alg_mumu = Algorithm.new(alg_name_mumu)
alg_mumu.set_header(["#{alg_name_mumu}Alg/#{alg_name_mumu}.h"])
        .set_constant({"ECMS" => [:double, 4.68192]})   # ECMS = 4.68192 GeV
        .set_alias({"std::vector<double>" => "Vdouble"})
        .note(:efficiency_curve, "the full analysis covers 19 c.m. energies from
          4395.38 to 4950.93 MeV; the BOSS selection is run on the representative
          4681.92 MeV data/inclusive-MC point and the signal MC is generated there,
          so the energy dependence of the selection efficiency is not described by
          the single-point BOSS job (handled with the scan samples in the later stage)")

selection_mumu = Selection.new
selection_mumu.select_track {        # Charged track selection
                cos_theta 0.93       # |cos(theta)| < 0.93
                Vz        10.0       # |Vz| < 10 cm
                Vr        1.0        # Vr < 1 cm
                nChrp     "==3"      # exactly 3 positive tracks
                nChrn     "==3"      # exactly 3 negative tracks
                nNet      "==0"      # net charge zero
              }
              .pid(method: :probability) {   # PID with the probability method
                prob_cut 0.001
                # high-momentum tracks (p > 1.0 GeV) are treated as leptons;
                # in the mu+mu- mode a lepton with EMC energy > 0.4 GeV is an electron
                identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                               treat_as_electron_if_energy_above: 0.4
                identify :kaon, against: [:pion, :proton]   # pi/K separation
                nkp "==1"      # one K+
                nkm "==0"      # zero K-
                nlp "==1"      # one l+
                nlm "==1"      # one l-
              }
              # remove the identified K+ and leptons, the rest are treated as pions
              .remove([:kp <= :chrgp, :lp <= :chrgp, :lm <= :chrgn])
              .assign({:chrgp => :pip, :chrgn => :pim})
              .secondary_vertex_fit([:pip, :pim]) {   # reconstruct K_S0 -> pi+ pi-
                build_virtual_particle(:K_S0).by_minimizing_mass_difference
                remove_used_particle_from_candidate_list
              }
              # 4C kinematic fit on K_S0 K+ pi- l+ l-
              .kinematic_fit([:K_S0, :kp, :pim, :lp, :lm]) {
                nominal
                constrain_four_momentum
                chi2_cut 200    # loose chi2 cut; tight cut applied in ROOT
              }
# NOTE: the additional MUC hit-depth sum > 75 cm requirement of the mu+mu- channel
# is a ROOT-stage cut applied after the kinematic fit, hence not part of this BOSS spec.

alg_mumu.with_decay_card(decay_card_mumu).apply(selection_mumu)
alg_mumu.execute_on([data_4680, incMC_4680, exMC_mumu])