### Dataset preparation ###
# X(3872) analysis at sqrt(s)=4.178-4.278 GeV, total 9.0 fb-1
# Multiple energy scan points under BOSS 703
data_4180 = DatasetManager.real_data.find("703_4180")
data_4230 = DatasetManager.real_data.find("703_4230")
data_4260 = DatasetManager.real_data.find("703_4260")

scan_points = [data_4180, data_4230, data_4260]

# Decay card: e+e- -> gamma X(3872) -> gamma J/psi, J/psi -> l+l-
decay_card_gJpsi = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma X_13872   VSP_PWAVE;
  Enddecay
  Decay X_13872
  1.0000 gamma J/psi      VSP_PWAVE;
  Enddecay
  Decay J/psi
  1.0000 e+ e-            PHOTOS VLL;
  Enddecay
  End
DECAYCARD

# Exclusive MC for gJpsi signal (multi-energy)
exMC_gJpsi = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_X3872_gJpsi"
  config.events        = 200_000
  config.decay_card    = decay_card_gJpsi
  config.cross_section = :default
end

### Algorithm 1: X(3872) -> gamma J/psi (l+l-) ###
alg_gJpsi = Algorithm.new("X3872_gJpsi")
alg_gJpsi.set_header(["X3872_gJpsiAlg/X3872_gJpsi.h"])
          .set_constant({"ECMS" => [:double, 4.260]})

sel_gJpsi = Selection.new
sel_gJpsi.select_track {
            cos_theta 0.93
            Vz 10.0
            Vr 1.0
            nChrp ">=1"
            nChrn ">=1"
            nNet "==0"
          }
          .select_photon {
            tdc_emc_start 0
            tdc_emc_end 14
            energyThreshold_b 0.025
            energyThreshold_e 0.050
            angle_to_track 10.0
            nGam ">=2"
          }
          .pid(method: :probability) {
            prob_cut 0.001
            identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                           treat_as_electron_if_energy_above: 0.6
            nlp ">=1"
            nlm ">=1"
          }
          # Remove leptons from charged lists; remainder assigned as pions
          .remove([:lp <= :chrgp])
          .remove([:lm <= :chrgn])
          .assign({chrgp: :pip, chrgn: :pim})
          # For_each: store highest-energy photon as gH, lower as gL
          .for_each(:gamma) {
            define(:egamma) { "energy" }
            best { maximize { "egamma" } }
            store(:gH_index)
          }
          # 4C kinematic fit: l+ l- gamma_H gamma_L
          .kinematic_fit([:lp, :lm, :gamma, :gamma]) {
            nominal
            constrain_four_momentum
            chi2_cut 40
          }

alg_gJpsi
  .note(:photon_energy_ordering, "Radiative photon with larger energy after kinematic fit denoted gamma_H, the other gamma_L; J/psi mass window |M(l+l-)-m_J/psi|<0.02 GeV/c2 applied; pi0/eta veto on M(gamma_L gamma_H); chi_c veto on |M(gamma_L J/psi)-m_chi_c|>0.02 GeV/c2")
  .note(:background_veto, "Bhabha suppression via |cos(theta_gamma)| in [-0.7,0.7] for J/psi->e+e- mode")
  .with_decay_card(decay_card_gJpsi)
  .apply(sel_gJpsi)

# Decay card: X(3872) -> D*0 D0bar, D*0 -> gamma D0, D0 -> K- pi+
decay_card_Dst0D0 = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma X_13872   VSP_PWAVE;
  Enddecay
  Decay X_13872
  1.0000 D*0 anti-D0      PHSP;
  Enddecay
  Decay D*0
  1.0000 gamma D0          VSP_PWAVE;
  Enddecay
  Decay D0
  1.0000 K- pi+           PHSP;
  Enddecay
  Decay anti-D0
  1.0000 K+ pi-           PHSP;
  Enddecay
  End
DECAYCARD

exMC_Dst0D0 = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_X3872_Dst0D0g"
  config.events        = 200_000
  config.decay_card    = decay_card_Dst0D0
  config.cross_section = :default
end

### Algorithm 2: X(3872) -> D*0 D0bar (D*0->gamma D0; D0->K-pi+) ###
alg_Dst0D0 = Algorithm.new("X3872_Dst0D0")
alg_Dst0D0.set_header(["X3872_Dst0D0Alg/X3872_Dst0D0.h"])
           .set_constant({"ECMS" => [:double, 4.260]})

sel_Dst0D0 = Selection.new
sel_Dst0D0.select_track {
             cos_theta 0.93
             Vz 10.0
             Vr 1.0
             nChrp ">=2"
             nChrn ">=2"
             nNet "==0"
           }
           .select_photon {
             tdc_emc_start 0
             tdc_emc_end 14
             energyThreshold_b 0.025
             energyThreshold_e 0.050
             angle_to_track 10.0
             nGam ">=2"
           }
           .pid(method: :probability) {
             prob_cut 0.001
             identify :kaon, against: [:pion]
             identify :pion, against: [:kaon]
             nkp ">=1"
             nkm ">=1"
             npip ">=1"
             npim ">=1"
           }
           # Reconstruct pi0 from photon pairs (for D*0->pi0 D0 mode)
           .kalman_kinematic_fit([:gamma, :gamma]) {
             invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
             chi2_cut 25
             npi0 ">=0"
           }
           # 4C kinematic fit: K+ pi- K- pi+ gamma gamma
           .kinematic_fit([:kp, :pim, :km, :pip, :gamma, :gamma]) {
             nominal
             constrain_four_momentum
             chi2_cut 60
           }

alg_Dst0D0
  .note(:background_veto, "pi0/eta veto on M(gamma_L gamma_H); D*0 mass window M(gamma_L D0) in [m_D*0-0.006, m_D*0+0.006] for gamma mode, M(pi0 D0) in [m_D*0-0.004, m_D*0+0.004] for pi0 mode; D meson mass window applied; kinematic fit combination with smallest chi2 selected among multi-photon ambiguities")
  .note(:tag_mode_unavailable, "D0 also reconstructed via K-pi+pi0 and K-pi+pi+pi- modes; D*0 also reconstructed via pi0 D0 mode; separate decay cards needed for each sub-mode combination")
  .with_decay_card(decay_card_Dst0D0)
  .apply(sel_Dst0D0)

# Decay card: X(3872) -> pi+ pi- J/psi (normalization mode)
decay_card_pipiJpsi = <<~DECAYCARD
  Decay psi(4260)
  1.0000 gamma X_13872   VSP_PWAVE;
  Enddecay
  Decay X_13872
  1.0000 pi+ pi- J/psi    JPIPI;
  Enddecay
  Decay J/psi
  1.0000 e+ e-            PHOTOS VLL;
  Enddecay
  End
DECAYCARD

exMC_pipiJpsi = DatasetManager.create_exclusive_mc_for(scan_points) do |config|
  config.sample_name   = "sig_X3872_pipiJpsi"
  config.events        = 200_000
  config.decay_card    = decay_card_pipiJpsi
  config.cross_section = :default
end

### Algorithm 3: X(3872) -> pi+ pi- J/psi (normalization) ###
alg_pipiJpsi = Algorithm.new("X3872_pipiJpsi")
alg_pipiJpsi.set_header(["X3872_pipiJpsiAlg/X3872_pipiJpsi.h"])
             .set_constant({"ECMS" => [:double, 4.260]})

sel_pipiJpsi = Selection.new
sel_pipiJpsi.select_track {
               cos_theta 0.93
               Vz 10.0
               Vr 1.0
               nChrp ">=2"
               nChrn ">=2"
               nNet "==0"
             }
             .select_photon {
               tdc_emc_start 0
               tdc_emc_end 14
               energyThreshold_b 0.025
               energyThreshold_e 0.050
               angle_to_track 10.0
               nGam ">=1"
             }
             .pid(method: :probability) {
               prob_cut 0.001
               identify_high_momentum_leptons treat_as_lepton_if_momentum_above: 1.0,
                                              treat_as_electron_if_energy_above: 0.6
               identify :pion, against: [:kaon]
               nlp ">=1"
               nlm ">=1"
               npip ">=1"
               npim ">=1"
             }
             .kinematic_fit([:pip, :pim, :lp, :lm, :gamma]) {
               nominal
               constrain_four_momentum
               chi2_cut 200
             }

alg_pipiJpsi
  .note(:background_veto, "J/psi mass window |M(l+l-)-m_J/psi|<0.02 GeV/c2; chi2 of kinematic fit < 40 for X(3872)->gJ/psi mode; pi0/eta veto applied")
  .with_decay_card(decay_card_pipiJpsi)
  .apply(sel_pipiJpsi)

# Execute all algorithms
alg_gJpsi.execute_on(scan_points)
alg_Dst0D0.execute_on(scan_points)
alg_pipiJpsi.execute_on(scan_points)