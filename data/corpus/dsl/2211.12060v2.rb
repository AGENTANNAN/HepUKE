# Paper: 2211.12060v2 — Search for hidden-charm tetraquark Z_cs'-
# e+e- → K+ D_s*- D*0 at √s = 4.661, 4.682, 4.699 GeV
# Uses partial-reconstruction technique with two tag methods:
#   Ds- tag: reconstruct K+ + Ds- → K+K-π- (or Ks0 K-), miss D*0
#   D*0 tag: reconstruct K+ + D*0 → D0π0, D0 → K-π+, miss D*s-

### Dataset preparation ###
data_4660 = DatasetManager.real_data.find("706_4660")
data_4680 = DatasetManager.real_data.find("706_4680")
data_4700 = DatasetManager.real_data.find("706_4700")
scan_data = [data_4660, data_4680, data_4700]

incMC_4660 = DatasetManager.inclusive_mc.find("706_4660")
incMC_4680 = DatasetManager.inclusive_mc.find("706_4680")
incMC_4700 = DatasetManager.inclusive_mc.find("706_4700")
scan_incMC = [incMC_4660, incMC_4680, incMC_4700]

# Decay card for non-resonant signal: e+e- → K+ D_s*- D*0
decay_card = <<~DECAYCARD
    Decay psi(4260)
    1.0 K+ D_s*- D*0 PHSP;
    Enddecay

    Decay D_s*-
    1.0 D_s- gamma PHSP;
    Enddecay

    Decay D_s-
    1.0 K+ K- pi- PHSP;
    Enddecay

    Decay D*0
    1.0 D0 pi0 PHSP;
    Enddecay

    Decay D0
    1.0 K- pi+ PHSP;
    Enddecay

    Decay pi0
    1.0 gamma gamma PHSP;
    Enddecay

    End
DECAYCARD

# Exclusive MC for energy scan points
exMC_signal = DatasetManager.create_exclusive_mc_for(scan_data) do |config|
  config.sample_name = "sig_K_Dsstar_Dstar0"
  config.events = 500000
  config.decay_card = decay_card
  config.cross_section = :default
end

# ============================================================
# Algorithm I: Ds- tag method
# Reconstruct K+ (bachelor) and Ds- → K+K-π-; miss D*0
# DecayCard recIDs:
#   0: psi(4260)  [skip]
#   1: K+         [bachelor — reconstruct]
#   2: D_s*-      [miss — gamma from D_s*- is soft]
#   3: D*0        [miss — inferred from recoil]
#   4: gamma      [from D_s*- — miss]
#   5: D_s-       [reconstruct]
#   6: K+         [from D_s- — reconstruct]
#   7: K-         [from D_s- — reconstruct]
#   8: pi-        [from D_s- — reconstruct]
#   9: D0         [miss]
#   10: pi0       [miss]
#   11: gamma     [miss]
#   12: gamma     [miss]
#   13: K-        [from D0 — miss]
#   14: pi+       [from D0 — miss]
# Reconstruct recIDs: 1 (K+ bachelor), 5 (D_s- expands to 6,7,8)
# ============================================================
alg_Ds_tag = Algorithm.new("DsTagDstar0")
alg_Ds_tag.set_header(["DsTagDstar0Alg/DsTagDstar0.h"])
            .set_constant({ "ECMS" => [:double, 4.682] })
            .note(:partial_reconstruction,
              "Ds- tag mode: reconstruct K+ and Ds- → K+K-pi- (or K_S0 K- via K_S0 reconstruction);
               D*0 inferred from RM(K+Ds-). Ds- mass window (1.955, 1.980) GeV.
               For Ds-→K+K-pi-: require M(K+K-) < 1.05 OR 0.850 < M(K+pi-) < 0.930 GeV.
               For Ks0 mode: K_S0 from pi+pi- with M in (0.485, 0.511) GeV, L/σ_L > 2.
               RM(K+Ds-) corrected by substituting reconstructed Ds- mass with known m(Ds-) from PDG;
               require RM(K+Ds-) > 2.11 GeV/c².")

sel_Ds_tag = Selection.new
sel_Ds_tag
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"
    nChrn ">=2"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 10
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon]
  }
  .remove([:pip <= :chrgp, :pim <= :chrgn])
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion]
    nkp ">=1"
    nkm ">=1"
  }
  # Reconstruct Ds- via partial_rec; K+ bachelor and Ds- daughters are in kp/km/pim lists
  .partial_rec([1, 5]) {
    best_combination_by_mass :D_s, 1.9683
    require_recoil_mass 2.11, 5.0
  }

alg_Ds_tag.with_decay_card(decay_card).apply(sel_Ds_tag)

# ============================================================
# Algorithm II: D*0 tag method
# Reconstruct K+ (bachelor) and D*0 → D0π0, D0 → K-π+; miss D*s-
# DecayCard recIDs (same card):
#   0: psi(4260)  [skip]
#   1: K+         [bachelor — reconstruct]
#   2: D_s*-      [miss — inferred from recoil]
#   3: D*0        [reconstruct]
#   4: gamma      [from D_s*- — miss]
#   5: D_s-       [miss]
#   6: K+         [miss]
#   7: K-         [miss]
#   8: pi-        [miss]
#   9: D0         [reconstruct]
#   10: pi0       [reconstruct via kalman_kinematic_fit]
#   11: gamma     [reconstruct as pi0 daughter]
#   12: gamma     [reconstruct as pi0 daughter]
#   13: K-        [from D0 — reconstruct]
#   14: pi+       [from D0 — reconstruct]
# Reconstruct recIDs: 1 (K+), 3 (D*0 expands to 9,10,13,14), with pi0 pre-built via kalman
# ============================================================
alg_Dstar_tag = Algorithm.new("DstarTagDsstar")
alg_Dstar_tag.set_header(["DstarTagDsstarAlg/DstarTagDsstar.h"])
               .set_constant({ "ECMS" => [:double, 4.682] })
               .note(:partial_reconstruction,
                 "D*0 tag mode: reconstruct K+ and D*0 → D0pi0 (D0 → K-pi+);
                  D_s*- inferred from RM(K+D*0). D0 mass window (1.850, 1.880) GeV.
                  D*0 mass window (2.000, 2.014) GeV; best candidate closest to known D*0 mass.
                  RM(K+D*0) corrected with known m(D*0) from PDG;
                  require RM(K+D*0) in (2.102, 2.122) GeV/c².")

sel_Dstar_tag = Selection.new
sel_Dstar_tag
  .select_track {
    cos_theta 0.93
    Vz 10.0
    Vr 1.0
    nChrp ">=2"
    nChrn ">=2"
  }
  .select_photon {
    tdc_emc_start 0
    tdc_emc_end 10
    energyThreshold_b 0.025
    energyThreshold_e 0.050
    angle_to_track 10.0
    nGam ">=2"
  }
  .pid(method: :probability) {
    prob_cut 0.001
    identify :pion, against: [:kaon]
  }
  .remove([:pip <= :chrgp, :pim <= :chrgn])
  .pid(method: :probability) {
    prob_cut 0.001
    identify :kaon, against: [:pion]
    nkp ">=1"
    nkm ">=1"
  }
  # Reconstruct π0 → γγ via Kalman kinematic fit
  .kalman_kinematic_fit([:gamma, :gamma]) {
    invariant_mass_of(:gamma, :gamma).constrain_to_nominal_mass_of(:pi0)
    chi2_cut 25
    npi0 ">=1"
  }
  # Reconstruct K+ bachelor and D*0 chain via partial_rec
  .partial_rec([1, 3, 9, 10, 13, 14]) {
    best_combination_by_mass :D0, 1.8648
    best_combination_by_mass :D_star0, 2.00685
    require_recoil_mass 2.102, 2.122
  }

alg_Dstar_tag.with_decay_card(decay_card).apply(sel_Dstar_tag)

# ============================================================
# Execute all algorithms on the scan datasets
# ============================================================
all_datasets = scan_data + scan_incMC + exMC_signal
alg_Ds_tag.execute_on(all_datasets)
alg_Dstar_tag.execute_on(all_datasets)